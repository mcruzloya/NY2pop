library('R2jags')
library('MCMCvis')
library('scales')
library('Hmisc')

# Relative R0 parametrization as in Eq. (A2) in Shocket et al. 
# (removing N and r).
R0 <- function(a, bc, lf, PDR, ER, pO, EV, pLA, MDR) {
  eps <- 1e-10 # Small constant to prevent denominators being zero 
  # Calculate as log of R0^2 for numerical stability.
  logR02 <- (3 * log(a) + log(bc) - (1 / (PDR * lf + eps)) + log(ER) + log(pO)
             + log(EV) + log(pLA) + log(MDR) + 3 * log(lf))
  return(sqrt(exp(logR02)))
}

# Briere TPC model.
briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax))
  return(result)
}

# Quadratic TPC model.
quad <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(result)
}

# Quadratic TPC model with upper limit at 1 for proportions.
quad_lim <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

# Linear TPC model with upper limit at 15C.
linear_lim <- function(T, m, b) {
  result <- (-m * T + b) * (T < b/m)
  result[T < 15] <- (-m * 15 + b) * (T < b/m) # Set maximum value at 15 degrees
  return(result)
}

## 1. Define common plot parameters and functions.

par(mar=c(5,6.5,4,1)+.1)
temps <- seq(0, 45, 0.1)
suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9"
  

plot_points_errbar <- function(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                               xlim=c(0, 45), ylim=c(0, 1), main="Trait",
                               xlab="Temperature [°C]", ylab="trait",
                               errbar.length=0.06, errbar.width=2,
                               cex.main=1.5, cex.lab=1.3, cex.axis=1.2, cex=1.5) {
  plot(c(22,25,28), data.mean.suf, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim, pch=20, col=suffolk.col, main=main, 
       cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, cex=cex)
  arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf+data.sderr.suf, col=suffolk.col,
         angle=90, length=errbar.length, lwd=errbar.width)
  arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf-data.sderr.suf, col=suffolk.col,
         angle=90, length=errbar.length, lwd=errbar.width)
  
  points(c(22,25,28), data.mean.alb, pch=20, col=albany.col,
         cex=cex)
  arrows(c(22,25,28), data.mean.alb, y1=data.mean.alb+data.sderr.alb, col=albany.col,
         angle=90, length=errbar.length, lwd=errbar.width)
  arrows(c(22,25,28), data.mean.alb, y1=data.mean.alb-data.sderr.alb, col=albany.col,
         angle=90, length=errbar.length, lwd=errbar.width)
}


plot_curves <- function(temps, curves.suf, curves.alb, 
                        suffolk.col="#CC79A7",
                        albany.col="#56B4E9",
                        mean.lwd=2, CI.lwd=1.5) {
  
  meancurve.suf <- apply(curves.suf, 1, mean)
  CI.suf <- apply(curves.suf, 1, quantile, c(0.025, 0.975))
  lines(temps, meancurve.suf, col=suffolk.col, lwd=mean.lwd)
  lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5), lwd=CI.lwd)
  lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5), lwd=CI.lwd)
  
  meancurve.alb <- apply(curves.alb, 1, mean)
  CI.alb <- apply(curves.alb, 1, quantile, c(0.025, 0.975))
  lines(temps, meancurve.alb, col=albany.col, lwd=mean.lwd)
  lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5), lwd=CI.lwd)
  lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5), lwd=CI.lwd)
}

plot_curves_shaded <- function(temps, curves.suf, curves.alb, 
                        suffolk.col="#CC79A7",
                        albany.col="#56B4E9",
                        mean.lwd=2.5, shade.alpha=0.2, scale.factor=1.0) {
  
  meancurve.suf <- scale.factor * apply(curves.suf, 1, mean)
  CI.suf <- scale.factor * apply(curves.suf, 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve.suf, col=suffolk.col, lwd=mean.lwd)
  polygon(c(temps, rev(temps)), c(CI.suf[1,], rev(CI.suf[2,])), 
          col=alpha(suffolk.col, shade.alpha), lty=0)
  
  meancurve.alb <- scale.factor * apply(curves.alb, 1, mean)
  CI.alb <- scale.factor * apply(curves.alb, 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve.alb, col=albany.col, lwd=mean.lwd)
  polygon(c(temps, rev(temps)), c(CI.alb[1,], rev(CI.alb[2,])), 
          col=alpha(albany.col, shade.alpha), lty=0)
}


### Biting rate

# Load data
data <- read.csv("biting_rate_lf.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

# Quadratic model

# Load MCMC chains
a.suf.chains <- readRDS("jagsout_a_suffolk.RDS")
a.alb.chains <- readRDS("jagsout_a_albany.RDS")

a.suf.chains
a.alb.chains

#View(a.suf.chains$BUGSoutput$summary[c(2,1,3,4), c(1,3,7,8,9)])
#View(a.alb.chains$BUGSoutput$summary[c(2,1,3,4), c(1,3,7,8,9)])

# Only keep chains for TPC model parameters.
a.suf.chains <- MCMCchains(a.suf.chains, params=c("Tmin", "Tmax", "c"))
a.alb.chains <- MCMCchains(a.alb.chains, params=c("Tmin", "Tmax", "c"))

# Calculate trait values for each sample at all temperatures.
a.suf.curves <-  apply(a.suf.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))
a.alb.curves <-  apply(a.alb.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))

# Calculate mean and std err.
data.mean.suf <- data.suffolk$biting.rate

data.sderr.suf <- sqrt((data.suffolk$biting.rate * (1 - data.suffolk$biting.rate))
                       / data.suffolk$mosquito.days)

data.mean.alb <- data.albany$biting.rate

data.sderr.alb <- sqrt((data.albany$biting.rate * (1 - data.albany$biting.rate))
                       / data.albany$mosquito.days)

# Plot data

plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 0.07), main="Biting rate",
                   ylab=expression(paste("Biting rate ", 
                   bgroup("[", frac("bites", "mosquito day"), "]"))))
plot_curves_shaded(temps, a.suf.curves, a.alb.curves)

### Lifespan

# Load data
data <- read.csv("lf_allpop.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

lf.suf.chains <- readRDS("jagsout_lf_suffolk.RDS")
lf.alb.chains <- readRDS("jagsout_lf_albany.RDS")

#View(lf.suf.chains$BUGSoutput$summary[c(5,1,6, 3), c(1,3,7,8,9)])
#View(lf.alb.chains$BUGSoutput$summary[c(5,1,6, 3), c(1,3,7,8,9)])

lf.suf.chains <- MCMCchains(lf.suf.chains, params=c("m", "b"))
lf.alb.chains <- MCMCchains(lf.alb.chains, params=c("m", "b"))

lf.suf.curves <- apply(lf.suf.chains, 1, function(x) linear_lim(temps, x[1], x[2]))
lf.alb.curves <- apply(lf.alb.chains, 1, function(x) linear_lim(temps, x[1], x[2]))

# Calculate mean and SD.
data.mean.suf <- c(mean(data.suffolk$lifespan[data.suffolk$temperature == 22]),
                   mean(data.suffolk$lifespan[data.suffolk$temperature == 25]),
                   mean(data.suffolk$lifespan[data.suffolk$temperature == 28]))
data.sd.suf <- c(sd(data.suffolk$lifespan[data.suffolk$temperature == 22]),
                 sd(data.suffolk$lifespan[data.suffolk$temperature == 25]),
                 sd(data.suffolk$lifespan[data.suffolk$temperature == 28]))
data.sderr.suf <- data.sd.suf / sqrt(c(length(data.suffolk$temperature == 22),
                                       length(data.suffolk$temperature == 25),
                                       length(data.suffolk$temperature == 28)))

data.mean.alb <- c(mean(data.albany$lifespan[data.albany$temperature == 22]),
                   mean(data.albany$lifespan[data.albany$temperature == 25]),
                   mean(data.albany$lifespan[data.albany$temperature == 28]))
data.sd.alb <- c(sd(data.albany$lifespan[data.albany$temperature == 22]),
                 sd(data.albany$lifespan[data.albany$temperature == 25]),
                 sd(data.albany$lifespan[data.albany$temperature == 28]))
data.sderr.alb <- data.sd.alb / sqrt(c(length(data.albany$temperature == 22),
                                       length(data.albany$temperature == 25),
                                       length(data.albany$temperature == 28)))


# Plot data
plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 50), main="Lifespan",
                   ylab=expression(paste("Lifespan ", 
                                         bgroup("[", "days", "]"))))
plot_curves_shaded(temps, lf.suf.curves, lf.alb.curves)


### Fecundity (eggs per raft)

# Load data
data <- read.csv("er_allpop.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

ER.suf.chains <- readRDS("jagsout_epr_suffolk.RDS")
ER.alb.chains <- readRDS("jagsout_epr_albany.RDS")

#View(ER.suf.chains$BUGSoutput$summary[c(2,1,3,4,6), c(1,3,7,8,9)])
#View(ER.alb.chains$BUGSoutput$summary[c(2,1,3,4,6), c(1,3,7,8,9)])

ER.suf.chains <- MCMCchains(ER.suf.chains, params=c("Tmin", "Tmax", "c"))
ER.alb.chains <- MCMCchains(ER.alb.chains, params=c("Tmin", "Tmax", "c"))

ER.suf.curves <- apply(ER.suf.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))
ER.alb.curves <- apply(ER.alb.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))

# Calculate mean and SD.
data.mean.suf <- c(mean(data.suffolk$er[data.suffolk$temperature == 22]),
                   mean(data.suffolk$er[data.suffolk$temperature == 25]),
                   mean(data.suffolk$er[data.suffolk$temperature == 28]))
data.sd.suf <- c(sd(data.suffolk$er[data.suffolk$temperature == 22]),
                 sd(data.suffolk$er[data.suffolk$temperature == 25]),
                 sd(data.suffolk$er[data.suffolk$temperature == 28]))
data.sderr.suf <- data.sd.suf / sqrt(c(length(data.suffolk$temperature == 22),
                                       length(data.suffolk$temperature == 25),
                                       length(data.suffolk$temperature == 28)))

data.mean.alb <- c(mean(data.albany$er[data.albany$temperature == 22]),
                   mean(data.albany$er[data.albany$temperature == 25]),
                   mean(data.albany$er[data.albany$temperature == 28]))
data.sd.alb <- c(sd(data.albany$er[data.albany$temperature == 22]),
                 sd(data.albany$er[data.albany$temperature == 25]),
                 sd(data.albany$er[data.albany$temperature == 28]))
data.sderr.alb <- data.sd.alb / sqrt(c(length(data.albany$temperature == 22),
                                       length(data.albany$temperature == 25),
                                       length(data.albany$temperature == 28)))


# Plot data
plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 200), main="Fecundity",
                   ylab=expression(paste("Fecundity ", 
                                         bgroup("[", frac("eggs", "raft"), "]"))))
plot_curves_shaded(temps, ER.suf.curves, ER.alb.curves)

### Proportion ovipositing

# Load data
data <- read.csv("pO_allpop.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

pO.suf.chains <- readRDS("jagsout_pO_suffolk.RDS")
pO.alb.chains <- readRDS("jagsout_pO_albany.RDS")

#View(pO.suf.chains$BUGSoutput$summary[c(2,1,3,4), c(1,3,7,8,9)])
#View(pO.alb.chains$BUGSoutput$summary[c(2,1,3,4), c(1,3,7,8,9)])

pO.suf.chains <- MCMCchains(pO.suf.chains, params=c("Tmin", "Tmax", "c"))
pO.alb.chains <- MCMCchains(pO.alb.chains, params=c("Tmin", "Tmax", "c"))

pO.suf.curves <- apply(pO.suf.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
pO.alb.curves <- apply(pO.alb.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

# Calculate mean and std err.
data.mean.suf <- data.suffolk$proportion
data.sderr.suf <- sqrt((data.suffolk$proportion * (1 - data.suffolk$proportion))
                       / data.suffolk$N)

data.mean.alb <- data.albany$proportion
data.sderr.alb <- sqrt((data.albany$proportion * (1 - data.albany$proportion))
                       / data.albany$N)


# Plot data
plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 1), main="Proportion ovipositing",
                   ylab="pO")
plot_curves_shaded(temps, pO.suf.curves, pO.alb.curves)


### Larval survival

# Load data
data <- read.csv("pLA_allpop.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

pLA.suf.chains <- readRDS("jagsout_pLA_suffolk.RDS")
pLA.alb.chains <- readRDS("jagsout_pLA_albany.RDS")

#View(pLA.suf.chains$BUGSoutput$summary[c(2,1,3,4), c(1,3,7,8,9)])
#View(pLA.alb.chains$BUGSoutput$summary[c(2,1,3,4), c(1,3,7,8,9)])

pLA.suf.chains <- MCMCchains(pLA.suf.chains, params=c("Tmin", "Tmax", "c"))
pLA.alb.chains <- MCMCchains(pLA.alb.chains, params=c("Tmin", "Tmax", "c"))

pLA.suf.curves <- apply(pLA.suf.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
pLA.alb.curves <- apply(pLA.alb.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

# Calculate mean and std err.
data.mean.suf <- data.suffolk$proportion
data.sderr.suf <- sqrt((data.suffolk$proportion * (1 - data.suffolk$proportion))
                       / data.suffolk$N)

data.mean.alb <- data.albany$proportion
data.sderr.alb <- sqrt((data.albany$proportion * (1 - data.albany$proportion))
                       / data.albany$N)


# Plot data
plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 1), main="Larval survival to adulthood",
                   ylab="pLA")
plot_curves_shaded(temps, pLA.suf.curves, pLA.alb.curves)


### Mosquito development rate

data <- read.csv("mdr_allpop.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

MDR.suf.chains <- readRDS("jagsout_MDR_suffolk.RDS")
MDR.alb.chains <- readRDS("jagsout_MDR_albany.RDS")

MDR.suf.chains
MDR.alb.chains
#View(MDR.suf.chains$BUGSoutput$summary[c(2,1,3,5,7,4), c(1,3,7,8,9)])
#View(MDR.alb.chains$BUGSoutput$summary[c(2,1,3,5,7,4), c(1,3,7,8,9)])

MDR.suf.chains <- MCMCchains(MDR.suf.chains, params=c("Tmin", "Tmax", "c"))
MDR.alb.chains <- MCMCchains(MDR.alb.chains, params=c("Tmin", "Tmax", "c"))

MDR.suf.curves <- apply(MDR.suf.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))
MDR.alb.curves <- apply(MDR.alb.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))

data.mean.suf <- c(mean(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 22]),
                   mean(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 25]),
                   mean(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 28]))
data.sd.suf <- c(sd(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 22]),
                 sd(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 25]),
                 sd(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 28]))
data.sderr.suf <- data.sd.suf / sqrt(c(length(data.suffolk$temperature == 22),
                                       length(data.suffolk$temperature == 25),
                                       length(data.suffolk$temperature == 28)))

data.mean.alb <- c(mean(1 / data.albany$days.to.pupation[data.albany$temperature == 22]),
                   mean(1 / data.albany$days.to.pupation[data.albany$temperature == 25]),
                   mean(1 / data.albany$days.to.pupation[data.albany$temperature == 28]))
data.sd.alb <- c(sd(1 / data.albany$days.to.pupation[data.albany$temperature == 22]),
                 sd(1 / data.albany$days.to.pupation[data.albany$temperature == 25]),
                 sd(1 / data.albany$days.to.pupation[data.albany$temperature == 28]))
data.sderr.alb <- data.sd.alb / sqrt(c(length(data.albany$temperature == 22),
                                       length(data.albany$temperature == 25),
                                       length(data.albany$temperature == 28)))

# Plot data
plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 0.2), main="Mosquito development rate",
                   ylab=expression(paste("MDR ", bgroup("[", frac(1, "day"), "]"))))
plot_curves_shaded(temps, MDR.suf.curves, MDR.alb.curves)


## R0 and mosquito surveillance data
temps <- seq(10, 40, 0.025)
suffolk.col <- "#CC79A7"
erie.col <- "#56B4E9" 

erie.weather <- read.csv("erie_weather.csv")
suffolk.weather <- read.csv("suffolk_weather.csv")

erie.weather.22 <- subset(erie.weather, erie.weather$year == 2022)
suffolk.weather.22 <- subset(suffolk.weather, suffolk.weather$year == 2022)

mosq.mle <- read.csv("mosquito_incidence.csv")
suffolk.mosq.22 <- subset(mosq.mle, (mosq.mle$County == "Suffolk") & (mosq.mle$year == 2022))
erie.mosq.22 <- subset(mosq.mle, (mosq.mle$County == "Erie") & (mosq.mle$year == 2022))

epi.week <- function(day, start.day=3) {
  return(1 + (day - start.day)  %/% 7)
}

# Start day of epi week 1.
start.day.2022 <- 2

erie.weather.22$week <- apply(erie.weather.22, 1, function(x) epi.week(as.integer(x[4]), start.day.2022))
suffolk.weather.22$week <- apply(suffolk.weather.22, 1, function(x) epi.week(as.integer(x[4]), start.day.2022))

suffolk.weather.22.weekly <- aggregate(suffolk.weather.22[,6:11], list(suffolk.weather.22$week), mean)
colnames(suffolk.weather.22.weekly)[1] <- "week"
erie.weather.22.weekly <- aggregate(erie.weather.22[,6:11], list(erie.weather.22$week), mean)
colnames(erie.weather.22.weekly)[1] <- "week"

suffolk.22 <- merge(suffolk.weather.22.weekly, suffolk.mosq.22, by=c("week"))
erie.22 <- merge(erie.weather.22.weekly, erie.mosq.22, by=c("week"))

suffolk.22$lagged.tmeanc <- Lag(suffolk.22$tmeanc, shift=2)
erie.22$lagged.tmeanc <- Lag(erie.22$tmeanc, shift=2)

plot(suffolk.22$lagged.tmeanc, suffolk.22$Infection.Rate, pch=15, col="darkred",
     xlab="Temperature [ºC]", ylab="WNV prevalence",  xlim = c(10, 35), ylim=c(0, 20),
     main="")


R0.suf.curves <- readRDS("R0_chains_suffolk.RDS")
R0.alb.curves <- readRDS("R0_chains_albany.RDS")
points(erie.22$lagged.tmeanc, erie.22$Infection.Rate, pch=15, col="steelblue")
plot_curves_shaded(temps, R0.suf.curves, R0.alb.curves, scale.factor=20, shade.alpha=0.2)


## Save mean and CIs for R0 for map plotting.
R0.suf.mean <- apply(R0.suf.curves, 1, mean)
R0.suf.CI <- apply(R0.suf.curves, 1, quantile, c(0.025, 0.975))

saveRDS(R0.suf.mean, "R0_mean_suffolk.RDS")
saveRDS(R0.suf.CI, "R0_CI_suffolk.RDS")

R0.alb.mean <- apply(R0.alb.curves, 1, mean)
R0.alb.CI <- apply(R0.alb.curves, 1, quantile, c(0.025, 0.975))

saveRDS(R0.alb.mean, "R0_mean_albany.RDS")
saveRDS(R0.alb.CI, "R0_CI_albany.RDS")
