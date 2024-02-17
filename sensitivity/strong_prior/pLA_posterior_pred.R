library('R2jags')
library('MCMCvis')
library('scales')


# Quadratic TPC model with upper limit at 1 for proportions.
quad_lim <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

# Shocket et al data
data.pLA <- read.csv("TraitData_pLA.csv") # Data from database for most traits (except below)
data.pLA.Cpip <- subset(data.pLA, host.code == "Cpip")
data.pLA.Cqui <- subset(data.pLA, host.code == "Cqui")

# Shocket et al fit.
load("./jagsout_pLA_Cpip_inf.Rdata")
pLA.Cpip.out.inf

pLA.chains.shocket
pLA.chains.shocket <- MCMCchains(pLA.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
pLA.curves.shocket <- apply(pLA.chains.shocket, 1, 
                            function(x) quad_lim(temps, x[1], x[2], x[3]))

# Read pLA data
data <- read.csv("pla_allpop.csv")
data

data.suffolk

data.suffolk <- subset(data, data$population == "suffolk")
data.albany <- subset(data, data$population == "albany")

pLA.alb.chains <- readRDS("jagsout_pLA_albany.RDS")
pLA.suf.chains <- readRDS("jagsout_pLA_suffolk.RDS")

data.suffolk
# Generate posterior samples

pLA.suf.chains <- MCMCchains(pLA.suf.chains, params=c("Tmin", "Tmax", "c"))
pLA.alb.chains <- MCMCchains(pLA.alb.chains, params=c("Tmin", "Tmax", "c"))

exp.temps <- c(22, 25, 28)

pred.survival.22 <- apply(pLA.suf.chains, 1, function(x) rbinom(1, 100, quad_lim(22, x[1], x[2], x[3]))) 
pred.survival.25 <- apply(pLA.suf.chains, 1, function(x) rbinom(1, 100, quad_lim(25, x[1], x[2], x[3]))) 
pred.survival.28 <- apply(pLA.suf.chains, 1, function(x) rbinom(1, 100, quad_lim(28, x[1], x[2], x[3]))) 



pLA.suf.curves <- apply(pLA.suf.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
dim(pLA.suf.curves)

pLA.suf.pred <- apply(pLA.suf.curves, 2, function(x) rbinom(length(temps), 100, x))
pLA.suf.pred <- pLA.suf.pred / 100
pLA.PI.suf <- apply(pLA.suf.pred, 1, quantile, c(0.025, 0.975))
pLA.PI.suf

dim(pLA.suf.curves)
dim(pLA.suf.pred)

length(temps)

pLA.alb.curves <- apply(pLA.alb.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

# Calculate mean and std err.
data.mean.suf <- data.suffolk$proportion
data.sderr.suf <- sqrt((data.suffolk$proportion * (1 - data.suffolk$proportion))
                       / data.suffolk$N)

data.mean.alb <- data.albany$proportion
data.sderr.alb <- sqrt((data.albany$proportion * (1 - data.albany$proportion))
                       / data.albany$N)

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

plot_point_errbar <- function(data.mean, data.stderr,
                               xlim=c(0, 45), ylim=c(0, 1), main="Trait",
                               xlab="Temperature [°C]", ylab="trait",
                               errbar.length=0.06, errbar.width=2,
                               cex.main=1.5, cex.lab=1.3, cex.axis=1.2, cex=1.5,
                               col=col) {
  plot(c(22,25,28), data.mean, xlab=xlab, ylab=ylab,
       xlim=xlim, ylim=ylim, pch=20, col=col, main=main, 
       cex.main=cex.main, cex.lab=cex.lab, cex.axis=cex.axis, cex=cex)
  arrows(c(22,25,28), data.mean, y1=data.mean+data.stderr, col=col,
         angle=90, length=errbar.length, lwd=errbar.width)
  arrows(c(22,25,28), data.mean, y1=data.mean-data.stderr, col=col,
         angle=90, length=errbar.length, lwd=errbar.width)
}

plot_curve_shaded <- function(temps, curves, 
                               col="#CC79A7",
                               mean.lwd=2.5, shade.alpha=0.2, scale.factor=1.0) {
  
  meancurve <- scale.factor * apply(curves, 1, mean)
  CI <- scale.factor * apply(curves, 1, quantile, c(0.025, 0.975))
  
  lines(temps, meancurve, col=col, lwd=mean.lwd)
  polygon(c(temps, rev(temps)), c(CI[1,], rev(CI[2,])), 
          col=alpha(col, shade.alpha), lty=0)
}


# Plot data for both populations
plot_points_errbar(data.mean.suf, data.stderr.suf, data.mean.alb, data.stderr.alb, 
                   xlim=c(0, 45), ylim=c(0, 1), main="Larval survival to adulthood",
                   ylab="pLA")
plot_curves_shaded(temps, pLA.suf.curves, pLA.alb.curves)


# Plot only Suffolk data.
plot_point_errbar(data.mean.suf, data.sderr.suf,
                   xlim=c(0, 45), ylim=c(0, 1), main="Larval survival to adulthood",
                   ylab="pLA", col=suffolk.col)
plot_curve_shaded(temps, pLA.suf.curves)

# Plot prediction interval
plot(c(22, 25, 28), data.mean.suf, col=suffolk.col, pch=20, ylim=c(0, 1), xlim=c(0, 45),
     xlab="Temperature [°C]", ylab="pLA", cex.main=1.5, cex.lab=1.3, cex.axis=1.2, cex=1.5,
     main="Larval survival prediction interval")
lines(temps, pLA.PI.suf[1,], col="gray", lwd=1.5)
lines(temps, pLA.PI.suf[2,], col="gray", lwd=1.5)

pLA.PI.suf[2,]

### Prior and posterior plot

par(mfrow=c(1,3))
plot(data.pLA.Cpip$T, data.pLA.Cpip$trait, pch=20, xlim=c(0, 45), ylim=c(0, 1),
     xlab="Temperature [°C]", ylab="pLA", main="Literature data and fit",
     cex.main=1.5, cex.lab=1.3, cex.axis=1.2, cex=1.5)
pLA.shocket.mean <- apply(pLA.curves.shocket, 1, mean)
pLA.shocket.CI <- apply(pLA.curves.shocket, 1, quantile, c(0.025, 0.975))

lines(temps, pLA.shocket.mean, lwd=2, col="darkgreen")
polygon(c(temps, rev(temps)), c(pLA.shocket.CI[1,], rev(pLA.shocket.CI[2,])), 
        col=alpha("darkgreen", 0.5), lty=0)

# Prior curves
mu_c <- 0.0036 
std_c <- 0.001 

prior.Tmin <- rnorm(10000, 7.8, 2) 
prior.Tmax <- rnorm(10000, 38.4, 2)
prior.c <- rgamma(10000, mu_c^2 /std_c^2, mu_c / std_c^2)

prior.chains <- cbind(prior.Tmin, prior.Tmax, prior.c)


pLA.curves.prior <- apply(prior.chains, 1, 
                          function(x) quad_lim(temps, x[1], x[2], x[3]))
pLA.prior.mean <- apply(pLA.curves.prior, 1, mean)
pLA.prior.CI <- apply(pLA.curves.prior, 1, quantile, c(0.025, 0.975))

plot(data.pLA.Cpip$T, data.pLA.Cpip$trait, pch=20, xlim=c(0, 45), ylim=c(0, 1),
     xlab="Temperature [°C]", ylab="pLA", main="Prior distribution",
     cex.main=1.5, cex.lab=1.3, cex.axis=1.2, cex=1.5)

lines(temps, pLA.prior.mean, lwd=2, col="orange")

# Sample ten prior curves randomly.
sample.idx <- sample(1:10000, 10)
for(i in 1:10) {
  lines(temps, pLA.curves.prior[, sample.idx[i]], lwd=1,
        col=alpha("grey", 0.8))
}




polygon(c(temps, rev(temps)), c(pLA.prior.CI[1,], rev(pLA.prior.CI[2,])), 
        col=alpha("orange", 0.5), lty=0)


pLA.suf.mean <- apply(pLA.suf.curves, 1, mean)
pLA.suf.CI <- apply(pLA.suf.curves, 1, quantile, c(0.025, 0.975))


plot(c(22,25,28), data.mean.suf, xlab="Temperature [°C]", 
     ylab="pLA",
     xlim=c(0, 45), ylim=c(0,1), pch=20, col=suffolk.col, 
     main="Posterior distribution", 
     cex.main=1.5, cex.lab=1.3, cex.axis=1.2, cex=1.5)

errbar.length=0.06
errbar.width=2
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf+data.sderr.suf, col=suffolk.col,
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf-data.sderr.suf, col=suffolk.col,
       angle=90, length=errbar.length, lwd=errbar.width)


lines(temps, pLA.suf.mean, lwd=2, col=suffolk.col)
#polygon(c(temps, rev(temps)), c(pLA.PI.suf[1,], rev(pLA.PI.suf[2,])), 
#        col=alpha("grey", 0.5), lty=0)
polygon(c(temps, rev(temps)), c(pLA.suf.CI[1,], rev(pLA.suf.CI[2,])), 
        col=alpha(suffolk.col, 0.5), lty=0)
lines(temps, pLA.PI.suf[1,], col="gray", lwd=1.5)
lines(temps, pLA.PI.suf[2,], col="gray", lwd=1.5)

#lines(temps, MDR.alb.mean, lwd=2, col=albany.col)
#polygon(c(temps, rev(temps)), c(MDR.alb.CI[1,], rev(MDR.alb.CI[2,])), 
#        col=alpha(albany.col, 0.5), lty=0)


#lines(temps, pLA.prior.mean, lwd=2, col="orange")
#polygon(c(temps, rev(temps)), c(pLA.prior.CI[1,], rev(pLA.prior.CI[2,])), 
#        col=alpha("orange", 0.5), lty=0)

