set.seed(42) # Set seed for reproducibility.

# We will need these packages. Please make sure you install them beforehand!
library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')
library('ggplot2')

briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax)  )
  return(result)
}

data <- read.csv("mdr_allpop.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

data.suffolk
plot(data.suffolk$temperature, 1 / data.suffolk$days.to.pupation, 
     xlim=c(10, 35), ylim=c(0, 0.2),
     pch=20, main="Suffolk")
plot(data.albany$temperature, 1 / data.albany$days.to.pupation, xlim=c(10, 35), ylim=c(0, 0.2),
     pch=20, main="Albany")

suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9" 


# Visualize data.
g <- (ggplot(data, aes(temperature, 1 / days.to.pupation)) +  
        geom_violin(aes(col=population, group=cut_width(temperature, 1))) + 
        geom_point(aes(col=population, alpha=I(0.3))) +
        facet_wrap(~population, nrow=1) +
        theme_classic() + scale_color_manual(values=c(albany.col, suffolk.col)))
g

# Visualize data.
g <- (ggplot(data, aes(temperature, days.to.pupation)) +  
        geom_violin(aes(col=population, group=cut_width(temperature, 1))) + 
        geom_point(aes(col=population, alpha=I(0.3))) +
        facet_wrap(~population, nrow=1) +
        theme_classic() + scale_color_manual(values=c(albany.col, suffolk.col)))
g

# MDR model with gamma likelihood and variable SD.
sink("briere_mdr.txt") 
cat("
    model{
    
    ## Priors
    mu_c_ <- 3.76 
    std_c_ <- 2 
    c_ ~ dgamma(mu_c_^2 /std_c_^2, mu_c_ / std_c_^2)
    c <- c_ * 10^-5
    
    Tmin ~ dnorm(5, 1/2.5^2)
    Tmax ~ dnorm(38.5, 1/1.5^2) 
    
    sigma10 ~ dunif(0, 20)
    beta ~ dnorm(0, 10)
    
    ## Likelihood
    for(i in 1:N.obs){
      trait.mu[i] <- min(1 / (c * temp[i] * (temp[i] - Tmin) * 
                     sqrt((Tmax - temp[i]) * (Tmax > temp[i])) * (Tmin < temp[i])),
                     1000)
      sigma[i] <- sigma10 * exp(beta*(trait.mu[i] - 10))
      trait[i] ~ dgamma(trait.mu[i]^2 / sigma[i]^2, trait.mu[i] / sigma[i]^2)
    }
    
    ## Derived Quantities
    Topt = (4*Tmax + 3*Tmin)/10 + sqrt(4 * Tmax^2 + (9/4) * Tmin^2 - 4*Tmax*Tmin) / 5
    
    } # close model
    ")
sink()

inits<-function(){list(
  c_ = runif(1, min=1, max=5),
  Tmax = runif(1, min=35, max=40),
  Tmin = runif(1, min=0, max=5),
  sigma10 = runif(1, min=1, max=10),
  beta = rnorm(1))}

##### Parameters to Estimate
parameters <- c("c", "Tmin", "Tmax","sigma10", "beta", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 1000000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
trait <- data.suffolk$days.to.pupation
N.obs <- length(trait)
temp <- data.suffolk$temperature

##### Bundle Data
jag.data<-list(trait = trait, N.obs = N.obs, temp = temp)
##### Run JAGS
suffolk.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="briere_mdr.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE,
                    working.directory=getwd())
suffolk.out
mcmcplot(suffolk.out)


chains.suffolk <- MCMCchains(suffolk.out, params=c("Tmin", "Tmax", "c"))

temps <- seq(0, 50, 0.1) 
curves <- apply(chains.suffolk, 1, function(x) briere(temps, x[1], x[2], x[3]))

# Find mean curve and credible intervals.
meancurve.suf <- apply(curves, 1, mean)
CI.suf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.suffolk$temperature, 1 / data.suffolk$days.to.pupation, 
     xlab="Temperature [°C]",
     ylab="MDR", xlim=c(0, 45), ylim=c(0, 0.2), pch=20, 
     col=alpha(suffolk.col, 0.3),
     main="Suffolk")
lines(temps, meancurve.suf, col="black")
lines(temps, CI.suf[1,], col="gray")
lines(temps, CI.suf[2,], col="gray")

## Albany analysis
trait <- data.albany$days.to.pupation
N.obs <- length(trait)
temp <- data.albany$temperature

##### Bundle Data
jag.data <- list(trait = trait, N.obs = N.obs, temp = temp)
##### Run JAGS
albany.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="briere_mdr.txt",
                   n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
albany.out

mcmcplot(albany.out)

chains.albany <- MCMCchains(albany.out, params=c("Tmin", "Tmax", "c"))


temps <- seq(0, 50, 0.1) 
curves <- apply(chains.albany, 1, function(x) briere(temps, x[1], x[2], x[3]))


# Find mean curve and credible intervals.
meancurve.alb <- apply(curves, 1, mean)
CI.alb <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.albany$temperature, 1 / data.albany$days.to.pupation, xlab="Temperature [°C]",
     ylab="MDR", xlim=c(0, 45), ylim=c(0, 0.2), pch=20, 
     col=alpha(albany.col, 0.3),
     main="Albany")
lines(temps, meancurve.alb, col="black")
lines(temps, CI.alb[1,], col="gray")
lines(temps, CI.alb[2,], col="gray")

# Plot both curves simultaneously.
plot(data.suffolk$temperature, 1 / data.suffolk$days.to.pupation, 
     xlab="Temperature [°C]",
     ylab="MDR [1 / day]", xlim=c(0, 45), ylim=c(0, 0.2), pch=20, 
     col=alpha(suffolk.col, 0.3),
     main="Mosquito Development Rate")
lines(temps, meancurve.suf, col=suffolk.col)
lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5))
lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5))

points(data.albany$temperature, 1 / data.albany$days.to.pupation, pch=20,
       col=alpha(albany.col, 0.3) )
lines(temps, meancurve.alb, col=albany.col)
lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5))
lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5))

data.mean.suf <- c(mean(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 22]),
                   mean(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 25]),
                   mean(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 28]))
data.sd.suf <- c(sd(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 22]),
                 sd(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 25]),
                 sd(1 / data.suffolk$days.to.pupation[data.suffolk$temperature == 28]))
data.sderr.suf <- data.sd.suf / sqrt(c(length(data.suffolk$temperature == 22),
                                       length(data.suffolk$temperature == 25),
                                       length(data.suffolk$temperature == 28)))

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


plot(c(22,25,28), data.mean.suf, xlab="Temperature [°C]",
     ylab="MDR [1 / day]", xlim=c(0, 45), ylim=c(0, 0.2), pch=20, col=suffolk.col,
     main="Mosquito Development Rate")
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf+data.sderr.suf, col=suffolk.col,
       angle=90, length=0.03)
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf-data.sderr.suf, col=suffolk.col,
       angle=90, length=0.03)

lines(temps, meancurve.suf, type="l", xlim=c(0, 45), ylim=c(0, 60),
      pch=20, col=suffolk.col)
lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5))
lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5))

points(c(22,25,28), data.mean.alb, pch=20, col=albany.col)
arrows(c(22,25,28), data.mean.alb, y1=data.mean.alb+data.sderr.alb, col=albany.col,
       angle=90, length=0.03)
arrows(c(22,25,28), data.mean.alb, y1=data.mean.alb-data.sderr.alb, col=albany.col,
       angle=90, length=0.03)
lines(temps, meancurve.alb, col=albany.col)
lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5))
lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5))

saveRDS(albany.out, file="jagsout_MDR_albany.RDS")
saveRDS(suffolk.out, file="jagsout_MDR_suffolk.RDS")

