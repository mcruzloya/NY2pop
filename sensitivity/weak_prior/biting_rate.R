set.seed(42) # Set seed for reproducibility.
setwd("/Users/cruzloya/git/NY2pop/sensitivity/weak_prior/")


library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax)  )
  return(result)
}

quad_model <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(result)
}

data <- read.csv("biting_rate_lf.csv")
data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany"))

data.suffolk

suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9" 

plot(data.suffolk$temperature, data.suffolk$biting.rate, xlim=c(10, 35), ylim=c(0, 0.1),
     pch=20, main="Suffolk")
plot(data.albany$temperature, data.albany$biting.rate, xlim=c(10, 35), ylim=c(0, 0.1),
     pch=20, main="Albany")

sink("quad_br_binom_wp.txt") 
cat("
    model{
    
    ## Priors
    c ~ dunif(0, 0.01) # Uniform prior since scale is different from prior biting rates.
    Tmin ~ dnorm(9.4, 1/(2*4)^2) # Based on fitted values from Shocket et. al
    Tmax ~ dnorm(39.6, 1/(2*4)^2) # Based on fitted values from Shocket et. al
    
    ## Likelihood
    for(i in 1:N.obs){
    
      p[i] <- min(c * (temp[i] - Tmin) * (Tmax - temp[i]) * (temp[i] < Tmax) * (Tmin < temp[i]),
              1)
      n[i] ~ dbin(p[i], N[i])
    }
    
    ## Derived Quantities
    Topt = (Tmin + Tmax) / 2
    
    } # close model
    ")
sink()


##### Calculate initial values for MCMC.
inits<-function(){list(
  c = runif(1, min=0.00001, 0.0001),
  Tmax = runif(1, min=30, max=40),
  Tmin = runif(1, min=5, max=10))}

##### Parameters to Estimate
parameters <- c("c", "Tmin", "Tmax", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
n.suf <- data.suffolk$feeding.events
N.suf <- data.suffolk$mosquito.days
N.obs.suf <- length(n.suf)
temp.suf <- data.suffolk$temperature

temp.suf

## Suffolk analysis

##### Bundle Data
jag.data<-list(N = N.suf, n = n.suf, N.obs = N.obs.suf, temp = temp.suf)
jag.data
##### Run JAGS
suffolk.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_br_binom_wp.txt",
                 n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
suffolk.out
View(suffolk.out$BUGSoutput$summary)
mcmcplot(suffolk.out)


chains.suffolk <- MCMCchains(suffolk.out, params=c("Tmin", "Tmax", "c"))


temps <- seq(0, 50, 0.1) 
curves <- apply(chains.suffolk, 1, function(x) quad_model(temps, x[1], x[2], x[3]))




# Find mean curve and credible intervals.
meancurve.suf <- apply(curves, 1, mean)
CI.suf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.suffolk$temperature, data.suffolk$biting.rate, xlab="Temperature [°C]",
     ylab="a (biting rate)", xlim=c(0, 50), ylim=c(0, 0.1), pch=20, col=suffolk.col,
     main="Suffolk")
lines(temps, meancurve.suf, col="black")
lines(temps, CI.suf[1,], col="gray")
lines(temps, CI.suf[2,], col="gray")

## Albany analysis
n.alb <- data.albany$feeding.events
N.alb <- data.albany$mosquito.days
N.obs.alb <- length(n.alb)
temp.alb <- data.albany$temperature


##### Bundle Data
jag.data<-list(N = N.alb, n = n.alb, N.obs = N.obs.alb, temp = temp.alb)
##### Run JAGS
albany.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, model.file="quad_br_binom_wp.txt",
                    n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
albany.out 
View(albany.out$BUGSoutput$summary)


mcmcplot(albany.out)

chains.albany <- MCMCchains(albany.out, params=c("Tmin", "Tmax", "c"))


temps <- seq(0, 50, 0.1) 
curves <- apply(chains.albany, 1, function(x) quad_model(temps, x[1], x[2], x[3]))


# Find mean curve and credible intervals.
meancurve.alb <- apply(curves, 1, mean)
CI.alb <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.albany$temperature, data.albany$biting.rate, xlab="Temperature [°C]",
     ylab="a (biting rate)", xlim=c(0, 50), ylim=c(0, 0.1), pch=20, col=albany.col,
     main="Albany")
lines(temps, meancurve.alb, col="black")
lines(temps, CI.alb[1,], col="gray")
lines(temps, CI.alb[2,], col="gray")

# Plot both curves simultaneously.
plot(data.suffolk$temperature, data.suffolk$biting.rate, xlab="Temperature [°C]",
     ylab="Biting rate [bites / mosquito * day]", xlim=c(0, 50), ylim=c(0, 0.1), pch=20, col=suffolk.col,
     main="Biting rate")
lines(temps, meancurve.suf, col=suffolk.col)
lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5))
lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5))

points(data.albany$temperature, data.albany$biting.rate, pch=20, col=albany.col)
lines(temps, meancurve.alb, col=albany.col)
lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5))
lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5))

saveRDS(albany.out, file="jagsout_a_albany_wp.RDS")
saveRDS(suffolk.out, file="jagsout_a_suffolk_wp.RDS")
