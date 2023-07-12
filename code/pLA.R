set.seed(42) # Set seed for reproducibility.

# We will need these packages. Please make sure you install them beforehand!
library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')


quad_model <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

temp <- seq(0, 50, 0.1)
plot(temp, quad_model(temp, 8.2, 33.2, 0.007), type="l")

data <- read.csv("pla_allpop.csv")
data

data.suffolk <- subset(data, (data$population == "suffolk") )
data.albany <- subset(data, (data$population == "albany") )

data.suffolk
data.albany

suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9" 

plot(data.suffolk$temperature, data.suffolk$proportion, xlim=c(10, 35), ylim=c(0, 1),
     pch=20, main="Suffolk")
plot(data.albany$temperature, data.albany$proportion, xlim=c(10, 35), ylim=c(0, 1),
     pch=20, main="Albany")


sink("quad_pLA_binom.txt")
cat("
    model{
    
    ## Priors
    c_mu <- 0.0036
    c_std <- 0.001
    c ~ dgamma(c_mu^2 / c_std^2, c_mu / c_std^2)
    
    Tmin ~ dnorm(7.8, 1/2^2) 
    Tmax ~ dnorm(38.4, 1/2^2 )
    
    ## Likelihood
    for(i in 1:N.obs){
    p[i] <- min(c * (temp[i] - Tmin) * (Tmax - temp[i]) * (temp[i] < Tmax) * (Tmin < temp[i]),
    1)
    n[i] ~ dbin(p[i], N[i])
    }
    
    ## Derived Quantities and Predictions
    Topt <- (Tmin + Tmax) / 2
    
    } # close model
    ",fill=TRUE)
sink()

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=4, max=10),
  Tmax = runif(1, min=30, max=40),
  c = runif(1, min=0.001, max=0.005))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "c", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 1000000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
N.suf <- data.suffolk$N
n.suf <- data.suffolk$alive
temp.suf <- data.suffolk$temperature
N.obs.suf <- length(N.suf)

N.alb <- data.albany$N
n.alb <- data.albany$alive
temp.alb <- data.albany$temperature
N.obs.alb <- length(N.alb)


##### Bundle Data
jag.data<-list(N=N.suf, n=n.suf, temp = temp.suf, N.obs=N.obs.suf)
jag.data
##### Run JAGS
suffolk.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_pLA_binom.txt", n.thin=nt, n.chains=nc, 
                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
suffolk.out
mcmcplot(suffolk.out)

chains.suf <- MCMCchains(suffolk.out, params=c("Tmin", "Tmax", "c"))
chains.suf
temps <- seq(0, 45, 0.1) 
curves <- apply(chains.suf, 1, 
                function(x) quad_model(temps, x[1], x[2], x[3]))

# Find mean curve and credible intervals.
meancurve.suf <- apply(curves, 1, mean)
CI.suf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.suffolk$temperature, data.suffolk$proportion, xlab="Temperature [°C]",
     ylab="pLA", xlim=c(0, 50), ylim=c(0, 1), pch=20, col=suffolk.col, 
     main="Suffolk (informative)")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.suf, col="black")
lines(temps, CI.suf[1,], col="gray")
lines(temps, CI.suf[2,], col="gray")

#lines(temps, quad_model(temps, 8.2, 33.2, 0.00445), type="l", col="green")

##### Bundle Data
jag.data<-list(N = N.alb, n = n.alb, temp = temp.alb, N.obs=N.obs.alb)
##### Run JAGS
albany.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_pLA_binom.txt", n.thin=nt, n.chains=nc, 
                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
albany.out
mcmcplot(albany.out)



chains.alb <- MCMCchains(albany.out, params=c("Tmin", "Tmax", "c"))
temps <- seq(0, 45, 0.1) 
curves <- apply(chains.alb, 1, 
                function(x) quad_model(temps, x[1], x[2], x[3]))
head(chains.alb)
# Find mean curve and credible intervals.
meancurve.alb <- apply(curves, 1, mean)
CI.alb <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.albany$temperature, data.albany$proportion, xlab="Temperature [°C]",
     ylab="pLA", xlim=c(0, 50), ylim=c(0, 1), pch=20, col=albany.col, 
     main="Albany (informative)")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.alb, col="black")
lines(temps, CI.alb[1,], col="gray")
lines(temps, CI.alb[2,], col="gray")


# Plot fit from Shocket et al.
#lines(temps, quad_model(temps, 7.8, 38.4, 0.0036), col="firebrick")

# Plot both curves simultaneously.
plot(data.suffolk$temperature, data.suffolk$proportion, xlab="Temperature [°C]",
     ylab="pLA", xlim=c(0, 45), ylim=c(0, 1), pch=20, col=suffolk.col,
     main="Larval survival to adulthood")
lines(temps, meancurve.suf, col=suffolk.col)
lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5))
lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5))

points(data.albany$temperature, data.albany$proportion, pch=20, col=albany.col)
lines(temps, meancurve.alb, col=albany.col)
lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5))
lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5))

saveRDS(albany.out, file="jagsout_pLA_albany.RDS")
saveRDS(suffolk.out, file="jagsout_pLA_suffolk.RDS")
