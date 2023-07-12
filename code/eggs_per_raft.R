set.seed(42) # Set seed for reproducibility.

library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

quad_model <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(result)
}

data <- read.csv("er_allpop.csv")
data

data.suffolk <- subset(data, (data$population == "suffolk"))
data.albany <- subset(data, (data$population == "albany") )

data.suffolk
data.albany

suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9" 
    

plot(data.suffolk$temperature, data.suffolk$er, xlim=c(10, 35), ylim=c(0, 200),
     pch=20, main="Suffolk")
plot(data.albany$temperature, data.albany$er, xlim=c(10, 35), ylim=c(0, 200),
     pch=20, main="Albany")

# Visualize data.
g <- (ggplot(data, aes(temperature, er)) +  
        geom_violin(aes(col=population, group=cut_width(temperature, 1))) + 
        geom_point(aes(col=population, alpha=I(0.3))) +
        facet_wrap(~population, nrow=1) +
        theme_classic() + scale_color_manual(values=c(albany.col, suffolk.col)))
g

sink("quad_inf_epr.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(5.3, 1 / 2.5^2) #  Shocket et al. estimate with inflated std
    Tmax ~ dnorm(38.9, 1 / 2.5^2) # Shocket. et al. estimate with inflated std
    
    c_mu <- 0.6
    c_std <- 0.2
    c ~ dgamma(c_mu^2 / c_std^2, c_mu / c_std^2)
    
    sigma_mu <- 25
    sigma_std <- 50
    sigma ~ dgamma(sigma_mu^2 / sigma_std^2, sigma_mu / sigma_std^2)
    tau <- 1 / (sigma*sigma)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- c * (temp[i] - Tmin) * (Tmax - temp[i]) * (temp[i] > Tmin) * (temp[i] < Tmax)
    trait[i] ~ dnorm(trait.mu[i], tau)
    }
    
    ## Derived quantities
    Topt = (Tmin + Tmax) / 2
    
    } # close model
    ",fill=TRUE)
sink()

##### Calculate initial values for MCMC.
inits<-function(){list(
  Tmin = runif(1, min=2, max=10),
  Tmax = runif(1, min=30, max=40),
  c = runif(1, min=0, max=0.7),
  sigma = runif(1, min=10, max=50))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "c", "sigma", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
trait.suf <- data.suffolk$er
temp.suf <- data.suffolk$temperature
N.obs.suf <- length(trait.suf)

trait.alb <- data.albany$er
temp.alb <- data.albany$temperature
N.obs.alb <- length(trait.alb)


##### Bundle Data
jag.data<-list(trait = trait.suf, N.obs = N.obs.suf, temp = temp.suf)
##### Run JAGS
suffolk.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_inf_epr.txt", n.thin=nt, n.chains=nc, 
                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
suffolk.out
mcmcplot(suffolk.out)

chains.suf <- MCMCchains(suffolk.out, params=c("Tmin", "Tmax", "c"))
temps <- seq(0, 45, 0.1) 
curves <- apply(chains.suf, 1, 
                function(x) quad_model(temps, x[1], x[2], x[3]))

# Find mean curve and credible intervals.
meancurve.suf <- apply(curves, 1, mean)
CI.suf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.suffolk$temperature, data.suffolk$er, xlab="Temperature [°C]",
     ylab="Fecundity (eggs/raft)", xlim=c(0, 50), ylim=c(0, 250), pch=20, 
     col=alpha(suffolk.col, 0.3),  main="Suffolk")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.suf, col="black")
lines(temps, CI.suf[1,], col="gray")
lines(temps, CI.suf[2,], col="gray")

##### Bundle Data
jag.data<-list(trait = trait.alb, N.obs = N.obs.alb, temp = temp.alb)
##### Run JAGS
albany.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="quad_inf_epr.txt", n.thin=nt, n.chains=nc, 
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

plot(data.albany$temperature, data.albany$er, xlab="Temperature [°C]",
     ylab="Fecundity (eggs/raft)", xlim=c(0, 50), ylim=c(0, 250), pch=20, 
     col=alpha(albany.col, 0.5), 
     main="Albany")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.alb, col="black")
lines(temps, CI.alb[1,], col="gray")
lines(temps, CI.alb[2,], col="gray")

# Plot both curves simultaneously.
plot(data.suffolk$temperature, data.suffolk$er, xlab="Temperature [°C]",
     ylab="fecundity [eggs/raft]", xlim=c(0, 50), ylim=c(0, 250), pch=20,
     col=alpha(suffolk.col, 0.3),
     main="Fecundity")
lines(temps, meancurve.suf, col=suffolk.col)
lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5))
lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5))

points(data.albany$temperature, data.albany$er, pch=20,
       col=alpha(albany.col, 0.3) )
lines(temps, meancurve.alb, col=albany.col)
lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5))
lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5))

## Plot means and stderr
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

data.mean.suf
plot(c(22,25,28), data.mean.suf, xlab="Temperature [°C]",
     ylab="Fecundity [eggs / raft]", xlim=c(0, 45), ylim=c(0, 250), pch=20, col=suffolk.col,
     main="Fecundity")
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf+data.sderr.suf, col=suffolk.col,
       angle=90, length=0.03)
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf-data.sderr.suf, col=suffolk.col,
       angle=90, length=0.03)

lines(temps, meancurve.suf, type="l", col=suffolk.col)
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


saveRDS(albany.out, file="jagsout_epr_albany.RDS")
saveRDS(suffolk.out, file="jagsout_epr_suffolk.RDS")
