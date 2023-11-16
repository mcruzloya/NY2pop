set.seed(42) # Set seed for reproducibility.

library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')
library('ggplot2')

linear <- function(T, m, b) {
  result <- (-m * T + b) * (T < b/m)
  result[T < 15] <- (-m * 15 + b) * (T < b/m) # Set maximum value at 15 degrees
  return(result)
}

data <- read.csv("lf_allpop.csv")
data

data.suffolk <- subset(data, (data$population == "suffolk") )
data.albany <- subset(data, (data$population == "albany")  )

suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9"  

hist(subset(data.suffolk, data$temperature == 28)$lifespan)

# Visualize data.
g <- (ggplot(data, aes(temperature, lifespan)) +  
        geom_violin(aes(col=population, group=cut_width(temperature, 1))) + 
        geom_point(aes(col=population, alpha=I(0.3))) +
        facet_wrap(~population, nrow=1) +
       theme_classic() + scale_color_manual(values=c(albany.col, suffolk.col)))
g

head(data.suffolk)
plot(data.suffolk$temperature, data.suffolk$lifespan, xlim=c(10, 35), ylim=c(0, 60),
     pch=20, main="Suffolk", col=alpha(suffolk.col, 0.3) )
plot(data.albany$temperature, data.albany$lifespan, xlim=c(10, 35), ylim=c(0, 60),
     pch=20, main="Albany", col=alpha(albany.col, 0.3) )

sink("linear_lf.txt")
cat("
    model{
    
    ## Priors
    m_mu <- 4.86
    m_sd <- 3.0
    m ~ dgamma(m_mu^2 / m_sd^2, m_mu / m_sd^2)
    
    Tmax ~ dnorm(34.9, 1 / 2.5^2)
    b <- m * Tmax
    
    mu_sigma <- 20
    sd_sigma <- 10
    
    # Standard deviation at 25 degrees Celsius
    sigma25 ~ dgamma(mu_sigma^2 / sd_sigma^2, mu_sigma / sd_sigma^2)
    
    beta ~ dnorm(0, 1)

    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <- (-1 * m * temp[i] + b)
    sigma[i] <- sigma25 * exp(beta * (temp[i] - 25))
    trait[i] ~ dnorm(trait.mu[i], 1 / sigma[i]^2)
    }
    
    } # close model
    ",fill=TRUE)
sink()

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 3 # number of chains

##### Organize Data for JAGS
trait.suf <- data.suffolk$lifespan
temp.suf <- data.suffolk$temperature
N.obs.suf <- length(trait.suf)

trait.alb <- data.albany$lifespan
temp.alb <- data.albany$temperature
N.obs.alb <- length(trait.alb)

##### Calculate initial values for MCMC.
inits<-function(){list(
  m = runif(1, min=0.5, max=5),
  Tmax = runif(1, min=30, max=40),
  sigma25 = runif(1, min=5, max=15),
  beta = rnorm(1, mean=-0.2, sd=0.1))}

##### Parameters to Estimate
parameters <- c("m", "b", "Tmax", "sigma25", "beta")

##### Bundle Data
jag.data<-list(trait = trait.suf, N.obs = N.obs.suf, temp = temp.suf)
##### Run JAGS
suffolk.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                    model.file="linear_lf.txt", n.thin=nt, n.chains=nc, 
                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
suffolk.out
mcmcplot(suffolk.out)

chains.suf <- MCMCchains(suffolk.out, params=c("m", "b"))
temps <- seq(0, 45, 0.1) 
curves <- apply(chains.suf, 1, function(x) linear(temps, x[1], x[2]))

# Find mean curve and credible intervals.
meancurve.suf <- apply(curves, 1, mean)
CI.suf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.suffolk$temperature, data.suffolk$lifespan, xlab="Temperature [°C]",
     ylab="Lifespan", xlim=c(0, 45), ylim=c(0, 80), pch=20, col=alpha(suffolk.col, 0.5), 
     main="Suffolk")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.suf, col="black")
lines(temps, CI.suf[1,], col="gray")
lines(temps, CI.suf[2,], col="gray")


##### Bundle Data
jag.data<-list(trait = trait.alb, N.obs = N.obs.alb, temp = temp.alb)
##### Run JAGS
albany.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
                   model.file="linear_lf.txt", n.thin=nt, n.chains=nc, 
                   n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
albany.out
mcmcplot(albany.out)

chains.alb <- MCMCchains(albany.out, params=c("m", "b"))
temps <- seq(0, 45, 0.1) 
curves <- apply(chains.alb, 1, function(x) linear(temps, x[1], x[2]))

# Find mean curve and credible intervals.
meancurve.alb <- apply(curves, 1, mean)
CI.alb <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.albany$temperature, data.albany$lifespan, xlab="Temperature [°C]",
     ylab="lf", xlim=c(0, 45), ylim=c(0, 80), pch=20, col=alpha(albany.col, 0.3), 
     main="Albany")

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.alb, col="black")
lines(temps, CI.alb[1,], col="gray")
lines(temps, CI.alb[2,], col="gray")

# Plot both curves simultaneously.


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


plot(c(22,25,28), data.mean.suf, xlab="Temperature [°C]",
     ylab="Lifespan [days]", xlim=c(0, 45), ylim=c(0, 60), pch=20, col=suffolk.col,
     main="Lifespan")
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf+data.sderr.suf, col=suffolk.col,
       angle=90, length=0.03)
arrows(c(22,25,28), data.mean.suf, y1=data.mean.suf-data.sderr.suf, col=suffolk.col,
       angle=90, length=0.03)

lines(temps, meancurve.suf, type="l", xlab="Temperature [°C]",
     ylab="Lifespan [days]", xlim=c(0, 45), ylim=c(0, 60), pch=20, col=suffolk.col,
     main="Lifespan (variable SD)")
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


saveRDS(albany.out, file="jagsout_lf_albany.RDS")
saveRDS(suffolk.out, file="jagsout_lf_suffolk.RDS")

