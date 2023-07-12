# R0 calculation

library('R2jags')
library('MCMCvis')
library('coda')

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

plot_curves <- function(temps, curves) {
  meancurve <- apply(curves, 1, mean)
  CI <- apply(curves, 1, quantile, c(0.025, 0.975))
  plot(temps, meancurve, type="l", col="black")
  lines(temps, CI[1,], col="gray")
  lines(temps, CI[2,], col="gray")
}

## 1. Load chains with posterior estimates of parameters and calculate trait
## values.

temps <- seq(10, 40, 0.025) # Test with small vector first to debug.

# Load data for biting rate
a.suf.chains <- readRDS("jagsout_a_suffolk.RDS")
a.alb.chains <- readRDS("jagsout_a_albany.RDS")

a.suf.chains
a.alb.chains
a.suf.chains$BUGSoutput$summary

# Only keep chains for TPC model parameters to save memory.
a.suf.chains <- MCMCchains(a.suf.chains, params=c("Tmin", "Tmax", "c"))
a.alb.chains <- MCMCchains(a.alb.chains, params=c("Tmin", "Tmax", "c"))

# Calculate trait values for each sample at all temperatures.
a.suf.curves <-  apply(a.suf.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))
a.alb.curves <-  apply(a.alb.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))

#plot_curves(temps, a.suf.curves)
#plot_curves(temps, a.alb.curves)

# Vector competence (from Shocket et al)
load("jagsout_bc_CpipWNV_inf.Rdata")
bc.CpipWNV.out.inf$BUGSoutput$summary
bc.chains <- MCMCchains(bc.CpipWNV.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
rm(bc.CpipWNV.out.inf)

bc.curves <- apply(bc.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))
head(bc.chains)

#plot_curves(temps, bc.curves)

# Lifespan
lf.suf.chains <- readRDS("jagsout_lf_suffolk.RDS")
lf.alb.chains <- readRDS("jagsout_lf_albany.RDS")

lf.suf.chains <- MCMCchains(lf.suf.chains, params=c("m", "b"))
lf.alb.chains <- MCMCchains(lf.alb.chains, params=c("m", "b"))

lf.suf.curves <- apply(lf.suf.chains, 1, function(x) linear_lim(temps, x[1], x[2]))
lf.alb.curves <- apply(lf.alb.chains, 1, function(x) linear_lim(temps, x[1], x[2]))

#plot_curves(temps, lf.suf.curves)
#plot_curves(temps, lf.alb.curves)

# Pathogen development rate (from Shocket et al)
load("jagsout_PDR_CpipWNV_inf.Rdata")

PDR.chains <- MCMCchains(PDR.CpipWNV.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
rm(PDR.CpipWNV.out.inf)

PDR.curves <- apply(PDR.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))

#plot_curves(temps, PDR.curves)

# Fecundity (eggs per raft)
ER.suf.chains <- readRDS("jagsout_epr_suffolk.RDS")
ER.alb.chains <- readRDS("jagsout_epr_albany.RDS")

ER.suf.chains <- MCMCchains(ER.suf.chains, params=c("Tmin", "Tmax", "c"))
ER.alb.chains <- MCMCchains(ER.alb.chains, params=c("Tmin", "Tmax", "c"))

ER.suf.curves <- apply(ER.suf.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))
ER.alb.curves <- apply(ER.alb.chains, 1, function(x) quad(temps, x[1], x[2], x[3]))

#plot_curves(temps, ER.suf.curves)
#plot_curves(temps, ER.alb.curves)


# Proportion ovipositing.
pO.suf.chains <- readRDS("jagsout_pO_suffolk.RDS")
pO.alb.chains <- readRDS("jagsout_pO_albany.RDS")

pO.suf.chains <- MCMCchains(pO.suf.chains, params=c("Tmin", "Tmax", "c"))
pO.alb.chains <- MCMCchains(pO.alb.chains, params=c("Tmin", "Tmax", "c"))

pO.suf.curves <- apply(pO.suf.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
pO.alb.curves <- apply(pO.alb.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

#plot_curves(temps, pO.suf.curves)
#plot_curves(temps, pO.alb.curves)


# Egg viability (from Shocket et al)
load("jagsout_EV_Cpip_inf.Rdata")
EV.chains <- MCMCchains(EV.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
rm(EV.Cpip.out.inf)

EV.curves <- apply(EV.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

#plot_curves(temps, EV.curves)


# Larval survival
pLA.suf.chains <- readRDS("jagsout_pLA_suffolk.RDS")
pLA.alb.chains <- readRDS("jagsout_pLA_albany.RDS")

pLA.suf.chains <- MCMCchains(pLA.suf.chains, params=c("Tmin", "Tmax", "c"))
pLA.alb.chains <- MCMCchains(pLA.alb.chains, params=c("Tmin", "Tmax", "c"))

# Sample random rows to make chains the same length as others.
pLA.suf.chains <- pLA.suf.chains[sample(nrow(pLA.suf.chains), 93750), ]
pLA.alb.chains <- pLA.alb.chains[sample(nrow(pLA.alb.chains), 93750), ]

pLA.suf.curves <- apply(pLA.suf.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))
pLA.alb.curves <- apply(pLA.alb.chains, 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

#plot_curves(temps, pLA.suf.curves)
#plot_curves(temps, pLA.alb.curves)

# Mosquito development rate
MDR.suf.chains <- readRDS("jagsout_MDR_suffolk.RDS")
MDR.alb.chains <- readRDS("jagsout_MDR_albany.RDS")

MDR.suf.chains <- MCMCchains(MDR.suf.chains, params=c("Tmin", "Tmax", "c"))
MDR.alb.chains <- MCMCchains(MDR.alb.chains, params=c("Tmin", "Tmax", "c"))

# Sample random rows to make chains the same length as others.
MDR.suf.chains <- MDR.suf.chains[sample(nrow(MDR.suf.chains), 93750), ]
MDR.alb.chains <- MDR.alb.chains[sample(nrow(MDR.alb.chains), 93750), ]


MDR.suf.curves <- apply(MDR.suf.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))
MDR.alb.curves <- apply(MDR.alb.chains, 1, function(x) briere(temps, x[1], x[2], x[3]))


#plot_curves(temps, MDR.suf.curves)
#plot_curves(temps, MDR.alb.curves)


R0.suf.curves <- R0(a.suf.curves, bc.curves, lf.suf.curves, 
                    PDR.curves, ER.suf.curves, pO.suf.curves, EV.curves, 
                    pLA.suf.curves, MDR.suf.curves)
R0.alb.curves <- R0(a.alb.curves, bc.curves, lf.alb.curves, 
                    PDR.curves, ER.alb.curves, pO.alb.curves, EV.curves, 
                    pLA.alb.curves, MDR.alb.curves)

suffolk.col <- "#CC79A7"
albany.col <- "#56B4E9" 

# Find mean curves and credible intervals.

meancurve.suf <- apply(R0.suf.curves, 1, mean)
CI.suf <- apply(R0.suf.curves, 1, quantile, c(0.025, 0.975))

meancurve.alb <- apply(R0.alb.curves, 1, mean)
CI.alb <- apply(R0.alb.curves, 1, quantile, c(0.025, 0.975))

plot(temps, meancurve.suf, col=suffolk.col, type="l", main="R0", 
     xlab="Temperature [°C]", ylab="Relative R0",
     xlim=c(12, 35), ylim=c(0, 1.0))
lines(temps, CI.suf[1,], col=alpha(suffolk.col, 0.5))
lines(temps, CI.suf[2,], col=alpha(suffolk.col, 0.5))

lines(temps, meancurve.alb, col=albany.col)
lines(temps, CI.alb[1,], col=alpha(albany.col, 0.5))
lines(temps, CI.alb[2,], col=alpha(albany.col, 0.5))


find_Topt <- function(temps, curves) {
  idx <- apply(curves, 2, which.max)
  return(temps[idx])
}

# Find optimal temperature.
R0.suf.Topt <- find_Topt(temps, R0.suf.curves)
R0.alb.Topt <- find_Topt(temps, R0.alb.curves)

print(mean(R0.suf.Topt))
print(quantile(R0.suf.Topt, c(0.025, 0.975)))

print(mean(R0.alb.Topt))
print(quantile(R0.alb.Topt, c(0.025, 0.975)))

saveRDS(R0.suf.curves, file="R0_chains_suffolk.RDS")
saveRDS(R0.alb.curves, file="R0_chains_albany.RDS")


