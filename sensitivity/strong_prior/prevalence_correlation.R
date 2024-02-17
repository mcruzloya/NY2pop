# Calculate correlation between mosquito data and R0.
library('scales')
library('Hmisc')
setwd("/Users/cruzloya/git/NY2pop/sensitivity/strong_prior/")

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

suffolk.22$lagged.tmeanc <- as.vector(Lag(suffolk.22$tmeanc, shift=2))
erie.22$lagged.tmeanc <- as.vector(Lag(erie.22$tmeanc, shift=2))

suffolk.22$tmeanc
suffolk.22$lagged.tmeanc


R0.suf.curves <- readRDS("R0_chains_suffolk_sp.RDS")
R0.alb.curves <- readRDS("R0_chains_albany_sp.RDS")

R0.suf.mean <- apply(R0.suf.curves, 1, mean)
R0.suf.CI <- apply(R0.suf.curves, 1, quantile, c(0.025, 0.975))

R0.alb.mean <- apply(R0.alb.curves, 1, mean)
R0.alb.CI <- apply(R0.alb.curves, 1, quantile, c(0.025, 0.975))


plot(suffolk.22$lagged.tmeanc, suffolk.22$Infection.Rate, pch=20, col=suffolk.col,
     xlab="Temperature [ºC]", ylab="WNV prevalence",  xlim = c(10, 35), ylim=c(0, 20),
     main="")
points(erie.22$lagged.tmeanc, erie.22$Infection.Rate, pch=20, col=erie.col)
lines(temps, 20 * R0.suf.mean, col=suffolk.col, lwd=2)
polygon(c(temps, rev(temps)), c(20 * R0.suf.CI[1,], rev(20 * R0.suf.CI[2,])), 
        col=alpha(suffolk.col, 0.5), lty=0)

lines(temps, 20 * R0.alb.mean, col=erie.col, lwd=2)
polygon(c(temps, rev(temps)), c(20 * R0.alb.CI[1,], rev(20 * R0.alb.CI[2,])), 
        col=alpha(erie.col, 0.5), lty=0)


correlation <- function(R0.temps, R0.curve, incidence.temps, incidence) {
  return(cor(approx(R0.temps, R0.curve, incidence.temps, yleft=0, yright=0)$y, 
             incidence, use="complete.obs", method="spearman"))
}
suf.cor <- apply(R0.suf.curves, 2, function(x) correlation(temps, x, suffolk.22$lagged.tmeanc,
                                                           suffolk.22$Infection.Rate))
erie.cor <- apply(R0.alb.curves, 2, function(x) correlation(temps, x, erie.22$lagged.tmeanc,
                                                            erie.22$Infection.Rate))
print(mean(suf.cor, na.rm=TRUE))
print(quantile(suf.cor, c(0.025, 0.975), na.rm=TRUE))

print(mean(erie.cor, na.rm=TRUE))
print(quantile(erie.cor, c(0.025, 0.975), na.rm=TRUE))



