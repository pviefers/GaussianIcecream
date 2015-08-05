# Markus Gesmann & Paul Viefers
library(arm) # for 'display' function only
icecream <- data.frame(
# http://www.statcrunch.com/5.0/viewreport.php?reportid=34965&groupid=1848
  temp=c(11.9, 14.2, 15.2, 16.4, 17.2, 18.1, 
         18.5, 19.4, 22.1, 22.6, 23.4, 25.1),
  units=c(185L, 215L, 332L, 325L, 408L, 421L, 
          406L, 412L, 522L, 445L, 544L, 614L)
  )
basicPlot <- function(...){
  plot(units ~ temp, data=icecream, bty="n", lwd=2,
       main="Number of ice creams sold", col="#00526D", 
       xlab="Temperature (Celsius)", 
       ylab="Units sold", ...)
  axis(side = 1, col="grey")
  axis(side = 2, col="grey")
}
basicPlot()

## ------------------------------------------------------------------------
basicPlot()
lsq.mod <- lsfit(icecream$temp, icecream$units)
abline(lsq.mod, col="orange", lwd=2)
legend(x="topleft", bty="n", lwd=c(2,2), lty=c(NA,1),
       legend=c("observation", "linear least square"),
         col=c("#00526D","orange"),  pch=c(1,NA))

## ------------------------------------------------------------------------
lin.mod <- glm(units ~ temp, data=icecream, 
              family=gaussian(link="identity"))
display(lin.mod)

## ------------------------------------------------------------------------
log.lin.mod <- glm(log(units) ~ temp, data=icecream, 
              family=gaussian(link="identity"))
display(log.lin.mod)
log.lin.sig <- summary(log.lin.mod)$dispersion
log.lin.pred <- exp(predict(log.lin.mod) + 0.5 * log.lin.sig)
basicPlot()
lines(icecream$temp, log.lin.pred, col="red", lwd=2)
legend(x="topleft", bty="n", lwd=c(2,2), lty=c(NA,1),
       legend=c("observation", "log-transformed LM"),
         col=c("#00526D","red"),  pch=c(1,NA))

## ------------------------------------------------------------------------
exp(coef(log.lin.mod)[1])

## ------------------------------------------------------------------------
pois.mod <- glm(units ~ temp, data=icecream, 
              family=poisson(link="log"))
display(pois.mod)
pois.pred <- predict(pois.mod, type="response")
basicPlot()
lines(icecream$temp, pois.pred, col="blue", lwd=2)
legend(x="topleft", bty="n", lwd=c(2,2), lty=c(NA,1),
       legend=c("observation", "Poisson (log) GLM"),
         col=c("#00526D","blue"),  pch=c(1,NA))

## ------------------------------------------------------------------------
predict(pois.mod, newdata=data.frame(temp=32), type="response")

## ------------------------------------------------------------------------
market.size <- 800
icecream$opportunity <- market.size - icecream$units
bin.glm <- glm(cbind(units, opportunity) ~ temp, data=icecream, 
    family=binomial(link = "logit"))
display(bin.glm)
bin.pred <- predict(bin.glm, type="response")*market.size
basicPlot()
lines(icecream$temp, bin.pred, col="purple", lwd=2)
legend(x="topleft", bty="n", lwd=c(2,2), lty=c(NA,1),
       legend=c("observation", "Binomial (logit) GLM"),
         col=c("#00526D","purple"),  pch=c(1,NA))

## ------------------------------------------------------------------------
# Sales at 0 Celsius
plogis(coef(bin.glm)[1])*market.size
# Sales at 35 Celsius
plogis(coef(bin.glm)[1] +  coef(bin.glm)[2]*35)*market.size

## ------------------------------------------------------------------------
# Paul Viefers
# Gaussian process with binomial link function in Stan
x1 <- seq(0, icecream$temp[1])
N1 <- length(x1)
idx1 <- 1:N1
x2 <- icecream$temp
N2 <- nrow(icecream)
idx2 <- (N1+1):(N1+N2)
x3 <- seq(icecream$temp[N2]+1.3, 40, length.out = 12)
N3 <- length(x3)
idx3 <- (N1+N2+1):(N1+N2+N3)
x <- c(x1, x2, round(x3, 1))
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

stan_data <- list(N1 = N1,
                  N2 = N2,
                  N3 = N3,
                  K=800,
                  x = x,
                  y2 = icecream$unit
)

ice.fit <- stan(file="gp-BIN-fit_v05.stan",   
     data=stan_data,
     iter=0, chains=0)

ice.res <- stan(fit = ice.fit, 
                data=stan_data,
                iter = 1000,
                warmup = 500,
                chains = 2,
                seed=123)

posterior.prediction <- extract(ice.res, pars = 'n_sold', permute = TRUE)
median.prediction <- apply(ys$n_sold, 2, median)
prediction.bands <- apply(posterior.prediction$n_sold, 2, quantile, probs = c(0.025,0.05,0.95, 0.975))
# Show the fit of the model against actual data
# Markus' basic plot with the Gaussian process against actual data
basicPlot()
lines(x[idx2], median.prediction[idx2],
      col = "red",
      lwd=2
      )
lines(x[idx2], prediction.bands[1, idx2], lty = 2, col = "grey")
lines(x[idx2], prediction.bands[2, idx2], lty = 2)
lines(x[idx2], prediction.bands[3, idx2], lty = 2)
lines(x[idx2], prediction.bands[4, idx2], lty = 2, col = "grey")

## ------------------------------------------------------------------------
temp <- 0:35
p.lm <- predict(lin.mod, data.frame(temp=temp), type="response")
p.log.lm <- exp(predict(log.lin.mod, data.frame(temp=0:35), type="response") + 
                  0.5 * summary(log.lin.mod)$dispersion)
p.pois <- predict(pois.mod, data.frame(temp=temp), type="response")
p.bin <- predict(bin.glm, data.frame(temp=temp), type="response")*market.size 
basicPlot(xlim=range(temp), ylim=c(-20,market.size))
lines(temp, p.lm, type="l", col="orange", lwd=2)
lines(temp, p.log.lm, type="l", col="red", lwd=2)
lines(temp, p.pois, type="l", col="blue", lwd=2)
lines(temp, p.bin, type="l", col="purple", lwd=2)
lines(x, median.prediction, col = "green", lwd=2)
legend(x="topleft", 
       legend=c("observation", 
                "linear model",
                "log-transformed LM",
                "Poisson (log) GLM",
                "Binomial (logit) GLM",
                "Binomial Gaussian Process"),
       col=c("#00526D","orange", "red", 
             "blue", "purple", "green"),  
       bty="n", lwd=rep(2,5), 
       lty=c(NA,rep(1,5)),
       pch=c(1,rep(NA,5)))

## ----Simulations---------------------------------------------------------
n <- nrow(icecream)
A <- model.matrix(units ~ temp, data=icecream)
set.seed(1234)
(rand.normal <- rnorm(n,
                     mean = A %*% coef(lin.mod),
                     sd = sqrt(summary(lin.mod)$dispersion)))
(rand.logtrans <- rlnorm(n,
                         meanlog = A %*% coef(log.lin.mod),
                         sdlog =  sqrt(summary(log.lin.mod)$dispersion)))
(rand.pois <- rpois(n,
                   lambda = exp(A %*% coef(pois.mod))))
(rand.bin <- rbinom(n,
                   size = market.size,
                   prob = plogis(A %*% coef(bin.glm))))
basicPlot(ylim=c(100,700))
cols <- adjustcolor(c("orange", "red", "blue", "purple"), 
                    alpha.f = 0.75)
points(icecream$temp, rand.normal, pch=19, col=cols[1])
points(icecream$temp, rand.logtrans, pch=19, col=cols[2])
points(icecream$temp, rand.pois, pch=19, col=cols[3])
points(icecream$temp, rand.bin, pch=19, col=cols[4])
legend(x="topleft",
       legend=c("observation",
                "linear model",
                "log-transformed LM",
                "Poisson (log) GLM",
                "Binomial (logit) GLM"),
       col=c("#00526D",cols), 
       lty=NA,
       bty="n", lwd=rep(2,5),
       pch=c(1,rep(19,4)))
## ----Posterior predictions using Gaussian process------------------------
# Plot over the full
basicPlot(xlim=c(0,max(x)), 
     ylim=c(-5, 800));
lines(x, median.prediction,
      col = "red",
      lwd=2);
abline(h=0, lwd=2)
lines(x, prediction.bands[1, ], lty = 2, col = "grey")
lines(x, prediction.bands[2, ], lty = 2)
lines(x, prediction.bands[3, ], lty = 2)
lines(x, prediction.bands[4, ], lty = 2, col = "grey")