################################################################################
# Setup
################################################################################

# Load required packages
library(magrittr)

# Bladder data
subject <- as.factor(rep(1:23, times = 8))
volume <- rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23) / 10
HD <- c(13.2, 11.1, 10.3, NA, 4.8, 7.7, NA, 5.9, 1.9, 6.5, 19.8,
        14.6, NA, NA, 9.7, 17.2, 10.6, 19.3, 8.5, 6.9, 8.1, 14.8, 13.7,
        27.4, 27.5, 15, 10, 18.6, 12.6, 24, 28.4, 12.5, 16.7, 29.6,
        27.1, 14, 18.7, 20.3, 35.8, 23.6, 37.4, 31.3, 23.7, 22, 34.3,
        28.5, 41.6, 58.1, 34.2, 28.8, 29.9, 31.4, 46.9, 44.4, 26.8,
        30.6, 51.7, 49.8, 19.1, 35.8, 38.9, 41.4, 49.9, 58.6, 54.8, 44,
        39.1, 58.5, 41.5, 60.1, 78.8, 49.4, 46.4, 39.4, 45.3, 50.4,
        70.7, 54.4, 41.8, 72.2, 67.5, 39.2, 49.6, 65.1, 69.7, 67.7,
        73.7, 78.3, 65.7, 44.7, 72.1, 59.8, 73.9, 91.5, 71.3, 54.8, NA,
        48, 67.8, 89.4, 63.1, 49.6, 81.9, 79.1, 48.7, 65.6, 65.1, 81.9,
        87.7, 79.4, 93, 80.3, 68.9, 90.9, 77.5, 85.5, 98.3, 81.3, 69.4,
        NA, 66.6, 81, 105.8, 83.5, 60.8, 95.1, 95.1, 67, 85.3, 86.9,
        96.6, 89.3, 102.6, NA, 93.6, 93.3, 105, 92.9, 95.6, 111.4, 94,
        73.9, NA, NA, 91.2, 113.5, 114.5, 80.1, 115.4, 109.8, 72.7,
        90.4, 98.6, 115, 108, 110.9, NA, 99.2, 102.4, 117.5, 99.4,
        107.4, 121, 104.3, NA, NA, NA, 99.8, 127.3, 124, 87.1, NA, NA,
        NA, NA, 107.2, 117, 114.8, 122.4, NA, 112.2, 104.7, 124.2, 113)
bladder <- data.frame(subject = subject, HD = HD, HDlin = HD ^ (3 / 2), 
                      volume = volume)
bladder <- na.omit(bladder)


################################################################################
# Graphics
################################################################################

# Load required packages
library(ggplot2)

# Spaghetti plots of original and transformed data
p1 <- ggplot(bladder, aes(x = volume, y = HD, group = subject)) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.7) +
  theme_light() +
  labs(x = "Volume (cl)") +
  scale_color_discrete(guide = FALSE)
p2 <- ggplot(bladder, aes(x = volume, y = HD ^ (3/2), group = subject)) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.7) +
  theme_light() +
  labs(x = "Volume (cl)", y = expression(HD^{3/2})) +
  scale_color_discrete(guide = FALSE)
gridExtra::grid.arrange(p1, p2, ncol = 2)


################################################################################
# Linear mixed-effects models
################################################################################

# Load required packages
library(investr)
library(lme4)
library(nlme)

# Random intercept and slope model for the transformed data
original.nlme <- lme(HD ~ volume + I(volume ^ 2), 
                     random = list(subject = pdDiag(~volume)),
                     data = bladder)
invest(original.nlme, y0 = 60)

# Random intercept and slope model for the transformed data
bladder.nlme <- lme(HD ^ (3/2) ~ volume, 
                    random = list(subject = pdDiag(~volume)),
                    data = bladder)
invest(bladder.nlme, y0 = 500)

# LMM for original (i.e., untransformed) data
original.lme <- lme(HD ~ volume + I(volume ^ 2), 
                    random = list(subject = pdDiag(~volume)),
                    data = bladder)
invest(original.lme, y0 = 80)


# Possible explanation for the bimodal posterior density
h <- function(x) ifelse(x > 0, x, 0)
fit <- lme(HD ~ volume + I(volume ^ (2/3)), # h(7.5 - volume) + h(volume - 7.5),
           random = list(subject = pdDiag(~volume)),
           data = bladder)
newx <- seq(from = 0, to = 19, length = 500)
pred <- predict(fit, newdata = data.frame(volume = newx), level = 0)
plot(HD ~ volume, data = bladder)
lines(newx, pred)

invest(fit, y0 = 80, upper = 50)
invest(fit, y0 = 80, interval = "Wald")


################################################################################
# Objective Bayes model using JAGS
################################################################################

library(rjags)

# Model file
model <- function() {
  
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau.y)
    mu[i] <- a[subject[i]] + b[subject[i]] * (x[i] - mean(x[])) +
      mu.c * pow(x[i] - mean(x[]), 2)
  }
  y0 ~ dnorm(mu.a0 + mu.b0 * (x0 - mean(x[])) + 
               mu.c * pow(x0 - mean(x[]), 2), tau.y)
  
  # Random effects
  for (j in 1:m) {
    a[j] ~ dnorm(mu.a, tau.a)
    b[j] ~ dnorm(mu.b, tau.b)
  }
  mu.a0 ~ dnorm(mu.a, tau.a)
  mu.b0 ~ dnorm(mu.b, tau.b)
  
  # Priors
  mu.a ~ dnorm(0, 0.0001)  
  mu.b ~ dnorm(0, 0.0001)  
  mu.c ~ dnorm(0, 0.0001)  
  tau.a <- pow(sigma.a, -2)  
  tau.b <- pow(sigma.b, -2)
  tau.y <- pow(sigma.y, -2)
  sigma.a ~ dunif(0, 100)
  sigma.b ~ dunif(0, 100)
  sigma.y ~ dunif(0, 100)
  
  # What's a good prior for x0?
  x0 ~ dnorm(0, 0.0001)%_%I(0, 10000)  # truncated to experimental range

}

# Write model to file
R2OpenBUGS::write.model(model, con = "bladder-bayes-original.txt")

# Data
data.list <- with(bladder, 
  list(y = HD, x = volume, subject = subject, y0 = 80, 
       n = length(HD), m = length(unique(subject)))
)

# Parameter estimates from model fit based on centered volume
bladder$volume.centered <- bladder$volume - mean(bladder$volume)
fit <- lmer(HD ~ volume.centered + I(volume.centered ^ (2/3)) + 
              (0 + 1 | subject) + (0 + volume | subject), data = bladder)
fe <- unname(fixef(fit))
vc <- as.data.frame(VarCorr(fit))$sdcor

## Initial values for chain 1
inits.list <- list(mu.a = fe[1], mu.b = fe[2], mu.c = fe[3], 
                   sigma.a = vc[1], sigma.b = vc[2], sigma.y = vc[3], 
                   x0 = invest(original.lme, y0 = 100, interval = "none"), 
                   .RNG.name = "base::Mersenne-Twister", .RNG.seed = 2)

# JAGS model
sim <- jags.model("bladder-bayes-original.txt", data = data.list, 
                  inits = inits.list)
update(sim, n.iter = 10000)  # burn-in

# Fixed-effects posterior
fe.post <- sim %>%
  coda.samples(c("mu.a", "mu.b", "mu.c"), n.iter = 100000, thin = 10) %T>%
  plot() %>%
  as.matrix()
plot(fe.coda)

# Variance components posterior
vc.post <- sim %>%
  coda.samples(c("tau.a", "tau.b", "tau.y"), n.iter = 100000, thin = 10) %T>%
  plot() %>%
  as.matrix()

# Unknown posterior
x0.post <- sim %>%
  coda.samples("x0", n.iter = 100000, thin = 10) %T>%
  plot() %>%
  as.matrix() %>%
  as.numeric()

# plot(fe.coda)
# plot(vc.coda)
# plot(x0.coda)

# Compare to bootstrap distribution
x0.boot <- na.omit(pb.bladder$t[, 1])
plot(density(x0.boot), lwd = 2, main = "", xlab = "")
lines(density(x0.post), col = "red2", lwd = 2)



## JAGS model for transformed data ---------------------------------------------

## Model file
model.transformed <- function() {
  
  # Likelihood
  for (i in 1:n) {
    y[i] ~ dnorm(mu[i], tau.y)
    mu[i] <- a[subject[i]] + b[subject[i]]*(x[i]-mean(x[]))
  }
  y0 ~ dnorm(mu.a0 + mu.b0*(x0-mean(x[])), tau.y)
  
  # Random effects
  for (j in 1:m) {
    a[j] ~ dnorm(mu.a, tau.a)
    b[j] ~ dnorm(mu.b, tau.b)
  }
  mu.a0 ~ dnorm(mu.a, tau.a)
  mu.b0 ~ dnorm(mu.b, tau.b)
  
  # Priors
  mu.a ~     dnorm(0, 0.0001)  
  mu.b ~     dnorm(0, 0.0001)  
  tau.a <-   pow(sigma.a, -2)  
  tau.b <-   pow(sigma.b, -2)
  tau.y <-   pow(sigma.y, -2)
  sigma.a ~  dunif(0, 100)
  sigma.b ~  dunif(0, 100)
  sigma.y ~  dunif(0, 100)
  
  # What's a good prior for x0?
  # x0 ~ dnorm(0, 0.0001)%_%I(1.0, 17.5)
  x0 ~ dnorm(0, 0.0001)%_%I(0, 10000)
  # x0 ~ dlnorm(0, 0.01)%_%I(0, 25)
  # x0 ~ dunif(1.0, 17.5)
  
}

# Write model to file
R2OpenBUGS::write.model(model.transformed, 
                        con = "bladder-bayes-transformed.txt")

## Data
data.list2 <- with(bladder, 
  list(y = HD ^ (3/2), x = volume, subject = subject, y0 = 80^(3/2), 
       n = length(HD), m = length(unique(subject)))
)

## Parameter estimates from model fit based on centered volume
fit2 <- lmer(HD ^ (3/2) ~ I(volume - mean(volume)) + (0 + 1 | subject) + 
               (0 + volume | subject), data = bladder)
fe2 <- unname(fixef(fit2))
vc2 <- as.data.frame(VarCorr(fit2))$sdcor

## Initial values for chain 1
inits.list2 <- list(mu.a = fe2[1], 
                    mu.b = fe2[2], 
                    sigma.a = vc2[1], 
                    sigma.b = vc2[2], 
                    sigma.y = vc2[3], 
                    x0 = 11.13502, 
                    .RNG.name = "base::Mersenne-Twister", 
                    .RNG.seed = 2)

## JAGS model
sim2 <- jags.model("bladder-bayes-transformed.txt", data = data.list2, 
                   inits = inits.list2)
update(sim2, n.iter = 10000)  # burn-in

# Fixed-effects posterior
fe.post <- sim %>%
  coda.samples(c("mu.a", "mu.b", "mu.c"), n.iter = 100000, thin = 10) %T>%
  plot() %>%
  as.matrix()
plot(fe.coda)

# Variance components posterior
vc.post <- sim2 %>%
  coda.samples(c("tau.a", "tau.b", "tau.y"), n.iter = 100000, thin = 10) %T>%
  plot() %>%
  as.matrix()

# Unknown posterior
x0.post <- sim2 %>%
  coda.samples("x0", n.iter = 100000, thin = 10) %T>%
  plot() %>%
  as.matrix() %>%
  as.numeric()

## Compare to bootstrap distribution
x0.boot2 <- na.omit(pb.bladder2$t[, 1])
plot(density(x0.boot2), lwd = 2, main = "", xlab = "")
lines(density(x0.post2), col = "red2", lwd = 2)

## Save results
post.bladder <- x0.post
post.bladder2 <- x0.post2
save(bladder, bladder2, pb.bladder, pb.bladder2, post.bladder,
     file = "/home/w108bmg/Desktop/Dropbox/bladder.RData")