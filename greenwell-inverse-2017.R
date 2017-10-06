################################################################################
# Setup
################################################################################

# Load required packages
library(boot)
library(investr)
library(ggplot2)
library(lme4)
library(magrittr)
library(nlme)
library(RColorBrewer); cols <- brewer.pal(9, "Set1")
library(rjags)

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
bladder <- data.frame(subject = subject, HD = HD, volume = volume)
bladder <- na.omit(bladder)  # omit missing values for this analysis


################################################################################
# Basic graphics
################################################################################

# Scatterplot of transformed data
p1 <- ggplot(bladder, aes(x = volume, y = HD ^ (3 / 2))) +
  geom_point(alpha = 0.7) +
  theme_light() +
  labs(x = "Volume (cl)", y = expression(HD^{3/2}))

# Spaghettiplot of transformed data
p2 <- ggplot(bladder, aes(x = volume, y = HD ^ (3 / 2), 
                          color = subject, group = subject)) +
  geom_line(alpha = 0.4) +
  geom_point(alpha = 0.7) +
  theme_light() +
  labs(x = "Volume (cl)", y = expression(HD^{3/2})) +
  scale_color_discrete(guide = FALSE)

# Layout plots on a grid
pdf("spaghettiplots.pdf", width = 8, height = 4)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()


################################################################################
# Linear mixed-effects models
################################################################################

# LMM for transformed data using nlme package
mod1 <- lme(HD ^ (3 / 2) ~ volume, random = list(subject = pdDiag( ~ volume)),
            data = bladder)

# knitr::kable(summary(transformed.lme)$tTable, format = "latex", digits = 3)

# LMM for transformed data using lme4 package
mod2 <- lmer(HD ^ (3 / 2) ~ volume + (0 + 1 | subject) + (0 + volume | subject), 
             data = bladder)


################################################################################
# Inverse estimation
################################################################################

# Point estimate
invest(mod1, y0 = 110 ^ (3/2), interval = "none")
# [1] 8.015521

# Wald-based interval
invest(mod1, y0 = 110 ^ (3/2), interval = "Wald")
# estimate     lower     upper        se 
# 8.015521  4.185376 11.845665  1.954191

# Estimated variance of Y0
z <- data.frame(volume = invest(mod1, y0 = 110 ^ (3/2), interval = "none"))
investr:::varY(mod1, newdata = z)

# Standard erros using the car package
params <- c(fixef(mod1), 110 ^ (3/2))
covmat <- diag(3)
covmat[1:2, 1:2] <- vcov(mod1)
covmat[3, 3] <-  investr:::varY(mod1, newdata = z)
names(params) <- c("b0", "b1", "y0")
car::deltaMethod(params, g = "(y0 - b0) / b1",  vcov. = covmat)$SE

# Inversion interval
invest(mod1, y0 = 110 ^ (3/2), upper = 30)
# estimate     lower     upper 
# 8.015521  4.227965 11.919182

# Inversion interval using different distribution
tvals <- qt(c(0.025, 0.975), df = nrow(bladder.na) - 1)
invest(mod1, y0 = 110 ^ (3/2), q1 = tvals[1L], q2 = tvals[2L], upper = 30)
# estimate     lower     upper 
# 8.015521  4.200185 11.948696 


################################################################################
# Parametric bootstrap
################################################################################

# Variance function
varY <- function(object, x) {
  vc <- as.data.frame(lme4::VarCorr(object))$vcov
  vc[1] + vc[2] * x ^ 2 + vc[3]
}

# Point estimate from lme4 fit
fe <- unname(fixef(mod2))  # fixed-effects
(x0.est <- (110 ^ (3/2) - fe[1]) / fe[2])

# fe3 <- unname(fixef(mod3))  # fixed-effects
# (x0.est <- ((110 - fe[1]) / fe[2]) ^ (3/2))


# Parametric bootstrap replicates
pboot <- bootMer(transformed.lmer, nsim = 999*2, seed = 105, FUN = function(.) {
  
  # Point estimate
  var.Y0.boot <- varY(., x = x0.est)
  fe.boot <- unname(fixef(.))
  if (all(getME(., "y") == bladder.na$HD)) {
    y0.boot <- 110 ^ (3/2)
  } else {
    y0.boot <- rnorm(1, mean = 110 ^ (3/2), sd = sqrt(var.Y0.boot))
  }
  (y0.boot - fe.boot[1L]) / fe.boot[2L]
})

# Save bootstrap replicates
# save(pboot, file = "pboot.RData")

plot(pboot)

# Summary of bootstrap replicates
summary(pboot)

boot.ci(pboot, type = c("norm", "perc", "stud"))
