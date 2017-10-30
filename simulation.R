################################################################################
# Setup
################################################################################

# Load packages
library(ggplot2)
library(investr)
library(lme4)
library(nlme)
library(plyr)
library(tidyr)

# Bladder volume by ultrasound data
bladder <- data.frame(
  subject = as.factor(rep(1:23, times = 8)),
  volume = rep(c(10, 25, 50, 75, 100, 125, 150, 175), each = 23),
  HD = c(13.2, 11.1, 10.3, NA, 4.8, 7.7, NA, 5.9, 1.9, 6.5, 19.8, 
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
         NA, NA, 107.2, 117, 114.8, 122.4, NA, 112.2, 104.7, 124.2, 113))
bladder <- na.omit(bladder)
bladder$volume <- bladder$volume/10  # convert to cL


################################################################################
# Spaghetti plots of both the original and transformed data
################################################################################

# Original data
p1 <- ggplot(bladder, aes(volume, HD, group = subject)) +
  geom_line(aes(color = subject), alpha = 0.75) +
  theme_light() +
  theme(legend.position = "none") +
  xlab("Volume (cl)") +
  ylab("Height (mm) times Depth (mm)")

# Transformed data
p2 <- ggplot(bladder, aes(volume, HD^(3/2), group = subject)) +
  geom_line(aes(color = subject), alpha = 0.75) +
  theme_light() +
  theme(legend.position = "none") +
  xlab("Volume (cl)") +
  ylab("Height (mm) times Depth (mm)")

# Display both plots side-by-side
gridExtra::grid.arrange(p1, p2, ncol = 2)


################################################################################
# Linear mixed-effects models
################################################################################

# Explore random effects structure
bladder.lmList <- lmList(HD^(3/2) ~ I(volume - mean(volume)) | subject, 
                         data = bladder)
plot(intervals(bladder.lmList))

# Fit random intercept and slope models with uncorrelated random effects
bladder.nlme <- lme(HD^(3/2) ~ volume, data = bladder,
                    random = list(subject = pdDiag(~volume)))
bladder.lme4 <- lmer(HD^(3/2) ~ volume + (0+1|subject) + (0+volume|subject), 
                     data = bladder)


################################################################################
# Helper functions
################################################################################

# Function to estimate response variance from an "lmerMod" object as a function 
# of volume
est_var <- function(object, x) {
  vc <- as.data.frame(VarCorr(object))$sdcor
  vc[1]^2 + vc[2]^2*x^2 + vc[3]^2
}

# Estimated standard deviation of HD^(3/2) when volume = 8.015521 cL
sd.y0 <- sqrt(est_var(bladder.lme4, x = 8.015521))


# Function to simulate data with different sample sizes
sim_data <- function(m, n) {
  volume <- rep(seq(from = 1, to = 17.5, length = n), times = m)
  subject <- rep(1:m, each = n)
  theta <- as.data.frame(VarCorr(bladder.lme4))$sdcor
  beta <- getME(bladder.lme4, "beta")
  alpha0 <- rnorm(m, mean = beta[1L], sd = theta[1L])
  alpha1 <- rnorm(m, mean = beta[2L], sd = theta[2L])
  HD <- rnorm(m*n, mean = alpha0[subject] + alpha1[subject] * volume, 
              sd = theta[3L])
  data.frame(subject = as.factor(subject), HD, volume)
}

# # Check simulated data
# ggplot(sim_data(50, 50), aes(volume, HD, group = subject)) +
#   geom_line(alpha = 0.5) +
#   theme_light() +
#   theme(legend.position = "none") +
#   xlab("Volume (cl)") +
#   ylab("Height (mm) times Depth (mm)")


# Function to generate a list of simulated data frames
make_data_list <- function(nsim, m, n) {
  set.seed(as.numeric(paste0(m, n)))  # for reproducibility
  rlply(nsim, sim_data(m, n))
}


# Function to fit LMM
fit_lmm <- function(x) {
  lme(HD ~ volume, data = x, random = list(subject = pdDiag(~volume)), 
      control = list(opt = "optim"))
}

# Function to compute asymptotic confidence interval
compute_ci <- function(x, interval = "Wald") {
  eta <- rnorm(1, mean = 500, sd = 132.5607)
  res <- invest(x, y0 = eta, interval = interval, mean.response = FALSE, 
                lower = -20, upper = 50)
  unlist(res[c("lower", "upper")])
}


# Function to simulate CIs
sim_cis <- function(nsim, m, n) {
  
  # List of nsim simulated data sets
  datasets <- make_data_list(nsim, m = m, n = n)
  
  # List of nsim fitted randon intercept and slope models
  lmms <- lapply(datasets, FUN = fit_lmm)
  
  # Array of Wald confidence intervals
  wald.cis <- do.call(rbind, lapply(lmms, FUN = function(x) {
    compute_ci(x, interval = "Wald")
  }))
  colnames(wald.cis) <- c("lwr", "upr")
  wald.cis <- data.frame("method" = "wald", wald.cis)
  
  # Array of inversion confidence intervals
  inversion.cis <- do.call(rbind, lapply(lmms, FUN = function(x) {
    compute_ci(x, interval = "Wald")
  }))
  colnames(inversion.cis) <- c("lwr", "upr")
  inversion.cis <- data.frame("method" = "inversion", inversion.cis)
  
  # Combine results and compute CI length and coverage
  res <- rbind(wald.cis, inversion.cis)
  res$len <- res$upr - res$lwr
  res$cov <- as.numeric(res$lwr < 8.015521 & res$upr > 8.015521)
  
  # Add coverage probability and mean length attributes
  attr(res, "coverage") <- tapply(res$cov, INDEX = res$method, FUN = mean)
  attr(res, "length") <- tapply(res$len, INDEX = res$method, FUN = mean)
  
  # Return results
  res
  
}

# Sample sizes
sample.sizes <- c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
mn.grid <- expand.grid(m = sample.sizes, n = sample.sizes)

# Run simulation
sim <- plyr::alply(as.matrix(mn.grid), .margins = 1, .fun = function(x) {
  sim_cis(nsim = 100, m = x[1L], n = x[2L])
}, .progress = "text")
save(sim, file = "sim.RData")

# Extract coverage probabilities
cov <- do.call(rbind, lapply(sim, FUN = function(x) attr(x, "coverage"))) %>%
  cbind(mn.grid, .) %>%
  as.data.frame() %>%
  gather(method, coverage, -m, -n)

# Plot results
ggplot(cov, aes(m, coverage, color = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ as.factor(n)) +
  theme_light()


len <- do.call(rbind, lapply(sim, FUN = function(x) attr(x, "length")))

cov <- do.call(rbind, lapply(sim, FUN = function(x) attr(x, "coverage"))) %>%
  as.data.frame() %>%
  gather(X, Y)


# m = 50, n = 50
data.50.50 <- make_data_list(nsim = 2, m = 50, n = 50)
lmm.50.50 <- lapply(data.50.50, FUN = fit_lmm)
wald.50.50 <- do.call(rbind, lapply(lmm.50.50, FUN = compute_ci))

# Function to summarize intervals
summarize.intervals <- function(intervals) {
  length.W <- apply(intervals[, 1:2], 1, diff)
  length.I <- apply(intervals[, 3:4], 1, diff)
  cover.W <- ifelse(intervals[, 1] < 8.015521 & intervals[, 2] > 8.015521, 1, 0)
  cover.I <- ifelse(intervals[, 3] < 8.015521 & intervals[, 4] > 8.015521, 1, 0)
  res = cbind(length.W = length.W, coverage.W = cover.W, 
              length.I = length.I, coverage.I = cover.I)
  list(results = res, summary = rbind(mean = apply(res, 2, mean),
                                        sd = apply(res, 2, sd)))
}



# Simulation function
simulation <- function(nsim = 1000, object = bladder.nlme, 
                       mean.response = TRUE, seed = 1234, ...) {
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Sample sizes
  m <- n <- sample.sizes

  # Initialize matrices to store results
  tab <- matrix(nrow = length(m), ncol = length(n))
  rownames(tab) <- colnames(tab) <- m
  tab.W.cov <- tab.W.len <- tab.W.std <- tab
  tab.I.cov <- tab.I.len <- tab.I.std <- tab  
  
  # Collect results
  pb <- txtProgressBar(min = 0, max = length(m)*length(n), style = 3)
  pb.counter <- 1
  for (i in seq_along(m)) {
    for (j in seq_along(n)) {
      print(paste("Simulation for sample sizes m =", m[i], "and n =", n[j]))
      z <- make.data.list(nsim, m = m[i], n = n[i])  # list of data frames
      intervals <- ldply(z, get.intervals, object = object, 
                         mean.response = mean.response, ...)  # get intervals
      res <- summarize.intervals(intervals)  # summarize the intervals
      tab.W.cov[i, j] <- res$summary[1, 2]  # Wald interval coverage
      tab.W.len[i, j] <- res$summary[1, 1]  # Wald interval length
      tab.W.std[i, j] <- res$summary[2, 1]  # Wald interval sd(length)
      tab.I.cov[i, j] <- res$summary[1, 4]  # inversion interval coverage
      tab.I.len[i, j] <- res$summary[1, 3]  # inversion interval length
      tab.I.std[i, j] <- res$summary[2, 3]  # inversion interval sd(length)
      setTxtProgressBar(pb, pb.counter)  # update progress bar
      pb.counter <- pb.counter + 1  # update counter
    }
  }
  close(pb)
  
  # Return results in a list
  list("Wald-coverage"      = tab.W.cov, 
       "inversion-coverage" = tab.I.cov, 
       "Wald-length"        = tab.W.len, 
       "inversion-length"   = tab.I.len,
       "Wald-std"           = tab.W.std, 
       "inversion-std"      = tab.I.std)
  
}

# Convert results to a single data frame
list2data <- function(x) {
  m <- n <- as.factor(sample.sizes)
  d <- expand.grid(m = m, n = n)
  d <- rbind(d, d)
  d$cp <- c(as.numeric(x[["Wald-coverage"]]), 
            as.numeric(x[["inversion-coverage"]]))
  d$length <- c(as.numeric(x[["Wald-length"]]), 
                as.numeric(x[["inversion-length"]]))
  d$std <- c(as.numeric(x[["Wald-std"]]), 
             as.numeric(x[["inversion-std"]]))
  d$method <- rep(c("Wald", "inversion"), each = 100)
  d
}

## Run simulations -------------------------------------------------------------
# library(doMC)
# registerDoMC(cores = 8)
# getDoParWorkers()  # check
# sim.nlme.reg <- simulation(nsim = 100, .progress = "text")
sim <- simulation(nsim = 100, .progress = "text")
save(sim, file = "simulation.RData")

## Plot coverage probabilities for nlme (regulation)
sim.cp <- list2data(sim)
levels(sim.cp$n) <- paste("n =", sample.sizes)
xyplot(cp ~ m|n, groups = method, data = sim.cp, type = "b", pch = 19, 
       xlab = "Number of subjects (m)", ylab = "Coverage probability",
       auto.key = list(columns = 2), 
       panel = function(x, y, ...) {
         panel.grid()
         panel.xyplot(x, y, ...)
         panel.abline(h = 0.95)
       })




###########################

sim_wald <- function(nsim = 100, m = 5, n = 5) {
  NULL  
}