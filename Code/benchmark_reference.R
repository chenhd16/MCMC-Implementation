# Same initial point: (10, -10)
# Same # of MCMC steps: 10000
# Same sigProp=0.05
# Samw burnin: 1000

# Install and load required packages
# install.packages("MCMCpack")
# install.packages("mvtnorm")
library(MCMCpack)
library(mvtnorm)

# Step 1: Define target distribution parameters
cov_mat <- matrix(c(5^2, 0.8*5*5, 0.8*5*5, 5^2), nrow = 2)  # Covariance matrix
true_mean <- c(1, 1)  # True mean of target distribution (for accuracy validation)


# Step 2: Run MCMC using MCMCpack (benchmark)
set.seed(42)
benchmark_mcmc <- MCMCmetrop1R(
  fun = function(theta) dmvnorm(theta, mean = true_mean, sigma = cov_mat, log = TRUE),
  theta.init = c(10, -10),
  burnin = 1000,
  mcmc = 10000,
  tune = 0.05,
  verbose = FALSE  # Disable redundant output
)

# Step 3: Re-run my Metropolis implementation (with burn-in)
myDenParams <- c(1, 1, 5, 5, 0.8)
den <- function(argmts, params) {
  x1 <- argmts[1]
  x2 <- argmts[2]
  u1 <- params[1]
  u2 <- params[2]
  sigma1 <- params[3]
  sigma2 <- params[4]
  ro <- params[5]
  
  z <- (x1 - u1)^2/(sigma1)^2 - 2*ro*(x1 - u1)*(x2 - u2)/(sigma1*sigma2) + (x2 - u2)^2/(sigma2)^2
  f <- 1/(2*pi*sigma1*sigma2*sqrt(1 - ro^2)) * exp(-z/(2*(1 - ro^2)))  # Fixed exp precision issue
  return(f)
}

oneStep <- function(denFun, initSt, denParams, sigProp) {
  prob.c <- denFun(initSt, denParams)  # Probability of current state
  theta.1 <- rnorm(1, initSt[1], sigProp)
  theta.2 <- rnorm(1, initSt[2], sigProp)
  theta <- c(theta.1, theta.2)
  prob.theta <- denFun(theta, denParams)
  
  acceptance <- prob.theta/prob.c
  if (acceptance > 1 | acceptance > runif(1)) {
    initSt <- theta
    jump <- 1 
    rej <- NA
  } else {
    jump <- 0
    rej <- theta
  }
  re <- c(initSt, jump, rej)
  return(re)
}

# Run my MCMC with burn-in
init <- c(10, -10)

for (i in 1:1000) {
  re <- oneStep(den, init, myDenParams, 0.05)
  init <- c(re[1], re[2])
}
# effective sampling: 10000
trajectory <- data.frame(x1 = init[1], x2 = init[2])
for (i in 1:10000) {
  re <- oneStep(den, init, myDenParams, 0.05)
  init <- c(re[1], re[2])
  if (re[3] > 0) {
    trajectory[nrow(trajectory)+1,] <- init
  }
}

# Step 4: Calculate core comparison metrics
# 4.1 Acceptance rate (only for effective sampling phase)
my_accept_rate <- nrow(trajectory)/10000
benchmark_trajectory <- as.data.frame(benchmark_mcmc)  
benchmark_accept_rate <- nrow(unique(benchmark_trajectory)) / nrow(benchmark_trajectory)


# 4.2 Sampling mean (compare with true mean)
my_mean <- colMeans(trajectory)
benchmark_sample_mean <- colMeans(benchmark_trajectory)

# 4.3 Runtime comparison (include burn-in + effective sampling)
my_time <- system.time({
  # Re-run my MCMC for time measurement (full pipeline)
  init <- c(10,-10)
  # Burn-in
  for (i in 1:1000) {
    re <- oneStep(den, init, myDenParams, 0.05)
    init <- c(re[1], re[2])
  }
  # Effective sampling
  trajectory <- data.frame(x1 = init[1], x2 = init[2])
  i <- 1
  while (i <= 10000) {
    re <- oneStep(den, init, myDenParams, 0.05)
    init <- c(re[1], re[2])
    if (re[3] > 0) trajectory[nrow(trajectory)+1,] <- init
    i <- i + 1
  }
})[3]

benchmark_time <- system.time({
  # Re-run MCMCpack for time measurement
  benchmark_mcmc <- MCMCmetrop1R(
    fun = function(theta) dmvnorm(theta, mean = true_mean, sigma = cov_mat, log = TRUE),
    theta.init = c(10, -10),
    burnin = 1000,
    mcmc = 10000,
    tune = 0.05,
    verbose = FALSE
  )
})[3]

# Print comparison results
cat(paste0("my implementation acceptance rate: ", round(my_accept_rate, 3), "\n"))
cat(paste0("MCMCpack acceptance rate: ", round(benchmark_accept_rate, 3), "\n"))
cat(paste0("my implementation sample mean: x1=", round(my_mean[1], 3), ", x2=", round(my_mean[2], 3), "\n"))
cat(paste0("MCMCpack sample mean: x1=", round(benchmark_sample_mean[1], 3), ", x2=", round(benchmark_sample_mean[2], 3), "\n"))
cat(paste0("my implementation runtime: ", round(my_time, 3), "s\n"))
cat(paste0("MCMCpack runtime: ", round(benchmark_time, 3), "s\n"))

