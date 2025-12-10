library(coda)

# 1. Bivariate Normal Density Function
myDenParams <- c(1, 1, 5, 5, 0.8)
den <- function(argmts, params) {
  x <- as.numeric(argmts)
  p <- as.numeric(params)
  
  if (length(x) < 2 || length(p) < 5) return(0)
  
  x1 <- x[1]
  x2 <- x[2]
  u1 <- p[1]
  u2 <- p[2]
  sigma1 <- p[3]
  sigma2 <- p[4]
  ro <- p[5]
  
  if (!is.numeric(x1) || !is.numeric(u1) || !is.numeric(sigma1) || 
      !is.numeric(x2) || !is.numeric(u2) || !is.numeric(sigma2) || 
      !is.numeric(ro)) return(0)
  
  z <- (x1 - u1)^2/(sigma1)^2 - 2*ro*(x1 - u1)*(x2 - u2)/(sigma1*sigma2) + (x2 - u2)^2/(sigma2)^2
  f <- 1/(2*pi*sigma1*sigma2*sqrt(1 - ro^2)) * exp(-z/(2*(1 - ro^2)))
  return(f)
}

# 2. Single Step of Metropolis-Hastings
oneStep <- function(denFun, initSt, denParams, sigProp) {
  initSt <- as.numeric(initSt) 
  
  prob.c <- denFun(initSt, denParams)
  theta.1 <- rnorm(1, initSt[1], sigProp)
  theta.2 <- rnorm(1, initSt[2], sigProp)
  theta <- c(theta.1, theta.2)
  prob.theta <- denFun(theta, denParams)
  acceptance <- prob.theta / prob.c
  
  if (acceptance > 1 | acceptance > runif(1)) {
    return(list(
      state = theta,   
      jump = 1,     
      rej = NA        
    ))
  } else {
    return(list(
      state = initSt,  
      jump = 0,     
      rej = theta 
    ))
  }
}

# 3. Metropolis Sampler Function (KEEP AS IS, using denFun from previous fix)
runMetropolis <- function(denFun, initSt, denParams, sigProp, n_steps) {
  # Ensure initial state is a numeric vector
  init <- as.numeric(initSt) 
  trajectory <- data.frame(x1 = numeric(n_steps), x2 = numeric(n_steps)) 
  rejections <- data.frame(x1 = numeric(0), x2 = numeric(0))            
  
  trajectory[1, ] <- init
  
  for (i in 2:n_steps) { 
    re <- oneStep(denFun, init, denParams, sigProp)
    
    init <- re$state
    
    if (re$jump == 1) {
      trajectory[i, ] <- init
    } else {
      if (!any(is.na(re$rej))) {
        rejections <- rbind(rejections, as.data.frame(t(re$rej)))
      }
      trajectory[i, ] <- init
    }
  }
  
  return(list(
    trajectory = trajectory,
    rejections = rejections
  ))
}

# 4. Running the MCMC Chains (FIXED: Robust state extraction)
n_chains <- 3
sigProp <- 0.05
burnin_steps <- 1000
effective_steps <- 10000
initial_points <- list(
  c(10, -10),
  c(-5, 8),
  c(7, -3)
)

chains_list <- list()

for (chain_id in 1:n_chains) {
  cat(paste0("Running chain ", chain_id, " with initial point: (", 
             initial_points[[chain_id]][1], ", ", initial_points[[chain_id]][2], ")\n"))
  
  burnin_result <- runMetropolis(
    denFun = den,
    initSt = initial_points[[chain_id]],
    denParams = myDenParams,
    sigProp = sigProp,
    n_steps = burnin_steps
  )
  
  init_burned <- unlist(burnin_result$trajectory[nrow(burnin_result$trajectory), ]) 
  
  effective_result <- runMetropolis(
    denFun = den,
    initSt = init_burned,
    denParams = myDenParams,
    sigProp = sigProp,
    n_steps = effective_steps
  )
  trajectory <- effective_result$trajectory
  
  chains_list[[chain_id]] <- as.mcmc(trajectory)
}

# 5. Convergence and Results
gelman_chains <- mcmc.list(chains_list)
gelman_result <- gelman.diag(gelman_chains, multivariate = FALSE)

cat("\n==================== Convergence Test Results ====================\n")
cat(paste0("PSRF for x1: ", round(gelman_result$psrf["x1", "Point est."], 3), "\n"))
cat(paste0("PSRF for x2: ", round(gelman_result$psrf["x2", "Point est."], 3), "\n"))

if (gelman_result$psrf["x1", "Point est."] < 1.1 & gelman_result$psrf["x2", "Point est."] < 1.1) {
  cat("\nConclusion: CONVERGED (PSRF < 1.1)\n")
} else {
  cat("\nConclusion: NOT CONVERGED (PSRF > 1.1)\n")
}

all_samples <- do.call(rbind, gelman_chains)
sample_mean <- colMeans(all_samples)
cat(paste0("\nSample Mean (x1, x2): (", round(sample_mean[1], 3), ", ", round(sample_mean[2], 3), ")\n"))
cat(paste0("True Mean: (1, 1)\n"))





calculate_acceptance_rate <- function(trajectory, initial_steps) {
  if (nrow(trajectory) <= 1) {
    return(0)
  }
  
  # Calculate differences between consecutive rows
  diffs <- trajectory[2:nrow(trajectory), ] - trajectory[1:(nrow(trajectory) - 1), ]
  
  # A 'jump' (acceptance) occurs if the state changes.
  # Sum the rows where at least one coordinate has changed (i.e., not a zero row vector)
  accepted_moves <- sum(apply(diffs, 1, function(row) any(row != 0)))
  
  # The total number of proposals is the total steps minus the initial step.
  total_proposals <- nrow(trajectory) - 1
  
  return(accepted_moves / total_proposals)
}

# 2. Adaptive Metropolis Sampler
runMetropolisAdaptive <- function(denFun, initSt, denParams, sigProp, n_steps, adapt_freq = 100, target_rate = 0.234) {
  init <- as.numeric(initSt)
  trajectory <- data.frame(x1 = numeric(n_steps), x2 = numeric(n_steps)) 
  rejections <- data.frame(x1 = numeric(0), x2 = numeric(0))
  current_sigProp <- sigProp
  
  trajectory[1, ] <- init
  
  for (i in 2:n_steps) { 
    # Use the current sigma proposal
    re <- oneStep(denFun, init, denParams, current_sigProp)
    
    init <- re$state
    
    if (re$jump == 1) {
      trajectory[i, ] <- init
    } else {
      if (!any(is.na(re$rej))) {
        rejections <- rbind(rejections, as.data.frame(t(re$rej)))
      }
      trajectory[i, ] <- init
    }
    
    # Adaptation logic: Adjust sigProp periodically
    if (i %% adapt_freq == 0 && i < n_steps) {
      # Calculate acceptance rate over the steps completed so far
      current_rate <- calculate_acceptance_rate(trajectory[1:i, ], i)
      
      # Simple multiplicative adjustment rule
      # If rate is too low, decrease step size (decrease sigProp)
      if (current_rate < target_rate) {
        current_sigProp <- current_sigProp * 0.99
      } 
      # If rate is too high, increase step size (increase sigProp)
      else {
        current_sigProp <- current_sigProp * 1.01
      }
    }
  }
  
  return(list(
    trajectory = trajectory,
    rejections = rejections,
    final_sigProp = current_sigProp
  ))
}

# 3. Run the Adaptive Sampler
init_adaptive <- c(10, -10)
burnin_adaptive <- 1000
n_steps_adaptive <- 10000
# Starting sigProp (will be tuned)
sigProp_start <- 0.05 
# Adaptation will happen every 50 steps
adapt_frequency <- 50 

cat("\n\n--- Running Adaptive Metropolis Sampler ---\n")

adaptive_result <- runMetropolisAdaptive(
  denFun = den, # Assuming 'den' is defined as in the previous solution
  initSt = init_adaptive,
  denParams = myDenParams, # Assuming 'myDenParams' is defined
  sigProp = sigProp_start,
  n_steps = n_steps_adaptive,
  adapt_freq = adapt_frequency
)

trajectory_adaptive <- adaptive_result$trajectory
rejections_adaptive <- adaptive_result$rejections

cat(paste0("Final Proposal Sigma (sigProp): ", round(adaptive_result$final_sigProp, 4), "\n"))

# Final Acceptance Rate
final_acc_rate <- calculate_acceptance_rate(trajectory_adaptive, n_steps_adaptive)
cat(paste0("Final Acceptance Rate: ", round(final_acc_rate, 4), "\n"))


# 4. Plot Adaptive Trajectory and Rejections
plot(trajectory_adaptive, 
     col = "blue", 
     ylim = c(-10, 10), 
     xlim = c(-5, 10), 
     xlab = "x1", 
     ylab = "x2", 
     main = "Adaptive MCMC Trajectory") # Main title revised
points(rejections_adaptive, col = "orange", pch = 16)

legend("bottomleft", 
       legend = c("trajectory", "rejections"), 
       col = c("blue", "orange"), 
       pch = c(1, 16))
