# Install and load visualization packages
# install.packages("ggplot2")
# install.packages("tidyr")
library(ggplot2)
library(tidyr)
library(coda)

oneStep <- function(denFun, initSt, denParams, sigProp) {
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

runMetropolis <- function(initSt, denParams, sigProp, n_steps) {
  init <- as.vector(initSt)
  trajectory <- data.frame(x1 = numeric(n_steps), x2 = numeric(n_steps)) 
  rejections <- data.frame(x1 = numeric(0), x2 = numeric(0))            
  
  trajectory[1, ] <- init
  
  for (i in 2:n_steps) { 
    re <- oneStep(den, init, denParams, sigProp)
    
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



# Test different sigProp values
sigProp_vals <- c(0.01, 0.05, 0.1, 0.2, 0.5)
results <- data.frame(
  sigProp = sigProp_vals,
  accept_rate = NA,
  mean_bias_x1 = NA,
  mean_bias_x2 = NA,
  ess_x1 = NA,
  ess_x2 = NA
)

# Iterate over sigProp values
for (i in 1:length(sigProp_vals)) {
  sig <- sigProp_vals[i]
  cat(paste0("Tuning sigProp = ", sig, "\n"))
  
  # Run MCMC using runMetropolis function (10000 steps total)
  init <- c(10, -10)
  n_steps <- 10000
  result <- runMetropolis(
    initSt = init,
    denParams = myDenParams,
    sigProp = sig,
    n_steps = n_steps
  )
  trajectory <- result$trajectory
  
  accept_count <- sum(diff(trajectory$x1) != 0 | diff(trajectory$x2) != 0)
  results$accept_rate[i] <- accept_count / (n_steps - 1)  # Total possible steps = n_steps - 1
  
  # Calculate mean bias (vs true mean (1,1))
  results$mean_bias_x1[i] <- abs(mean(trajectory$x1) - 1)
  results$mean_bias_x2[i] <- abs(mean(trajectory$x2) - 1)
  
  # Calculate effective sample size (ESS)
  mcmc_chain <- as.mcmc(trajectory)
  ess <- effectiveSize(mcmc_chain)
  results$ess_x1[i] <- ess[1]
  results$ess_x2[i] <- ess[2]
}


# Reshape data for visualization
results_long <- pivot_longer(results, cols = -sigProp, names_to = "metric", values_to = "value")

# Plot hyperparameter tuning results
ggplot(results_long, aes(x = sigProp, y = value, color = metric)) +
  geom_line(linewidth = 1) +
  facet_wrap(~metric, scales = "free_y") +
  labs(title = "MCMC Performance vs Proposal Standard Deviation (sigProp)",
       x = "Proposal Standard Deviation (sigProp)",
       y = "Metric Value",
       color = "Evaluation Metric") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))