# Load necessary libraries
library(parallel)
library(VGAM)
library(fMultivar)

# Define the number of cores to use
num_cores <- detectCores()
k = 10000
# Initialize a parallel cluster
cl <- makeCluster(num_cores)

# Load necessary libraries and export functions/objects to workers
clusterEvalQ(cl, {
  library(VGAM)
  library(fMultivar)
  
  # Define the scorr function
  scorr <- function(arr) {
    cor(arr[,1], arr[,2], method = 'spearman')
  }
  
  # Define the cut_off function
  cut_off <- function(n, prob, alt) {
    k <- 10000
    if(n <= 10) {
      rho <- replicate(k, scorr(rbinorm(n)))
      cp <- quantile(rho, prob)
    } else {
      cp <- qnorm(prob, sd = 1/sqrt(n-1))
    }
    names(cp) <- NULL
    return(cp)
  }
})

# Define power function
# k =10000
n <- 7
alpha <- 0.05
delta_upper <- c(seq(0.01,0.5,0.01))

power <- function(n, alpha, alt, func, delta) {
  k = 10000
  rho = replicate(k, scorr(func(n, delta)))
  if(alt == "upper") {
    cp = cut_off(n, 1-alpha, alt)
    power = mean(rho > cp)
  }
  
  else if(alt == "lower") {
    cp = cut_off(n, alpha, alt)
    power = mean(rho < cp)
  }
  
  else {
    cp = cut_off(n, 1-alpha/2, alt)
    power = mean((rho > cp) | (rho < cp))
  }
  
  return(power)
}

# Export necessary objects and functions to the workers
clusterExport(cl, c("n", "alpha", "rbinorm", "power"))


# Run parSapply in parallel
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rnorm2d, delta))

# Stop the parallel cluster
stopCluster(cl)

# Print the result
print(powers)

plot(delta_upper, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")
