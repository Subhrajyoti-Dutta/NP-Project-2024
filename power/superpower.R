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
  cut_off <- function(n, prob) {
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
  
  genGamma = function(n, rho){
    v1 = 5
    v2 = v1 * rho^2 / (1-rho^2)
    x = rgamma(n, v2)
    y = rgamma(n, v1)
    cbind(x, x+y)
  }
  
  genCusExpo = function(n, rho){
    l = 5
    c = (l+1)^2/(l*(l+2))
    a = c * (l+1)^2 * rho^2
    b = c * rho^2
    mu = max(((2*b+2) + sqrt((2*b+2)^2 - 4*(b-1)*a))/(2*(b-1)),
             ((2*b+2) - sqrt((2*b+2)^2 - 4*(b-1)*a))/(2*(b-1)))
    z1 = rexp(n,l)
    z2 = rexp(n,mu)
    cbind(z1,exp(z2-z1))
  }
})

# Define power function
# k =10000
n <- 20
alpha <- 0.05

delta_upper = seq(0,0.5,0.01)
delta_lower = seq(-0.5,0,0.01)
delta_both = seq(-0.5,0.5,0.01)

power <- function(n, alpha, alt, func, delta) {
  k = 10000
  rho = replicate(k, scorr(func(n, delta)))
  if(alt == "upper") {
    cp = cut_off(n, 1-alpha)
    power = mean(rho > cp)
  }
  
  else if(alt == "lower") {
    cp = cut_off(n, alpha)
    power = mean(rho < cp)
  }
  
  else {
    cp = cut_off(n, 1-alpha/2)
    # print(cp)
    power = mean(abs(rho) > cp)
  }
  
  return(power)
}

# Export necessary objects and functions to the workers
clusterExport(cl, c("n", "alpha", "rbinorm", "power"))

#========================================================================================================================
#                                                          rho > 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rnorm2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_upper, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")

#========================================================================================================================

#2. BVG

# Run parSapply in parallel
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", genGamma, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_upper, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Gamma")

#========================================================================================================================

#3. BVT

# Run parSapply in parallel
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rt2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_upper, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Student T")

#========================================================================================================================
#                                                          rho < 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rnorm2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")

#========================================================================================================================

#2. BVT

# Run parSapply in parallel
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rt2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Student T")

#========================================================================================================================

#3. CBVE

# Run parSapply in parallel
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", genCusExpo, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Custom Bivar Exponential")

#========================================================================================================================
#                                                          rho != 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rnorm2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

plot(delta_both, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")

#========================================================================================================================

#2. BVT

# Run parSapply in parallel
powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rt2d, delta))

# Stop the parallel cluster
stopCluster(cl)

plot(delta_both, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Student T")
