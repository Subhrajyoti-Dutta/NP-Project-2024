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
  
  set.seed(56)
  
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
  
  genGammaN = function(n, N){
    M = 5
    x = rgamma(n, M)
    y = rgamma(n, N)
    cbind(x, x+y)
  }
  
  genGammaM = function(n, M){
    N = 5
    x = rgamma(n, M)
    y = rgamma(n, N)
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
#                                                          vary n
#========================================================================================================================

#1. Upper
# Run parSapply in parallel
nvals = 2:30
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "upper", rnorm2d, 0.3))

# Stop the parallel cluster
# stopCluster(cl)

# plot(nvals[1:, powers[1:19], type = 'p', col = 'darkred', xlab = "n values", ylab = "Power", main = "Power vs n")
png(file = paste(".\\var_with_n\\powervsn.png",sep=""), width = 960, height = 480)
plot(nvals, powers, type = 'p', col = 'darkred', xlab = "n values", ylab = "Power", main = "Upper-Tailed")
title(sub = paste("BVN (rho = ",0.3,")",sep=""), line = -14)

#========================================================================================================================

#2. Lower
# Run parSapply in parallel
nvals = 2:30
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "lower", rt2d, -0.3))

# Stop the parallel cluster
# stopCluster(cl)

# plot(nvals[1:, powers[1:19], type = 'p', col = 'darkred', xlab = "n values", ylab = "Power", main = "Power vs n")
plot(nvals, powers, type = 'p', col = 'darkred', xlab = "n values", ylab = "Power", main = "Lower-Tailed")
title(sub = paste("BVT (rho = ",-0.3,")",sep=""), line = -14)

#========================================================================================================================

#2. Both
# Run parSapply in parallel
nvals = 2:30
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "both", rbiamhcop, 0.7))

# Stop the parallel cluster
stopCluster(cl)

# plot(nvals[1:, powers[1:19], type = 'p', col = 'darkred', xlab = "n values", ylab = "Power", main = "Power vs n")

plot(nvals, powers, type = 'p', col = 'darkred', xlab = "n values", ylab = "Power", main = "Two-Tailed")
title(sub = paste("AMH (asso = ",0.7,")",sep=""), line = -14)
#at asso = 0.7, scorr = 0.28 and pcorr = 0.29

#========================================================================================================================
#                                                          rho > 0
#========================================================================================================================
# 
# #1. BVN
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rnorm2d, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_upper, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")
# title(sub = paste("n =",n), line = -14)

# #========================================================================================================================
# 
# #2a. BVG
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, 0:50, function(delta) power(n, alpha, "upper", genGammaN, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(0:50, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Gamma")
# title(sub = paste("n =",n), line = -14)
# 
# #2b. BVG
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, 0:50, function(delta) power(n, alpha, "upper", genGammaM, delta))
# 
# # Stop the parallel cluster
# stopCluster(cl)
# 
# plot(0:50, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Gamma")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# 
# #3. BVT
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rt2d, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_upper, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Student T")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# #                                                          rho < 0
# #========================================================================================================================
# 
# #1. BVN
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rnorm2d, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# 
# #2. BVT
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rt2d, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Student T")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# 
# #3. CBVE
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", genCusExpo, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Prop. Bivar Life Dist")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# #                                                          rho != 0
# #========================================================================================================================
# 
# #1. BVN
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rnorm2d, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_both, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# 
# #2. BVT
# 
# # Run parSapply in parallel
# powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rt2d, delta))
# 
# # Stop the parallel cluster
# # stopCluster(cl)
# 
# plot(delta_both, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Student T")
# # plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Prop. Bivar Life Dist")
# title(sub = paste("n =",n), line = -14)
# 
# #========================================================================================================================
# 
#3. AMH
# 
# alpha_both = seq(-0.7,0.7,0.02)
# # Run parSapply in parallel
# powers <- parSapply(cl, alpha_both, function(delta) power(n, alpha, "both", rbiamhcop, delta))
# 
# # Stop the parallel cluster
# stopCluster(cl)
# #spearman's rho for the above range is (-0.2,0.3)
# 
# plot(alpha_both, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "AMH (Uniform Marginals)")
# # plot(delta_lower, powers, type = 'p', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Prop. Bivar Life Dist")
# title(sub = paste("n =",n), line = -14)
# 
