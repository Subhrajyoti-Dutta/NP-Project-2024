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
width_img = 480
height_img = 480
type = 'p'
color = "red"
cex_sub = 1.5
cex_main = 2
main_line = 2.2
sub_line = -25.5
lwd = 5

#========================================================================================================================
#                                                          vary n
#========================================================================================================================

#1. Upper
# Run parSapply in parallel
nvals = 2:30
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "upper", rnorm2d, 0.3))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Upper-Tailed.png", width = width_img, height = height_img)
plot(nvals, powers, xlab = "n values", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Upper-Tailed", line=main_line, cex.main=cex_main)
title(sub = paste("BVN (rho = ",0.3,")", sep=""), line = sub_line, cex.sub=cex_sub)
dev.off()
#========================================================================================================================

#2. Lower
# Run parSapply in parallel
nvals = 2:30
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "lower", rt2d, -0.3))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Lower-Tailed.png", width = width_img, height = height_img)
plot(nvals, powers, xlab = "n values", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Lower-Tailed", line=main_line, cex.main=cex_main)
title(sub = paste("BVT (rho = ",-0.3,")",sep=""), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. Both
# Run parSapply in parallel
nvals = 2:30
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "both", rbiamhcop, 0.7))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Both-Tailed.png", width = width_img, height = height_img)
plot(nvals, powers, xlab = "n values", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Both-Tailed", line=main_line, cex.main=cex_main)
title(sub = paste("AMH (asso = ",0.7,")",sep=""), line = sub_line, cex.sub=cex_sub)
dev.off()

# at asso = 0.7, scorr = 0.28 and pcorr = 0.29

#========================================================================================================================
#                                                          rho > 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rnorm2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Normal_Upper.png", width = width_img, height = height_img)
plot(delta_upper, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate Normal", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#2a. BVG

# Run parSapply in parallel
powers <- parSapply(cl, 0:50, function(delta) power(n, alpha, "upper", genGammaN, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Gamma1_Upper.png", width = width_img, height = height_img)
plot(0:50, powers, xlab = "Association Parameter N", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate Gamma", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#2b. BVG

# Run parSapply in parallel
powers <- parSapply(cl, 0:50, function(delta) power(n, alpha, "upper", genGammaM, delta))

# Stop the parallel cluster
# stopCluster(cl)
png(file = ".\\power\\Gamma2_Upper.png", width = width_img, height = height_img)
plot(0:50, powers, xlab = "Association Parameter M", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate Gamma", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. BVT

# Run parSapply in parallel
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rt2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\T_Upper.png", width = width_img, height = height_img)
plot(delta_upper, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate T", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================
#                                                          rho < 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rnorm2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Normal_Lower.png", width = width_img, height = height_img)
plot(delta_lower, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate Normal", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#2. BVT

# Run parSapply in parallel
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rt2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\T_Lower.png", width = width_img, height = height_img)
plot(delta_lower, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate T", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. CBVE

# Run parSapply in parallel
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", genCusExpo, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\BVLD_Lower.png", width = width_img, height = height_img)
plot(delta_lower, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Prop. Bivar Life Dist", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================
#                                                          rho != 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rnorm2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\Normal_Both.png", width = width_img, height = height_img)
plot(delta_both, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate Normal", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#2. BVT

# Run parSapply in parallel
powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rt2d, delta))

# Stop the parallel cluster
# stopCluster(cl)

png(file = ".\\power\\T_Both.png", width = width_img, height = height_img)
plot(delta_both, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "Bivariate Student T", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. AMH

alpha_both = seq(-0.7,0.7,0.02)
# Run parSapply in parallel
powers <- parSapply(cl, alpha_both, function(delta) power(n, alpha, "both", rbiamhcop, delta))

# Stop the parallel cluster
stopCluster(cl)
#spearman's rho for the above range is (-0.2,0.3)

png(file = ".\\power\\AMH_Both.png", width = width_img, height = height_img)
plot(alpha_both, powers, xlab = "Association Parameter", ylab = "Power", lwd=lwd, type = type, col = color, cex.lab=cex_lab) 
title(main = "AMH (Uniform Marginals)", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()