# Comparison with U statistics
library("combinat")
library("parallel")

compute_kernel_sym <- function(arr) {
  count <- 0
  index <- permn(1:3)
  # print(arr)
  for (i in 1:6) {
    indi <- index[[i]][1]
    indj <- index[[i]][2]
    indk <- index[[i]][3]
    count <- count + (arr[indi,, 1] < arr[indj,, 1] & arr[indi,, 2] < arr[indk,,2])
  }
  return(count / 6 - 1/4)
}

compute_kernel_sym_2 <- function(arr) {
  count <- 0
  index <- permn(1:3)
  # print(arr)
  for (i in 1:6) {
    indi <- index[[i]][1]
    indj <- index[[i]][2]
    indk <- index[[i]][3]
    count <- count + (arr[,indi,, 1] < arr[,indj,, 1] & arr[,indi,, 2] < arr[,indk,,2])
  }
  return(count / 6 - 1/4)
}

compute_kernel_sym_0 <- function(arr) {
  count <- 0
  index <- permn(1:3)
  # print(arr)
  for (i in 1:6) {
    indi <- index[[i]][1]
    indj <- index[[i]][2]
    indk <- index[[i]][3]
    count <- count + (arr[indi, 1] < arr[indj, 1] & arr[indi, 2] < arr[indk,2])
  }
  return(count / 6 - 1/4)
}

compute_u_stat <- function(arr) {
  n <- nrow(arr)
  indices <- combn(1:n, 3)
  # print(indices)
  vals = arr[indices,]
  # print(length(indices))
  dim(vals) = c(3,n*(n-1)*(n-2)/6,2)
  mean(compute_kernel_sym(vals))
  # mean(compute_kernel_sym_2(vals))
  # return(vals)
  # mean(apply(indices, 2, function(ind) compute_kernel_sym_t(arr[ind, ])))
}

compute_u_stat_super <- function(arr,k) {
  n <- nrow(arr)/k
  indices <- combn(1:n, 3)
  # print(indices)
  dim(arr) = c(k,n,2)
  vals = arr[,indices,]
  # print(vals)
  # print(length(indices))
  dim(vals) = c(k,3,n*(n-1)*(n-2)/6,2)
  # dim(compute_kernel_sym(vals))
  arr = compute_kernel_sym_2(vals)
  apply(arr,1,mean)
  # return(vals)
  # mean(apply(indices, 2, function(ind) compute_kernel_sym_t(arr[ind, ])))
}

# c = compute_u_stat(t)
# c
# 
# # Power Comparison
U.stat.power.super <- function(n, alpha, alt, func, delta) {

  k <- 1000
  U.stat = compute_u_stat_super(func(n*k, delta),k)

  if (alt == "upper") {
    powers <- mean(sqrt(n) * U.stat > qnorm(1 - alpha, sd = 1/12))
  } else if (alt == "lower") {
    powers <- mean(sqrt(n) * U.stat < qnorm(alpha, sd = 1/12))
  } else {
    cp <- qnorm(1 - alpha / 2, sd = 1/12)
    powers <- mean((sqrt(n) * U.stat > cp) | (sqrt(n) * U.stat < -cp))
  }

  return(powers)
}

U.stat.power <- function(n, alpha, alt, func, delta) {
  num_cores <- detectCores()
  cl <- makeCluster(num_cores)

  clusterEvalQ(cl, {
    library("fMultivar")
    library("VGAM")
    library("combinat")
  })

  k <- 1000
  clusterExport(cl, c("compute_kernel_sym", "compute_u_stat",'n',"genGammaN","genGammaM","genCusExpo"))
  U.stat <- parSapply(cl, 1:k, function(i) compute_u_stat(func(n, delta)))
  stopCluster(cl)

  if (alt == "upper") {
    powers <- mean(sqrt(n) * U.stat > qnorm(1 - alpha, sd = 1/12))
  } else if (alt == "lower") {
    powers <- mean(sqrt(n) * U.stat < qnorm(alpha, sd = 1/12))
  } else {
    cp <- qnorm(1 - alpha / 2, sd = 1/12)
    powers <- mean((sqrt(n) * U.stat > cp) | (sqrt(n) * U.stat < -cp))
  }

  return(powers)
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

scorr <- function(arr) {
  cor(arr[,1], arr[,2], method = 'spearman')
}

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

power <- function(n, alpha, alt, func, delta) {
  k = 10000
  # rho = replicate(k, scorr())
  arr = func(n*k, delta)
  dim(arr) = c(k,n,2)
  rho = apply(arr,1,scorr)
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

# Define power credentials
n <- 15
alpha <- 0.05

delta_upper = seq(0.025,0.5,0.025)
delta_lower = seq(-0.5,-0.025,0.025)
delta_both = seq(-0.5,0.5,0.05)
nvals = 10:30

# Export necessary objects and functions to the workers
num_cores <- detectCores()
cl <- makeCluster(num_cores)
cldata = clusterEvalQ(cl, {
  library(VGAM)
  library(fMultivar)
})

clusterExport(cl, c("n", "alpha","scorr","cut_off","genGammaN","genGammaM","genCusExpo","power"))

#========================================================================================================================

width_img = 480
height_img = 480
type = 'b'
color = "blue"
color2 = "red"
cex_sub = 1.5
cex_lab = 1.5
cex_main = 2
main_line = 2.2
sub_line = -25.5
lwd = 1

#========================================================================================================================
#                                                          vary n
#========================================================================================================================

#1. Upper
powersU = sapply(nvals, function(n) {U.stat.power.super(n, alpha, "upper", rnorm2d, 0.3)})
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "upper", rnorm2d, 0.3))
maxm = max(powersU,powers)
minm = min(powersU,powers)

#
png(file = ".\\Ustatpower\\Upper-Tailed.png", width = width_img, height = height_img)
plot(nvals, powersU, xlab = "n values", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(nvals, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Upper-Tailed", line=main_line, cex.main=cex_main)
title(sub = paste("BVN (rho = ",0.3,")", sep=""), line = sub_line, cex.sub=cex_sub)
dev.off()
#========================================================================================================================

#2. Lower
powersU = sapply(nvals, function(n) {U.stat.power.super(n, alpha, "lower", rt2d, -0.3)})
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "lower", rt2d, -0.3))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Lower-Tailed.png", width = width_img, height = height_img)
plot(nvals, powersU, xlab = "n values", ylab = "Power", lwd=lwd, ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(nvals, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Lower-Tailed", line=main_line, cex.main=cex_main)
title(sub = paste("BVT (rho = ",-0.3,")",sep=""), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. Both
powersU = sapply(nvals, function(n) {U.stat.power.super(n, alpha, "both", rbiamhcop, 0.7)})
powers <- parSapply(cl, nvals, function(n) power(n, alpha, "both", rbiamhcop, 0.7))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Both-Tailed.png", width = width_img, height = height_img)
plot(nvals, powersU, xlab = "n values", ylab = "Power", lwd=lwd, ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(nvals, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Two-Tailed", line=main_line, cex.main=cex_main)
title(sub = paste("AMH (asso = ",0.7,")",sep=""), line = sub_line, cex.sub=cex_sub)
dev.off()

# at asso = 0.7, scorr = 0.28 and pcorr = 0.29

#========================================================================================================================
#                                                          rho > 0
#========================================================================================================================

#1. BVN

# Run parSapply in parallel
powersU = sapply(delta_upper, function(delta) {U.stat.power.super(n, alpha, "upper", rnorm2d, delta)})
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rnorm2d, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Normal_Upper.png", width = width_img, height = height_img)
plot(delta_upper, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(delta_upper, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate Normal", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#2a. BVG

# Run parSapply in parallel
powersU = sapply(0:50, function(delta) {U.stat.power.super(n, alpha, "upper", genGammaN, delta)})
powers <- parSapply(cl, 0:50, function(delta) power(n, alpha, "upper", genGammaN, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Gamma1_Upper.png", width = width_img, height = height_img)
plot(0:50, powersU, xlab = "Association Parameter N", ylab = "Power", lwd=lwd, ylim=c(minm,maxm),type = type, col = color, cex.lab=cex_lab)
lines(0:50, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate Gamma", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#2b. BVG

# Run parSapply in parallel
powersU <- sapply(seq(0.2,10,0.2), function(delta) {U.stat.power.super(n, alpha, "upper", genGammaM, delta)})
powers <- parSapply(cl, seq(0.2,10,0.2), function(delta) power(n, alpha, "upper", genGammaM, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Gamma2_Upper.png", width = width_img, height = height_img)
plot(seq(0.2,10,0.2), powersU, xlab = "Association Parameter M", ylab = "Power", ylim=c(minm,maxm), lwd=lwd, type = type, col = color, cex.lab=cex_lab)
lines(seq(0.2,10,0.2), powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate Gamma", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. BVT

# Run parSapply in parallel
powersU <- sapply(delta_upper, function(delta) {U.stat.power.super(n, alpha, "upper", rt2d, delta)})
powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rt2d, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

# powers <- parSapply(cl, delta_upper, function(delta) power(n, alpha, "upper", rt2d, delta))

png(file = ".\\Ustatpower\\T_Upper.png", width = width_img, height = height_img)
plot(delta_upper, powersU, xlab = "Association Parameter", ylab = "Power",ylim=c(minm,maxm), lwd=lwd, type = type, col = color, cex.lab=cex_lab)
lines(delta_upper, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate T", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================
#                                                          rho < 0
#========================================================================================================================

#1. BVN

powersU <- sapply(delta_lower, function(delta) {U.stat.power.super(n, alpha, "lower", rnorm2d, delta)})
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rnorm2d, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Normal_Lower.png", width = width_img, height = height_img)
plot(delta_lower, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(delta_lower, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate Normal", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#2. BVT

powersU <- sapply(delta_lower, function(delta) {U.stat.power.super(n, alpha, "lower", rt2d, delta)})
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", rt2d, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\T_Lower.png", width = width_img, height = height_img)
plot(delta_lower, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(delta_lower, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate T", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. CBVE

powersU <- sapply(delta_lower, function(delta) {U.stat.power.super(n, alpha, "lower", genCusExpo, delta)})
powers <- parSapply(cl, delta_lower, function(delta) power(n, alpha, "lower", genCusExpo, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\BVLD_Lower.png", width = width_img, height = height_img)
plot(delta_lower, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(delta_lower, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Prop. Bivar Life Dist", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================
#                                                          rho != 0
#========================================================================================================================

#1. BVN

powersU <- sapply(delta_both, function(delta) {U.stat.power.super(n, alpha, "both", rnorm2d, delta)})
powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rnorm2d, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\Normal_Both.png", width = width_img, height = height_img)
plot(delta_both, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(delta_both, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate Normal", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#2. BVT

powersU <- sapply(delta_both, function(delta) {U.stat.power.super(n, alpha, "both", rt2d, delta)})
powers <- parSapply(cl, delta_both, function(delta) power(n, alpha, "both", rt2d, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\T_Both.png", width = width_img, height = height_img)
plot(delta_both, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd,ylim=c(minm,maxm), type = type, col = color, cex.lab=cex_lab)
lines(delta_both, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "Bivariate Student T", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()

#========================================================================================================================

#3. AMH

alpha_both = seq(-0.7,0.7,0.02)
powersU <- sapply(alpha_both, function(delta) {U.stat.power.super(n, alpha, "both", rbiamhcop, delta)})
powers <- parSapply(cl, alpha_both, function(delta) power(n, alpha, "both", rbiamhcop, delta))
maxm = max(powersU,powers)
minm = min(powersU,powers)

png(file = ".\\Ustatpower\\AMH_Both.png", width = width_img, height = height_img)
plot(alpha_both, powersU, xlab = "Association Parameter", ylab = "Power", lwd=lwd, ylim=c(minm,maxm),type = type, col = color, cex.lab=cex_lab)
lines(alpha_both, powers, lwd=lwd, type = type, col = color2, cex.lab=cex_lab)
title(main = "AMH (Uniform Marginals)", line=main_line, cex.main=cex_main)
title(sub = paste("n =",n), line = sub_line, cex.sub=cex_sub)
dev.off()
#
# powers <- sapply(nvals, function(n) {
#   U.stat.power(n, 0.05, "upper", rnorm2d, 0.5)
# })
#
# print(powers)
#
# plot(nvals,powers)
#Right Tailed
#bivariate normal
#mckay bivariate gamma
#bivariate t

#Left Tailed Alternative
#bivariate normal
#custom exponential
#bivariate t


#Two-Sided Alternative
#bivariate normal
#bivariate t
#bivariate logistic
