# install.packages("VGAM")
# install.packages("fMultivar")
library("VGAM")
library("fMultivar")

k = 10000

# computing spearman's correlation coefficient over an array of random number pairs
scorr = function(arr) {
  cor(arr[,1], arr[,2], method = 'spearman')
}

# obtaining cut-off points for rho 
# by simulating null distribution of rho for small n, 
# and using asymptotic normality for large n

cut_off = function(n, prob, alt) {

  if(n <=10) {
    rho = replicate(k, scorr(rbinorm(n)))
    cp = quantile(rho, prob)
  }
  
  else 
    cp = qnorm(prob, sd = 1/sqrt(n-1))
  
  names(cp) = NULL
  return(cp)
}

power = function(n, alpha, alt, func, delta) {
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
    power = sum((rho > cp) | (rho < cp))/k
  }
  
  return(power)
}

# Level of significance fixed as 0.05
alpha = 0.05

# Power as a function of n

nvals = c(2:10,seq(12,20,2),seq(30,100,10),seq(200,1000,100))
powers = lapply(nvals, function(n) power(n, alpha, "upper", rbifgmcop, 0.5))#,mc.cores=4)
plot(nvals,powers,log="x")


# Power as a function of association parameter (Fix n)
n = 7

delta_upper = c(seq(0.01,0.5,0.01))
delta_lower = c(seq(-0.5,-0.01,0.01))
delta_two.tailed = c(delta_lower, delta_upper)

powers = sapply(delta_upper, function(delta) power(n, alpha, "upper", rnorm2d, delta))
plot(delta_upper, powers, type = 'l', col = 'darkred', xlab = "Association Parameter", ylab = "Power", main = "Bivariate Normal", sub = "n = 7")
 
powers = sapply(delta_lower, function(delta) power(n, alpha, "lower", rnorm2d, delta))
plot(delta_lower,powers)
