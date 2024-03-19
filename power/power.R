# install.packages("VGAM")
library("VGAM")

# computing spearman's correlation coefficient over an array of random number pairs
scorr = function(arr) {
  cor(arr[,1], arr[,2], method = 'spearman')
  }

# obtaining cut-off points for rho 
# by simulating null distribution of rho for small n, 
# and using asymptotic normality for large n

cut_off = function(n, alpha, alt) {
 
  if (alt =='upper')
    prob = 1 - alpha
  
  else if (alt =='lower')
    prob = alpha
  
  else
    prob = 1 - (alpha/2)    
    # assuming null distribution of rho is symmetric about it's mean 0 (check)
    # make this remark when presenting the histograms showing variation with n
 
  k = 100000
  
  if(n <=10) {
    
    data = replicate(k, rnorm2d(n))
    rho = apply(data, 3, scorr)
    cp = quantile(rho, prob)
    
  }
  
  else 
    cp = qnorm(prob, sd = 1/sqrt(n-1))
  
  names(cp) = NULL
  
  return(cp)
}

power = function(n, alpha, alt, func, delta) {
  
  cp = cut_off(n, alpha, alt)
  
  k = 10000
  data = replicate(k, func(n, delta))
  rho = apply(data, 3, scorr)
  
  if(alt == "upper") {
    power = sum(rho >= cp)/k
  }
  
  else if(alt == "lower") {
    power = sum(rho >= cp)/k
  }
  
  else {
    power = sum((rho >= cp) | (rho <= cp))/k
  }
  
  return(power)
}


# Set-Ups
# rbifgmcop() - uniform; -1/3 to 1/3 srho
# rbiamhcop() - uniform; -0.271 to 0.478 rho
# rbilogis() 
# rbinorm() - normal; -1 to 1 rho
# rbiclaytoncop() - uniform; 
# x_bar & S

