# install.packages("VGAM")
library("VGAM")

# computing spearman's correlation coefficient over an array of random number pairs
scorr = function(arr) {
  cor(arr[,1], arr[,2], method = 'spearman')
  }

# obtaining cut-off points for rho 
# by simulating null distribution of rho for small n, 
# and using asymptotic normality for large n

cut_off = function(n, prob, alt) {
  # 
  # if (alt =='upper')
  #   prob = 1 - alpha
  # 
  # else if (alt =='lower')
  #   prob = alpha
  # 
  # else
  #   prob = 1 - (alpha/2)    
    # assuming null distribution of rho is symmetric about it's mean 0 (check)
    # make this remark when presenting the histograms showing variation with n
 
  k = 10000
  
  if(n <=10) {
    rho = replicate(k, scorr(rbinorm(n)))
    # rho = apply(data, 3, scorr)
    # print(mean(rho))
    cp = quantile(rho, prob)
    # hist(rho)
    # print((ecdf(rho)(0.8)))
    # plot(ecdf(rho))
  }
  
  else 
    cp = qnorm(prob, sd = 1/sqrt(n-1))
  
  names(cp) = NULL
  # print(cp)
  return(cp)
}

power = function(n, alpha, alt, func, delta) {
  
  
  
  k = 10000
  print("hell")
  rho = replicate(k, scorr(func(n, delta)))
  # rho = apply(data, 3, scorr)
  # print(c(mean(rho),var(rho)))
  # print()
  print("hello")
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
  print("helllloooo")
  
  return(power)
}


# Set-Ups
# rbifgmcop() - uniform; -1/3 to 1/3 srho
# rbiamhcop() - uniform; -0.271 to 0.478 rho
# rbilogis() 
# rbinorm() - normal; -1 to 1 rho
# rbiclaytoncop() - uniform; 
# x_bar & S
nvals = c(2:10,seq(12,20,2),seq(30,100,10))


powers = sapply(nvals,function(n) power(n, 0.05, "upper", rbifgmcop, 0.5))
plot(nvals,powers)