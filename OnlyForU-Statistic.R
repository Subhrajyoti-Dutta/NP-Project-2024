#Comparison with U statistics
library('copula')
library("combinat")
library("parallel")

compute_kernel_sym = function(arr){
  count = 0
  index = permn(1:3)
  for (i in 1:6){
    indi = index[[i]][1]
    indj = index[[i]][2]
    indk = index[[i]][3]
    count = count + (arr[indi,1] < arr[indj,1] && arr[indi,2] < arr[indk,2])
  }
  return (count/6 - 1/4)
}

compute_u_stat = function(arr){
  n = nrow(arr)
  indices = combn(1:n,3)
  mean(apply(indices, 2, function(ind){ compute_kernel_sym(arr[ind,]) }))
}

gen.sample <- function(n, param.cop, margins, param.margins){
     cop <- fgmCopula(param.cop, dim = 2)
     mv.BivariateD <- mvdc(cop, margins, param.margins)
     rMvdc(n, mv.BivariateD)
  }

# 
# #calculate cut-offs
U.stat.null.cutoff <- function(k,n,margins, param.margins){
  cl <- makeCluster(num_cores)
  
  clusterEvalQ(cl, {
    library('copula')
    library("combinat")
  })
  
  clusterExport(cl, c("compute_kernel_sym","compute_u_stat","gen.sample"))
  
  stat = parSapply(cl, 1:k, function(i){compute_u_stat(gen.sample(n,0,margins, param.margins))})
  hist(stat, freq = FALSE)
  lines(density(stat), col = "darkred")
}
# 


#Power Comparison
U.stat.power <- function(n, alpha, alt, func, delta){
  U.stat <- replicate(k, compute_u_stat(func(n,delta)))
  if(alt == "upper") {
    #use asymptotic normality of U statistic
    cp = cut_off(n, 1-alpha, alt)
    power = mean(U.stat > cp)
  }
  
  else if(alt == "lower") {
    cp = cut_off(n, alpha, alt)
    power = mean(U.stat < cp)
  }
  
  else {
    cp = cut_off(n, 1-alpha/2, alt)
    power = mean((U.stat > cp) | (U.stat < -cp))
  }
  
  return(power)

}

# Define the number of cores to use
num_cores <- detectCores()
k = 10000

# #testrun
# U.stat.null.cutoff(10000, 15, c("exp","exp"), list(list(rate = 2),list(rate = 4)))
U.stat.null.cutoff(10000, 15, c("norm","norm"), list(list(0,1),list(2,9)))


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
