#Null Distribution of Spearman's Rho 
library('copula')
library('ggplot2')
set.seed(22)

#gives Spearman's rho for a given 2D array
scorr <- function(arr){
  cor(arr[,1],arr[,2], method = "spearman")
}

#generate random sample of a given size with fixed marginals
gen.sample <- function(n, param.cop, margins, param.margins){
  cop <- fgmCopula(param.cop, dim = 2)
  mv.BivariateD <- mvdc(cop, margins, param.margins)
  rMvdc(n, mv.BivariateD)
}

#computes k values of Spearman's rho to simulate its distribution
nullD <- function(k, n, margins, param.margins){
  samples <- replicate(k, gen.sample(n,0, margins, param.margins)) 
  apply(samples, 3, scorr)
}

#histograms with density

## Bivariate Exponential
rho.spearman.exp <- nullD(5000, 15, c("exp","exp"), list(list(rate = 2),list(rate = 4)))
hist(rho.spearman.exp, freq = FALSE , col = "blue", xlab = "Spearman's rho" , main = "Distribution of Spearman's Rho Under Ho")
lines(density(rho.spearman.exp), col = "darkred")
cut.off <- quantile(rho.spearman.exp, probs = c(0.95)) #0.436
##Bivariate Normal
rho.spearman.norm <- nullD(5000, 15, c("norm","norm"), list(list(0,1),list(2,9)))
hist(rho.spearman.norm, freq = FALSE , col = "blue", xlab = "Spearman's rho" , main = "Distribution of Spearman's Rho Under Ho")
lines(density(rho.spearman.norm), col = "darkred")
cut.off <- quantile(rho.spearman.norm, probs = c(0.95)) #0.432
## Bivariate Gamma
rho.spearman.gamma <- nullD(5000, 15, c("gamma","gamma"), list(list(shape = 4, rate = 2),list(shape = 6, rate = 4)))
hist(rho.spearman.gamma, freq = FALSE , col = "blue", xlab = "Spearman's rho" , main = "Distribution of Spearman's Rho Under Ho")
lines(density(rho.spearman.gamma), col = "darkred")
cut.off <- quantile(rho.spearman.gamma, probs = c(0.95)) #0.436
## Bivariate Cauchy
rho.spearman.cauchy <- nullD(5000, 15, c("cauchy","cauchy"), list(list(location = 0, scale = 1),list(location = 0, scale = 1)))
hist(rho.spearman.cauchy, freq = FALSE , col = "blue", xlab = "Spearman's rho" , main = "Distribution of Spearman's Rho Under Ho")
lines(density(rho.spearman.cauchy), col = "darkred")
cut.off <- quantile(rho.spearman.cauchy, probs = c(0.95)) #0.45
## Bivariate Gumbel
cut.off <- quantile(rho.spearman, probs = c(0.95))
## Bivariate Logistic



#qqplots

##Exponential vs Gamma
qqplot(rho.spearman.exp,rho.spearman.gamma)

##Gamma vs Normal
qqplot(rho.spearman.gamma,rho.spearman.norm)

##Normal vs Cauchy
qqplot(rho.spearman.norm,rho.spearman.cauchy)
##Cauchy vs Gumbel

##Gumbel vs Logistic


#Limiting Distribution under Ho

lim.distbn <- function(n){
  lim.dist.form <- sqrt(n-1)*nullD(10000, n, c("exp","exp"), list(list(rate = 2),list(rate = 4)))
  par(mfrow = c(1,2))
  hist(lim.dist.form, freq = FALSE, main = paste("Histogram of rho when n = ",n))
  lines(density(lim.dist.form), col = "darkred")
  #combine histogram and qqplot in a grid
  qqnorm(lim.dist.form)
}

lim.distbn(25)

