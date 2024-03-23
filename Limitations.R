#Limitation 
library('copula')
set.seed(22)

##Discrete

scorr <- function(arr){
  cor(arr[,1],arr[,2], method = "spearman")
}

gen.sample <- function(n, param.cop, margins, param.margins){
  cop <- fgmCopula(param.cop, dim = 2)
  mv.BivariateD <- mvdc(cop, margins, param.margins)
  rMvdc(n, mv.BivariateD)
}

nullD <- function(k, n, margins, param.margins){
  samples <- replicate(k, gen.sample(n,0, margins, param.margins)) 
  apply(samples, 3, scorr)
}

#null distribution
hist(nullD(10000, 7, c("pois","pois"),list(list(2),list(4))), freq = FALSE)
hist(nullD(10000,7, c("geom","geom"),list(list(0.2),list(0.5))), freq = FALSE) #getting warnings for small values of n

#power analysis



##Non independent Samples
spearman.rho <- c()
for (k in 1:10000){
  x <- rnorm(10)
  y <- c(runif(1))
  for(i in 2:10){
    y[i] <- y[1]*i
  }
  spearman.rho[k] <- cor(x,y+x,method = "spearman")
}





##Zero Correlation is not equivalent to independence

nvals <- seq(5,100, by=1)
cor_vals <- c()
i <- 1
for (n in nvals){
  X <- rnorm(n)
  Y <- X^2
  cor_vals[i] <- cor(X,Y, method = "spearman")
  i <- i+1
}
hist(cor_vals, freq = FALSE)
plot(nvals, cor_vals, main = "Spearman's Correlation between X and Y")
abline(h = 0, col = "red", lty = 2)
abline(h = c(0.5,-0.5), col = "blue", lty = 2)
mean(cor_vals)

#as the values of n increase, the mean value of cor gets closer and closer to 0
