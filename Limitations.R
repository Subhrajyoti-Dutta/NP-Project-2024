library("fMultivar")

#set.seed(423)

# Violation of Assumptions 
# Discrete Case

scorr <- function(arr){
  t = cor(arr[,1],arr[,2], method = "spearman")
  return(t)
}

Normal = function(n){
  cbind(rnorm(n),rnorm(n))
}

Binom = function(n){
  cbind(rbinom(n,5,0.5), rbinom(n,5,0.5))
}

Geom = function(n){
  cbind(rgeom(n,0.5),rgeom(n,0.5))
}

Pois = function(n){
  cbind(rpois(n,1),rpois(n,1))
}

k = 10000
dists = c(Normal, Binom, Geom, Pois)
names = c("BVN", "BVB", "BVG", "BVP")

for (n in c(7,20)){
  png(file=paste(".\\viol_lim\\together_discrete",n,".png",sep=""),width=1200,height=800)
  par(mfrow=c(2,2))
  for (i in 1:4){
    rho <- replicate(k, scorr(dists[[i]](n)))
    hist(rho,xlim = c(-1,1),main = names[i],xlab="S. Correlation (rho)",cex.lab=1.5,cex.main=2)
  }
  dev.off()
}

# non_identical sample

non_iden.1 = function(n) {
  x = rep(1,n)
  y = rep(1,n)
  
  for (i in i:n) {
    x[i] = rnorm(1) 
    y[i] = rnorm(1)
  }
  cbind(x,y)
}

non_iden.2 = function(n) {
  x = rep(1,n)
  y = rep(1,n)
  
  for (i in i:n) {
    x[i] = rnorm(1, i, i) 
    y[i] = rnorm(1, i, i)
  }
  cbind(x,y)
}

non_iden.3 = function(n) {
  x = rep(1,n)
  y = rep(1,n)
  
  for (i in i:n) {
    x[i] = rexp(1, 1/i)
    y[i] = rexp(1, 1/i)
  }
  cbind(x,y)
}

non_iden.4 = function(n) {
  x = rep(1,n)
  y = rep(1,n)
  
  for (i in i:n) {
    x[i] = rexp(1, i)
    y[i] = rexp(1, i)
  }
  cbind(x,y)
}


dists = c(non_iden.1, non_iden.2, non_iden.3, non_iden.4)
names = c("BVN1", "BVN2", "BVE1", "BVE2")

for (n in c(7,20)){
  png(file=paste(".\\viol_lim\\together_non_iden",n,".png",sep=""),width=1200,height=800)
  par(mfrow=c(2,2))
  for (i in 1:4){
    rho <- replicate(k, scorr(dists[[i]](n)))
    hist(rho,xlim = c(-1,1),main = names[i],xlab="S. Correlation (rho)",cex.lab=1.5,cex.main=2)
  }

  dev.off()
}

# Non-independent Samples: Power Computation

non_ind = function(n){
  x <- rnorm(n)
  y[1] <- c(runif(1))
  for(i in 2:n){
    y[i] <- i*y[1] + x[i]
  }
  cbind(x,y)
}

rho = replicate(k, scorr(non_ind.1(n)))
power = mean(rho > cut_off(n, 0.95))
power

# Only Monotonic Association can be captured

n <- 7

# eg 1.1
x <- rnorm(n)
y <- x^2

genlimdista1 = function(n){
  x = rnorm(n)
  y = x^2
  cbind(x, y)
}

# eg 1.2
x <- rexp(n)
y <- x^2

genlimdista2 = function(n){
  x = rexp(n) #rexp(n) for showing monotone
  y = x^2
  cbind(x, y)
}

# eg 2.1
x <- runif(n, -1, 1) #runif(n) for showing monotone
y <- cos(x)

genlimdistb1 = function(n){
  x = runif(n,-1,1)
  y = cos(x)
  cbind(x, y)
}

# eg 2.2
x <- runif(n, -1, 1) #runif(n) for showing monotone
y <- cos(x)

genlimdistb2 = function(n){
  x = runif(n,0,1)
  y = cos(x)
  cbind(x, y)
}

# not sure whether to expect positive or negative association, so will carry out two-tailed test, alpha = 0.05
cor.test(x, y, alternative = "two.sided", method="spearman")

# Power Estimation

critic <- cut_off(n, 0.975)

limdists = c(genlimdista1, genlimdista2, genlimdistb1, genlimdistb2)

for (i in 1:4) {
  rho = replicate(k, scorr(limdists[[i]](n)))
  power = mean((rho < -critic) | (rho > critic))
  
  print(paste("Power of this test is ", power))
}
