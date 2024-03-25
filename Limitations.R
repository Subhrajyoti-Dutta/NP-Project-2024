library("fMultivar")

#set.seed(423)

# Violation of Assumptions 
# Case of Discrete Bivariate Distributions

scorr <- function(arr){
  t = cor(arr[,1],arr[,2], method = "spearman")
  return(t)
}

Normal = function(n){
  cbind(rnorm(n),rnorm(n))
}

Tdist = function(n){
  cbind(rt(n,4), rt(n,4))
}

Geom = function(n){
  cbind(rgeom(n,0.5),rgeom(n,0.5))
}

Pois = function(n){
  cbind(rpois(n,1),rpois(n,1))
}

k = 10000
dists = c(Normal, Geom, Tdist, Pois)
names = c("BVN", "BVG", "BVT", "BVP")

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

non_iden = function(n) {
  x = rep(1,n)
  y = rep(1,n)
  
  for (i in i:n) {
    x[i] = rnorm(1, i/2, sqrt(i/4))
    y[i] = rnorm(1, i, sqrt(i/4))
  }
  cbind(x,y)
}

dists = c(non_iden, Normal, Tdist)
names = c("Non-Identical", "BVN", "BVT")

for (n in c(7,20)){
  png(file=paste(".\\viol_lim\\together_non_iden",n,".png",sep=""),width=1200,height=600)
  par(mfrow=c(1,3))
  for (i in 1:3){
    rho <- replicate(k, scorr(dists[[i]](n)))
    hist(rho,xlim = c(-1,1),main = names[i],xlab="S. Correlation (rho)",cex.lab=1.5,cex.main=2)
  }

  dev.off()
}

# Non-independent Samples

non_ind = function(n){
  x <- rnorm(n)
  y[1] <- c(runif(1))
  for(i in 2:n){
    y[i] <- i*y[1] + x[i]
  }
  cbind(x,y)
}

rho = replicate(k, scorr(non_ind(n)))
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
