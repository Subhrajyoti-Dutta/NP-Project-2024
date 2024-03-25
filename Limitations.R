library("fMultivar")

set.seed(42)

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
    rho <- replicate(500, scorr(dists[[i]](n)))
    hist(rho,xlim = c(-1,1),main = names[i],xlab="S. Correlation (rho)",cex.lab=1.5,cex.main=2)
  }
  dev.off()
}


# Non-independent Samples

non_iid = function(n){
  x <- rnorm(7)
  y[1] <- rbeta(1,3,4) + x[1]
  for(i in 2:n){
    y[i] <- i*y[i-1] + x[i]
  }
  cbind(x,y)
}

k = 10000
dists = c(non_iid, Normal, Tdist)
names = c("Non-IID Model", "BVN", "BVT")

for (n in c(7,20)){
  png(file=paste(".\\viol_lim\\together_non_iid",n,".png",sep=""),width=900,height=600)
  par(mfrow=c(1,3))
  for (i in 1:3){
    rho <- replicate(500, scorr(dists[[i]](n)))
    hist(rho,xlim = c(-1,1),main = names[i],xlab="S. Correlation (rho)",cex.lab=1.5,cex.main=2)
  }
  dev.off()
}

# Only Monotonic Association can be captured
n <- 7
# eg 1
x <- rnorm(n)
y <- x^2
rho <- cor(x,y, method = 'spearman')
print(paste("rho =", rho))

# expecting positive association, so will carry out upper-tailed test, alpha = 0.05
critic <- cut_off(n, 0.95)
print(paste("critical value =", critic))

if(rho < critic) {
  print("Hence, we fail to Reject H_0, so X and Y are independent.")
}

genlimdista = function(n){
  x = rexp(n) #rexp(n) for showing monotone
  y = x^2
  cbind(x, y)
}

k = 10000
rho = replicate(k, scorr(genlimdista(n)))
cp = cut_off(n, 0.95)
power = mean(rho > cp)

print(paste("Power of this test is ", power))

# eg 2
x <- runif(n, -1, 1) #runif(n) for showing monotone
y <- cos(x)
rho <- cor(x,y, method = 'spearman')
print(paste("rho =", rho))

# not sure whether to expect positive or negative association, so will carry out two-tailed test, alpha = 0.05
critic <- cut_off(n, 0.975)
print(paste("critical value =", critic))

if((-critic < rho) & (rho < critic)) {
  print("Hence, we fail to Reject H_0, so X and Y are independent.")
}

genlimdistb = function(n){
  x = runif(n,0,1)
  y = cos(x)
  cbind(x, y)
}

rho = replicate(k, scorr(genlimdistb(n)))
cp = cut_off(n, 0.975)
power = mean((rho < -cp) | (rho > cp))

print(paste("Power of this test is ", power))

