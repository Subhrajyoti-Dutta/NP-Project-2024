library("fMultivar")

scorr = function(arr) {
  cor(arr[,1], arr[,2], method = 'spearman')
}

Normal = function(n){
  cbind(rnorm(n),rnorm(n))
}

Tdist = function(n){
  cbind(rt(n,4), rt(n,4))
}

Gamma = function(n){
  cbind(rgamma(n, 5), rgamma(n, 12))
}

Exponential = function(n) {
  cbind(rexp(n,7), rexp(n,3))
}

AMH= function(n) {
  rbiamhcop(n,0)
}

ChiSqr = function(n) {
  x1 = rchisq(n,3)
  x2 = rchisq(n,6)
  cbind(x1/(x1+x2), x1+x2)
}

set.seed(42)

k = 10000
dists = c(Normal, Tdist, Gamma, Exponential, AMH, ChiSqr)
names = c("BVN", "BVT", "BVG", "BVE", "AMH", "BVCS")

for (n in c(7,20)){
  png(file=paste(".\\dist_free\\together",n,".png",sep=""),width=1200,height=800)
  par(mfrow=c(2,3))
  for (i in 1:length(dists)){
    rho <- replicate(k, scorr(dists[[i]](n)))
    hist(rho,xlim = c(-1,1),main = names[i],xlab="Correlation (rho)",cex.lab=1.5,cex.main=2)
  }
  dev.off()
}