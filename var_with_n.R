library("fMultivar")

gen_dist_and_scorr = function(n){
  x = rnorm(n)     #dist1
  y = rnorm(n)     #dist2
  cor(x,y, method = "spearman")
}

set.seed(42)

nvals = c(2:10,seq(12,20,2),seq(30,100,10))
k = 10000
for (n in nvals){
  rho <- replicate(k, sqrt(n-1) * gen_dist_and_scorr(n))
  png(file = paste(".\\var_with_n\\dist_n=",n,".png",sep=""), width = 960, height = 480)
  par(mfrow=c(1,2))
  hist(rho,xlim = c(-3,3),main = "",xlab="Correlation (rho)")
  plot(ecdf(rho),xlim=c(-3,3),main = "",xlab="Correlation (rho)")
  lines(seq(-3,3,1/1000),pnorm(seq(-3,3,1/1000)),col="red")
  mtext(paste("Distribution at n =",n), side = 3, line = - 2, outer = TRUE)
  dev.off()
}
