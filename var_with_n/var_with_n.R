library("fMultivar")

set.seed(3.14)

nvals = c(2:10,seq(12,20,2),seq(30,100,10))
k = 10000
for (n in nvals){
  gen_dist = function(){
    x = rnorm(n)     #dist1
    y = rnorm(n)     #dist2
    cbind(x,y)
  }
  scorr = function(arr) {cor(arr[,1], arr[,2], method = "spearman")}
  dist1 <- replicate(k, gen_dist())
  rho <- apply(dist1, 3 , scorr)
  png(file = paste("dist_n=",n,".png",sep=""), width = 960, height = 480)
  par(mfrow=c(1,2))
  hist(rho,xlim = c(-3,3),main = "",xlab="Correlation (rho)")
  plot(ecdf(rho),xlim=c(-3,3),main = "",xlab="Correlation (rho)")
  mtext(paste("Distribution at n =",n), side = 3, line = - 2, outer = TRUE)
  dev.off()
}