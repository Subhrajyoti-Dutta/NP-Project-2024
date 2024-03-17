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
  png(file = paste("dist_n=",n,".png",sep=""))
  hist(rho,xlim = c(-1,1),main = "Normal Distribution",xlab="Correlation (rho)")
  dev.off()
}