library("fMultivar")

gen_dist_and_scorr = function(n){
  x = rnorm(n)     #dist1
  y = rnorm(n)     #dist2
  cor(x,y, method = "spearman")
}

set.seed(42)

k = 10000
n = 12
dists = c(rexp, rnorm, rlogis, function(n){rpois(n,2)})
names = c("Exponential(mean=1)","Normal(0,1)","Logistic(0,1)","Poisson(2)")

for (i in 1:length(dists)){
  rho <- replicate(k, gen_dist_and_scorr(n))
  png(file = paste(".\\dist_free\\",names[i],".png",sep=""), width = 480, height = 480)
  hist(rho,xlim = c(-1,1),main = names[i],xlab="Correlation (rho)")
  dev.off()
}