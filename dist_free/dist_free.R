library("fMultivar")

set.seed(3.14)

k = 10000
n = 12
gen_dist = function(){
  x = rlogis(n)     #dist1
  y = rlogis(n)     #dist2
  cbind(x,y)
}
scorr = function(arr) {cor(arr[,1], arr[,2], method = "spearman")}
dist1 <- replicate(k, gen_dist())
rho <- apply(dist1, 3 , scorr)
hist(rho,xlim = c(-1,1),main = "Logistic Distribution",xlab="Correlation (rho)")

