library("fMultivar")

k = 10000
n = 20
gen_dist = function(dist1, dist2, n){
  x = dist1(n)
  y = dist2(n)
  cbind(x,y)
}
scorr = function(arr) {cor(arr[,1], arr[,2], method = "spearman")}
dist1 <- replicate(k, gen_dist(rexp,rexp,n))
iter <- apply(dist1, 3 , scorr)
hist(iter)
