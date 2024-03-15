library("fMultivar")

k = 10000
n = 12
gen_dist = function(){
  x = rexp(n)     #dist1
  y = rexp(n)     #dist2
  cbind(x,y)
}
scorr = function(arr) {cor(arr[,1], arr[,2], method = "spearman")}
dist1 <- replicate(k, gen_dist())
iter <- apply(dist1, 3 , scorr)
hist(iter)

# 