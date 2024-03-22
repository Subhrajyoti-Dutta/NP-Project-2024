library("fMultivar")

gen_dist_and_scorr = function(n){
  x = rnorm(n)     #dist1
  y = rnorm(n)     #dist2
  cor(x,y, method = "spearman")
}

Exponential = function(n) {
  rexp(n,2)
}

Normal = function(n){
  rnorm(n,1,2)
}

Logistic = function(n){
  rlogis(n,1,2)
}

Gamma = function(n){
  rgamma(n,shape=3,scale=4)
}

Cauchy = function(n){
  rcauchy(n,2,3)
}

Chisqr = function(n){
  rchisq(n,3)
}

set.seed(42)

k = 10000
n = 12
dists = c(Exponential, Normal, Logistic, Gamma, Cauchy, Chisqr)
names = c("Exponential(Mean=1)","Normal(1,2)","Logistic(1,2)","Gamma(Shape=3,Scale=4)","Cauchy(2,3)","Chisqr(3)")

png(file=".\\dist_free\\together.png",width=1200,height=800)
par(mfrow=c(2,3))
for (i in 1:length(dists)){
  rho <- replicate(k, gen_dist_and_scorr(n))
  hist(rho,xlim = c(-1,1),main = names[i],xlab="Correlation (rho)")
}
dev.off()