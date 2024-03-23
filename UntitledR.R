rho = c()
rang = seq(-1,1,0.02)
for (i in rang){
  # print(cor(rbiamhcop(10000, i)))
  rho = c(rho,cor(rbiamhcop(10000, i),method="spearman")[1,2])
}
mat = cbind(rang,rang^2,rang^3,rang^4,rang^5,rang^6,rang^7,rang^8,rang^9,rang^10)
plot(rang, rho, col="red")
model = lm(rho ~ mat-1)
lines(rang, mat %*% unname(model$coeff))
model$coeff
cbind(rang,rho)