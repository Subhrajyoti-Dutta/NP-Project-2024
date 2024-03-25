
a=1
b=10
for (i in 1:20){
  for (j in seq(a,b,10^(-i))){
    rangx = seq(-5,5,0.001)
    rangy = exp(rangx)-rangx*j
    
    print(c(j,cor(rangx, rangy), cor(rangx, rangy, method='spearman')))
    if (cor(rangx, rangy) < 0) {
      b = j
      a = b - 10^(-i)
      break
    }
  }
}

rangx = seq(-5,5,0.00001)
rangy = exp(rangx)-rangx*t

print(c(cor(rangx, rangy),
        cor(rangx, rangy, method='spearman')))
