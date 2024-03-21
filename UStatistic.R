#Comparison with U statistics


x <- cbind(seq(1:5),seq(1,5))
combn(seq(1:5),2)

compute.kernel.sy <- function(arr){
  count <- 0
  index <- combn(seq(1:3),2)
  for (i in 1:3){
    index.first <- index[1,i] #1
    index.second <- index[2,i] #3
    index.third <- setdiff(1:3, c(index.first,index.second)) #2
    if(arr[index.first,1]< arr[index.second,1]){
      if(arr[index.first,2] < arr[index.third,2]){
        count =count + 1
      }
    }
    else{
      if(arr[index.second,2] < arr[index.third,2]){
        count = count+1
      }
    }  
  }
 
 count/6-(1/4)
}

#n obs

compute.u.statistic <- function(arr){
  n = nrow(arr)
  indices <- t(combn(1:n,3))
  val <- 0               
  for (i in 1:nrow(indices)){
    val <- val + compute.kernel.sy(arr[indices[i,],])
  }
  6*val/(n*(n-1)*(n-2))

}
