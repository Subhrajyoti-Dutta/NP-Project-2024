calc = function(rho,l){
  c = (l+1)^2/(l*(l+2))
  a = c * (l+1)^2 * rho^2
  b = c * rho^2
  res = max(((2*b+2) + sqrt((2*b+2)^2 - 4*(b-1)*a))/(2*(b-1)), ((2*b+2) - sqrt((2*b+2)^2 - 4*(b-1)*a))/(2*(b-1)))
  return(res)
}

lambda = 5
mu = calc(-0.7, lambda)
print(mu)

check = -sqrt(lambda*(lambda+2))/(lambda+1)*sqrt((mu*(mu+2))/((lambda+1)^2 + mu*(mu-2)))
print(check)