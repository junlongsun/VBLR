bsxfun4times <- function(X, y)
{
  temp = dim(X)
  n = temp[1]
  d = temp[2]
  output = array(0,c(n,d))
  for (i in 1:d){
  output[,i] = X[,i] * y
  }
  return(output)
}
