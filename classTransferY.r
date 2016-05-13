classTransferY <- function(y, class)
{
  temp = dim(y)
  n = temp[1]
  newY = array(-1, c(n,1) )
  for (i in 1:n){
    if(y[i]==class){
        newY[i] = 1
    }
  }
  return(newY)
}
