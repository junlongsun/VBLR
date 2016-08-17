VB_Multinominal_Logit_Predict01 <- function(X, beta)
{
  temp = dim(X)
  n = temp[1]

  temp1 = dim(beta)
  typeNum = temp[2]

  py_class_each = array(0,c(n,typeNum))
  py_class = array(0,c(n,1))

  for(i in 1:n){
    for(k in 1:typeNum){
      py_class_each[i,k] = exp( t(beta[,k]) %*% t(X_test[i,]) )
    }
    py_class_each[i,] = py_class_each[i,] / sum(py_class_each[i,])

    py_class[i] = which.max( py_class_each[i,] )
  }
  return(py_class)
}
