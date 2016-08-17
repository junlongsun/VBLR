VB_Multinominal_Logit_Fit01 <- function(X, y, typeNum, loop = 10)
{
  temp = dim(X)
  n = temp[1]
  d = temp[2]


  muHatVB = array(0, c(d, typeNum) )
  sigmaHatVB = array(0, c(d, d, typeNum) )
  for (i in 1:typeNum){
    output = VB_Binary_Logit_Fit01(X, classTransferY(y,i) );
    muHatVB[,i] = output$mu
    sigmaHatVB[,,i] = output$Sigma
  }

  #iterative
  muVB = muHatVB;
  sigmaVB = sigmaHatVB;
  for(t in 1:loop){
      #fix 2~20, predict 1
      for(k in 1:typeNum){
          yhat = array(0, c(n,1) )
          for(i in 1:n){
              yhat[i] = sum( muHatVB[,k] * X[i,] )
              temp = 0
              for(j in 1:typeNum){
                  if(j!=k){
                      temp = temp + exp( sum(muVB[,j] * X[i,])  + 0.5 * (t(X[i,]) %*% sigmaVB[,,j] %*% X[i,] ) )
                  }
              }
              yhat[i] = yhat[i] + log(temp)
          }

          muVB[,k] = ginv( t(X) %*% X ) %*% t(X) %*% yhat
      }
  }

  return(muVB)
}
