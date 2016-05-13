VB_Multinominal_Logit_Fit <- function(X, Y, typeNum, loop = 10)
{
  temp = dim(X)
  n = temp[1]
  d = temp[2]

  
  muHatVB = array(0, c(d, typeNum) )
  sigmaHatVB = array(0, c(d, d, typeNum) )
  for (i in 1:typeNum){
    output = VB_Binary_Logit_Fit(X, classTransferY(y,i) );
    muHatVB[,i] = output$w
    sigmaHatVB[,,i] = output$V
  }

  #iterative
  muVB = muHatVB;
  sigmaVB = sigmaHatVB;
  for(t in 1:loop){
      #fix 2~20, predict 1
      for(k in 1:typeNum){
          yhat = array(0, c(n,1) )
          for(i in 1:n){
              yhat[i] = ( t(muHatVB[,k]) %*% t(X[i,]) )
              temp = 0
              for(j in 1:typeNum){
                  if(j!=k){
                      temp = temp + exp( t(muVB[,j]) * t(X[i,])  + 0.5 * (X[i,] * sigmaVB[,,j]) %*% t(X[i,]) )
                  }
              }
              yhat[i] = yhat[i] + log(temp)
          }

          muVB[,k] = ginv( t(X) %*% X ) %*% X %*% yhat
      }
  }

  return(muVB)
}
