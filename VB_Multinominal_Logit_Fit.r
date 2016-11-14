VB_Multinominal_Logit_Fit <- function(X, y, typeNum, loop = 10, startWithAlpha = 0)
{
  temp = dim(X)
  n = temp[1]
  d = temp[2]


  muHatVB = array(0, c(d, typeNum) )
  sigmaHatVB = array(0, c(d, d, typeNum) )
  for (i in 1:typeNum){
    if (startWithAlpha){
      output = VB_Binary_Logit_Fit_startAlpha(X, classTransferY(y,i) );
    } else {
      output = VB_Binary_Logit_Fit(X, classTransferY(y,i) );
    }
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
