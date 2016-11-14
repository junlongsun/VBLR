VB_Binary_Logit_Fit_startAlpha <- function(X, y, a0=1e-2, b0=1e-4)
{
  y[y==-1]=0

  library(MASS)

  temp = dim(X)
  N = temp[1]
  D = temp[2]

  max_iter = 500

  E_a = a0 / b0
  aN = a0 + 0.5 * D

  lambda_epsilon = array(1,c(N,1)) / 8
  invSigmaN = E_a * diag(1, D, D) + 2 * t(X) %*% bsxfun4times(X, lambda_epsilon)
  SigmaN = ginv(invSigmaN)
  muN = SigmaN %*% colSums( bsxfun4times(X, y-0.5) )

  bN = b0 + 0.5 * (t(muN) %*% muN + sum(diag(SigmaN)) )

  L_last = - N * log(2) + 0.5 * (t(muN) %*% invSigmaN %*% muN - log(det(invSigmaN)) ) - aN / bN * b0 - aN * log(bN) + lgamma(aN) + aN

  ## update xi, bn, (V, w) iteratively
  for (i in 1:max_iter){
    # update xi by EM-algorithm
    epsilon2 = rowSums( X * (X %*% ( SigmaN + muN %*% t(muN) ) ) )
    epsilon = sqrt( epsilon2 )
    lambda_epsilon = lambdaF(epsilon)

    # update posterior parameters of a based on xi
    bN = b0 + 0.5 * (t(muN) %*% muN + sum(diag(SigmaN)) )
    E_a = aN / bN

    # recompute posterior parameters of w
    invSigmaN = E_a[1,1] * diag(1, D, D) + 2 * t(X) %*% bsxfun4times(X, lambda_epsilon)
    SigmaN = ginv(invSigmaN)
    logdetSigmaN = - log(det(invSigmaN))
    muN = SigmaN %*% colSums( bsxfun4times(X, y-0.5) )

    # variational bound, ingnoring constant terms for now
    L = 0.5*t(muN) %*% invSigmaN %*% muN + 0.5*logdetSigmaN - sum(log(1 + exp(- epsilon))) - 0.5*sum(epsilon)  + sum(lambda_epsilon * epsilon ^ 2) - E_a * b0 - aN * log(bN) + lgamma(aN) + aN
    # either stop if variational bound grows or change is < 0.001%

    if ((L_last > L) || (abs(L_last - L) < abs(0.00001 * L)))
        break
    end

    L_last = L
  }

  if(i == max_iter)
      warning('Bayesian logistic regression reached maximum number of iterations.')

  ## add constant terms to variational bound
  L = L - lgamma(a0) + a0 * log(b0)

  output = list(mu=muN, Sigma=SigmaN, invSigma = invSigmaN, logdetSigma = logdetSigmaN, E_a = E_a, L = L)
  return(output)
}
