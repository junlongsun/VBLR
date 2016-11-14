VB_Binary_Logit_Fit <- function(X, y, a0=1e-2, b0=1e-4)
{
  y[y==0]=-1

  library(MASS)

  temp = dim(X)
  N = temp[1]
  D = temp[2]

  max_iter = 500

  an = a0 + 0.5 * D
  gammaln_an_an = lgamma(an) + an

  t_w = 0.5 * t(colSums( bsxfun4times(X, y) ) )

  lam_xi = array(1,c(N,1)) / 8
  E_a = a0 / b0
  invV = E_a * diag(1, D, D) + 2 * t(X) %*% bsxfun4times(X, lam_xi)
  V = ginv(invV)
  w = V %*% t(t_w)
  bn = b0 + 0.5 * (t(w) %*% w + sum(diag(V)) )
  L_last = - N * log(2) + 0.5 * (t(w) %*% invV %*% w - log(det(invV)) ) - an / bn * b0 - an * log(bn) + gammaln_an_an


  ## update xi, bn, (V, w) iteratively
  for (i in 1:max_iter){
    # update xi by EM-algorithm
    xi = sqrt( rowSums( X * (X %*% ( V + w %*% t(w) ) ) ) )
    lam_xi = lam(xi)

    # update posterior parameters of a based on xi
    bn = b0 + 0.5 * (t(w) %*% w + sum(diag(V)) )
    E_a = an / bn

    # recompute posterior parameters of w
    invV = E_a[1,1] * diag(1, D, D) + 2 * t(X) %*% bsxfun4times(X, lam_xi)
    V = ginv(invV)
    logdetV = - log(det(invV))
    w = V %*% t(t_w)

    # variational bound, ingnoring constant terms for now
    L = - sum(log(1 + exp(- xi))) + sum(lam_xi * xi ^ 2) + 0.5 * (t(w) %*% invV %*% w + logdetV - sum(xi))   - E_a * b0 - an * log(bn) + gammaln_an_an

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

  output = list(w=w,V=V, invV = invV, logdetV = logdetV, E_a = E_a, L = L)
  return(output)
}

lam <- function(xi)
{
  out = tanh(xi / 2) / (4 * xi)
  out[is.nan(out)] = 1/8;

  return(out)
}
