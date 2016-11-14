VB_Binary_Logit_Predict01 <- function(X, w, V, invV)
{
  mu = w
  Sigma = V
  invSigma = invV

  max_iter = 500

  temp = dim(X)
  N = temp[1]
  Dx = temp[2]

  # precompute some constants
  mu_t = array(1,c(N,1)) %*% (t(mu) %*% invSigma) + 0.5 * X
  # w_t = V^-1 w + x / 2 as rows
  SigmaX = X %*% Sigma
  # W x as rows
  VxxVwt = SigmaX * (rowSums(mu_t * SigmaX) %*% array(1,c(1,Dx)) )
  # V x x^T V^T w_t as rows
  Vwt = mu_t %*% Sigma
  # V w_t as rows
  xVx = rowSums(SigmaX * X)
  # x^T V x as rows
  xVx2 = xVx ^ 2
  # x^T V x x^T V x as rows

  # start first iteration with xi = 0, lam_xi = 1/8
  xi = array(0,c(N,1))
  lam_xi = lam(xi)
  a_xi = 1 / (4 + xVx)
  w_xi = Vwt - (a_xi %*% array(1,c(1,Dx))) * VxxVwt
  logdetV_xi = - log(1 + xVx / 4)
  wVw_xi = rowSums(w_xi * (w_xi %*% invSigma) ) + rowSums(w_xi * X) ^ 2 / 4
  L_last = 0.5 * (sum(logdetV_xi) + sum(wVw_xi)) - N * log(2)

  # iterate to from xi's that maximise variational bound
  for (i in 1:max_iter){
    # update xi by EM algorithm
    xi = sqrt(xVx - a_xi * xVx2 + rowSums(w_xi * X) ^ 2)
    lam_xi = lam(xi)

    # Sherman Morrison formula and Matrix determinant lemma
    a_xi = 2 * lam_xi / (1 + 2 * lam_xi * xVx)
    w_xi = Vwt - (a_xi %*% array(1,c(1,Dx))) * VxxVwt
    logdetV_xi = - log(1 + 2 * lam_xi * xVx)

    # variational bound, omitting constant terms
    wVw_xi = rowSums(w_xi * (w_xi %*% invSigma)) + 2 * lam_xi * (rowSums(w_xi * X) ^ 2)
    L = sum(0.5 * (logdetV_xi + wVw_xi - xi) - log(1 + exp(-xi)) + lam_xi * xi ^ 2)

    # variational bound must grow!
    if(L_last > L)
        stop('Variational bound should not reduce')

    # stop if change in variation bound is < 0.001%
    if(abs(L_last - L) < abs(0.00001 * L))
        break

    L_last = L
  }

  # posterior from optimal xi's
  output = 1 / (1 + exp(-xi)) / sqrt(1 + 2 * lam_xi * xVx) * exp(0.5 * (-xi - t(mu) %*% invSigma %*% mu + wVw_xi) + lam_xi * xi ^ 2)

  return(output)
}
