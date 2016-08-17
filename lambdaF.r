lambdaF <- function(x)
{
  sigmoid = 1/(1+exp(-x))
  out = (1/(2*x)) * (sigmoid-0.5)
  out[is.nan(out)] = 1/8

  return(out)
}
