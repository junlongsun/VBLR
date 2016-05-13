lam <- function(xi)
{

  # returns 1 / (4 * xi) * tanh(xi / 2)
  #divby0_w = warning('query', 'MATLAB:divideByZero');

  #warning('off', 'MATLAB:divideByZero');

  out = tanh(xi / 2) / (4 * xi)
  #warning(divby0_w.state, 'MATLAB:divideByZero');
  # fix values where xi = 0
  out[is.nan(out)] = 1/8;

  return(out)

}
