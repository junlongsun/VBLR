#demo

beta = matrix(c(0.4, 0.8, 1.3, -0.5, 1.1, 1.2, 1.5, 0.9, 0.2, -0.5), nrow = 2, ncol = 5)

set.seed(999)
n=50000;
n_train = n*0.8;
K = 3;

x1 = rnorm(n, mean = 0, sd = 1)
x2 = rnorm(n, mean = 0, sd = 1)
x3 = rnorm(n, mean = 0, sd = 1)
x4 = rnorm(n, mean = 0, sd = 1)