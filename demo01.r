library(R.matlab)
data=readMat('simulation.mat')
X_train = data$X.train
X_test = data$X.test
Y_train = data$Y.train
Y_test = data$Y.test

output = VB_Multinominal_Logit_Fit01(X_train, Y_train, typeNum=3, loop = 10)


library(R.matlab)
data=readMat('simulation01.mat')
X=data$X.ext
y=data$y01
output = VB_Binary_Logit_Fit01(X, y, a0=1e-2, b0=1e-4)
y2=data$y
output1 = VB_Binary_Logit_Fit(X, y2, a0=1e-2, b0=1e-4)

output2 = VB_Binary_Logit_Predict01(X, output$mu, output$Sigma, output$invSigma)
