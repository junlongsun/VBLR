library(R.matlab)
data=readMat('example-for-multinominal-logistic.mat')
X_train = data$X.train
X_test = data$X.test
Y_train = data$Y.train
Y_test = data$Y.test

output1 = VB_Multinominal_Logit_Fit(X_train, Y_train, typeNum=3, loop = 10)
output2 = VB_Multinominal_Logit_Fit(X_train, Y_train, typeNum=3, loop = 10, startWithAlpha = 1)
