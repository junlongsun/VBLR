library(R.matlab)
data=readMat('example-for-logistic.mat')
X=data$X.ext
y=data$y
output1 = VB_Binary_Logit_Fit(X, y, a0=1e-2, b0=1e-4)
output2 = VB_Binary_Logit_Fit_startAlpha(X, y, a0=1e-2, b0=1e-4)

output3 = VB_Binary_Logit_Predict(X, output1$w, output1$V, output1$invV)
output4 = VB_Binary_Logit_Predict(X, output2$w, output2$V, output2$invV)
