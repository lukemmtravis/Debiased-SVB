rm(list=ls())
setwd("/Users/lmt15/Documents/phd/Variational Inference/paper_codes")
source('SVB.R')

n = 100
p = 200
s0 = 10
dat = make_data(n, p, s0, beta_0_1 = log(n), signal_size = log(n), k = 1)

X = dat$X
Y = dat$Y
beta_0 = dat$beta_0

isvb = isvb.fit(X, Y)
jm = jm.fit(X, Y)

# All fitting methods now have the same structure
cat('I-SVB Estimate: ', isvb$beta_hat, '\n')
cat('I-SVB Credible Interval: ', isvb$CI, '\n')
cat('I-SVB Fit Time: ', isvb$fit_time, 'seconds \n')

cat('JM Estimate: ', jm$beta_hat, '\n')
cat('JM Credible Interval: ', jm$CI, '\n')
cat('JM Fit Time: ', jm$fit_time, 'seconds \n')

## However if we go into higher dimensions, we return a covariance matrix
# instead of a credible interval
dat = make_data(n, p, s0, beta_0_1 = log(n), signal_size = log(n), k = 2)
isvb = isvb.fit(dat$X, dat$Y, k = 2)
jm = jm.fit(dat$X, dat$Y, k = 2)

isvb$cov_hat
jm$cov_hat
