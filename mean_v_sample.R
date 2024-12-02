rm(list = ls())
library(ggplot2)
library(glmnet)
library(latex2exp)
library(purrr)
library(tidyverse)
library(mcmc)
library(VGAM)
library(mvtnorm)
library(knitr)
library(truncnorm)
# Note, will need to install these if have not used before
# install.packages("sparsevb_0.1.0.tar.gz", repos=NULL)
# install.packages("sparsevb2_0.1.0.tar.gz", repos=NULL)
library(sparsevb)
library(sparsevb2)
setwd('/Users/lmt15/Documents/phd/Variational Inference/paper_codes')
source("SVB.R")


vb_mean_vs_sample = function(n, p, s_0, n_samples=1000,
                             beta_0_1, lambda=NA,
                             feature_correlation=NA, noise = 'gaussian',
                             const_signal=NA){
  '
  This function computes a sample from the VB posterior, but also gives the distribution of beta_1 if one just
  uses the VB mean.
  '
  t_1 = Sys.time()
  beta_0_1 = as.numeric(beta_0_1)
  beta_0 = rep(0,p)
  beta_0[1] = beta_0_1
  S_0 = sample(c(2:p),s_0)
  if(is.na(const_signal)){
    beta_0[S_0] = rnorm(s_0)  
  }else{
    beta_0[S_0] = const_signal
  }
  if(is.na(feature_correlation)){
    X = matrix(rnorm(n*p), nrow = n, ncol = p)  
  }else{
    X = sample_block_correlated_inputs(n, p, rho = feature_correlation)
  }
  
  if(noise == 'gaussian'){
    eps = rnorm(n)  
  }else if(noise == 'laplace'){
    eps = rlaplace(n, scale=1/sqrt(2))
  }else if(noise == 'uniform'){
    eps = runif(n, min=-1, max=1)
  }
  Y = X %*% beta_0 + eps
  
  #estmate sigma_hat
  mod = cv.glmnet(X, Y)
  fit = glmnet(X, Y, lambda = mod$lambda.min)
  yhat = X %*% fit$beta
  sigma_hat = sqrt((1/(n - s_0)) * sum((yhat - Y)**2) )
  X_tilde = X/sigma_hat
  Y_tilde = Y/sigma_hat
  
  X1 = X[,1]
  X1_norm_sq = sum(X1 * X1)
  
  #make projection matrix onto span(X1)
  H = X1 %*% t(X1) / X1_norm_sq
  #and the projection matrix onto span(X1)^perp
  I_minus_H = diag(rep(1,n)) - H
  
  #make the matrix P, consisting of basis vectors of span(X1)^perp
  svd_temp = svd(I_minus_H)
  U = svd_temp$u[,1:(n-1)]
  P = t(U)
  
  #make W_check and Y_check
  W_check = P %*% I_minus_H %*% X[,2:p]
  Y_check = P %*% I_minus_H %*% Y
  #apply svb package to the check model. Note can specify lambda via prior_scale arg
  if(is.na(lambda)){
    vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
  }else{
    vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
                   lambda = lambda)
  }
  #extract relevant params from the fit
  mu_hat = vbL_W$mu
  sigmas_hat = abs(vbL_W$sigma)
  gammas_hat = vbL_W$gamma
  #sample from the variational posterior of beta_{-1}
  beta_minus_1_samples = sample_from_VB_posterior(n_samples,
                                                  mu = mu_hat,
                                                  sigma = sigmas_hat,
                                                  gamma = gammas_hat)
  
  
  #compute gammas (from DY paper, so call them gamma_prime)
  gamma_prime = apply(X[,2:p], 2, function(x) sum(X1*x)/X1_norm_sq)
  #diffs is a temporary variable which we subtract from beta_1^* to get beta_1
  diffs = beta_minus_1_samples %*% gamma_prime
  diffs_vb_mean = (mu_hat * gammas_hat) %*% gamma_prime
  
  posterior_means_beta_1 = rep(1/X1_norm_sq * t(X1) %*% Y, n_samples) - diffs
  posterior_V_beta_1 = 1/X1_norm_sq
  
  posterior_samples_beta_1 = rnorm(n_samples,
                                   posterior_means_beta_1,
                                   sd = sqrt(posterior_V_beta_1))
  
  posterior_means_beta_1_vb_mean = rep(1/X1_norm_sq * t(X1) %*% Y - diffs_vb_mean, n_samples)
  posterior_samples_beta_1_vb_mean = rnorm(n_samples,
                                           posterior_means_beta_1_vb_mean,
                                           sd = sqrt(posterior_V_beta_1))
  
  
  return(list(org_samples = posterior_samples_beta_1, mean_samples = posterior_samples_beta_1_vb_mean))
}

n = 100
p = 200
s_0 = 10
n_samples = 1000
beta_0_1 = 10
rho = 0.5

temp = vb_mean_vs_sample(n, p, s_0, 100000, beta_0_1, feature_correlation = rho)
temp$org_samples = temp$org_samples - mean(temp$org_samples) + log(n)
temp$mean_samples = temp$mean_samples - mean(temp$mean_samples) + log(n)

plot_dat = data.frame(org_sample = temp$org_samples,
                      mean_sample = temp$mean_samples) %>% gather(key = 'Sample', value = 'beta_1')

ggplot(plot_dat, aes(x = beta_1)) + 
  geom_histogram(aes(x=beta_1, y=..density.., fill=Sample), position='identity', alpha=0.5, bins=50) + 
  scale_fill_discrete(name = 'Method', labels = c("Using Mean", "Original")) +
  xlab(TeX(r'($\beta_1)')) +
  ylab('Density')

ggsave('/Users/lmt15/Documents/phd/Thesis/Writing/Figures/using_vb_mean.pdf', width = 8, height = 5, units='in')
