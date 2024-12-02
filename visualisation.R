rm(list=ls())
setwd("/Users/lmt15/Documents/phd/Variational Inference/paper_codes")
source('SVB.R')
library(latex2exp)
Rcpp::sourceCpp("SASGibbsSampler.cpp")

MS.posterior.fit = function(X, Y, k = 1, n_samples = 5000, burnin=1000){
  t_1 = Sys.time()
  n = dim(X)[1]
  p = dim(X)[2]
  z_init = rep(0, p)
  beta_init = rep(0, p)
  mcmc_chain = sample_mc_cpp(burnin + n_samples, beta_init, z_init, X, Y, verbose=TRUE)
  beta_samples = mcmc_chain$beta[-c(1:burnin),]
  z_samples = mcmc_chain$z[-c(1:burnin),]
  posterior_samples = beta_samples*z_samples
  if(k == 1){
    posterior_samples_beta_1 = posterior_samples[,1]
    beta_hat = mean(posterior_samples_beta_1)
    credible_interval = quantile(posterior_samples_beta_1, probs = c(0.025, 0.975))
    t_2 = Sys.time()
    out = list(beta_hat=beta_hat,
               CI=as.numeric(credible_interval),
               fit_time=as.numeric(difftime(t_2, t_1, units = 'secs')))
    return(out)
  }else{
    posterior_samples_beta_k = posterior_samples[,c(1:k)]
    beta_hat = colMeans(posterior_samples_beta_k)
    cov_hat = empirical_covariance(posterior_samples_beta_k)
    t_2 = Sys.time()
    out = list(beta_hat=beta_hat,
               cov_hat=cov_hat,
               fit_time=as.numeric(difftime(t_2, t_1, units = 'secs')))
    return(out)
  }
  
}
i.posterior.fit = function(X, Y, k = 1, n_samples = 5000, burnin=1000){
  # Fit the I-SVB method on the first k coordinates.
  t_1 = Sys.time()
  n = dim(X)[1]
  p = dim(X)[2]
  if(k == 1){
    # Treat k == 1 case separately because here we will use quantiles for cred interval
    X1 = X[,1]
    X1_norm_sq = sum(X1 * X1)
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
    # if(is.na(lambda)){
    #   vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
    # }else{
    #   vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
    #                  lambda = lambda)
    # }
    # mu_hat = vbL_W$mu
    # sigmas_hat = abs(vbL_W$sigma)
    # gammas_hat = vbL_W$gamma
    # #sample from the variational posterior of beta_{-1}
    # beta_minus_1_samples = sample_from_VB_posterior(n_samples,
    #                                                 mu = mu_hat,
    #                                                 sigma = sigmas_hat,
    #                                                 gamma = gammas_hat)
    z_init = rep(0, p-1)
    beta_init = rep(0, p-1)
    mcmc_chain = sample_mc_cpp(burnin + n_samples, beta_init, z_init, W_check, Y_check, verbose=TRUE)
    beta_samples = mcmc_chain$beta[-c(1:burnin),]
    z_samples = mcmc_chain$z[-c(1:burnin),]
    beta_minus_1_samples = beta_samples*z_samples
    
    
    
    #compute gammas (from DY paper, so call them gamma_prime)
    gamma_prime = apply(X[,2:p], 2, function(x) sum(X1*x)/X1_norm_sq)
    #diffs is a temporary variable which we subtract from beta_1^* to get beta_1
    diffs = beta_minus_1_samples %*% gamma_prime
    
    # improper specific part
    posterior_means_beta_1 = rep(1/X1_norm_sq * t(X1) %*% Y, n_samples) - diffs
    posterior_V_beta_1 = 1/X1_norm_sq
    posterior_samples_beta_1 = rnorm(n_samples, posterior_means_beta_1, 
                                     sd=sqrt(posterior_V_beta_1))
    # end improper specific part
    
    beta_hat = mean(posterior_samples_beta_1)
    credible_interval = quantile(posterior_samples_beta_1, probs = c(0.025, 0.975))
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
                CI=as.numeric(credible_interval),
                fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
  }else{
    A_k = X[,1:k]
    L = chol(t(A_k)%*% A_k)
    H = A_k %*% solve(t(A_k)%*%A_k, t(A_k))
    I_minus_H = diag(rep(1,n)) - H
    #make the matrix P, consisting of basis vectors of span(X1)^perp
    svd_temp = svd(I_minus_H)
    U = svd_temp$u[,1:(n-k)]
    P = t(U)
    W_check = P %*% I_minus_H %*% X[,(k+1):p]
    Y_check = P %*% I_minus_H %*% Y
    #apply svb package to the check model. Note can specify lambda via prior_scale arg
    # if(is.na(lambda)){
    #   vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
    # }else{
    #   vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
    #                  lambda = lambda)
    # }
    # #extract relevant params from the fit
    # mu_hat = vbL_W$mu
    # sigmas_hat = abs(vbL_W$sigma)
    # gammas_hat = vbL_W$gamma
    # beta_minus_k_samples = sample_from_VB_posterior(n_samples,
    #                                                 mu = mu_hat,
    #                                                 sigma = sigmas_hat,
    #                                                 gamma = gammas_hat)
    
    z_init = rep(0, p-k)
    beta_init = rep(0, p-k)
    mcmc_chain = sample_mc_cpp(burnin + n_samples, beta_init, z_init, W_check, Y_check, verbose=TRUE)
    beta_samples = mcmc_chain$beta[-c(1:burnin),]
    z_samples = mcmc_chain$z[-c(1:burnin),]
    beta_minus_k_samples = beta_samples*z_samples
    
    
    Gamma_mat = solve(t(A_k)%*%A_k, t(A_k)) %*% X[,(k+1):p]
    
    diffs = beta_minus_k_samples %*% t(Gamma_mat)
    
    # Improper specific part
    Sigma_p = solve(t(A_k) %*% A_k, diag(k))  
    Mu_p = Sigma_p %*% t(A_k) %*% Y
    beta_star_k_samples = rmvnorm(n_samples, mean = Mu_p, sigma = Sigma_p)
    beta_k_samples = beta_star_k_samples - diffs
    # Improper specific part end
    
    beta_hat = apply(beta_k_samples, 2, mean)
    cov_hat = empirical_covariance(beta_k_samples)
    cov_hat_beta_star = Sigma_p 
    cov_hat_diffs = empirical_covariance(diffs)
    
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
                cov_hat=cov_hat,
                cov_hat_beta_star = cov_hat_beta_star,
                cov_hat_nuisance = cov_hat_diffs,
                fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
  }
  
}

set.seed(1)
n = 200
p = 400
s0 = 10
AR = TRUE
rand_coords=TRUE
rho = 0.95
k = 2

dat = make_data(n, p, s0, feature_correlation = rho,
                AR=AR, randomized_coords = rand_coords,
                beta_0_1 = 5, k = k)

X = dat$X
Y = dat$Y
beta_0 = dat$beta_0

S0 = which(beta_0 != 0)
Sigma = matrix(rho, p, p)
diag(Sigma) = 1
if(s0 > k){
  cond_on = beta_0[S0][-c(1:k)]
}else{
  cond_on = NA
}

fits = list(
  'JM' = jm.fit(X, Y, k=k),
  'I-SVB' = isvb.fit(X, Y, k=k),
  'MF' = mf.fit(X, Y, k=k),
  'Oracle' = oracle.fit(X, Y, S0, k),
  'MS.Posterior' = MS.posterior.fit(X, Y, k=k),
  'I.Posterior' = i.posterior.fit(X, Y, k=k)
)



level_sets = lapply(c(1:length(fits)), function(i) compute_level_sets(fit=fits[[i]],
                                                                      name=names(fits)[i]))
names(level_sets) = names(fits)
level_sets = do.call(rbind, level_sets)
centerings = data.frame(do.call(rbind, lapply(fits, function(fit) fit$beta_hat)))
centerings['Method'] = names(fits)

colnames(centerings) = c('beta1', 'beta2', 'Method')
truth = data.frame(beta1=beta_0[1],beta2=beta_0[2],Method='Truth')
data = list(level_sets=level_sets, centerings=centerings, truth=truth)

data = lapply(data, function(df){
  df$Scenario=rho
  return(df)
})

# fits_to_show = c("I-SVB", "MS.Posterior", "I.Posterior", "Oracle")
fits_to_show = c("I-SVB", "MS.Posterior", "I.Posterior", "Oracle")

ggplot(data$level_sets[(data$level_sets$Method %in% fits_to_show),], aes(x=beta1, y=beta2, color=Method)) + geom_path() +
  geom_point(data=centerings[centerings$Method %in% fits_to_show,], shape='cross', size=4) +
  geom_point(data=truth) +
  facet_wrap(~Scenario, nrow=1, scales='free') +
  theme(legend.position = 'bottom') +
  # geom_blank(data=dummy_dat, aes(x=beta1,y=beta2, shape=Method)) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))


#### Difference between MS and New Method ####
methods = c("I-SVB", "MS.Posterior", "I.Posterior", "MF")
filtered = level_sets %>% filter(Method %in% methods) %>% 
  mutate(Prior = ifelse(Method %in% c('MS.Posterior', 'MF'), 'MS', 'New Method')) %>% 
  mutate(Posterior = ifelse(Method %in% c('MS.Posterior', 'I.Posterior'), 'Exact', 'VB Approximation'))

filtered_centerings = centerings %>% filter(Method %in% methods) %>% 
  mutate(Prior = ifelse(Method %in% c('MS.Posterior', 'MF'), 'MS', 'New Method')) %>% 
  mutate(Posterior = ifelse(Method %in% c('MS.Posterior', 'I.Posterior'), 'Exact', 'VB Approximation'))

filtered_truth = rbind(truth, truth)
filtered_truth['Prior'] = c('MS', 'New Method')
filtered_truth['Posterior'] = 'Truth'


ggplot(filtered, aes(x=beta1, y=beta2, color=Posterior)) + 
  # geom_path(aes(linetype=Posterior)) +
  geom_path() +
  # geom_point(data=filtered_centerings, aes(color=Posterior), size=1) +
  geom_point(data=filtered_centerings, size=2 ) +
  geom_point(data=filtered_truth, size = 2, color='magenta') +
  facet_wrap(~Prior, nrow=1) +
  theme(legend.position = 'bottom') +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)')) +
  theme(
    # panel.background = element_rect(fill='transparent') 
    plot.background = element_rect(fill='transparent', color=NA)
  )

# ggsave('results/posterior_samples.pdf', units = 'in', width = 9, height = 5)
ggsave('/Users/lmt15/Documents/phd/Thesis/Writing/Figures/old_vs_new_hdlr.pdf', units = 'in', width = 8, height = 6)


# prior='New Method'
prior='MS'
ggplot(filtered %>% filter(Prior == prior), aes(x=beta1, y=beta2)) + geom_path(aes(linetype=Posterior, color=Posterior)) +
  geom_point(data=filtered_centerings %>% filter(Prior == prior), aes(color=Posterior), size=1) +
  geom_point(data=filtered_truth %>% filter(Prior == prior), shape='cross', size = 2) +
  facet_wrap(~Prior, nrow=1) +
  theme(legend.position = 'bottom') +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)')) +
  theme(
    # panel.background = element_rect(fill='transparent') 
    plot.background = element_rect(fill='transparent', color=NA)
  )
# ggsave('results/posterior_samples_new.pdf', units = 'in', width = 5, height = 5)

ggplot(filtered %>% filter(Prior == 'MS' & Posterior=='Exact'), aes(x=beta1, y=beta2)) + geom_path(aes(linetype=Posterior, color=Posterior)) +
  geom_point(data=filtered_centerings %>% filter(Prior == 'MS' & Posterior=='Exact'), aes(color=Posterior), size=1) +
  geom_point(data=filtered_truth %>% filter(Prior == 'MS'), shape='cross', size = 2) +
  facet_wrap(~Prior, nrow=1) +
  theme(legend.position = 'bottom') +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)')) +
  theme(
    # panel.background = element_rect(fill='transparent') 
    plot.background = element_rect(fill='transparent', color=NA)
  )
ggsave('results/posterior_samples_MS_no_MF.pdf', units = 'in', width = 5, height = 5)


#### Plot Shapes of Various Methods ####

####
centering=c(0, 0)

MS_approx_cov = fits[['MF']]$cov_hat
MS_cov = fits[['MS.Posterior']]$cov_hat

I_approx_cov = fits[['I-SVB']]$cov_hat
I_cov = fits[['I.Posterior']]$cov_hat

mats = list('MS'= MS_cov, 'MS.Approx' = MS_approx_cov, 'Decoupled MS' = I_cov, 'Decoupled.Approx' = I_approx_cov)
# mats = list(fits[[1]]$cov_hat, fits[[1]]$cov_hat_beta_star, fits[[1]]$cov_hat_nuisance, fits[[2]]$cov_hat)
# names=c('Abeta', 'Bbeta*', 'Diffs', 'COracle')
level_sets = lapply(c(1:length(mats)), function(i) compute_level_sets(cent=centering,
                                                                      mat=mats[[i]],
                                                                      name=names(mats)[i]))
level_sets = do.call(rbind, level_sets) %>%
  mutate(Prior = ifelse(Method %in% c('MS', 'MS.Approx'), 'MS', 'New Method')) %>% 
  mutate(Posterior = ifelse(Method %in% c('MS', 'Decoupled MS'), 'Exact', 'VB Approximation'))


ggplot(level_sets %>% filter(Prior=='MS'), aes(x=beta1, y=beta2, linetype=Posterior)) + 
  geom_path() +
  facet_wrap(~Prior, ncol=2) +
  theme(legend.position = 'bottom') +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))
ggsave('results/MS_Posterior_Sample.pdf', units = 'in', width = 7, height = 5)

ggplot(level_sets, aes(x=beta1, y=beta2, linetype=Posterior)) + 
  geom_path() +
  facet_wrap(~Prior, ncol=2) +
  theme(legend.position = 'bottom') +
  theme(
    # panel.background = element_rect(fill='transparent') 
    plot.background = element_rect(fill='transparent', color=NA)
    ) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))

# ggsave('results/Decoupled_vs_MS.pdf', units = 'in', width = 10, height = 5)
ggsave('results/Decoupled_vs_MS_Posterior_Sample.pdf', units = 'in', width = 7, height = 5)

lapply(fits, function(fit) fit$fit_time)
####



level_sets = lapply(c(1:length(fits)), function(i) compute_level_sets(fit=fits[[i]],
                                                                      name=names(fits)[i]))
names(level_sets) = names(fits)
level_sets = do.call(rbind, level_sets)
centerings = data.frame(do.call(rbind, lapply(fits, function(fit) fit$beta_hat)))
centerings['Method'] = names(fits)

colnames(centerings) = c('beta1', 'beta2', 'Method')
truth = data.frame(beta1=beta_0[1],beta2=beta_0[2],Method='Truth')
data = list(level_sets=level_sets, centerings=centerings, truth=truth)

data = lapply(data, function(df){
  df$Scenario=rho
  return(df)
})

# fits_to_show = c("I-SVB", "MS.Posterior", "I.Posterior", "Oracle")
fits_to_show = c("I-SVB", "MS.Posterior", "I.Posterior", "Oracle")

ggplot(data$level_sets[(data$level_sets$Method %in% fits_to_show),], aes(x=beta1, y=beta2, color=Method)) + geom_path() +
  geom_point(data=centerings[centerings$Method %in% fits_to_show,], shape='cross', size=4) +
  geom_point(data=truth) +
  facet_wrap(~Scenario, nrow=1, scales='free') +
  theme(legend.position = 'bottom') +
  # geom_blank(data=dummy_dat, aes(x=beta1,y=beta2, shape=Method)) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))


#### Plots all centeres at 0 ####
level_sets = lapply(c(1:length(fits)), function(i) compute_level_sets(fit=fits[[i]],
                                                                      name=names(fits)[i], recenter = TRUE))
names(level_sets) = names(fits)
level_sets = do.call(rbind, level_sets)
centerings = data.frame(do.call(rbind, lapply(fits, function(fit) fit$beta_hat)))
centerings['Method'] = names(fits)

colnames(centerings) = c('beta1', 'beta2', 'Method')
truth = data.frame(beta1=beta_0[1],beta2=beta_0[2],Method='Truth')
data = list(level_sets=level_sets, centerings=centerings, truth=truth)

data = lapply(data, function(df){
  df$Scenario=rho
  return(df)
})

#### Difference between MS and New Method ####
methods = c("I-SVB", "MS.Posterior", "I.Posterior", "MF")
filtered = level_sets %>% filter(Method %in% methods) %>% 
  mutate(Prior = ifelse(Method %in% c('MS.Posterior', 'MF'), 'MFVI', 'New Method')) %>% 
  mutate(Posterior = ifelse(Method %in% c('MS.Posterior', 'I.Posterior'), 'Exact', 'VB Approximation'))

filtered_centerings = centerings %>% filter(Method %in% methods) %>% 
  mutate(Prior = ifelse(Method %in% c('MS.Posterior', 'MF'), 'MS', 'New Method')) %>% 
  mutate(Posterior = ifelse(Method %in% c('MS.Posterior', 'I.Posterior'), 'Exact', 'VB Approximation'))

filtered_truth = rbind(truth, truth)
filtered_truth['Prior'] = c('MFVI', 'New Method')
filtered_truth['Posterior'] = 'Truth'


ggplot(filtered, aes(x=beta1, y=beta2)) + geom_path(aes(linetype=Posterior)) +
  # geom_point(data=filtered_centerings, aes(color=Posterior), size=1) +
  # geom_point(data=filtered_centerings, size=1) +
  # geom_point(data=filtered_truth, shape='cross', size = 2) +
  facet_wrap(~Prior, nrow=1) +
  theme(legend.position = 'bottom') +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)')) +
  theme(
    # panel.background = element_rect(fill='transparent') 
    plot.background = element_rect(fill='transparent', color=NA)
  )

ggsave('/Users/lmt15/Documents/phd/Thesis/Writing/Figures/MF_vs_Full.pdf', 
       units = 'in', width = 8, height = 5)
