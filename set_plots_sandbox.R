rm(list=ls())
setwd("/Users/lmt15/Documents/phd/Variational Inference/paper_codes")
source('SVB.R')
library(latex2exp)
simulate_set = function(n, p, s0, rho=0.0, seed=1, scenario=1, AR=FALSE, rand_coords=TRUE, dataset='generated'){
  k = 2
  set.seed(seed)
  dat = make_data(n, p, s0, beta_0_1=5,
                  k=2, feature_correlation=rho, AR=AR, randomized_coords=rand_coords,
                  dataset=dataset)
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  
  # Compute Oracle quantities
  S0 = which(beta_0 != 0)
  Sigma = matrix(rho, p, p)
  diag(Sigma) = 1
  if(s0 > k){
    cond_on = beta_0[S0][-c(1:k)]
  }else{
    cond_on = NA
  }
  
  fits = list(
    # 'JM' = fit.jmo(X, Y, Sigma=Sigma, k=k),
    'JM' = jm.fit(X, Y, k=k),
    # 'JM18' = jm18.fit(X, Y, Sigma, k=k),
    'I-SVB' = isvb.fit(X, Y, k=k),
    'MF' = mf.fit(X, Y, k=k),
    # 'Freq' = freq.fit(X, Y, k=k),
    'Oracle' = oracle.fit(X, Y, S0, k)
    # 'Conditional Oracle' = oracle.fit.cond(X, Y, S0=S0, k=k, cond_on=cond_on)
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
  
  return(data)
}

visualise_covariances = function(n, p, s0, rho=0.0, seed=1, scenario=1, AR=FALSE, rand_coords=TRUE){
  k = 2
  set.seed(seed)
  dat = make_data(n, p, s0, beta_0_1=5, k=2, feature_correlation=rho, AR=AR, randomized_coords=rand_coords)
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  
  # Compute Oracle quantities
  S0 = which(beta_0 != 0)
  Sigma = matrix(rho, p, p)
  diag(Sigma) = 1
  if(s0 > k){
    cond_on = beta_0[S0][-c(1:k)]
  }else{
    cond_on = NA
  }
  fits = list(
    'I-SVB' = isvb.fit(X, Y, k=k),
    'Oracle' = oracle.fit(X, Y, S0, k)
  )
  centering=c(0, 0)
  mats = list(fits[[1]]$cov_hat, fits[[1]]$cov_hat_beta_star, fits[[1]]$cov_hat_nuisance, fits[[2]]$cov_hat)
  names=c('Abeta', 'Bbeta*', 'Diffs', 'COracle')
  length(mats)
  level_sets = lapply(c(1:length(mats)), function(i) compute_level_sets(cent=centering,
                                                                        mat=mats[[i]],
                                                                        name=names[i]))
  names(level_sets) = names
  level_sets = do.call(rbind, level_sets)
  
  level_sets$rho=rho
  level_sets$Scenario=scenario
  level_sets$s0 = s0
  if(!AR){
    level_sets$DataGen = 'Fully Correlated' 
  }else{
    if(rand_coords){
      level_sets$DataGen = 'AR Random Signal Placement'   
    }else{
      level_sets$DataGen = 'AR'   
    }
  }
  
  
  return(level_sets)
}

# p1 = list(n = 200, p = 400, s0=10, rho = 0, seed = 1, scenario = 1)
# p2 = list(n = 400, p = 800, s0=10, rho = 0.2, seed = 2, scenario = 2)
# p3 = list(n = 600, p = 1600, s0=15, rho = 0.5, seed = 3, scenario = 3)
# # p4 = list(n = 600, p = 400, s0=2, rho = 0.9, seed = 3, scenario = 4)
# params = list(p1 = p1, p2 = p2, p3 = p3)

p1 = list(n = 200, p = 400, s0=10, rho = 0.0, seed = 1, scenario = 1, AR=TRUE)
p2 = list(n = 200, p = 400, s0=10, rho = 0.5, seed = 1, scenario = 2, AR=TRUE)
p3 = list(n = 200, p = 400, s0=10, rho = 0.95, seed = 1, scenario = 3, AR=TRUE)
params = list(p1 = p1, p3 = p3)
# params = list(p3 = p3)


data = lapply(params, function(par) 
  simulate_set(n=par$n, p=par$p, s0=par$s0, rho=par$rho,
               seed=par$seed, scenario=par$scenario, AR = par$AR))

level_sets = do.call(rbind, lapply(data, function(x) x$level_sets))
centerings = do.call(rbind, lapply(data, function(x) x$centerings))
truth = do.call(rbind, lapply(data, function(x) x$truth))

# If one wants to define axis limits individually, 
# use the below dummy dataframe and uncomment the geom_blank line below
# dummy_p1 = data.frame(beta1 = c(4.5, 5.2), beta2 = c(4.5, 5.2), Scenario=1, Method='Z')
# dummy_p2 = data.frame(beta1 = c(4.7, 5.2), beta2 = c(4.7, 5.2), Scenario=2, Method='Z')
# dummy_p3 = data.frame(beta1 = c(3.9, 5.2), beta2 = c(3.9, 5.2), Scenario=3, Method='Z')
# dummy_dat = rbind(dummy_p1, dummy_p2, dummy_p3)


ggplot(level_sets, aes(x=beta1, y=beta2, color=Method)) + geom_path() +
  geom_point(data=centerings) +
  geom_point(data=truth) +
  facet_wrap(~Scenario, nrow=1, scales='free') +
  theme(legend.position = 'bottom') +
  # geom_blank(data=dummy_dat, aes(x=beta1,y=beta2, shape=Method)) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))

# ggsave('results/pres_credible_regions.pdf', units = 'in', width = 10, height = 5)
ggsave('/Users/lmt15/Documents/phd/Thesis/Writing/Figures/AR_credible_regions.pdf', units = 'in', width = 6, height = 4)


## Compare Covariances
n = 200
p = 400
s0 = 10
p1 = list(n = 200, p = 400, s0=s0, rho = 0.0, seed = 1, scenario = 1, AR=TRUE, RS=TRUE)
p2 = list(n = 200, p = 400, s0=s0, rho = 0.5, seed = 1, scenario = 2, AR=TRUE, RS=TRUE)
p3 = list(n = 200, p = 400, s0=s0, rho = 0.95, seed = 1, scenario = 3, AR=TRUE, RS=TRUE)

p4 = list(n = 200, p = 400, s0=s0, rho = 0.0, seed = 1, scenario = 1, AR=TRUE, RS=FALSE)
p5 = list(n = 200, p = 400, s0=s0, rho = 0.5, seed = 1, scenario = 2, AR=TRUE, RS=FALSE)
p6 = list(n = 200, p = 400, s0=s0, rho = 0.95, seed = 1, scenario = 3, AR=TRUE, RS=FALSE)

p7 = list(n = 200, p = 400, s0=s0, rho = 0.0, seed = 1, scenario = 1, AR=FALSE, RS=TRUE)
p8 = list(n = 200, p = 400, s0=s0, rho = 0.5, seed = 1, scenario = 2, AR=FALSE, RS=TRUE)
p9 = list(n = 200, p = 400, s0=s0, rho = 0.95, seed = 1, scenario = 3, AR=FALSE, RS=TRUE)

params = list(
  p1 = p1, p3 = p3
  # p4 = p4, p5 = p5, p6 = p6
  # p7 = p7, p8 = p8, p9 = p9
)
# params = list(p1 = p1, p2 = p2, p3 = p3)
# params = list(p2 = p2)


data = lapply(params, function(par) 
  visualise_covariances(n=par$n, p=par$p, s0=par$s0, rho=par$rho,
                        seed=par$seed, scenario=par$scenario, AR = par$AR, rand_coords = par$RS))

level_sets = do.call(rbind, data)

ggplot(level_sets, aes(x=beta1, y=beta2, color=Method)) + geom_path() +
  facet_wrap(~rho, nrow=1) +
  theme(legend.position = 'bottom') +
  scale_color_discrete(name='Covariance Matrix', labels=c(TeX(r'($\beta$)'), TeX(r'($\beta^*$)'), 'Oracle', 'Debiasing Quantity')) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)')) +
  theme(panel.spacing.x = unit(6, "mm"))
# +
#   ggtitle(paste('n: ', n, 'p: ', p, 's0: ', s0))

# ggsave('results/pres_cov_structure.pdf', units = 'in', width = 10, height = 5)
ggsave('/Users/lmt15/Documents/phd/Thesis/Writing/Figures/cov_structure_AR.pdf', units = 'in', width = 6, height = 5)



## Compare s0 behaviour
n = 200
p = 400
rho=0.5
p1 = list(n = n , p = p , s0=2  , rho = rho , seed = 1 , scenario = 1 , AR=TRUE  , RS=TRUE)
p2 = list(n = n , p = p , s0=8  , rho = rho , seed = 1 , scenario = 2 , AR=TRUE  , RS=TRUE)
p3 = list(n = n , p = p , s0=32 , rho = rho , seed = 1 , scenario = 3 , AR=TRUE  , RS=TRUE)
p4 = list(n = n , p = p , s0=64 , rho = rho , seed = 1 , scenario = 4 , AR=TRUE  , RS=TRUE)
p5 = list(n = n , p = p , s0=2  , rho = rho , seed = 1 , scenario = 1 , AR=TRUE  , RS=FALSE)
p6 = list(n = n , p = p , s0=8  , rho = rho , seed = 1 , scenario = 2 , AR=TRUE  , RS=FALSE)
p7 = list(n = n , p = p , s0=32 , rho = rho , seed = 1 , scenario = 4 , AR=TRUE  , RS=FALSE)
p8 = list(n = n , p = p , s0=128 , rho = rho , seed = 1 , scenario = 3 , AR=TRUE  , RS=FALSE)
p9 = list(n = n , p = p , s0=2  , rho = rho , seed = 1 , scenario = 1 , AR=FALSE , RS=TRUE)
p10 = list(n = n , p = p , s0=8  , rho = rho , seed = 1 , scenario = 2 , AR=FALSE , RS=TRUE)
p11 = list(n = n , p = p , s0=32 , rho = rho , seed = 1 , scenario = 3 , AR=FALSE , RS=TRUE)
p12 = list(n = n , p = p , s0=128 , rho = rho , seed = 1 , scenario = 4 , AR=FALSE , RS=TRUE)

params = list(p1 = p1, p2 = p2, p3 = p3,
              p4 = p4, p5 = p5, p6 = p6,
              p7 = p7, p8 = p8, p9 = p9,
              p10 = p10, p11 = p11, p12 = p12)
# params = list(p1 = p1, p2 = p2, p3 = p3)
# params = list(p2 = p2)


data = lapply(params, function(par) 
  visualise_covariances(n=par$n, p=par$p, s0=par$s0, rho=par$rho,
                        seed=par$seed, scenario=par$scenario, AR = par$AR, rand_coords = par$RS))

level_sets = do.call(rbind, data)

ggplot(level_sets, aes(x=beta1, y=beta2, color=Method)) + geom_path() +
  facet_wrap(DataGen~s0, nrow=3) +
  theme(legend.position = 'bottom') +
  scale_color_discrete(name='Covariance Matrix', labels=c(TeX(r'($\beta$)'), TeX(r'($\beta^*$)'), 'Oracle', 'Debiasing Quantity')) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)')) 
# +
#   ggtitle(paste('n: ', n, 'p: ', p, 's0: ', s0))

ggsave('results/s0_cov_structure.pdf', units = 'in', width = 10, height = 5)



#### Look at riboflavin dataset
simulate_set_riboflavin = function(s0, indices=c(1:2), seed=1){
  k = 2
  set.seed(seed)
  n = 71
  p = 4088
  dat = make_data(n, p, s0, beta_0_1=5, index=indices,
                  k=2, dataset='riboflavin')
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  n = dim(X)[1]
  p = dim(X)[2]
  
  # Compute Oracle quantities
  S0 = which(beta_0 != 0)
  
  fits = list(
    'JM' = jm.fit(X, Y, k=k),
    'I-SVB' = isvb.fit(X, Y, k=k),
    'MF' = mf.fit(X, Y, k=k),
    'Oracle' = oracle.fit(X, Y, S0, k)
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
    df$Scenario=paste(indices, collapse=', ')
    return(df)
  })
  
  return(data)
}

n = 71
p = 4088
num_plots = 10

p1 = list(s0 = 10, indices = sample(c(1:p), 2), seed = 1)
p2 = list(s0 = 10, indices = sample(c(1:p), 2), seed = 2)
p3 = list(s0 = 10, indices = sample(c(1:p), 2), seed = 3)
params = list(p1 = p1, p2 = p2, p3 = p3)
# params = list(p2 = p2)


data = lapply(params, function(par) 
  simulate_set_riboflavin(s0=par$s0, indices=par$indices, seed=par$seed))

level_sets = do.call(rbind, lapply(data, function(x) x$level_sets))
centerings = do.call(rbind, lapply(data, function(x) x$centerings))
truth = do.call(rbind, lapply(data, function(x) x$truth))

# If one wants to define axis limits individually, 
# use the below dummy dataframe and uncomment the geom_blank line below
# dummy_p1 = data.frame(beta1 = c(4.5, 5.2), beta2 = c(4.5, 5.2), Scenario=1, Method='Z')
# dummy_p2 = data.frame(beta1 = c(4.7, 5.2), beta2 = c(4.7, 5.2), Scenario=2, Method='Z')
# dummy_p3 = data.frame(beta1 = c(3.9, 5.2), beta2 = c(3.9, 5.2), Scenario=3, Method='Z')
# dummy_dat = rbind(dummy_p1, dummy_p2, dummy_p3)


ggplot(level_sets, aes(x=beta1, y=beta2, color=Method)) + geom_path() +
  geom_point(data=centerings, shape='cross', size=4) +
  geom_point(data=truth) +
  facet_wrap(~Scenario, nrow=1, scales='free') +
  theme(legend.position = 'bottom') +
  # geom_blank(data=dummy_dat, aes(x=beta1,y=beta2, shape=Method)) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))

ggsave('results/riboflavin_credible_regions.pdf', units = 'in', width = 10, height = 5)

n = 200
p = 400
rho = 0.95
AR = TRUE
rand_coords = TRUE
k = 3
set.seed(1)
dat = make_data(n, p, s0, beta_0_1=5,
                k=k, feature_correlation=rho, AR=AR, randomized_coords=rand_coords,
                dataset='generated')
X = dat$X
Y = dat$Y
beta_0 = dat$beta_0

# Compute Oracle quantities
S0 = which(beta_0 != 0)
Sigma = matrix(rho, p, p)
diag(Sigma) = 1

fits = list(
  'JM' = jm.fit(X, Y, k=k),
  'I-SVB' = isvb.fit(X, Y, k=k),
  'MF' = mf.fit(X, Y, k=k),
  'Oracle' = oracle.fit(X, Y, S0, k)
)

fits$'I-SVB'$cov_hat
