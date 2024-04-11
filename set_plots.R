rm(list=ls())
setwd("/Users/lmt15/Documents/phd/Variational Inference/paper_codes")
source('SVB.R')
library(latex2exp)
simulate_set = function(n, p, s0, rho=0.0, seed=1, scenario=1){
  k = 2
  set.seed(seed)
  dat = make_data(n, p, s0, beta_0_1=5, k=2, feature_correlation=rho)
  X = dat$X
  Y = dat$Y
  beta_0 = dat$beta_0
  
  # Compute Oracle quantities
  S0 = which(beta_0 != 0)
  Sigma = matrix(rho, p, p)
  diag(Sigma) = 1
  
  fits = list(
    # 'JM' = fit.jmo(X, Y, Sigma=Sigma, k=k),
              # 'JM' = fit.jm(X, Y, k=k),
              'JM' = jm.fit(X, Y, k=k),
              'I-SVB' = isvb.fit(X, Y, k=k),
              'MF' = mf.fit(X, Y, k=k),
              'Oracle' = oracle.fit(X, Y, S0, k)
              # 'OracleCond' = fit.oracle.cond(X, Y, S0, k)
              )
  level_sets = lapply(c(1:length(fits)), function(i) compute_level_sets(fits[[i]], names(fits)[i]))
  names(level_sets) = names(fits)
  level_sets = do.call(rbind, level_sets)
  centerings = data.frame(do.call(rbind, lapply(fits, function(fit) fit$beta_hat)))
  centerings['method'] = names(fits)
  colnames(centerings) = c('beta1', 'beta2', 'method')
  truth = data.frame(beta1=beta_0[1],beta2=beta_0[2],method='truth')
  data = list(level_sets=level_sets, centerings=centerings, truth=truth)
  
  data = lapply(data, function(df){
    df$Scenario=scenario
    return(df)
  })
  
  return(data)
}

p1 = list(n = 200, p = 400, s0=10, rho = 0, seed = 1, scenario = 1)
p2 = list(n = 600, p = 800, s0=10, rho = 0.1, seed = 1, scenario = 2)
p3 = list(n = 600, p = 1600, s0=15, rho = 0.5, seed = 3, scenario = 3)
params = list(p1 = p1, p2 = p2, p3 = p3)
# params = list(p1 = p1)

data = lapply(params, function(par) 
  simulate_set(n=par$n, p=par$p, s0=par$s0, rho=par$rho, seed=par$seed, scenario=par$scenario))

level_sets = do.call(rbind, lapply(data, function(x) x$level_sets))
centerings = do.call(rbind, lapply(data, function(x) x$centerings))
truth = do.call(rbind, lapply(data, function(x) x$truth))

ggplot(level_sets, aes(x=beta1, y=beta2, color=method)) + geom_path() +
  geom_point(data=centerings, shape='cross') +
  geom_point(data=truth) +
  # geom_path(data=adj_vecs)+
  facet_wrap(~Scenario) +
  xlab(TeX(r'($\beta_1$)')) + 
  ylab(TeX(r'($\beta_2$)'))
# ggsave('results/3_credible_regions.pdf', units = 'in', width = 12, height = 6)



