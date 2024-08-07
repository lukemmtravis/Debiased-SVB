rm(list=ls())

# Specify working directory (where SVB.R is)#
setwd("/path/to/SVB.R")
source('SVB.R')
library(kableExtra)

# Specify the number of cores you wish to use here
# If you need to debug, it is best to set this to 1.
CORES = 1


# Initialise parameters for experiment and combine into list.
# All param sets must have the same arguments, but if you want to leave one out it is
# okay to set parameter = NA, for instance if you want beta_0_1 to be generated randomly
# then set beta_0_1 = NA.
# Below are for Table 1 in the paper.
param_set_1 = list(n = 100, p = 1000,
                    s0 = 3, beta_0_1 = log(100), signal_size = log(100), signal_scheme=NA,
                    feature_correlation = 0.0, noise_var=1)
param_set_2 = list(n = 100, p = 1000,
                   s0 = 3, beta_0_1 = log(100), signal_size = log(100), signal_scheme=NA,
                   feature_correlation = 0.5, noise_var=16)
param_set_3 = list(n = 400, p = 1500,
                   s0 = 32, beta_0_1 = NA, signal_size = NA, signal_scheme = 'normal',
                   feature_correlation = 0.25, noise_var=1)
param_set_4 = list(n = 200, p = 800,
                   s0 = 10, beta_0_1 = log(200), signal_size= log(200), signal_scheme=NA,
                   feature_correlation = 0.9, noise_var = 1)
param_set_5 = list(n = 200, p = 1000,
                   s0 = 5, beta_0_1 = 0, signal_size = log(100), signal_scheme=NA,
                   feature_correlation = 0.5, noise_var = 1)
param_set_6 = list(n = 500, p = 1000,
                   s0 = 10, beta_0_1 = NA, signal_size = NA, signal_scheme = 'uniform',
                   feature_correlation = 0.5, noise_var = 1)

param_list = list(p1 = param_set_1, p2 = param_set_2,
                  p3 = param_set_3, p4 = param_set_4,
                  p5 = param_set_5, p6 = param_set_6)

# Define methods that we want to compare, this should be a list where you write
# `name you want to have` = relevant fitting function from SVB.R

# For comparing various SVB methods
# fits = list(isvb = isvb.fit, 
#             gsvb.1 = function(X, Y) gsvb.fit(X, Y, prior_sd=1), 
#             gsvb.4 = function(X, Y) gsvb.fit(X, Y, prior_sd=4),
#             lsvb.1 = function(X, Y) lsvb.fit(X, Y, prior_sd=1), 
#             lsvb.4 = function(X, Y) lsvb.fit(X, Y, prior_sd=4))

# For comparison to frequentist methods
fits = list(isvb = isvb.fit, 
            mf = mf.fit,
            zz = zz.fit,
            jm = jm.fit)


# Define number of replicates that should be used for estimation
# Try a small number first to check everything works as expected, but set to 500
# to match what we do in the paper
n_replicates = 500

t1 = Sys.time()
# Estimation Function, if you add extra parameters to the param_sets above,
# then remember to add these. For instance if you set a new k then add k=par$k
# as an argument. 
results = lapply(param_list, function(par) estimate_stats(n=par$n,
                                                          p=par$p,
                                                          s0=par$s0,
                                                          beta_0_1=par$beta_0_1,
                                                          signal_size=par$signal_size,
                                                          feature_correlation=par$feature_correlation,
                                                          signal_scheme=par$signal_scheme,
                                                          noise_var=par$noise_var,
                                                          fits=fits,
                                                          n_replicates=n_replicates,
                                                          mc.cores=CORES))


t2 = Sys.time()
cat('Experiment Time: ', difftime(t2, t1, units='secs'), '\n')

# Below formats the results into a latex table, but if you just want to look at 
# the results in R then you can look at the formatted_results object that is produced.
if('k' %in% names(param_list[[1]])){
  if(param_list[[1]]$k > 1){
    k = 2
  }else{
    k = 1
  }
}else{
  k = 1
}
round = 3
formatted_results = lapply(results, function(p) format_results(p, round=round, methods=names(fits), k=k))
formatted_results = do.call(rbind, formatted_results)
colnames(formatted_results) = c('Method', 'Cov. ', 'MAE', 'Length', 'Time')
res_tex = kable(formatted_results, format = 'latex', booktabs = TRUE,
             align = 'rcccc',
             row.names=TRUE,
             # digits = c(1,1,3,3),
             linesep = c('', '', '', '\\hline'))
writeLines(as.character(res_tex), "results/experiment_results.tex")

