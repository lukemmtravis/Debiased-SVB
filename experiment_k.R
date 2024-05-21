rm(list=ls())
setwd("/Users/lmt15/Documents/phd/Variational Inference/paper_codes")
source('SVB.R')
library(kableExtra)

CORES = 5

# Initialise parameters for experiment and combine into list
param_set_1 = list(n = 200, p = 400, k = 2,
                    s0 = 10, beta_0_1 = log(200), signal_size = log(200),
                    feature_correlation = 0.0)
param_set_2 = list(n = 200, p = 400, k = 2,
                   s0 = 10, beta_0_1 = log(200), signal_size = log(200),
                   feature_correlation = 0.5)
param_set_3 = list(n = 400, p = 800, k = 4,
                   s0 = 10, beta_0_1 = log(400), signal_size = log(400),
                   feature_correlation = 0.5)
param_set_4 = list(n = 400, p = 800, k = 4,
                   s0 = 20, beta_0_1 = c(0, 0, log(400), log(400)), signal_size=log(400), signal_scheme=NA,
                   feature_correlation = 0.25)
param_set_5 = list(n = 400, p = 1000, k = 6,
                   s0 = 10, beta_0_1 = log(400), signal_size = log(400),
                   feature_correlation = 0.5)


param_set_6 = list(n = 400, p = 1000, k = 6,
                   s0 = 10, beta_0_1 = 0, signal_size = NA, signal_scheme='uniform',
                   feature_correlation = 0.9)

param_list = list(p1 = param_set_1, p2 = param_set_2,
                  p3 = param_set_3, p4 = param_set_4,
                  p5 = param_set_5, p6 = param_set_6)

# Define methods that we want to compare
fits = list(isvb = isvb.fit,
            mf = mf.fit,
            jm = jm.fit,
            oracle=oracle.fit)

# Define number of replicates that should be used for estimation
n_replicates = 5

t1 = Sys.time()
results = lapply(param_list, function(par) estimate_stats(n=par$n,
                                                          p=par$p,
                                                          s0=par$s0,
                                                          beta_0_1=par$beta_0_1,
                                                          signal_size=par$signal_size,
                                                          signal_scheme=par$signal_scheme,
                                                          feature_correlation=par$feature_correlation,
                                                          k=par$k,
                                                          fits=fits,
                                                          n_replicates=n_replicates,
                                                          mc.cores=CORES))


t2 = Sys.time()
cat('Experiment Time: ', difftime(t2, t1, units='secs'), '\n')
 
round = 3
k = 2
formatted_results = lapply(results, function(p) format_results(p, round=round, methods=names(fits), k=k))
formatted_results = do.call(rbind, formatted_results)
colnames(formatted_results) = c('Method', 'Cov. ', 'L2', 'Volume', 'Time')
res_tex = kable(formatted_results, format = 'latex', booktabs = TRUE,
             align = 'rcccc',
             row.names=TRUE,
             # digits = c(1,1,3,3),
             linesep = c('', '', '', '\\hline'))
writeLines(as.character(res_tex), "results/experiments_k.tex")
