rm(list=ls())

# Specify working directory (where SVB.R is)#
setwd("/path/to/SVB.R")
source('SVB.R')
library(kableExtra)

# Specify the number of cores you wish to use here
# If you need to debug, it is best to set this to 1.
CORES = 1

NUM_COLUMNS = 100

base_params = list(s0 = 5, beta_0_1 = log(71), signal_size = log(71), signal_scheme=NA)
indexed_params = lapply(c(1:NUM_COLUMNS), function(ind){
  ret <- base_params
  ret['index'] = ind
  return(ret)
})
param_list = indexed_params

lambda_hat=sqrt(71*log(4088)/6)

fits = list(
            isvb = function(X, Y, k) isvb.fit(X, Y, lambda=lambda_hat, k=k)
            mf = mf.fit,
            zz = zz.fit,
            # br = br23.fit,
            jm = jm.fit,
            oracle = oracle.fit
            )

# Define number of replicates that should be used for estimation
# Try a small number first to check everything works as expected.
n_replicates = 100

t1 = Sys.time()
# Estimation Function, if you add extra parameters to the param_sets above,
# then remember to add these. For instance if you set a new k then add k=par$k
# as an argument. 
results = lapply(param_list, function(par) estimate_stats(n=71,
                                                          p=4088,
                                                          s0=par$s0,
                                                          beta_0_1=par$beta_0_1,
                                                          signal_size=par$signal_size,
                                                          signal_scheme=par$signal_scheme,
                                                          dataset='riboflavin',
                                                          index=par$index,
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
# formatted_results = lapply(results, function(p) format_results(p, round=round, methods=names(fits), k=k))
# formatted_results = do.call(rbind, formatted_results)
formatted_results = format_riboflavin_results(results, round=round, methods=names(fits),k=1)

colnames(formatted_results) = c('Method', 'Cov. ', 'MAE', 'Length', 'Time')
res_tex = kable(formatted_results, format = 'latex', booktabs = TRUE,
             align = 'rcccc',
             row.names=TRUE,
             # digits = c(1,1,3,3),
             linesep = c('', '','', '', '', '', '\\hline'))
writeLines(as.character(res_tex), "results/riboflavin_experiment.tex")












