library(ggplot2)
library(glmnet)
library(latex2exp)
library(purrr)
library(tidyverse)
library(mcmc)
library(VGAM)
library(mvtnorm)
library(truncnorm)
library(pbapply)
library(mcreplicate)

library(sparsevb)
library(sparsevb2)
source('lasso_inference.r')

#### Data Generation Functions ####

make_Sigma = function(p, rho, block_size=NA, AR=FALSE){
	if(!AR){
		if(is.na(block_size)){
      block_size = p
    }
	  if(p %% block_size != 0){
	      stop('Block size does not divide p')
	    }
	  B = matrix(rho, nrow = block_size, ncol = block_size)
	  diag(B) = rep(1,block_size)
	    
	  Sigma = matrix(0, nrow = p, ncol = p)
	  for(block_num in 1:(p/block_size) ){
	    start_ix = (block_num - 1)*block_size + 1
	    end_ix = block_num*block_size
	    Sigma[start_ix:end_ix, start_ix:end_ix] = B
	  }
  	return(Sigma)	
	}else{
		Sigma = matrix(nrow=p, ncol=p)
		for(i in c(1:p)){
			for(j in c(1:p)){
				Sigma[i, j] = rho**abs(i - j)
			}
		}
		return(Sigma)
	}
}

sample_block_correlated_inputs = function(n, p, rho, block_size=NA, AR=FALSE){
    Sigma = make_Sigma(p, rho, block_size, AR)
    U = chol(Sigma)
    normal_sample = matrix(rnorm(n*p), nrow = n, ncol = p)
    X = normal_sample%*%U
    return(X)
}

empirical_covariance = function(samples){
  n = dim(samples)[1]
  p = dim(samples)[2]
  mean_samples = apply(samples, 2, mean)
  demeaned_samples = t(apply(samples, 1, function(v) v - mean_samples))
  sigma_hat = matrix(apply(apply(demeaned_samples, 1, function(v) v %*% t(v)), 1, mean), nrow=p)
  return(sigma_hat)
}


#### Experiment Functions ####
make_beta = function(p, s0, k = 1, beta_0_1 = NA, signal_size = NA, signal_scheme='normal'){
	# Produce beta_0 with given structure
    beta_0 = rep(0,p)
    S_0 = c(1:k, sample(c((k+1):p), s0-k))
    if(is.na(signal_size)){
    	if(signal_scheme == 'normal'){
    		beta_0[S_0] = rnorm(s0)  
			}else if(signal_scheme=='uniform'){
				beta_0[S_0] = runif(s0, -5, 5)
			}else if(signal_scheme=='very_weak'){
				if(s0 < 10){
					stop('For very_weak signal_scheme, s0 must be larger than 10')
				}
				weak_indices = S_0[1:5]
				strong_indices = S_0[6:s0]
				beta_0[weak_indices] = 0.01
				beta_0[strong_indices] = log(p)/s0 
				return(beta_0)
			}else{
				stop('signal_scheme, must be one in {normal, uniform, very_weak}.')
			}
    }else{
      beta_0[S_0] = signal_size
    }
    if(any(!is.na(beta_0_1))){
    	beta_0[1:k] = beta_0_1
    }
    return(beta_0)
}

make_X = function(n, p, 
				feature_correlation = 0, block_size = NA, rescale_first_column = FALSE, AR = FALSE){
	# Produce X with given structure.
	if(is.na(feature_correlation)){
      X = matrix(rnorm(n*p), nrow = n, ncol = p)  
      if(rescale_first_column){
        X[,1] = X[,1]/sqrt(sum(X[,1]*X[,1]))
      }
    }else{
      X = sample_block_correlated_inputs(n, p, rho = feature_correlation, block_size = block_size, AR = AR)
      if(rescale_first_column){
        X[,1] = X[,1]/sqrt(sum(X[,1]*X[,1]))
      }
    }
    return(X)
}

make_data = function(n, p, s0, noise_var=1, noise = 'gaussian',
					 beta_0_1 = NA, signal_size = NA, signal_scheme = 'normal',
					 feature_correlation = 0, block_size = NA, rescale_first_column = FALSE, AR = FALSE,
					 k = 1, dataset = 'generated', index=1){
	# Produce X, Y and beta_0 with given structure for simulations.

	# make design
	if(dataset=='generated'){
		X = make_X(n, p, feature_correlation=feature_correlation, block_size=block_size,
										 rescale_first_column=FALSE, AR=AR)	
		beta_0 = make_beta(p, s0, k, beta_0_1=beta_0_1, signal_size=signal_size, signal_scheme=signal_scheme)
	}else if(dataset=='riboflavin'){
		X = as.matrix(read.csv('riboflavin_normalized.csv'))
		# reorder columns so that index is first.
		X = X[,c(index, c(1:p)[-index])]
		n = dim(X)[1]
		p = dim(X)[2]
		beta_0 = rep(0,p)
    if(is.na(signal_size)){
    	if(signal_scheme == 'normal'){
    		beta_0[1:s0] = rnorm(s0)  
			}else if(signal_scheme=='uniform'){
				beta_0[1:s0] = runif(s0, -5, 5)
			}else{
				stop('signal_scheme, must be one in {normal, uniform}.')
			}
    }else{
      beta_0[1:s0] = signal_size
    }
    if(!is.na(beta_0_1)){
    	beta_0[1:k] = beta_0_1
    }
    # beta_0 = make_beta(p, s0, k, beta_0_1=beta_0_1, signal_size=signal_size, signal_scheme=signal_scheme)
	}else{
		stop('dataset must be one of `generated` or `riboflavin`.')
	}

	# make beta_0
	# beta_0 = make_beta(p, s0, k, beta_0_1=beta_0_1, signal_size=signal_size, signal_scheme=signal_scheme)
	

	# make eps
	if(noise == 'gaussian'){
      eps = rnorm(n, sd = sqrt(noise_var))  
    }else if(noise == 'laplace'){
      eps = rlaplace(n, scale=sqrt(noise_var/2))
    }else if(noise == 'uniform'){
      eps = runif(n, min=-sqrt(3*noise_var), max=sqrt(3*noise_var))
    }else{
    	stop('noise must be in {gaussian, laplace, uniform}.')
    }

  # make Y
  Y = X%*%beta_0 + eps

  return(list(Y=Y, X=X, beta_0=beta_0)) 
}

#### SVB Functions ####
sample_from_VB_posterior = function(n, mu, sigma, gamma){
	# Produce n samples from the VB posterior with variational parameters mu, sigma, gamma
	p = length(mu)
	samples = sapply(c(1:p), function(i) (runif(n) <= gamma[i])*rnorm(n,
	                                                                  mean = mu[i],
	                                                                  sd = sigma[i]))
	return(samples)
}

isvb.fit = function(X, Y, n_samples = 1000, lambda = NA, k = 1){
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
		if(is.na(lambda)){
			vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
		}else{
			vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
		               	   lambda = lambda)
		}
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
		beta_minus_k_samples = sample_from_VB_posterior(n_samples,
		                                                mu = mu_hat,
		                                                sigma = sigmas_hat,
		                                                gamma = gammas_hat)
		
		
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
		# cov_hat = Sigma_p + Gamma_mat %*% diag(sigmas_hat**2) %*% t(Gamma_mat)/n_samples
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					cov_hat=cov_hat,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
	}
	
}

isvb.qr.fit = function(X, Y, n_samples = 1000, lambda = NA, k = 1){
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
    
    w = X1 + sign(X1[1]) * sqrt(X1_norm_sq)*c(1, rep(0, n-1))
    Q = diag(n) - 2 * (w %*% t(w))/sum(w*w)
    
    P = t(Q[,2:n])
    
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
    w = X1 + sign(X1[1]) * sqrt(X1_norm_sq)*c(1, rep(0, n-1))
    Q = diag(n) - 2 * (w %*% t(w))/sum(w*w)
    P = t(Q[,2:n])

    W_check = P %*% I_minus_H %*% X[,(k+1):p]
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
    beta_minus_k_samples = sample_from_VB_posterior(n_samples,
                                                    mu = mu_hat,
                                                    sigma = sigmas_hat,
                                                    gamma = gammas_hat)
    
    
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
    # cov_hat = Sigma_p + Gamma_mat %*% diag(sigmas_hat**2) %*% t(Gamma_mat)/n_samples
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
                cov_hat=cov_hat,
                fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
  } 
}

lsvb.fit = function(X, Y, n_samples = 1000, lambda = NA, prior_sd=1/sqrt(2), k = 1){
	# Fit the L-SVB method with prior sd given by prior_sd. 
	if(k != 1){
		stop('LSVB Method only implemented for k = 1 at present.')
	}
	n = dim(X)[1]
	p = dim(X)[2]
	t_1 = Sys.time()
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
	if(is.na(lambda)){
		vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
	}else{
		vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
	               	   lambda = lambda)
	}
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

	# laplace specific part
	tau = 1/prior_sd
	X1_norm = sqrt(sum(X1*X1))
	# Compute parameters for truncated normal distributions
	mu_plus = (sum(X1 * Y) + tau)/X1_norm_sq
	mu_minus = (sum(X1 * Y) - tau)/X1_norm_sq
	V = 1/X1_norm_sq
	# Compute weights, note lower.tail=FALSE gives Phi tilde, lower.tail=TRUE gives Phi
	w1 = exp(tau* sum(X1*Y)/X1_norm_sq) * pnorm(mu_plus/sqrt(V), lower.tail = FALSE)
	w2 = exp(-tau* sum(X1*Y)/X1_norm_sq) * pnorm(mu_minus/sqrt(V), lower.tail = TRUE)

	d1_samples = rtruncnorm(n_samples, a=-Inf, b=0, mean=mu_plus, sd=sqrt(V))
	d2_samples = rtruncnorm(n_samples, a=0, b=Inf , mean=mu_minus,sd=sqrt(V))

	prob_1 = w1/(w1+w2)
	mask = (runif(n_samples) <= prob_1)*1
	posterior_samples_beta_1_star = mask*d1_samples + (1-mask)*d2_samples
	posterior_samples_beta_1 = posterior_samples_beta_1_star - diffs
	# end laplace specific part
    beta_hat = mean(posterior_samples_beta_1)
    credible_interval = quantile(posterior_samples_beta_1, probs = c(0.025, 0.975))
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
    			CI=as.numeric(credible_interval),
    			fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
}

gsvb.fit = function(X, Y, n_samples = 1000, lambda = NA, prior_sd = 1, k = 1){
	# Fit the G-SVB method with prior sd given by prior_sd.
	if(k != 1){
		stop('GSVB Method only implemented for k = 1 at present.')
	}
	n = dim(X)[1]
	p = dim(X)[2]
	t_1 = Sys.time()
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
	if(is.na(lambda)){
		vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear")  
	}else{
		vbL_W<-svb.fit(W_check,Y_check,tol = 10e-5, max_iter = 10^4, family ="linear",
	               	   lambda = lambda)
	}
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
	
	# gaussian specific part
	posterior_means_beta_1 = rep(prior_sd^2/(1+prior_sd^2*X1_norm_sq) * t(X1) %*% Y, n_samples) - diffs
    posterior_V_beta_1 = prior_sd^2/(1+X1_norm_sq*prior_sd^2)
    posterior_samples_beta_1 = rnorm(n_samples, posterior_means_beta_1, 
    								sd=sqrt(posterior_V_beta_1))
    # end gaussian specific part

    beta_hat = mean(posterior_samples_beta_1)
    credible_interval = quantile(posterior_samples_beta_1, probs = c(0.025, 0.975))
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
    			CI=as.numeric(credible_interval),
    			fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
}

mf.fit = function(X, Y, n_samples = 1000, lambda = NA, k = 1){
	# Fit the mean-field approach.
	n = dim(X)[1]
	p = dim(X)[2]
	t_1 = Sys.time()   
	if(k == 1){
		if(is.na(lambda)){
			vbL<-sparsevb::svb.fit(X,Y,tol = 10e-5, max_iter = 10^4, family ="linear")  
		}else{
			vbL<- sparsevb::svb.fit(X,Y,tol = 10e-5, max_iter = 10^4, family ="linear",
		                         lambda = lambda)
		}

		#extract relevant params from the fit
		mu_hat = vbL$mu
		sigmas_hat = abs(vbL$sigma)
		gammas_hat = vbL$gamma

		beta_hat = mu_hat[1]*gammas_hat[1]
		credible_interval = beta_hat + c(qnorm(0.025, mean = 0, sd = sigmas_hat[1]),
		                                        qnorm(0.975, mean = 0, sd = sigmas_hat[1]))
	    t_2 = Sys.time()
	    return(list(beta_hat=beta_hat,
	    			CI=as.numeric(credible_interval),
	    			fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
    }else{
    	if(is.na(lambda)){
        	vbL<-sparsevb2::svb.fit(X,Y,tol = 10e-5, max_iter = 10^4, family ="linear")  
      	}else{
        	vbL<- sparsevb2::svb.fit(X,Y,tol = 10e-5, max_iter = 10^4, family ="linear",
                                 lambda = lambda)
      	}
		#extract relevant params from the fit
		mu_hat = vbL$mu
		sigmas_hat = abs(vbL$sigma)
		gammas_hat = vbL$gamma

		vb_posterior_mean = mu_hat[1:k]
		vb_posterior_var = diag(sigmas_hat[1:k]**2)
		vb_posterior_inc = gammas_hat[1:k]

		beta_hat = vb_posterior_mean*vb_posterior_inc
		cov_hat = vb_posterior_var
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					cov_hat=cov_hat,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
	}
}

mf2.fit = function(X, Y, n_samples = 1000, lambda = NA){
	# Fit the mean field method, but with a slab (no spike) on the first coordinate.
	n = dim(X)[1]
	p = dim(X)[2]
	t_1 = Sys.time()   
	if(is.na(lambda)){
		vbL<-sparsevb2::svb.fit(X,Y,tol = 10e-5, max_iter = 10^4, family ="linear")  
	}else{
		vbL<- sparsevb2::svb.fit(X,Y,tol = 10e-5, max_iter = 10^4, family ="linear",
	                         lambda = lambda)
	}

	#extract relevant params from the fit
	mu_hat = vbL$mu
	sigmas_hat = abs(vbL$sigma)
	gammas_hat = vbL$gamma

	beta_hat = mu_hat[1]*gammas_hat[1]
	credible_interval = beta_hat + c(qnorm(0.025, mean = 0, sd = sigmas_hat[1]),
	                                        qnorm(0.975, mean = 0, sd = sigmas_hat[1]))
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
    			CI=as.numeric(credible_interval),
    			fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
}

zz.fit = function(X, Y, k = 1){
	# Fit the ZZ method
	if(k != 1){
		stop('Only implemented for k = 1')
	}
	n = dim(X)[1]
	p = dim(X)[2]
	t_1 = Sys.time()
	gamma = 0.95
	cv_fit = cv.glmnet(X, Y)
	beta_lasso = glmnet(X, Y,lambda =cv_fit$lambda.min)$beta
	X_minus_1 = X[,-1]
	cv_fit_gamma = cv.glmnet(x = X_minus_1, y = X[,1])
	gamma_hat = glmnet(x = X_minus_1, y = X[,1], lambda = cv_fit_gamma$lambda.min)$beta
	
	z_1 = as.matrix(X[,1] - X_minus_1 %*% gamma_hat)
	beta_hat = as.numeric(beta_lasso[1] + t(z_1)%*%(Y - X%*%beta_lasso)/(t(z_1)%*%X[,1]))
	CI = c(beta_hat - qnorm((1+gamma)/2)/sqrt(t(z_1) %*% z_1),
	     beta_hat + qnorm((1+gamma)/2)/sqrt(t(z_1) %*% z_1))
	t_2 = Sys.time()
	return(list(beta_hat=beta_hat,
    			CI=CI,
    			fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
}

jm.fit = function(X, Y, k = 1){
	# Fit JM method using their code. Note, also returns a vector add_length which quantifies 
	# additional uncertainty about the centering.
	if(k == 1){
		n = dim(X)[1]
		p = dim(X)[2]
		t_1 = Sys.time()
		jm_fit = SSLassoFirstCoordAdjusted(X, Y)
		beta_hat = jm_fit$unb.coef
		CI = c(jm_fit$low.lim, jm_fit$up.lim)
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					CI=CI,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}else{
		n = dim(X)[1]
		p = dim(X)[2]
		t_1 = Sys.time()
		jm_fit = SSLassok(X, Y, k)
		beta_hat = jm_fit$beta_hat
		cov_hat = jm_fit$cov_hat
		addlength = jm_fit$addlength
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					cov_hat=cov_hat,
					add_length=addlength,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}
}

jm18.fit = function(X, Y, Sigma, k = 1){
	# Fit JM method using their code. Note, also returns a vector add_length which quantifies 
	# additional uncertainty about the centering.
	n = dim(X)[1]
	p = dim(X)[2]
	if(k == 1){
		t_1 = Sys.time()
    # gamma = 0.95
    # cv_fit = cv.glmnet(X, Y)
    # beta_lasso = glmnet(X, Y,lambda =cv_fit$lambda.min)$beta
    # M = solve(Sigma, diag(rep(1, p)))
    # beta_hat = (beta_lasso + (1/n)* M %*% t(X) %*% (Y - X%*%beta_lasso))[1]
    # Sigma_hat = t(X)%*%X / n
    # C = M %*% Sigma_hat  %*% M
    # CI = c(beta_hat - qnorm((1+gamma)/2)*sqrt(C[1,1]/n),
    #        beta_hat + qnorm((1+gamma)/2)*sqrt(C[1,1]/n))
    jm18_fit = SSLassoKnownSigmaK(X, Y, Sigma, k)
 		CI = c(jm18_fit$low.lim, jm18_fit$up.lim)
		t_2 = Sys.time()
		return(list(beta_hat=jm18_fit$beta_hat,
					CI=CI,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}else{
		t_1 = Sys.time()
    # gamma = 0.95
    # cv_fit = cv.glmnet(X, Y)
    # beta_lasso = glmnet(X, Y,lambda =cv_fit$lambda.min)$beta
    # M = solve(Sigma, diag(rep(1, p)))
    # beta_hat = (beta_lasso + (1/n)* M %*% t(X) %*% (Y - X%*%beta_lasso))[1:k]
    # Sigma_hat = t(X)%*%X / n
    # cov_hat = (M %*% Sigma_hat  %*% M/n)[1:k,1:k]
    # # explore adding length later
    # add_length = rep(0, k)
    jm18_fit = SSLassoKnownSigmaK(X, Y, Sigma, k)
		t_2 = Sys.time()
		return(list(beta_hat=jm18_fit$beta_hat,
					cov_hat=jm18_fit$cov_hat,
					add_length=jm18_fit$add_length,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}
}

jm.org.fit = function(X, Y, k = 1){
	# Fit the JM method without optimising their code for just k coordinates (much slower).
	if(k == 1){
		n = dim(X)[1]
		p = dim(X)[2]
		t_1 = Sys.time()
		jm_fit = SSLassoFirstCoord(X, Y)
		beta_hat = jm_fit$unb.coef
		CI = c(jm_fit$low.lim, jm_fit$up.lim)
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					CI=CI,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}else{
		n = dim(X)[1]
		p = dim(X)[2]
		t_1 = Sys.time()
		jm_fit = SSLasso(X, Y)
		beta_hat = jm_fit$unb.coef[1:k]
		cov_hat = jm_fit$M_est[1:k, 1:k]
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					cov_hat=cov_hat,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}
}

jmo.fit = function(X, Y, Sigma, k = 1){
	# Fitting function for the JM method with known Sigma (the oracle).
	t_1 = Sys.time()
	n = dim(X)[1]
	p = dim(X)[2]
	
	gamma = 0.95
	cv_fit = cv.glmnet(X, Y)
	beta_lasso = glmnet(X, Y,lambda =cv_fit$lambda.min)$beta
	M = solve(Sigma, diag(rep(1, p)))
	beta_hat = beta_lasso + (1/n)* M %*% t(X) %*% (Y - X%*%beta_lasso)
	beta_hat = beta_hat[1:k]
	Sigma_hat = t(X)%*%X / n
	C =  M %*% Sigma_hat  %*% M/n
	if(k == 1){
		t_2 = Sys.time()
		CI = c(beta_hat - qnorm((1+gamma)/2)*sqrt(C[1,1]),
            	beta_hat + qnorm((1+gamma)/2)*sqrt(C[1,1]))
		return(list(beta_hat=beta_hat,
					CI=CI,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
	}else{
		t_2 = Sys.time()
		return(list(beta_hat = beta_hat,
				cov_hat = C[1:k, 1:k],
				fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
	}
}

br23.fit = function(X, Y, k = 1, delta = 1){
	n = dim(X)[1]
	p = dim(X)[2]

	t_1 = Sys.time()
	tau=1
	gamma=0.95
	if(k > 1){
		stop('Only implemented for k = 1.')
	}
	Xv <- X[,1]
  Xnotv <- X[, -1, drop = FALSE]
  M <- Xnotv %*% t(Xnotv)
  qOpt <- solve( delta*diag(n) + M, Xv)
  scale <- t(qOpt) %*% Xv
  variance <- tau * (t(qOpt) %*% qOpt) / (t(qOpt) %*% Xv)^2
  beta_hat <- t(qOpt) %*% Y / t(qOpt) %*% Xv

  cv = qnorm((1+gamma)/2)
  upperlimit <- beta_hat + cv * sqrt(variance)
  lowerlimit <- beta_hat - cv * sqrt(variance)
  CI = c(lowerlimit, upperlimit)

  t_2 = Sys.time()

  return(list(beta_hat=beta_hat,
					CI=CI,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
}

freq.fit = function(X, Y, k = 1){
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
		mod = cv.glmnet(W_check, Y_check)
	  fit = glmnet(W_check, Y_check, lambda = mod$lambda.min)
	  beta_minus_1_hat = fit$beta

	  beta_hat = as.numeric(t(X1) %*% (Y - X[,-1] %*% beta_minus_1_hat) / X1_norm_sq)

	  z_95 = as.numeric(qnorm(0.975))
    credible_interval = c(
    	beta_hat - z_95 / sqrt(X1_norm_sq),
    	beta_hat + z_95 / sqrt(X1_norm_sq)
    	)
    t_2 = Sys.time()
    return(list(beta_hat=beta_hat,
    			CI=credible_interval,
    			fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))	
	}else{
		
		A_k = X[,1:k]
    L = chol(t(A_k)%*% A_k)
    H = A_k %*% solve(t(A_k)%*%A_k, t(A_k))
		I_minus_H = diag(rep(1,n)) - H
		# #make the matrix P, consisting of basis vectors of span(X1)^perp
		svd_temp = svd(I_minus_H)
		U = svd_temp$u[,1:(n-k)]
		P = t(U)
		W_check = P %*% I_minus_H %*% X[,(k+1):p]
		Y_check = P %*% I_minus_H %*% Y
		mod = cv.glmnet(W_check, Y_check)
	  fit = glmnet(W_check, Y_check, lambda = mod$lambda.min)
	  beta_minus_k_hat = fit$beta

	  beta_hat = as.numeric(solve(t(A_k)%*% A_k ) %*% t(A_k) %*% (Y - X[,-c(1:k)] %*% beta_minus_k_hat))
		cov_hat = solve(t(A_k) %*% A_k )

		# beta_hat = apply(beta_k_samples, 2, mean)
		# cov_hat = empirical_covariance(beta_k_samples)
		# # cov_hat = Sigma_p + Gamma_mat %*% diag(sigmas_hat**2) %*% t(Gamma_mat)/n_samples
		t_2 = Sys.time()
		return(list(beta_hat=beta_hat,
					cov_hat=cov_hat,
					fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
	}
	
}

oracle.fit = function(X, Y, S0, k = 1){
	if(k == 1){
		gamma=0.95
		t_1 = Sys.time()
		X_S0 = X[,S0]
		oracle_centering = (solve(t(X_S0) %*% X_S0) %*% t(X_S0) %*% Y)[1]
		oracle_matrix = solve(t(X_S0) %*% X_S0)

		CI = c(oracle_centering - qnorm((1+gamma)/2)*sqrt(oracle_matrix[1,1]),
					 oracle_centering + qnorm((1+gamma)/2)*sqrt(oracle_matrix[1,1]))
		t_2 = Sys.time()
	  return(list(beta_hat = oracle_centering, 
	  						CI = CI,
	  						fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
	}else{
		t_1 = Sys.time()
		if(length(intersect(S0, c(1:k))) == 0){
			t_2 = Sys.time()
			return(list(beta_hat = rep(0, k), 
	  						cov_hat = matrix(0, nrow=k, ncol=k),
	  						fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
		}else if(length(intersect(S0, c(1:k))) == k){
			X_S0 = X[,S0]
			oracle_centering = (solve(t(X_S0) %*% X_S0) %*% t(X_S0) %*% Y)[1:k]
			oracle_matrix = solve(t(X_S0) %*% X_S0)[1:k, 1:k]
			t_2 = Sys.time()
		  return(list(beta_hat = oracle_centering, 
	  						cov_hat = oracle_matrix,
	  						fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
		}else{
			X_S0 = X[,S0]
			full_oracle_centering = (solve(t(X_S0) %*% X_S0) %*% t(X_S0) %*% Y)
			full_oracle_matrix = solve(t(X_S0) %*% X_S0)
			# Assume that beta_0 consists of zero coordinates followed by non-zero coordinates
			nnz = length(intersect(S0, c(1:k)))
			nz = k - nnz
			beta_hat = c(rep(0, nz), full_oracle_centering[1:nnz])
			cov_hat = matrix(0, nrow=k, ncol=k)
			cov_hat[(nz+1):k, (nz+1):k] = full_oracle_matrix[1:nnz, 1:nnz]
			t_2 = Sys.time()
		  return(list(beta_hat = beta_hat, 
	  						cov_hat = cov_hat,
	  						fit_time=as.numeric(difftime(t_2, t_1, units = 'secs'))))
		}
	}
}

oracle.fit.cond = function(X, Y, S0, cond_on = NA, k = 1){
	if(length(S0) == k){
		return(oracle.fit(X, Y, S0, k))
	}
	X_S0 = X[,S0]

	oracle_centering = (solve(t(X_S0) %*% X_S0) %*% t(X_S0) %*% Y)
	oracle_matrix = solve(t(X_S0) %*% X_S0)

	mu_k = oracle_centering[1:k]
	mu_minus_k = oracle_centering[-c(1:k)]

	M_k = oracle_matrix[1:k, 1:k]
	M_k_minus_k = oracle_matrix[1:k, -c(1:k)]
	M_minus_k = oracle_matrix[-c(1:k), -c(1:k)]

	beta_hat = as.numeric(mu_k + M_k_minus_k %*% solve(M_minus_k) %*% (cond_on - mu_minus_k))
	cov_hat = M_k - M_k_minus_k %*% solve(M_minus_k) %*% t(M_k_minus_k)

	return(list(beta_hat = beta_hat, cov_hat = cov_hat))
}

#### Compute and Evaluate Fits ####

# All functions in this section take as an argument `fit`,
# which should be a list containing {beta_hat, CI, fit_time} in the 1D case and
# {beta_hat, cov_hat, fit_time} in the kD case, and an argument `truth' which represents
# the true value of the whole parameter beta_0.
#

# 1 dimensional metrics
MAE = function(fit, truth){
	est = fit$beta_hat
	return(abs(est-truth))
}

CI_hit = function(fit, truth){
	ind = 1*(fit$CI[1] <= truth && truth <= fit$CI[2])
	return(ind)
}

int_length = function(fit, truth){
	return(fit$CI[2] - fit$CI[1])
}

# k dimensional metrics
L2 = function(fit, truth){
	est = fit$beta_hat
	k = length(est)
	beta_0_k = truth[1:k]
	return( sqrt(sum( (est-beta_0_k)**2) ) ) 
}

kHit = function(fit, truth){
	# Note, can cope with first diagonal elements of the covariance being 0, then just checks beta_hat
	# in those coordinates is equal to the truth, and performs the routine on the submatrix
	est = fit$beta_hat
	k = length(est)
	beta_0_k = truth[1:k]
	chi_zval = qchisq(0.95, df = k)
	mat = fit$cov_hat

	if('add_length' %in% names(fit)){
		add_length = fit$add_length
  	cov=cov.adj(mat, add_length, chi_zval)
  }else{
  	cov=mat
  }

  zero_inds = which(diag(mat) == 0)
	non_zero_inds = which(diag(mat) != 0)

	if(length(non_zero_inds) == k){
		test_val = t(beta_0_k - est) %*% solve(cov) %*% (beta_0_k - est) 
		hit = test_val <= chi_zval
	}else if(length(zero_inds) == k){
	  hit = all(est == beta_0_k)
	}else{
	  zero_equality = all(est[zero_inds] == beta_0_k[zero_inds])
	  sub_cov = cov[non_zero_inds, non_zero_inds]
	  sub_est = est[non_zero_inds]
	  sub_beta_0 = beta_0_k[non_zero_inds]
	  test_val = t(sub_beta_0 - sub_est) %*% solve(sub_cov) %*% (sub_beta_0 - sub_est) 
	  sub_hit = test_val <= chi_zval
	  hit = sub_hit && zero_equality
	}
	return(hit)
}

volume = function(fit){
	k = length(fit$beta_hat)
	mat = fit$cov_hat
	chi_zval = qchisq(0.95, df = k)
	if('add_length' %in% names(fit)){
		add_length = fit$add_length
  	cov=cov.adj(mat, add_length, chi_zval)
  }else{
  	cov=mat
  }
	svd_cov = svd(cov)
	vol = sqrt(prod(svd_cov$d))
	return(vol)
}

#### Utility Functions ####

cov.adj = function(mat, add_length, chi_zval){
  svd_mat = svd(mat)
  evals = svd_mat$d
  V = svd_mat$u
 

  if(min(svd(V**2)$d) < 1e-15){
  	# happens with very small probability, in this case cannot add length correctly in
  	# e1 and e2 directions, as both dimension are lengthened simultaneously
  	lambda_tilde = evals + add_length
  }else{
  	# below adds length in directions of canonical basis (correct)
	  M_inv = solve(mat)
	  ys = 1/diag(M_inv) + 2*add_length*sqrt(1/(chi_zval * diag(M_inv))) + add_length**2/chi_zval
  	# cat('mat:  \n', mat,  '\n-------\n')
	  # cat('V**2: \n', V**2, '\n-------\n')
	  # cat('1/ys: \n', 1/ys, '\n-------\n')
	  lambda_tilde_inverse= solve(V**2, 1/ys)
	  lambda_tilde = abs(1/lambda_tilde_inverse)	
  }
  
  adj_mat = V %*% diag(lambda_tilde) %*% t(V)
  return(adj_mat)
}

compute_level_sets = function(fit, name){
  #only works for k = 2
  cent = fit$beta_hat
  mat = fit$cov_hat
  if('add_length' %in% names(fit)){
  	addLength=TRUE
  	add_length = fit$add_length
  	# print(add_length)
  }else{
  	addLength=FALSE
  }

  k = length(cent)
  if(k != 2){
  	stop('Only works for k = 2')
  }
  chi_zval = qchisq(0.95, k)
  if(addLength){
  	new_mat = cov.adj(mat, add_length, chi_zval)
  	L = chol(solve(new_mat))	
  	#uncomment below if want to check dimension
  	# cat('Old Dimensions: ', sqrt(chi_zval/diag(solve(mat))), '\n')
  	# cat('New Dimensions: ', sqrt(chi_zval/diag(solve(new_mat))), '\n')
  }else{
  	L = chol(solve(mat))	
  }
  P = solve(L)
  r = sqrt(chi_zval)
  
  
  thetas = seq(0, 2*pi, length.out = 100)
  xs = r*cos(thetas)
  ys = r*sin(thetas)
  circle_points = cbind(xs, ys)
  
  levels = as.data.frame(t(apply(circle_points, 1, function(v) P %*% v + cent)))
  colnames(levels) = c('beta1', 'beta2')
  
  if(is.na(name)){
    return(levels)    
  }else{
    levels['Method'] = name
    return(levels)
  }
  
}

#### Experiment Functions ####

sample_fits = function(n, p, s0, noise_var=1, noise='gaussian',
					   fits = list(isvb=isvb.fit, zz=zz.fit, jm=jm.fit),
					   beta_0_1 = NA, signal_size = NA, signal_scheme = 'normal',
					   feature_correlation = 0, block_size = NA, rescale_first_column = FALSE, AR = FALSE,
					   k=1,
					   dataset='generated', index=1){
	'
	Sample data according to params given in the function. 
	Fit each method supplied in `fits` to X and Y
	fits should be a list of functions which take X and Y and return ... TODO: Add specification of function returns.
	'
	# Generate data according to input parameters
	dat = make_data(n, p, s0, noise_var, noise=noise,
					beta_0_1=beta_0_1, signal_size=signal_size, signal_scheme=signal_scheme,
					feature_correlation=feature_correlation, block_size=block_size, rescale_first_column=rescale_first_column, AR=AR,
					k=k, dataset=dataset, index=index)
	X = dat$X
	Y = dat$Y
	beta_0 = dat$beta_0	

	# Noise estimation and normalisation the same for all appraoches
	if(noise_var != 1){
		mod = cv.glmnet(X, Y)
	  fit = glmnet(X, Y, lambda = mod$lambda.min)
	  yhat = X %*% fit$beta
	  sigma_hat = sqrt((1/(n - s0)) * sum((yhat - Y)**2) )
	  X_scaled = X/sigma_hat
	  Y_scaled = Y/sigma_hat	
	}else{
		X_scaled = X
		Y_scaled = Y
	}
	
	# use anonymous functions to place 'oracle' quantities as necessary
	if('oracle' %in% names(fits)){
		S0 = which(beta_0!=0)
		fits[['oracle']] = function(x, y, k) oracle.fit(x, y, S0=S0, k=k)
	}
	if('jm18' %in% names(fits)){
		Sigma = make_Sigma(p=p, rho=feature_correlation, block_size=block_size)
		fits[['jm18']] = function(x, y, k) jm18.fit(x, y, Sigma, k=k)	
	}

  fitted = lapply(fits, function(fn) fn(X_scaled, Y_scaled, k=k))	
  
	# Compute Losses
	if(k==1){
		fit_maes = data.frame(lapply(fitted, function(fit) MAE(fit, beta_0[1])))
		fit_hits = data.frame(lapply(fitted, function(fit) CI_hit(fit, beta_0[1])))
		fit_lengths = data.frame(lapply(fitted, function(fit) int_length(fit, beta_0[1])))
		fit_times = data.frame(lapply(fitted, function(fit) fit$fit_time))	
		return(list('hit'=fit_hits, 'mae'=fit_maes, 'length'=fit_lengths, 'time' = fit_times))
	}else{
		fit_l2s = data.frame(lapply(fitted, function(fit) L2(fit, beta_0)))
		fit_hits = data.frame(lapply(fitted, function(fit) kHit(fit, beta_0)))
		fit_vols = data.frame(lapply(fitted, function(fit) volume(fit)))
		fit_times = data.frame(lapply(fitted, function(fit) fit$fit_time))	
		return(list('hit'=fit_hits, 'l2'=fit_l2s, 'volume'=fit_vols, 'time' = fit_times))
	}	
}

estimate_stats = function(n, p, s0, beta_0_1, noise_var=1, noise='gaussian',
					   fits = list(isvb=isvb.fit, zz=zz.fit, jm=jm.fit),
					   signal_size = NA, signal_scheme = 'normal',
					   feature_correlation = 0, block_size = NA, rescale_first_column = FALSE, AR = FALSE,
					   n_replicates = 100, mc.cores = 1,
					   k = 1, dataset = 'generated', index = 1, override_core_warning=FALSE){
	'
	Sample data according to params given in the function. 
	Fit each method supplied in fits to X and Y
	fits should be a list of functions which take X and Y and return.
	Returns mean and sd of metrics over n_replicates realisations.
	'
	cat('\n---------------------------',
		'\nRunning experiment with:',
		'\n index: ', index,
		'\nn: ',n,
		'\np: ', p,'\ns0: ',s0,'\nbeta_0_1: ',beta_0_1,
		'\nnoise_var: ',noise_var,'\nnoise: ',noise,
		'\nsignal_size: ',signal_size,'\nsignal_scheme: ',signal_scheme,
		'\nfeature_correlation: ', feature_correlation,'\nn_replicates: ',n_replicates,
		'\ncores: ',mc.cores,
		'\n----------------------------\n')
	if((mc.cores > 8) && !override_core_warning){
		stop('Are you sure you want to use more than 8 cores?
		 If so, call this funciton with override_core_warning=TRUE')
	}
	reps = mc_replicate(n_replicates,
	 sample_fits(n, p, s0, noise_var, noise,
					   fits, beta_0_1, signal_size, signal_scheme,
					   feature_correlation, block_size, rescale_first_column, AR, k, dataset, index),
	 mc.cores=mc.cores
	 )
	print('Finished replicates.')
	if(k == 1){
		metrics = 		c('hit', 'mae', 'length','time')
		compute_sd = 	c(FALSE,  TRUE, TRUE,  TRUE)
	}else{
		metrics = 		c('hit', 'l2', 'volume','time')
		compute_sd = 	c(FALSE,  TRUE, TRUE,  TRUE)
	}
	results=list()
	results_sds=list()
	for(i in c(1:length(metrics))){
		metric=metrics[i]
		met_mean = apply(do.call(rbind, reps[metric,]), 2, mean)
		results = append(results, list(met_mean))
		
		if(compute_sd[i]){
			met_sd = apply(do.call(rbind, reps[metric,]), 2, sd)
			results_sds = append(results_sds, list(met_sd))
		}
	}
	names(results)=metrics
	names(results_sds) = metrics[compute_sd]
	return(list('mean'=results, 'sd'=results_sds))
}

format_results = function(res, round=3, methods=names(fits), k=1){
  mean_df = res$mean
  sd_df = res$sd
  
  if(k == 1){
    metrics = 		c('hit', 'mae', 'length','time')
    metrics_sd = 	c(FALSE,  TRUE, TRUE,  TRUE)
  }else{
    metrics = 		c('hit', 'l2', 'volume','time')
    metrics_sd = 	c(FALSE,  TRUE, TRUE,  TRUE)
  }
  formatted_results = methods
  for(i in c(1:length(metrics))){
    metric = metrics[i]
    metric_sd = metrics_sd[i]
    if(metric == 'volume'){
      # Then normalise
      if('oracle' %in% methods){
      		normalisation = mean_df[[metric]]['oracle']
      		if(normalisation == 0){normalisation = mean_df[[metric]]['isvb']}
      		met_col = paste(format(round(mean_df[[metric]]/normalisation,round), round),
                      format(round(sd_df[[metric]]/normalisation,round), round), sep = ' ± ')		
      	}else{
      		met_col = paste(format(round(mean_df[[metric]]/mean_df[[metric]]['isvb'],round), round),
                      format(round(sd_df[[metric]]/mean_df[[metric]]['isvb'],round), round), sep = ' ± ')		
      	}
    }else if(metric_sd){
      met_col = paste(format(round(mean_df[[metric]],round), round),
                      format(round(sd_df[[metric]],round), round), sep = ' ± ')
    }else{
      met_col = format(round(mean_df[[metric]],round), round)
    }
    formatted_results = cbind(formatted_results, met_col)
  }
  formatted_results = as.data.frame(formatted_results)
  colnames(formatted_results) = c('method', metrics)
  rownames(formatted_results) = c()
  return(formatted_results)
}

format_riboflavin_results = function(results, round=3, methods=names(fits), k=1){
	if(k != 1){
		stop('format_riboflavin_results() only implemented with k = 1.')
	}
  mean_data <- lapply(results, function(x) x$mean)
	sd_data <- lapply(results, function(x) x$sd)
	metrics = c('hit', 'mae', 'length','time')
	metrics_sd = 	c(FALSE,  TRUE, TRUE,  TRUE)
	mean_results = list()
	sd_results = list()
	for(i in c(1:length(metrics))){
	  metric = metrics[i]
	  metric_sd  = metrics_sd[i]
	  metric_list = lapply(mean_data, function(x) x[[metric]])
	  metric_mean = colMeans(do.call(rbind, metric_list))  
	  mean_results[[metric]] = metric_mean
	  if(metric_sd){
	    sd_metric_list = lapply(sd_data, function(x) x[[metric]])
	    sd_metric_mean = colMeans(do.call(rbind, sd_metric_list))  
	    sd_results[[metric]] = sd_metric_mean
	  }
	}

	return(format_results(list('mean'=mean_results, 'sd'=sd_results)))
}
