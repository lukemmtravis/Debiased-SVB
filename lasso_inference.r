#
# Authors: Hamid Javadi, Adel Javanmard, Andrea Montanari, Sven Schmit
# Date: April 24th, 2014
#
# Confidence intervals for high-dimensional regression using the method of
# A. Javanmard, A. Montanari
# "Confidence intervals and hypothesis testing for high-dimensional regression"
# 2013, arXiv:1306.3171
#
library(Matrix);
library(glmnet);
library(expm);
library(flare);

SoftThreshold <- function( x, lambda ) {
#
# Standard soft thresholding
#
  if (x>lambda){
    return (x-lambda);}
  else {
    if (x< (-lambda)){
      return (x+lambda);}
    else {
      return (0); }
  }
}

InverseLinftyOneRow <- function ( sigma, i, mu, maxiter=50, threshold=1e-2 ) {
  p <- nrow(sigma);
  rho <- max(abs(sigma[i,-i])) / sigma[i,i];
  mu0 <- rho/(1+rho);
  beta <- rep(0,p);
    
  if (mu >= mu0){
    beta[i] <- (1-mu0)/sigma[i,i];
    returnlist <- list("optsol" = beta, "iter" = 0);
    return(returnlist);
  }
  
  diff.norm2 <- 1;
  last.norm2 <- 1;
  iter <- 1;
  iter.old <- 1;
  beta[i] <- (1-mu0)/sigma[i,i];
  beta.old <- beta;
  sigma.tilde <- sigma;
  diag(sigma.tilde) <- 0;
  vs <- -sigma.tilde%*%beta;
  
  while ((iter <= maxiter) && (diff.norm2 >= threshold*last.norm2)){    
    
    for (j in 1:p){
      oldval <- beta[j];
      v <- vs[j];
      if (j==i)
        v <- v+1;    
      beta[j] <- SoftThreshold(v,mu)/sigma[j,j];
      if (oldval != beta[j]){
      	vs <- vs + (oldval-beta[j])*sigma.tilde[,j];
      }
    }
    
    iter <- iter + 1;
    if (iter==2*iter.old){
      d <- beta - beta.old;
      diff.norm2 <- sqrt(sum(d*d));
      last.norm2 <-sqrt(sum(beta*beta));
      iter.old <- iter;
      beta.old <- beta;
      if (iter>10)
         vs <- -sigma.tilde%*%beta;
    }
  }

  returnlist <- list("optsol" = beta, "iter" = iter)
  return(returnlist)
}

InverseLinfty <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
  	isgiven <- 0;
  }

  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:p) {
        if ((i %% xp)==0){
          xperc = xperc+10;
          if (verbose) {
            print(paste(xperc,"% done",sep="")); }
        }
  	if (isgiven==0){
  	  mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
  	}
  	mu.stop <- 0;
  	try.no <- 1;
  	incr <- 0;
  	while ((mu.stop != 1)&&(try.no<10)){
  	  last.beta <- beta
  	  output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
  	  beta <- output$optsol
  	  iter <- output$iter
  	  if (isgiven==1){
  	  	mu.stop <- 1
  	  }
  	  else{
            if (try.no==1){
              if (iter == (maxiter+1)){
                incr <- 1;
                mu <- mu*resol;
              } else {
                incr <- 0;
                mu <- mu/resol;
              }
            }
            if (try.no > 1){
              if ((incr == 1)&&(iter == (maxiter+1))){
                mu <- mu*resol;
              }
              if ((incr == 1)&&(iter < (maxiter+1))){
                mu.stop <- 1;
              }
              if ((incr == 0)&&(iter < (maxiter+1))){
                mu <- mu/resol;
              }
              if ((incr == 0)&&(iter == (maxiter+1))){
                mu <- mu*resol;
                beta <- last.beta;
                mu.stop <- 1;
              }                        
            }
          }
  	  try.no <- try.no+1
  	}
  	M[i,] <- beta;
  }
  return(M)
}

InverseLinftyFirstCoord <- function(sigma, n, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }

  p <- nrow(sigma);
  M <- matrix(0, p, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:1) {
        if ((i %% xp)==0){
          xperc = xperc+10;
          if (verbose) {
            print(paste(xperc,"% done",sep="")); }
        }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
            if (try.no==1){
              if (iter == (maxiter+1)){
                incr <- 1;
                mu <- mu*resol;
              } else {
                incr <- 0;
                mu <- mu/resol;
              }
            }
            if (try.no > 1){
              if ((incr == 1)&&(iter == (maxiter+1))){
                mu <- mu*resol;
              }
              if ((incr == 1)&&(iter < (maxiter+1))){
                mu.stop <- 1;
              }
              if ((incr == 0)&&(iter < (maxiter+1))){
                mu <- mu/resol;
              }
              if ((incr == 0)&&(iter == (maxiter+1))){
                mu <- mu*resol;
                beta <- last.beta;
                mu.stop <- 1;
              }                        
            }
          }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M[1,])
}

InverseLinftyKCoords <- function(sigma, n, k, resol=1.5, mu=NULL, maxiter=50, threshold=1e-2, verbose = TRUE) {
  # Compute M[1:k,]
  isgiven <- 1;
  if (is.null(mu)){
    isgiven <- 0;
  }

  p <- nrow(sigma);
  M <- matrix(0, k, p);
  xperc = 0;
  xp = round(p/10);
  for (i in 1:k) {
        if ((i %% xp)==0){
          xperc = xperc+10;
          if (verbose) {
            print(paste(xperc,"% done",sep="")); }
        }
    if (isgiven==0){
      mu <- (1/sqrt(n)) * qnorm(1-(0.1/(p^2)));
    }
    mu.stop <- 0;
    try.no <- 1;
    incr <- 0;
    while ((mu.stop != 1)&&(try.no<10)){
      last.beta <- beta
      output <- InverseLinftyOneRow(sigma, i, mu, maxiter=maxiter, threshold=threshold)
      beta <- output$optsol
      iter <- output$iter
      if (isgiven==1){
        mu.stop <- 1
      }
      else{
            if (try.no==1){
              if (iter == (maxiter+1)){
                incr <- 1;
                mu <- mu*resol;
              } else {
                incr <- 0;
                mu <- mu/resol;
              }
            }
            if (try.no > 1){
              if ((incr == 1)&&(iter == (maxiter+1))){
                mu <- mu*resol;
              }
              if ((incr == 1)&&(iter < (maxiter+1))){
                mu.stop <- 1;
              }
              if ((incr == 0)&&(iter < (maxiter+1))){
                mu <- mu/resol;
              }
              if ((incr == 0)&&(iter == (maxiter+1))){
                mu <- mu*resol;
                beta <- last.beta;
                mu.stop <- 1;
              }                        
            }
          }
      try.no <- try.no+1
    }
    M[i,] <- beta;
  }
  return(M)
}

NoiseSd <- function( yh, A, n ){
  ynorm <- sqrt(n)*(yh/sqrt(diag(A)));
  sd.hat0 <- mad(ynorm);

  zeros <- (abs(ynorm)<3*sd.hat0);
  y2norm <- sum(yh[zeros]^2);
  Atrace <- sum(diag(A)[zeros]);
  sd.hat1 <- sqrt(n*y2norm/Atrace);

  ratio <- sd.hat0/sd.hat1;
  if (max(ratio,1/ratio)>2)
    print("Warning: Noise estimate problematic");

  s0 <- sum(zeros==FALSE);
  return (list( "sd" = sd.hat1, "nz" = s0));
}


Lasso <- function( X, y, lambda = NULL, intercept = TRUE){
#
# Compute the Lasso estimator:
# - If lambda is given, use glmnet and standard Lasso
# - If lambda is not given, use square root Lasso
#
  p <- ncol(X);
  n <- nrow(X);
   
  if  (is.null(lambda)){
    lambda <- sqrt(qnorm(1-(0.1/p))/n);
    outLas <- slim(X,y,lambda=c(lambda),method="lq",q=2,verbose=FALSE);
                                        # Objective : sqrt(RSS/n) +lambda *penalty
    if (intercept==TRUE) {
      return (c(as.vector(outLas$intercept),as.vector(outLas$beta)))
    }  else {
      return (as.vector(outLas$beta));
    }
  } else {
    outLas <- glmnet(X, y, family = c("gaussian"), alpha =1, intercept = intercept );
                                        # Objective :1/2 RSS/n +lambda *penalty
    if (intercept==TRUE){
      return (as.vector(coef(Las,s=lambda)));
    } else {
      return (as.vector(coef(Las,s=lambda))[2:(p+1)]);
    }
  }
}

SSLasso <- function (X, y, alpha=0.05, lambda = NULL, mu = NULL, intercept = TRUE, 
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = TRUE) {
#
# Compute confidence intervals and p-values.
#
# Args:
#   X     :  design matrix
#   y     :  response
#   alpha :  significance level
#   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
#   mu    :  Linfty constraint on M (if null, searches)
#   resol :  step parameter for the function that computes M
#   maxiter: iteration parameter for computing M
#   threshold : tolerance criterion for computing M
#   verbose : verbose?
#
# Returns:
#   noise.sd: Estimate of the noise standard deviation
#   norm0   : Estimate of the number of 'significant' coefficients
#   coef    : Lasso estimated coefficients
#   unb.coef: Unbiased coefficient estimates
#   low.lim : Lower limits of confidence intervals
#   up.lim  : upper limit of confidence intervals
#   pvals   : p-values for the coefficients						 
#
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);

  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);

  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  if ((n>=2*p)){
	 tmp <- eigen(sigma.hat)
	 tmp <- min(tmp$values)/max(tmp$values)
  }else{
	tmp <- 0
  }
  if ((n>=2*p)&&(tmp>=1e-4)){
	  M <- solve(sigma.hat)
  }else{
    M <- InverseLinfty(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  }
  
  unbiased.Lasso <- as.numeric(htheta + (M%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  A <- M %*% sigma.hat %*% t(M);
  noise <- NoiseSd(unbiased.Lasso, A, n );
  s.hat <- noise$sd;
  
  interval.sizes <- qnorm(1-(alpha/2))*s.hat*sqrt(diag(A))/(sqrt(n));

  if  (is.null(lambda)){
    lambda <- s.hat*sqrt(qnorm(1-(0.1/p))/n);
  }
  
  addlength <- rep(0,pp);
  MM <- M%*%sigma.hat - diag(pp);
  for (i in 1:pp){
  	effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
	  effectivemuvec <- effectivemuvec[0:(noise$nz-1)];
	  addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  }  
 
  htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm;
  interval.sizes <- interval.sizes*col.norm;
  addlength <- addlength*col.norm;

  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    interval.sizes <- interval.sizes[2:pp];
    addlength <- addlength[2:pp];
  }  
 p.vals <- 2*(1-pnorm(sqrt(n)*abs(unbiased.Lasso)/(s.hat*col.norm[(pp-p+1):pp]*sqrt(diag(A[(pp-p+1):pp,(pp-p+1):pp])))))
						 
  returnList <- list("noise.sd" = s.hat,
                     "norm0" = noise$nz,
                     "coef" = htheta,
                     "unb.coef" = unbiased.Lasso,
                     "low.lim" = unbiased.Lasso - interval.sizes - addlength,
                     "up.lim" = unbiased.Lasso + interval.sizes + addlength,
					 "pvals" = p.vals,
           "M_est" = A/n
					 )
  return(returnList)
}

SSLassoFirstCoord <- function (X, y, alpha=0.05, lambda = NULL, mu = NULL, intercept = FALSE, 
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = TRUE) {
#
# Compute confidence intervals and p-values.
#
# Args:
#   X     :  design matrix
#   y     :  response
#   alpha :  significance level
#   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
#   mu    :  Linfty constraint on M (if null, searches)
#   resol :  step parameter for the function that computes M
#   maxiter: iteration parameter for computing M
#   threshold : tolerance criterion for computing M
#   verbose : verbose?
#
# Returns:
#   noise.sd: Estimate of the noise standard deviation
#   norm0   : Estimate of the number of 'significant' coefficients
#   coef    : Lasso estimated coefficients
#   unb.coef: Unbiased coefficient estimates
#   low.lim : Lower limits of confidence intervals
#   up.lim  : upper limit of confidence intervals
#   pvals   : p-values for the coefficients            
#
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);

  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);

  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  if ((n>=2*p)){
  tmp <- eigen(sigma.hat)
  tmp <- min(tmp$values)/max(tmp$values)
  }else{
  tmp <- 0
  }
  # if ((n>=2*p)&&(tmp>=1e-4)){
  #   M <- solve(sigma.hat)
  # }else{
  #   M <- InverseLinfty(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  # }
  M1 = InverseLinftyFirstCoord(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose)
  
  unbiased.Lasso <- as.numeric(htheta[1] + (M1%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  
  A <- t(M1) %*% sigma.hat %*% M1;
  # noise <- NoiseSd(unbiased.Lasso, A, n );
  s.hat <- 1 ;
  
  interval.sizes <- qnorm(1-(alpha/2))*s.hat*sqrt(diag(A))/(sqrt(n));

  # if  (is.null(lambda)){
  #   lambda <- s.hat*sqrt(qnorm(1-(0.1/p))/n);
  # }
  
  addlength <- rep(0,1);
  MM <- M1%*%sigma.hat - 1;
  # for (i in 1:1){
  #   effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
  #   effectivemuvec <- effectivemuvec[0:(noise$nz-1)];
  #   addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  # }  
 
  # htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm[1];
  interval.sizes <- interval.sizes*col.norm[1];
  # addlength <- addlength*col.norm;

  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    interval.sizes <- interval.sizes[2:pp];
    addlength <- addlength[2:pp];
  }  
 # p.vals <- 2*(1-pnorm(sqrt(n)*abs(unbiased.Lasso)/(s.hat*col.norm[(pp-p+1):pp]*sqrt(diag(A[(pp-p+1):pp,(pp-p+1):pp])))))
             
  returnList <- list("unb.coef" = unbiased.Lasso,
                     "low.lim" = unbiased.Lasso - interval.sizes,
                     "up.lim" = unbiased.Lasso + interval.sizes
           )
  return(returnList)
}


SSLassok <- function (X, y, k, alpha=0.05, lambda = NULL, mu = NULL, intercept = FALSE, 
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = FALSE) {
#
# Compute confidence intervals and p-values.
#
# Args:
#   X     :  design matrix
#   y     :  response
#   k     :  number of parameters to compute for
#   alpha :  significance level
#   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
#   mu    :  Linfty constraint on M (if null, searches)
#   resol :  step parameter for the function that computes M
#   maxiter: iteration parameter for computing M
#   threshold : tolerance criterion for computing M
#   verbose : verbose?
#
# Returns:
#   noise.sd: Estimate of the noise standard deviation
#   norm0   : Estimate of the number of 'significant' coefficients
#   coef    : Lasso estimated coefficients
#   unb.coef: Unbiased coefficient estimates
#   low.lim : Lower limits of confidence intervals
#   up.lim  : upper limit of confidence intervals
#   pvals   : p-values for the coefficients            
#
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);

  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);

  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  if ((n>=2*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  # if ((n>=2*p)&&(tmp>=1e-4)){
  #   M <- solve(sigma.hat)
  # }else{
  #   M <- InverseLinfty(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  # }
  Mk = InverseLinftyKCoords(sigma.hat, n, k, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose)

  unbiased.Lasso <- as.numeric(htheta[1:k] + (Mk%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  
  A <- Mk %*% sigma.hat %*% t(Mk);
  # noise <- NoiseSd(unbiased.Lasso, A, n );
  s.hat <- 1;
  
  interval.sizes <- qnorm(1-(alpha/2))*s.hat*sqrt(diag(A))/(sqrt(n));

  if  (is.null(lambda)){
    lambda <- s.hat*sqrt(qnorm(1-(0.1/p))/n);
  }
  
  addlength <- rep(0,k);
  MM <- Mk%*%sigma.hat - diag(pp)[1:k,];
  nz = sum(htheta != 0)
  for (i in 1:k){
    effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
    effectivemuvec <- effectivemuvec[0:(nz-1)];
    addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  }  
 
  # htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm[1:k];
  interval.sizes <- interval.sizes*col.norm[1:k];
  addlength <- addlength*col.norm[1:k];

  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    interval.sizes <- interval.sizes[2:pp];
    addlength <- addlength[2:pp];
  }  
             
  returnList <- list("beta_hat" = unbiased.Lasso,
                     "cov_hat" = A/n,
                     'addlength' = addlength
           )
  return(returnList)
}

SSLassoFirstCoordAdjusted <- function (X, y, alpha=0.05, lambda = NULL, mu = NULL, intercept = FALSE, 
                     resol=1.3, maxiter=50, threshold=1e-2, verbose = TRUE) {
#
# Compute confidence intervals and p-values.
#
# Args:
#   X     :  design matrix
#   y     :  response
#   alpha :  significance level
#   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
#   mu    :  Linfty constraint on M (if null, searches)
#   resol :  step parameter for the function that computes M
#   maxiter: iteration parameter for computing M
#   threshold : tolerance criterion for computing M
#   verbose : verbose?
#
# Returns:
#   noise.sd: Estimate of the noise standard deviation
#   norm0   : Estimate of the number of 'significant' coefficients
#   coef    : Lasso estimated coefficients
#   unb.coef: Unbiased coefficient estimates
#   low.lim : Lower limits of confidence intervals
#   up.lim  : upper limit of confidence intervals
#   pvals   : p-values for the coefficients            
#
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);

  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);

  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  if ((n>=2*p)){
    tmp <- eigen(sigma.hat)
    tmp <- min(tmp$values)/max(tmp$values)
  }else{
    tmp <- 0
  }
  # if ((n>=2*p)&&(tmp>=1e-4)){
  #   M1 <- solve(sigma.hat)[1,]
  # }else{
  #   M1 <- InverseLinftyFirstCoord(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose);
  # }
  M1 = InverseLinftyFirstCoord(sigma.hat, n, resol=resol, mu=mu, maxiter=maxiter, threshold=threshold, verbose=verbose)
  
  unbiased.Lasso <- as.numeric(htheta[1] + (M1%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  A <- t(M1) %*% sigma.hat %*% M1;
  # noise <- NoiseSd(unbiased.Lasso, A, n );
  # s.hat <- noise$sd;
  s.hat = 1
  nz = sum(htheta != 0)
  
  interval.sizes <- qnorm(1-(alpha/2))*s.hat*sqrt(diag(A))/(sqrt(n));

  if  (is.null(lambda)){
    lambda <- s.hat*sqrt(qnorm(1-(0.1/p))/n);
  }
  
  addlength <- rep(0,pp);
  MM <- M1 %*% sigma.hat - c(1, rep(0, pp - 1));

  effectivemuvec <- sort(abs(MM),decreasing=TRUE);
  effectivemuvec <- effectivemuvec[0:(nz-1)];
  addlength <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
 
  # print(addlength)
  htheta <- htheta*col.norm[1];
  unbiased.Lasso <- unbiased.Lasso*col.norm[1];
  interval.sizes <- interval.sizes*col.norm[1];
  addlength <- addlength*col.norm[1];

  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    interval.sizes <- interval.sizes[2:pp];
    addlength <- addlength[2:pp];
  }  
 # p.vals <- 2*(1-pnorm(sqrt(n)*abs(unbiased.Lasso)/(s.hat*col.norm[(pp-p+1):pp]*sqrt(diag(A[(pp-p+1):pp,(pp-p+1):pp])))))
             
  returnList <- list("unb.coef" = unbiased.Lasso,
                     "low.lim" = unbiased.Lasso - interval.sizes - addlength,
                     "up.lim" = unbiased.Lasso + interval.sizes + addlength
           )
  return(returnList)
}

SSLassoKnownSigmaK <- function (X, y, Sigma, k, alpha=0.05, lambda = NULL, mu = NULL, intercept = FALSE, 
                     resol=1.3, maxiter=50, threshold=1e-2) {
#
# Compute confidence intervals and p-values.
#
# Args:
#   X     :  design matrix
#   y     :  response
#   alpha :  significance level
#   lambda:  Lasso regularization parameter (if null, fixed by sqrt lasso)
#   mu    :  Linfty constraint on M (if null, searches)
#   resol :  step parameter for the function that computes M
#   maxiter: iteration parameter for computing M
#   threshold : tolerance criterion for computing M
#   verbose : verbose?
#
# Returns:
#   noise.sd: Estimate of the noise standard deviation
#   norm0   : Estimate of the number of 'significant' coefficients
#   coef    : Lasso estimated coefficients
#   unb.coef: Unbiased coefficient estimates
#   low.lim : Lower limits of confidence intervals
#   up.lim  : upper limit of confidence intervals
#   pvals   : p-values for the coefficients            
#
  p <- ncol(X);
  n <- nrow(X);
  pp <- p;
  col.norm <- 1/sqrt((1/n)*diag(t(X)%*%X));
  X <- X %*% diag(col.norm);

  htheta <- Lasso (X,y,lambda=lambda,intercept=intercept);

  if (intercept==TRUE){
    Xb <- cbind(rep(1,n),X);
    col.norm <- c(1,col.norm);
    pp <- (p+1);
  } else {
    Xb <- X;
  }
  sigma.hat <- (1/n)*(t(Xb)%*%Xb);
  
  if ((n>=2*p)){
   tmp <- eigen(sigma.hat)
   tmp <- min(tmp$values)/max(tmp$values)
  }else{
  tmp <- 0
  }
  M = solve(Sigma)
  
  unbiased.Lasso <- as.numeric(htheta + (M%*%t(Xb)%*%(y - Xb %*% htheta))/n);
  A <- M %*% sigma.hat %*% t(M);
  noise <- NoiseSd(unbiased.Lasso, A, n );
  s.hat <- noise$sd;
  
  interval.sizes <- qnorm(1-(alpha/2))*s.hat*sqrt(diag(A))/(sqrt(n));

  if  (is.null(lambda)){
    lambda <- s.hat*sqrt(qnorm(1-(0.1/p))/n);
  }
  
  addlength <- rep(0,pp);
  MM <- M%*%sigma.hat - diag(pp);
  for (i in 1:pp){
    effectivemuvec <- sort(abs(MM[i,]),decreasing=TRUE);
    effectivemuvec <- effectivemuvec[0:(noise$nz-1)];
    addlength[i] <- sqrt(sum(effectivemuvec*effectivemuvec))*lambda;
  }  
 
  htheta <- htheta*col.norm;
  unbiased.Lasso <- unbiased.Lasso*col.norm;
  interval.sizes <- interval.sizes*col.norm;
  addlength <- addlength*col.norm;

  if (intercept==TRUE){
    htheta <- htheta[2:pp];
    unbiased.Lasso <- unbiased.Lasso[2:pp];
    interval.sizes <- interval.sizes[2:pp];
    addlength <- addlength[2:pp];
  }  
  p.vals <- 2*(1-pnorm(sqrt(n)*abs(unbiased.Lasso)/(s.hat*col.norm[(pp-p+1):pp]*sqrt(diag(A[(pp-p+1):pp,(pp-p+1):pp])))))
  if(k > 1){
    returnList <- list("beta_hat" = unbiased.Lasso[1:k],
                     "cov_hat" = A[1:k, 1:k]/n,
                      "add_length" = addlength[1:k]
           )
  }else{
    returnList <- list("beta_hat" = unbiased.Lasso[1],
                      "low.lim" = unbiased.Lasso[1] - interval.sizes[1] - addlength[1],
                      "up.lim" = unbiased.Lasso[1] + interval.sizes[1] + addlength[1]
           )
  }            
  
  return(returnList)
}