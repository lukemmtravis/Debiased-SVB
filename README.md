## Debiased SVB ##
This repo contains functions to fit the methods explored in our paper.

* The main bulk of the functions are in `SVB.R`, and can be fitted to an $n \times p$ matrix $X$ accompanied by an $n-$ vector Y. For example,

  * `isvb.fit(X,Y)` returns an estimate for $\beta_1$, associated 95% confidence interval and fit time.
  * `isvb.fit(X,Y, k=2)` returns an estimate for $\beta_{1:2}$, associated covariance matrix for a 2-dimensional credible interval and fit time.

* `experiment.R` contains an example of how to run our experiments.

* `set_plots.R` contains an example of how to plot the 2-dimensional credible regions.
