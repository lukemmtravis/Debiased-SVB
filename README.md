## Debiased SVB ##
This repo contains functions to fit the methods explored in our paper.

* If one has not used the `sparsevb` package before, install it using the following command `install.packages('sparsevb_0.1.0.tar.gz', repos = NULL, type="source")`.

* The main bulk of the functions are in `SVB.R`, and can be fitted to an $n \times p$ matrix $X$ accompanied by an $n-$ vector Y. For example,

  * `isvb.fit(X,Y)` returns an estimate for $\beta_1$, associated 95% confidence interval and fit time.
  * `isvb.fit(X,Y, k=2)` returns an estimate for $\beta_{1:2}$, associated covariance matrix for a 2-dimensional credible interval and fit time.

* `example.R` contains an example of how to make simulated data and fit various methods to it.

* `experiment.R` contains the code to run our 1D experiments. Running this with `n_replicates` set to 500 will produce Table 1 of our paper.

* `experiment_k.R` contains the code to run our higher dimensional experiments. Running this with `n_replicates` set to 500 will produce Table 2 of our paper.

* `riboflavin_experiment.R` contains the code to run our riboflavin experiments. Running this with `n_replicates` set to 100 will produce Table 3 of our paper.

* `set_plots.R` contains an example of how to plot the 2-dimensional credible regions.


