## Debiased SVB ##
This repo contains functions to fit the methods explored in our paper.

All of the functions are in `SVB.R`, and can be fitted to an $n \times p$ matrix $X$ accompanied by an $n-$vector Y.

E.g. `isvb.fit(X,Y)` returns an estimate for $\beta_1$, associated $95\%-$confidence interval and fit time.
