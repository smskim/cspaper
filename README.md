This repository contains an R package with all code necessary to reproduce
simulation results and all plots from

Howard, S. R., Ramdas, A., McAuliffe, J, and Sekhon, J. (2018), [Time-uniform,
nonparametric, nonasymptotic confidence
sequences](https://arxiv.org/abs/1810.08240), preprint, arXiv:1810.08240.

For implementations of the uniform boundaries themselves, look at the [confseq
package](https://github.com/gostevehoward/confseq).

You can install the package like so:

```R
install.packages('devtools')
devtools::install_github('gostevehoward/cspaper')
```

Building the simulation code requires a C++ compiler with C++14 support. You
will need to have the `BH` and `confseq` R packages installed, since we link
against C++ headers from those libraries.

To run simulations for Figure 1 and generate the plots:
```R
intro_simulations(save=TRUE)
plot_intro(save=TRUE)
```

To run simulations for Figure 7 and generate the plots:
```R
bounded_simulations(save=TRUE)
save_all_plots()
```

To generate the rest of the plots:
```R
# Figure 3:
poly_stitching_plot(save=TRUE)
# Figure 4:
plot_normalized_boundaries(save=TRUE)
# Figure 5:
bernoulli_ate_plot(save=TRUE)
# Figure 6:
save_ellipse_plots()
# Figure S1:
discrete_mixture_plot(save=TRUE)
# Figure S2:
make_finite_lil_plot_simplified(save=TRUE)
# Figure S3:
make_finite_lil_plot(save=TRUE)
```
