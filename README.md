
### Current Status

[![build status](https://travis-ci.org/ibartomeus/cxr.svg?branch=master)](https://travis-ci.org/ibartomeus/cxr)
[![CRAN status](https://www.r-pkg.org/badges/version/cxr)](https://www.r-pkg.org/badges/version/cxr)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/cxr)](https://cran.r-project.org/package=cxr)
[![DOI](https://zenodo.org/badge/115796966.svg)](https://zenodo.org/badge/latestdoi/115796966)

# cxr 0.1.0

CXR provides a complete toolbox for modelling interactions between species, calculate fitness and niche differences, and project species abundances. The functions are flexible and can include covariates, use different optimization algorithms, or accept user-defined mathematical population models as starting points. 

### Installation

The package can't be installed from CRAN yet, but when available will be at:

```R
install.packages("cxr")
library("cxr")
```

The development version can be installed via the `devtools` package:

```R
devtools::install_github("ibartomeus/cxr")
library("cxr")
```

### Setup requirements and use

The best way to start is to follow our [vignettes](https://github.com/ibartomeus/cxr/tree/master/vignettes).
Note that for installing the vignettes alongside the package, you may use following options at install:

```R
devtools::install_github("ibartomeus/cxr", build_opts = c("--no-resave-data", "--no-manual"))
```

### Key features

The package has several key functions:

- `cxr_pm_fit()`
- `cxr_pm_multifit()`
- `cxr_er_fit()` 
- `niche_overlap()`
- `species_fitness()`
- `avg_fitness_diff()` 
- `competitive_ability()` 
- `abundance_projection()`

And a set of internal functions, data and models.

### Prefix information

The functions used to numerically fit model parameters have the prefix `cxr_`. The acronym `pm` means 'population models', and `er` means 'effect and response', for specifying the type of models these functions parameterize. Other functions to calculate coexistence metrics have no prefix. We have implemented four model families to use with the package. These are 'Beverton-Holt' (`BH`), 'Lotka-Volterra' (`LV`), 'Law-Watkinson' (`LW`), and 'Ricker' (`RK`). Models from each family are named starting with the model family prefix, e.g. 'BH_', etc.

### Citation information

When citing, please refer to both the [package citation](https://github.com/ibartomeus/cxr/blob/master/inst/CITATION) and the release paper (in prep).  

## Bug reports and contributions.  

We expect to update the package with new developments. We welcome contributions (e.g. via pull request) and [bug reports](https://github.com/ibartomeus/cxr/issues).

## Code of Conduct  

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ibartomeus/cxr/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

