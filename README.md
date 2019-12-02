
### Current Status

[![build status](https://travis-ci.org/ibartomeus/cxr.svg?branch=master)](https://travis-ci.org/ibartomeus/cxr)
[![CRAN status](https://www.r-pkg.org/badges/version/cxr)](https://www.r-pkg.org/badges/version/cxr)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/cxr)](https://cran.r-project.org/package=cxr)
[![DOI](https://zenodo.org/badge/115796966.svg)](https://zenodo.org/badge/latestdoi/115796966)


# cxr 0.1.0

CXR provides a complete toolbox for modelling competitive effects between species, calculate fitness and niche differences, and calculate and predict coexistence regions. The functions are flexible and can include covariates, use different optimization algorithms, or accept user-defined mathematical population models as starting points. 

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

- `pm_optim()`
- `er_optim()`
- `NicheOverlap()` 
- `SpeciesFitness()`
- `AvgFitnessRatio()` 
- `PredictAbundances()` 
- `GenerateTestData()`

And a set of internal functions, data and models.

### Prefix information

The package functions have three different prefixes:

- `model_`: Population models.
- `pm_`: Population model parametrization.
- `er_`: Effect-Response model parametrization.
- `cxr_`: Internal functions.

Other functions have no prefix.

### Citation information

When citing, please refer to both the [package citation](https://github.com/ibartomeus/cxr/blob/master/inst/CITATION) and the release paper (in prep).  

## Bug reports and contributions.  

We expect to update the package with new developments. We welcome contributions (e.g. via pull request) and [bug reports](https://github.com/ibartomeus/cxr/issues).

## Code of Conduct  

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ibartomeus/cxr/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

