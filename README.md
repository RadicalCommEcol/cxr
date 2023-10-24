
### Current Status

[![R-CMD-check](https://github.com/RadicalCommEcol/cxr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RadicalCommEcol/cxr/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/cxr)](https://www.r-pkg.org/badges/version/cxr)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/cxr)](https://cran.r-project.org/package=cxr)
[![DOI](https://zenodo.org/badge/115796966.svg)](https://zenodo.org/badge/latestdoi/115796966)

# cxr 1.1.1

cxr provides a complete toolbox for modelling interactions between species, calculate coexistence metrics (e.g. niche and fitness differences), and project species abundances. The functions are flexible and can include covariates, use different optimization algorithms, or accept user-defined mathematical population models as starting points. It furthermore includes a series of functions to model metapopulation dynamics of stage-structured populations.

### Installation

The package can be installed from CRAN:

```R
install.packages("cxr")
```

Furthermore, the development version can be installed via the `remotes` package:

```R
remotes::install_github("RadicalCommEcol/cxr")
library("cxr")
```

### Setup requirements and use

The best way to start is to follow our [vignettes](https://github.com/RadicalCommEcol/cxr/tree/master/vignettes).
Note that for installing the vignettes alongside the package, you need to use following options at install:

```R
remotes::install_github("RadicalCommEcol/cxr", 
                         build_opts = c("--no-resave-data", "--no-manual"), 
                         build_vignettes = TRUE)
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

And a set of internal functions, data and models. The functions with prefix `cxr` estimate model parameters from observational data, and functions without prefix calculate metrics related to species coexistence from model parameters. This basic workflow of the package is explained in different vignettes:

- Vignette 1 "Getting started" shows how to fit observational data to population dynamics models, in order to obtain vital rates and interaction coefficients.
- Vignette 2 "Data formats" explains the data structures accepted by these fitting functions, using as an example the dataset included with the package.
- Vignette 3 "Coexistence metrics" details the complete workflow, from estimating model parameters to calculate the coexistence metrics available.
- Vignette 4 "Using your own models" explains how to extend the package using user-defined population models.
- Vignette 5 "Projecting species abundance" shows how to use model fits to project species abundances in time.

The complementary functionality on metapopulation dynamics is described in another vignette:

- Vignette 6 "Metapopulation projections", that showcases the functions and methodologies for modelling metapopulations of stage-structured species across a number of sites, potentially including environmental effects on their vital rates.

Once the package is installed, vignettes can be accessed in the standard way:

```R
vignette("V1_Getting_started",package = "cxr")
vignette("V2_Data_formats",package = "cxr")
vignette("V3_Coexistence_metrics",package = "cxr")
vignette("V4_Models",package = "cxr")
vignette("V5_Abundance_projections",package = "cxr")
vignette("V6_Metapopulation_projections",package = "cxr")
```

### Citation information

When citing, please refer to both the [package citation](https://github.com/RadicalCommEcol/cxr/blob/master/inst/CITATION) and the release paper [García-Callejas, Godoy, and Bartomeus, 2020](https://doi.org/10.1111/2041-210X.13443).  

### Future developments

`cxr` is in continuous development. This is a partial list of features we aim to implement in future releases:

- integrate optimizer arguments, importantly for allowing scale normalization (parscale).
- uncertainty estimation: propagate standard error calculation to coexistence metrics and abundance projections.
- diagnostic functions: guidance on whether model fits/coexistence metrics are meaningful for a certain data and model.
- update the interface to allow users to pass other arguments to their models, 
e.g. to implement properly the annual plant model with 'g' and 's' constants.
- package design and style guide: provide a complete rationale and set of recommendations to contribute new features.
- package website with pkgdown.
- continue the integration of stage-specific demography, species interactions, and ecological forecasts (in collaboration with Maria Paniw).

## Bug reports and contributions.  

We welcome contributions (e.g. via pull request) and [bug reports](https://github.com/RadicalCommEcol/cxr/issues).

## Code of Conduct  

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/RadicalCommEcol/cxr/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

