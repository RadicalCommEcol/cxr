## Resubmission
This is a new release of a package currently on CRAN, in which I add new functionality and fix a couple of bugs. Specifically, in this version I have:

* added and documented the following functions: `build_params`,`densities_to_df`,`fill_demography_matrix`,
`fill_dispersal_matrix`,`fill_transition_matrix`, `generate_vital_rate_coefs`,`get_densities`,`vital_rate`.
* added the vignette `V6_Metapopulation_projections`.
* added automated tests to the new functionality, in the file `test-metapopulations.R`.
* fixed small bugs in the files `cxr_sort_params.R`, `cxr_pm_fit.R`, and `cxr_pm_multifit.R`.
* improved the function `cxr_er_fit.R`.
* removed the dependency on the package `hydroPSO`, which has been discarded from CRAN. This does not alter the core functionality of `cxr`.
* updated the `DESCRIPTION.md`, `NEWS.md`, `README.md` files accordingly.

## Test environments
* local ubuntu 20.04 install, R 4.3.1
* ubuntu 20.04 (on github actions), R 4.3.1
* mac OS (on github actions), R 4.3.1
* windows (on github actions), R 4.3.1
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

## Downstream dependencies

None.
