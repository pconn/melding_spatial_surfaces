# melding_spatial_surfaces

Summary
=============

This repository includes code to estimate a single relative abundance map when multiple
relative abundance maps (and associated covariance matrices) are available.  It is intended
to accompany the paper

"A GLMM approach for combining multiple relative abundance surfaces"
by Paul B. Conn, Jay M. Ver Hoef, Devin S. Johnson, Brett T. McClintock, and Brian Brost
(currently submitted to Methods in Ecology and Evolution)

Our implementation relies on the TMB R package to calculate likelihoods, and INLA to set up
spatial autocorrelation structures (see below for installation instructions for these packages).
Likelihoods were programmed in templated C++ files (located in the "src"" directory), while scripts used to execute particular analyses are located in the "inst" directory, and Steller sea lion data are located in the "data" directory. 

Data
=====

Steller sea lion estimated telemetry location data are provided in "./data_integration/data/SSL_crawl.RData"
This workspace includes one object, "clipcrwData" which is standard
output from fitting the *crawl* R package to Steller Sea Lion location data.  The object includes interpolated Steller sea lion positions at regular 20 minute intervals for 36 telemetered sea lions.

Steller sea lion habitat data are provided in "./data_integration/data/ssl_habitat.RData."  This R workspace contains
two object. "hab_cov" is a rasterBrick including ocean depth (depth; in meters), distance from nearest rookery or haul-out location (dist_site; in meters), distance from land (dist_land; in meters), and sea floor slope (slope).  "Grid_sf" is an sf polygons object containing the analysis grid reported in the paper.

A log-scale Platforms-of-Opportunity (POP) relative abundance surface ("logPOP"), together with a diagonalized version of the inverse variance-covariance matrix ("VCinv_logPOP") are provided in "./data_integration/data/POPests.RData".  The order of these estimates matches the cells in our analysis grid.  Note that raw POP data cannot be provided because of privacy concerns related to fishing vessel locations.


R scripts (inst directory)
=========

*fitLangevin.R*:  This script fits a Langevin diffusion model to interpolated Steller sea lion locations to produce a utilization distribution (UD) that relies on habitat variables.  The output of this script is *fitLangevin.RData*

*prep_SSL_inputs_rr.R*: This script prepares the UD and POP surfaces for our 
GLMM analysis; the output is *SSL_TMB_data.RData*

*integrate_SSL_red_rank.R*: This script calls TMB to combine the UD and POP SSL surfaces, using two different models (one with process errors estimated, and one with them turned off). It outputs two R workspaces with results ("output_rr_bias.RData" and "output_rr_nobias.RData" respectively)

*plot_SSL_ests.R* This script produces SSL relative abundance maps from the previous script's output

*run_TMB_sims.R*  This script runs our simulation study.  

*calc_plot_sim_results.R*  This script makes summary plots from simulation results.

*run_TMB_misalign.Rmd*  This R markdown file includes code to replicate Online Appendix S1 - where simulated data occurred on a variety of different spatial resolutions.


TMB source code (src directory)
===============================

*fit_1_surface.cpp*  This templated .cpp file contains TMB code to
calculate the negative log likelihood for a simple spatial model fitted to  count data from a single survey.  This code is
called by *run_TMB_sims.R* and *run_TMB_misalign.Rmd* to produce estimates of individual surfaces and associated covariance matrices.   

*fit_multiple_surfaces.cpp*  This templated .cpp file contains TMB code to
calculate the negative log likelihood for our integrated model when the 
spatial resolution of each fitted surface is the same.  This was used in the primary simulation study reported in the manuscript, and is called in 
*run_TMB_sims.R*.   

*fit_multiple_surfaces_misalign.cpp*  This templated .cpp file contains TMB code to calculate the negative log likelihood for our integrated model when the spatial resolution of fitted and estimated surfaces differs.  It is used in Online Appendix S1 and is called from *run_TMB_misalign.Rmd*.

*fit_multiple_surfaces_rr.cpp*  This templated .cpp file contains TMB code to calculate the negative log likelihood for a version of our integrated model where a reduced rank model is included for log-relative abundance (i.e. $\boldsymbol{\mu}$) and process errors (i.e. $\boldsymbol{\xi}_i$).  This version is useful when the number of pixels is large (e.g., $>10,000$) which can make log-likelihood calculations too computationally demanding to be practical.  It was used in our Steller sea lion example, being called from *integrate_SSL_red_rank.R*.



Installation Instructions for TMB and INLA
==========================================
The TMB/INLA implementation depends on R version >=4.1.0 and a variety of other tools.

First, install the "devtools" package from CRAN

    # Install and load devtools package
    install.packages("devtools")
    library("devtools")

Second, please install the following:
* TMB (Template Model Builder): https://github.com/kaskr/adcomp
* INLA (integrated nested Laplace approximations): http://www.r-inla.org/download

Note: at the moment, TMB and INLA can be installed using the commands 

    # devtools command to get TMB from GitHub
    install_github("kaskr/adcomp/TMB") 
    # source script to get INLA from the web
    source("http://www.math.ntnu.no/inla/givemeINLA.R")  

