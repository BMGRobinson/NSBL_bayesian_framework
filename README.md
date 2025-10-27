# Project Description
This repository contains a general Bayesian computational framework for the concurrent estimation of time-varying and time-invariant parameters in nonlinear stochastic dynamical systsems.

The code used to generate the results presented in the manuscript *Semi-analytical hierarchical Bayesian inference of nonlinear
model structure in stochastic dynamics: Applied to compartmental models of infectious diseases* can be found in the `/framework/Examples` directory.

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [License](#license)


# Installation
The code is tested on Mint 18.2, Ubuntu 22.04, and various machines in the Digital Reseach Alliance of Canada consortium (currently for the machine 'Graham'). README files and scripts relevant to the compilation of the framework on these systems can be found at:  
`/framework/readme_mint182`  
`/framework/readme_ubuntu2204`  
`/framework/readme_graham`  

# Usage

Existing results can be plotted, and cases can be modified or re-run from the following subdirectories:  
`/framework/Examples` 

Within these directories, the parameter estimation and state estimation results can be plotted by running either of the following ipython notebooks:  
`/framework/Examples/Case_PLOS/PlotPDF.ipynb`  
`/framework/Examples/Case_PLOS/PlotSE.ipynb`  

State space models can be altered by modifying the file. Existing state spaces can be edited, or new ones may be defined following the template.  
`/framework/Examples/functions.cpp`  

Parameter estimation is performed by toggling "run = True/False" for the 4 models:  
`/framework/Examples/caseSIRS.cfg`

The MCMC samples (prior pdf and posterior pdf) and the MAP estimate are stored in the directory  
`/framework/Examples/Case_PLOS/chains`

Following the parameter estimation, state estimation is performed using  
`/framework/Examples/caseSIRSse.cfg`

This will by default perform a single instance of state estimation, performed at the MAP estimate of the paramters, you can alternatively run state estimation without performing parameter estimation, by defning an input vector of time-invariant parameters in the .cfg file.

NSBL is run after performing MCMC using 
`/framework/Examples/Case_PLOS/NSBL/NSBL.ipynb`

