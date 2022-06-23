## mEXNEX model for borrowing information between baskets in a basket clinical trial. 

This repository is the code used in the paper: 'Bayesian Information Borrowing Methods in Basket Trials: A Modified Exchangeability-Nonexchangeability Approach' authored by Libby Daniells, Pavel Mozgunov, Alun Bedding and Thomas Jaki.

The files are as follows:
* JAGS Code for each of the borrowing models. 
* Delta-Calibration.R - calibration of the cut-off value for final analysis to control type I error rate under the null at 10%. 
* Simulation.R - computation of operating characteristics for the simulation study within the paper. 
* c_e-analysis.R - pre-trial simulation study to pick a cut-off, c_e, for removal of 'extreme' baskets in the mEXNEX method. 
