# The 'LCTMC' package

This R package allows users to fit a **L**atent class **C**ontinuous-**T**ime **M**arkov **C**hain model to longitudinal categorical data. 

## Overview

This README is a work in progress at the moment. 

For now, please see the ['LCTMC.simulate'](https://github.com/j-kuo/LCTMC.simulate) repository for more info on the latent class CTMC model.

## Work in Progress

1. Update on `@seealso` so that functions can reference each other in the documentation
2. Update on prediction function 
   * is the function working as intended?
   * model seems to over predict disease state 1 when actual state is "death"
   * work out probability by hand and compare to those from prediction function
   * is there an improvement from by add the "W" covariates when there are only few observations?
3. Allow model to alternative between init 02 & EM until some criteria is met
4. Can model run smoothly for all 3 scenarios?
   * 2x2 model (only "X" covariate)
   * 3x3 model (only "X" covariate)
   * 3x3 model ("X" and "W" covariate)
   * Try for 10 sim, and how about 50?
5. Add function to quickly compute various statistics from the model
   * **Q** matrix
   * **P** matrix
   * Hazard ratios
   * Sojourn Time
   * Other stats for between cluster comparisons?