# The 'LCTMC' package

This R package allows users to fit a **L**atent class **C**ontinuous-**T**ime **M**arkov **C**hain model to longitudinal categorical data. 

## Overview

This README is a work in progress at the moment. 

For now, please see the ['LCTMC.simulate'](https://github.com/j-kuo/LCTMC.simulate) repository for more info on the latent class CTMC model.

## Work in Progress

1. Update on prediction function 
   * ~~is the function working as intended?~~ (yes)
   * ~~model seems to over predict disease state 1 when actual state is "death"~~ (it was issue with input data, not function itself)
   * ~~work out probability by hand and compare to those from prediction function~~ (done)
   * is there an improvement from by adding the "W" covariates when there are only few observations?
1. Can model run smoothly for all 3 scenarios?
   * 2x2 model (only "X" covariate)
   * 3x3 model (only "X" covariate)
   * 3x3 model ("X" and "W" covariate)
   * Try for 10 sim, and how about 50?
1. Add function to quickly compute various statistics from the model
   * ~~**Q** matrix~~
   * ~~**P** matrix~~
   * Transition Rate ratio
   * ~~Sojourn Time~~
   * Other stats for between cluster comparisons?

## Authors

* **Jacky Kuo** - _author_, _maintainer_
* **Wenyaw Chan**, PhD - _dissertation advisor_