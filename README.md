# The 'LCTMC' package

This R package allows users to fit a **L**atent class **C**ontinuous-**T**ime **M**arkov **C**hain model to longitudinal categorical data. 

## Overview

This README is a work in progress at the moment. 

For now, please see the ['LCTMC.simulate'](https://github.com/j-kuo/LCTMC.simulate) repository for more info on the latent class CTMC model.

## Work in Progress

1. Can model run smoothly for all 3 scenarios?
   * 2x2 model (only "X" covariate)
   * 3x3 model (only "X" covariate)
   * 3x3 model ("X" and "W" covariate)
   * Try for 10 sim, and how about 50?
1. ~~Add Examples~~
   * ~~lctmc, init01, init02, EM, SE, rescale bundled together (use \dontrun{} and code to generate example data, `LCTMC.simulate`)~~
   * ~~add output of model fit on example data~~
   * ~~use fitted model to create example for compute functions~~
1. Miscelleneous
   * add badge
   * ~~license~~
   * is there an improvement from by adding the "W" covariates when there are only few observations

## Authors

* **Jacky Kuo** - _author_, _maintainer_
* **Wenyaw Chan**, PhD - _dissertation advisor_