# The 'LCTMC' package

<!-- badges: start -->
[![R-CMD-check](https://github.com/j-kuo/LCTMC/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/j-kuo/LCTMC/actions/workflows/R-CMD-check.yaml)
[![](https://img.shields.io/badge/R%20version-4.2.1-steelblue.svg)](https://cran.r-project.org/bin/windows/base/old/4.2.1)
<!-- badges: end -->

This R package allows users to fit a **L**atent class **C**ontinuous-**T**ime **M**arkov **C**hain model to longitudinal categorical data. 

## Overview

This package is a work in progress at the moment. 

For now, please see the ['LCTMC.simulate'](https://github.com/j-kuo/LCTMC.simulate) repository for more info on the latent class CTMC model.

## Work in Progress

1. Can model run smoothly for all 3 scenarios?
   1. 2x2 model (only "X" covariate)
   1. 3x3 model (only "X" covariate)
   1. 3x3 model ("X" and "W" covariate)
1. Try for 10 sim, and how about 50?
   * How's the coverage & bias estimate
   * Is there an improvement from by adding the "W" covariates when there are only few observations

## Authors

* **Jacky Kuo** - _author_, _maintainer_
* **Wenyaw Chan**, PhD - _dissertation advisor_