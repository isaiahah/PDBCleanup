# PDBCleanup

<!-- badges: start -->
<!-- badges: end -->

## Description
PDBCleanup is an R package developed as part of the University of Toronto Fall 2023 course BCB410. It provides utilities for improving the quality of PDB protein structure files by selecting high resolution regions and improving relative domain positions by superimposing individual regions against corresponding regions of another reference protein structure. It further allows visualization of a quality score along an input protein's sequence, visualization of the 3D protein structure produced after selecting high quality regions and performing superimposition, and visualization of the structural difference between a pair of protein structures.

This package is beneficial as a pre-processing step before performing analysis using protein structural data. It is compatible with other common structural biology packages but provides new highly customizable functions to streamline the first steps of structural data preparation.

PDBCleanup was developed in an environment using `R version 4.3.1` on an `aarch64-apple-darwin20 (64-bit)` platform running `macOS Sonoma 14.0`.

## Installation

You can install the development version of PDBCleanup from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library("devtools")
devtools::install_github("isaiahah/PDBCleanup", buildVignettes = TRUE)
library("PDBCleanup")
```

## Overview
```
ls("package:PDBCleanup")
data(package = "PDBCleanup") 
browseVignettes("PDBCleanup")
```

There are X functions in `PDBCleanup`:

## Contributions
This package was created by Isaiah Hazelwood.

## Acknowledgements
This package was developed as part of an assessment for 2023 BCB410H: Applied Bioinformatics course at the University of Toronto, Toronto, Canada. PDBCleanup welcomes issues, enhancement requests, and other contributions. To submit an issue, use the GitHub issues.


## References
Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.

Jumper, J., Evans, R., Pritzel, A. et al. Highly accurate protein structure prediction with AlphaFold. Nature 596, 583–589 (2021). https://doi.org/10.1038/s41586-021-03819-2

Zhoutong Sun, Qian Liu, Ge Qu, Yan Feng, and Manfred T. Reetz. Utility of B-Factors in Protein Science: Interpreting Rigidity, Flexibility, and Internal Motion and Engineering Thermostability. Chemical Reviews 2019 119 (3), 1626-1665. https://doi.org/10.1021/acs.chemrev.8b00290
