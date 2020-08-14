[![](https://travis-ci.com/js2264/VplotR.svg?branch=master)](https://travis-ci.com/js2264/VplotR)
[![](https://codecov.io/gh/js2264/VplotR/branch/master/graph/badge.svg)](https://codecov.io/github/js2264/VplotR?branch=master)
[![](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](https://img.shields.io/github/languages/code-size/js2264/VplotR.svg)](https://github.com/js2264/VplotR)
[![](https://img.shields.io/badge/license-GPL--3-orange.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# VplotR <img src="inst/figures/logo.png" align="right" alt="" />

## Introduction

This R package makes the process of generating fragment density plots 
(also known as "V-plots") straightforward. V-plots have been introduced 
[for the first time](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215028/) 
by the Henikoff lab in 2011. Recently, V-plots have proven to be very 
instructive to understand the molecular organization of the chromatin. 
For instance, the 
[nucleoATAC](https://genome.cshlp.org/content/early/2015/08/27/gr.192294.115)
package relies on cross-correlation of ATAC-seq fragment density plots to 
accurately map nucleosome occupancy along the genome.  
VplotR aim is to streamline the process of generating V-plots. It contains 
wrapping functions to import paired-end sequencing bam files and generate 
V-plots around genomic loci of interest.  
VplotR is designed around [ggplot2](https://ggplot2.tidyverse.org/) and 
makes full use of its potential. As such, it is easy to generate V-plots 
in batches and combine them with other plots to make 
publication-ready figures.  

## Installation

VplotR can be installed from Github as follows:

```r
install.packages("devtools")
devtools::install_github("js2264/VplotR")
```

## Main functions 

The main user-level functions of VplotR are `getFragmentsDistribution()`, 
`plotVmat()`, `plotFootprint()` and `plotProfile()`.

* `getFragmentsDistribution()` computes the distribution of fragment sizes
  over sets of genomic ranges;
* `plotVmat()` is used to compute fragment density and generate V-plots;
* `plotFootprint()` generates the MNase-seq or ATAC-seq footprint at a 
  set of genomic ranges.
* `plotProfile()` is used to plot the distribution of paired-end fragments 
  at a single locus of interest.

See the package vignettes for full description of each function 
and examples of use cases.

## Contributions
Code contributions, bug reports, fixes and feature requests are most welcome.
Please make any pull requests against the master branch at 
https://github.com/js2264/VplotR
and file issues at https://github.com/js2264/VplotR/issues

## License 
**VplotR** is licensed under the GPL-3 license.
