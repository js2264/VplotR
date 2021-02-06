[![](https://img.shields.io/badge/release%20version-1.1.0-orange.svg)](https://www.bioconductor.org/packages/VplotR)
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

## System requirements 

VplotR has been tested on Mac OS (`>= 10.13`), Ubuntu (`>= 18.04`) and Windows Server (`2012`).  
VplotR requires `R >= 4.0` and is available from Bioconductor (release `> 3.12`). 

## Installation

VplotR and all its dependencies can be installed from Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VplotR")
library("VplotR")
```

VplotR developpment version can be installed from Github as follows:

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

See the [package vignette](https://jserizay.com/VplotR/articles/VplotR.html) 
for an in-depth introduction to the package functionalities and 
the [reference](https://jserizay.com/VplotR/reference/index.html) 
for a full description of each function, a list of the available data 
provided in the package and examples of use cases.

## Get started with VplotR 

### Plotting a V-plot

V-plots can be generated using the `plotVmat()` function: 

```r
data(MNase_sacCer3_Henikoff2011)
data(ABF1_sacCer3)
p <- plotVmat(
    x = MNase_sacCer3_Henikoff2011, 
    granges = ABF1_sacCer3
)
```

### Importing sequencing fragments are `GRanges`

Paired-end `.bam` files are imported using the `importPEBamFiles()` 
function as follows:

```r
library(VplotR)
bamfile <- system.file("extdata", "ex1.bam", package = "Rsamtools")
fragments <- importPEBamFiles(
    bamfile,
    shift_ATAC_fragments = TRUE
)
```

### Plotting fragment size distribution 

The distribution of fragment sizes can be computed with 
the `getFragmentsDistribution()` function: 

```r
library(VplotR)
data(MNase_sacCer3_Henikoff2011)
data(ABF1_sacCer3)
df <- getFragmentsDistribution(
    MNase_sacCer3_Henikoff2011,
    ABF1_sacCer3
)
ggplot(df, aes(x = x, y = y)) + 
    geom_line() + 
    theme_ggplot2()
```

### Datasets provided in VplotR

Several dataset samples are provided as `GRanges` 
object in this package:

- A set of tissue-specific ATAC-seq experiments in
young adult _C. elegans_ (Serizay et al., 2020):

```r
library(VplotR)
data(ATAC_ce11_Serizay2020)
```

- An MNase-seq experiment in yeast (Henikoff et al., 2011) 
and associated ABF1 binding sites:

```r
library(VplotR)
data(MNase_sacCer3_Henikoff2011)
```

- A `.bam` file sample containing 100,000 reads mapped to _C. elegans_:

```r
library(VplotR)
data(bam_test)
```

For a full list of the data available in this package, please go to the 
[reference page](https://jserizay.com/VplotR/reference/index.html#section-data). 

## Contributions

Code contributions, bug reports, fixes and feature requests are most welcome.
Please make any pull requests against the master branch at 
https://github.com/js2264/VplotR
and file issues at https://github.com/js2264/VplotR/issues.  
Feel free to reach out to [J. Serizay](mailto:jacquesserizay@gmail.com) for any 
query. 

## License 

**VplotR** is licensed under the GPL-3 license.

## Code of Conduct

Please note that `VplotR` is released with a 
(Contributor Code of Conduct)[https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html].  
By contributing to this project, you agree to abide by its terms.
