# VplotR

**VplotR is still in Alpha and has not been thorouhly tested yet.**

![VplotR](examples/png/Comparison_tissue-specific-normalized-Vmats.png)

## Introduction

This R package makes the process of generating "V-plots" straighforward. 
V-plots have been introduced 
[for the first time](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215028/) 
by the Henikoff lab in 2011.  

Recently, V-plots have proven to be very instructive to understand the molecular 
organization of the chromatin. For instance, the 
[nucleoATAC]((https://genome.cshlp.org/content/early/2015/08/27/gr.192294.115)) 
package relies on cross-correlation of ATAC-seq-derived V-plots to accurately map 
nucleosome occupancy along the genome.

VplotR aim is to streamline the process of generating V-plots. 
It contains wrapping functions to import paired-end sequencing bam files and 
generate V-plots around genomic loci of interest.

VplotR is designed around [ggplot2](https://ggplot2.tidyverse.org/) 
and makes full use of its potential. As such, it is easy to generate V-plots 
in batches and combine them with other plots to make publication-ready figures.

## Installation

VplotR can be ran installed from Github as follow:

```r
install.packages("devtools")
devtools::install_github("js2264/VplotR")
library(VplotR)
```

## Overview

Firstly BAM files are read using the `importPEBamFiles()` function and loci of
interest from a BED file, for instance.

```r
bam_files <- list.files(
    path = 'path/to/bam/files/', 
    pattern = '*.bam', 
    full.names = TRUE
)
granges <- rtracklayer::import('loci_of_interest.bed')
bam_list <- importPEBamFiles(
    bam_files, 
    where = GenomicRanges::resize(granges, width = 2000, fix = 'center'), 
    shift_ATAC_fragments = TRUE
)
```

*Note: to allow for a background normalization, the `where` argument should be 
omitted.*

Then V-plots of the bam files over the set of loci of interest (`granges`) 
are generated using the `plotVmat()` function:

```r
Vplot <- plotVmat(bam_list, granges)
```

The generation of multiple V-plots can be parallelized as follow:

```r
list_params <- list(
    list("bam" = bam1, "granges" = granges1), 
    list("bam" = bam2, "granges" = granges2), 
    ..., 
    list("bam" = bamN, "granges" = grangesN)
)
plots <- plotVmat(
    list_params, 
    cores = length(list_params)
) + facet_wrap(~Cond.)
```

Finally, the `nucleosomeEnrichment()` function is useful to statistically quantify 
and compare nucleosome enrichment (e.g. flanking nucleosome at promoters). To do so:

```r
nucleosomeEnrichment(
    bam_granges = bam_files[[1]], 
    granges = granges
)
```

For more details and additional functions, read the 
[Introduction vignette](vignettes/Introduction.md).


---

```r
sessionInfo()

R version 3.5.2 (2018-12-20)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.2 LTS

Matrix products: default
BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
[1] VplotR_0.3.0         ggplot2_3.1.1        GenomicRanges_1.34.0
[4] GenomeInfoDb_1.18.2  IRanges_2.16.0       S4Vectors_0.20.1
[7] BiocGenerics_0.28.0  magrittr_1.5

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1                        lattice_0.20-38
 [3] prettyunits_1.0.2                 Rsamtools_1.34.1
 [5] ps_1.3.0                          Biostrings_2.50.2
 [7] zoo_1.8-5                         assertthat_0.2.1
 [9] rprojroot_1.3-2                   digest_0.6.18
[11] R6_2.4.0                          plyr_1.8.4
[13] backports_1.1.3                   pillar_1.3.1
[15] zlibbioc_1.28.0                   rlang_0.3.2
[17] lazyeval_0.2.2                    curl_3.3
[19] callr_3.2.0                       Matrix_1.2-17
[21] labeling_0.3                      desc_1.2.0
[23] devtools_2.0.1                    BiocParallel_1.16.6
[25] stringr_1.4.0                     RCurl_1.95-4.12
[27] munsell_0.5.0                     DelayedArray_0.8.0
[29] rtracklayer_1.42.2                compiler_3.5.2
[31] pkgconfig_2.0.2                   pkgbuild_1.0.3
[33] tcltk_3.5.2                       tidyselect_0.2.5
[35] SummarizedExperiment_1.12.0       tibble_2.1.1
[37] GenomeInfoDbData_1.2.0            matrixStats_0.54.0
[39] XML_3.98-1.19                     crayon_1.3.4
[41] dplyr_0.8.0.1                     withr_2.1.2
[43] GenomicAlignments_1.18.1          bitops_1.0-6
[45] grid_3.5.2                        gtable_0.3.0
[47] scales_1.0.0                      cli_1.1.0
[49] stringi_1.3.1                     XVector_0.22.0
[51] reshape2_1.4.3                    fs_1.2.7
[53] remotes_2.0.2                     testthat_2.0.1
[55] cowplot_0.9.4                     RColorBrewer_1.1-2
[57] tools_3.5.2                       BSgenome_1.50.0
[59] Biobase_2.42.0                    glue_1.3.1
[61] purrr_0.3.2                       processx_3.3.0
[63] pkgload_1.0.2                     colorspace_1.4-1
[65] sessioninfo_1.1.1                 BSgenome.Celegans.UCSC.ce11_1.4.2
[67] memoise_1.1.0                     usethis_1.4.0
```