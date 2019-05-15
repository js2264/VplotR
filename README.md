# VplotR

**VplotR is still in pre-Alpha and has not been thorouhly tested yet.  
Documentation will come soon.**

![VplotR](https://raw.githubusercontent.com/js2264/VplotR/master/examples/pdf/Vplot.png)

## Introduction

This R package makes the process of generating "V-plots" straighforward. 
V-plots have been introduced for the first time in 
[a Chip-seq paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3215028/) from the 
Henikoff lab, in 2011.  

Recently, V-plots have proven to be very instructive to understand the molecular 
organization of the chromatin. For instance, the [nucleoATAC]((https://genome.cshlp.org/content/early/2015/08/27/gr.192294.115)) package allows the 
user to map nucleosome occupancy from ATAC-seq data. 
However, nucleoATAC is written in Python and unfortunately, its [port
to R](https://github.com/GreenleafLab/NucleoATACR) only allows to generate V-plots
once the Python pipeline has been ran. More generally, it does not allow the user 
to dynamically explore nucleosome patterns at specific locations in the genome.  

VplotR aim is to streamline the process of generating V-plots without having to 
leave R. It contains functions to quickly feed in .bam and .bw files and methods
designed to generate V-plots from a variety of different inputs.  

Finally, VplotR is designed around [ggplot2](https://ggplot2.tidyverse.org/) 
and makes full use of its potential. As such, it is easy to make "grobs" 
(graphic objects) and combine them to make publication-ready 
figures.  

## Installation

VplotR can be ran installed from Github as follow:

```{r}
    install.packages("devtools")
    devtools::install_github("js2264/VplotR")
```

## Overview

Firstly a .bam file is read using the `importBam()` function:
```{r}
    bam <- importBam("path/to/bam/file.bam")
```

Then a V-plot of the .bam file over a set of GRanges of interest (`granges`) 
is generated using the `plotVmat()` function:

```{r}
    vmat <- plotVmat(bam, granges, estimate.background = TRUE)
    ggplot2::ggsave('vmat.pdf', vmat)
```

The generation of multiple V-plots can be parallelized as follow:

```{r}
    params <- list(
        list("bam" = bam1, "granges" = granges1), 
        list("bam" = bam2, "granges" = granges2), 
        ..., 
        list("bam" = bamN, "granges" = grangesN)
    )
    plots <- parallel::mclapply(seq_along(params), function(K) {
        plotVmat(
            params[[K]]$bam, 
            params[[K]]$granges, 
            estimate.background = TRUE
        )
    }, mc.cores = length(params))
    plots <- cowplot::plot_grid(plotlist = plots)
    ggplot2::ggsave('mulit-plots.pdf', plots)
```

Finally, the `nucleosomeEnrichment()` function is useful to statistically quantify 
and compare nucleosome enrichment (e.g. flanking nucleosome at promoters). To do so:

```{r}
    nucleosomeEnrichment(
        bam.granges = bam, 
        granges = granges, 
        estimate.background = TRUE
    )
```

For more details and additional functions, read the 
[Introduction vignette](vignettes/Introduction.Rmd) (Not available yet).

