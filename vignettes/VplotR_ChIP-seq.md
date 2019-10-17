---
title: "Advanced use of VplotR package to look at ChIP-seq data"
author: "Jacques Serizay"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

**WARNING:** This section is still under development. 

### 3.4. Make Vplots and footprint plots @ TFs binding sites (using ATAC or MNase)

VplotR has been developed primarily to work with paired-end fragments to look at 
patterns of accessibility (either from ATAC-seq, MNase-seq or DNase-seq) at
promoters. However, it can also be used to estimate the binding pattern 
of a transcription factor at its binding sites.

```r
.sourceDirs('~/shared/bin/R/packages/memeR/R')
pfms <- list.files('~/shared/bin/R/packages/memeR/data/pfms/')
TFs_granges <- mclapply(seq_along(pfms), function(K) {
    FILE <- pfms[K]
    TF <- gsub('.pfm', '', FILE)
    motif <- importMotif(pfm_file = paste0('~/shared/bin/R/packages/memeR/data/pfms/', FILE))
    g_motif <- scanGenome(motif, genome = 'ce11')
    message('.. ', TF, ' imported.')
    return(g_motif)
}, mc.cores = 41) %>% setNames(gsub('.pfm', '', pfms))
# ATAC Vmats
list_params <- lapply(seq_along(TFs_granges), function(K) {
    TF <- names(TFs_granges)[K]
    g_motif <- TFs_granges[[TF]]
    g_motif <- g_motif[order(g_motif$relScore, decreasing = TRUE)]
    g <- g_motif[g_motif$relScore > 0.85 & strand(g_motif) == '+']
    if (length(g) > 2000) g <- g[order(g$relScore, decreasing = TRUE)][1:2000]
    l <- list(
        bam_list[['mixed']],
        resize(g, 300, fix = 'center')
    )
    return(l)
}) %>% setNames(names(TFs_granges))
plots <- plotVmat(
    list_params,
    estimate_background = FALSE, 
    ylim = c(50, 450),
    xlim = c(-300, 300), 
    cores = length(TFs_granges)
) + 
    coord_fixed(1) +
    facet_wrap(~Cond.) + 
    theme(legend.position = 'bottom') +
    theme(panel.spacing = unit(1, "lines"))
ggsave('examples/ATAC-Vmats_TFs-binding-sites.pdf', height = 15, width = 15)
# ATAC footprint profils
ATAC_cut_ce <- coverage(c(
    resize(bam_list[['mixed']], fix = 'start', width = 1),
    resize(bam_list[['mixed']], fix = 'end', width = 1)
))
plots <- plotAggregateCoverage(
    ATAC_cut_ce, 
    lapply(list_params, function(L) {L[[2]] %>% '['(strand(.) == '+') %>% resize(300, fix = 'center')}),
    BIN = 1
)
p <- plots + 
    facet_wrap(~grange) + 
    theme(legend.position = 'bottom') +
    theme(panel.spacing = unit(1, "lines")) + 
    xlim(c(-150, 150)) + 
    ylim(c(0, 1.5))
ggsave('examples/ATAC-footprints_TFs-binding-sites.pdf', height = 15, width = 15)
# MNase Vmats
pool_mnase <- unlist(GRangesList(MNase_mixed[1:4]))
list_params <- lapply(seq_along(TFs_granges), function(K) {
    TF <- names(TFs_granges)[K]
    g_motif <- TFs_granges[[TF]]
    g_motif <- g_motif[order(g_motif$relScore, decreasing = TRUE)]
    g <- g_motif[g_motif$relScore > 0.85 & strand(g_motif) == '+']
    if (length(g) > 2000) g <- g[order(g$relScore, decreasing = TRUE)][1:2000]
    l <- list(
        pool_mnase,
        resize(g, 300, fix = 'center')
    )
    return(l)
}) %>% setNames(names(TFs_granges))
plots <- plotVmat(
    list_params,
    estimate_background = FALSE, 
    ylim = c(50, 190),
    xlim = c(-300, 300), 
    cores = length(TFs_granges)
) + 
    coord_fixed(3) +
    facet_wrap(~Cond.) + 
    theme(legend.position = 'bottom') +
    theme(panel.spacing = unit(1, "lines"))
ggsave('examples/MNase-Vmats_TFs-binding-sites.pdf', height = 15, width = 15)
# MNase footprint profil
MNase_cut <- coverage(c(
    resize(c(MNase_mixed[[1]], MNase_mixed[[2]], MNase_mixed[[3]], MNase_mixed[[4]]), fix = 'start', width = 1),
    resize(c(MNase_mixed[[1]], MNase_mixed[[2]], MNase_mixed[[3]], MNase_mixed[[4]]), fix = 'end', width = 1)
))
plots <- plotAggregateCoverage(
    MNase_cut, 
    lapply(list_params, function(L) {L[[2]] %>% '['(strand(.) == '+') %>% resize(300, fix = 'center')}),
    BIN = 1
) + 
    facet_wrap(~grange) + 
    theme(legend.position = 'bottom') +
    theme(panel.spacing = unit(1, "lines")) + 
    xlim(c(-40, 40)) + 
    ylim(c(0, 2))
ggsave('examples/MNase-footprints_TFs-binding-sites.pdf', height = 15, width = 15)
```

### 3.5. Look at CTCF binding sites motifs in human ATAC

VplotR can be used to study TF binding pattern / footprint in any organism.
Here we focus on CTCF binding in human.

```r
load('~/20190730_ATAC_hPGCs_Walfred/.bam.list.RData')
ATAC_cut <- coverage(c(
    resize(bam.list$ATAC_hPGC, fix = 'start', width = 1),
    resize(bam.list$ATAC_hPGC, fix = 'end', width = 1)
))
CTCF <- importMotif(pfm_file = '~/shared/bin/R/packages/memeR/data/MA0139.1.pfm')
g_CTCF <- scanGenome(CTCF, genome = 'hg38')
p <- plotVmat(
    bam.list$ATAC_hPGC, 
    resize(g_CTCF[g_CTCF$relScore > 0.90], 300, fix = 'center'),
    estimate_background = FALSE, 
    ylim = c(0, 350),
    xlim = c(-200, 200),
    stranded = TRUE,
    main = paste0('ATAC-seq (hPGCs)\nfragments over CTCF binding sites')
) + coord_fixed(1)
ggsave('examples/CTCF-binding-motif_Vplot.pdf')
p <- plotAggregateCoverage(
    ATAC_cut, 
    resize(g_CTCF[g_CTCF$relScore > 0.90 & seqnames(g_CTCF) %in% seqlevels(ATAC_cut) & strand(g_CTCF) == '+'], 200, fix = 'center'),
    BIN = 1
)
ggsave('examples/CTCF-binding-motif_footprint.pdf')
```
