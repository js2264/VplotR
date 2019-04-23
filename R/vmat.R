computeVmat <- function(reads.granges, target.granges, XAXIS.CUTOFF = 350, YAXIS.CUTOFF = 300) {
    `%over%` <- IRanges::`%over%`
    center.targets <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(target.granges), 
        IRanges::IRanges(
            start = IRanges::start(target.granges) + floor((IRanges::end(target.granges) - IRanges::start(target.granges)) / 2), 
            width = 1
        )
    )
    extended.targets <- center.targets + 2000
    center.reads <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(reads.granges), 
        IRanges::IRanges(
            start = IRanges::start(reads.granges) + floor((IRanges::end(reads.granges) - IRanges::start(reads.granges)) / 2),
             width = 1
         )
     )
    granges <- reads.granges[center.reads %over% extended.targets]
    center.reads.subset <- center.reads[center.reads %over% extended.targets]
    center.target.per.read <- center.targets[S4Vectors::subjectHits(GenomicRanges::distanceToNearest(center.reads.subset, center.targets))]
    dists <- factor(as.character(IRanges::start(center.target.per.read) - IRanges::start(center.reads.subset)), levels = -2000:2000)
    widths <- factor(as.character(IRanges::width(granges)), levels = 0:YAXIS.CUTOFF)
    tab <- table(dists, widths)
    tab <- tab[seq(round(nrow(tab) / 2) - XAXIS.CUTOFF, round(nrow(tab) / 2) + XAXIS.CUTOFF, 1), 1:YAXIS.CUTOFF]
    return(tab)
}

scaleVmat <- function(Vmat, FUN = 'pctmax') {
    if (FUN == 'pctmedian') {
        Vmat <- Vmat/median(Vmat)*100
    } else if (FUN == 'pctsum') {
        Vmat <- Vmat/sum(Vmat)*1000000
    } else if (FUN == 'pctmax') {
        Vmat <- Vmat/max(Vmat)*100
    } else if (FUN == 'colZscore') {
        Vmat <- apply(Vmat, 2, scale)
    } else if (FUN == 'rowZscore') {
        Vmat <- t(apply(Vmat, 1, scale))
    }
    return(Vmat)
}

normalizeVmat <- function(Vmat1, background = NULL, scale = TRUE, normFUN = 'pctsum', roll = 1) {
    if (is.null(background))
        background <- 0
    if (scale == T) {
        Vmat <- scaleVmat(Vmat1 - background, normFUN)
    } else {
        Vmat <- Vmat1 - background
    }
    if (roll > 1) {
        Vmat <- zoo::rollmean(Vmat, roll)
    }
    return(Vmat)
}

# --------- PlotVmat function ---------

plotVmat <- function(x, ...) {
    UseMethod("plotVmat")
}

plotVmat.default <- function(
    Vmat, 
    pdf = NULL, 
    HM.COLOR.CUTOFF = 90, 
    colors = c(
        colorRampPalette(rev(RColorBrewer::brewer.pal('Spectral', n = 10))[1:5])(30), 
        colorRampPalette(rev(RColorBrewer::brewer.pal('Spectral', n = 10))[6:10])(30)
    ), 
    breaks = NULL, 
    xlim = c(-250, 250),
    ylim = c(50, 300),
    main = '', 
    xlab = 'Distance from center of elements',
    ylab = 'Fragment length',
    ...
    ) {
    # Replace NA / inf values by 0
    Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
    # Extract the Vmat to the specified limits
    Vmat <- Vmat[(round(nrow(Vmat)/2) + xlim[1]) : (round(nrow(Vmat)/2) + xlim[2]), ylim[1] : ylim[2]]
    # Define breaks and clamp matrix within breaks
    if (is.null(breaks)) {
        Vmat[Vmat > quantile(c(Vmat), probs = seq(0, 1, length.out = 101))[HM.COLOR.CUTOFF+(100-HM.COLOR.CUTOFF)/2]] <- quantile(c(Vmat), probs = seq(0, 1, length.out = 101))[HM.COLOR.CUTOFF+(100-HM.COLOR.CUTOFF)/2]
        Vmat[Vmat < quantile(c(Vmat), probs = seq(0, 1, length.out = 101))[(100-HM.COLOR.CUTOFF)/2]] <- quantile(c(Vmat), probs = seq(0, 1, length.out = 101))[(100-HM.COLOR.CUTOFF)/2]
        breaks <- seq(min(Vmat), max(Vmat), length.out = length(colors)+1)
    } else {
        Vmat[Vmat < min(breaks)] <- min(breaks)
        Vmat[Vmat > max(breaks)] <- max(breaks)
    }
    # Plot
    p <- ggplot2::ggplot(reshape2::melt(Vmat), ggplot2::aes(Var1, Var2, fill = value)) +
        ggplot2::geom_raster() +
        ggplot2::scale_fill_gradientn(
            breaks = c(min(breaks), max(breaks)), 
            colors = colors
        ) +
        ggplot2::theme(
            plot.background = ggplot2::element_blank(),
            legend.background = ggplot2::element_blank(), 
            panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1)
        ) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::labs(
            title = main,
            x = xlab,
            y = ylab, 
            fill = 'Score'
        )
    if (!is.null(pdf))
        ggplot2::ggsave(pdf, p)
    invisible(p)
}

plotVmat.GRanges <- function(
    bam.granges, 
    granges,
    XAXIS.CUTOFF = 350, 
    YAXIS.CUTOFF = 300, 
    normalize = TRUE,
    Vmat2 = NULL,
    estimate.background = FALSE,
    normFun = 'pctsum',
    roll = 3,
    return_Vmat = FALSE,
    ...) {
    # Calculate Vmat
    Vmat <- computeVmat(bam.granges, granges, XAXIS.CUTOFF, YAXIS.CUTOFF)
    # Normalize Vmat using Vmat2
    if (normalize) {
        if (estimate.background == TRUE & is.null(Vmat2)) {
            random.granges <- shuffleGRanges(granges)
            Vmat2 <- computeVmat(bam.granges, random.granges, XAXIS.CUTOFF, YAXIS.CUTOFF)
        }
        Vmat <- normalizeVmat(Vmat, background = Vmat2, scale = TRUE, normFUN = normFun, roll = roll)
    }
    if (return_Vmat == TRUE) {
        return(Vmat)
    } else {
        plotVmat(Vmat, ...)
    }
}

# ---------- nucleosomeEnrichment function -----------

shuffleVmat <- function(Vmat, SEED = 999) {
    shuffled.Vmat <- c(Vmat)
    set.seed(SEED)
    shuffled.Vmat <- shuffled.Vmat[sample(x = 1:length(shuffled.Vmat), size = length(shuffled.Vmat), replace = F)]
    shuffled.Vmat <- matrix(shuffled.Vmat, nrow = nrow(Vmat), ncol = ncol(Vmat))
    return(shuffled.Vmat)
}

getVec <- function(
    Vmat, 
    background = NULL, 
    nuc1 = list(c(100, 150), c(120, 200)), 
    nuc2 = list(c(350, 400), c(120, 200)), 
    nuc1_neg = list(c(100, 150), c(0, 100)), 
    nuc2_neg = list(c(350, 400), c(0, 100))
    ) {
    if (is.null(background)) {
        background <- shuffleVmat(Vmat)
    }
    vec <- c(
        sum(Vmat[nuc1[[1]][1] : nuc1[[1]][2], nuc1[[2]][1] : nuc1[[2]][2]]) + sum(Vmat[nuc2[[1]][1] : nuc2[[1]][2], nuc2[[2]][1] : nuc2[[2]][2]]), 
        sum(Vmat[nuc1_neg[[1]][1] : nuc1_neg[[1]][2], nuc1_neg[[2]][1] : nuc1_neg[[2]][2]]) + sum(Vmat[nuc2_neg[[1]][1] : nuc2_neg[[1]][2], nuc2_neg[[2]][1] : nuc2_neg[[2]][2]]),
        sum(background[nuc1[[1]][1] : nuc1[[1]][2], nuc1[[2]][1] : nuc1[[2]][2]]) + sum(background[nuc2[[1]][1] : nuc2[[1]][2], nuc2[[2]][1] : nuc2[[2]][2]]), 
        sum(background[nuc1_neg[[1]][1] : nuc1_neg[[1]][2], nuc1_neg[[2]][1] : nuc1_neg[[2]][2]]) + sum(background[nuc2_neg[[1]][1] : nuc2_neg[[1]][2], nuc2_neg[[2]][1] : nuc2_neg[[2]][2]])
    )
    return(vec)
}

nucleosomeEnrichment <- function(x, ...) {
    UseMethod("nucleosomeEnrichment")
}

nucleosomeEnrichment.default <- function(Vmat, background, subset = list(c(100:600), c(51:300)), ...) {
    mat <- matrix(getVec(
        Vmat = Vmat[subset[[1]], subset[[2]]], 
        background = background[subset[[1]], subset[[2]]]
    ), ncol = 2)
    attr(Vmat, "fisher.test") <- fisher.test(mat)
    return(Vmat)
}

nucleosomeEnrichment.GRanges <- function(bam.granges, granges, estimate.background = TRUE, subset = list(c(100:600), c(51:300)), verbose = TRUE, ...) {
    if (verbose) 
        message("Computing Vmat...")
    Vmat <- computeVmat(
        bam.granges, 
        granges
    )
    if (verbose) 
        message("Computing background...")
    if (estimate.background) {
        background <- computeVmat(
            bam.granges, 
            shuffleGRanges(granges)
        )
    } else {
        background = NULL
    }
    if (verbose) 
        message("Computing enrichment...")
    mat <- matrix(getVec(
        Vmat = Vmat[subset[[1]], subset[[2]]], 
        background = background[subset[[1]], subset[[2]]]
    ), ncol = 2)
    attr(Vmat, "fisher.test") <- fisher.test(mat)
    return(Vmat)
}