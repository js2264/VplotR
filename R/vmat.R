#### ---- Accessory functions ---- ####

computeVmat <- function(reads.granges, target.granges, XAXIS.CUTOFF = 1000, YAXIS.CUTOFF = 1000, stranded = TRUE) {
    `%over%` <- IRanges::`%over%`
    center.targets <- GenomicRanges::GRanges(
        seqnames = GenomicRanges::seqnames(target.granges), 
        IRanges::IRanges(
            start = IRanges::start(target.granges) + floor((IRanges::end(target.granges) - IRanges::start(target.granges)) / 2), 
            width = 1
        ), 
        strand = GenomicRanges::strand(target.granges)
    )
    if (stranded) center.targets <- deconvolveBidirectionalPromoters(center.targets)
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
    dists <- c(IRanges::start(center.target.per.read) - IRanges::start(center.reads.subset))
    if (stranded) dists <- ifelse(strand(center.target.per.read) == '-', dists, -dists)
    dists %<>% factor(levels = c(-2000:2000))
    widths <- c(IRanges::width(granges)) %>% 
        factor(levels = 0:YAXIS.CUTOFF)
    tab <- table(dists, widths)
    tab <- tab[seq(floor(nrow(tab) / 2) - XAXIS.CUTOFF, floor(nrow(tab) / 2) + XAXIS.CUTOFF, 1)[1:(XAXIS.CUTOFF*2)], 1:YAXIS.CUTOFF]
    return(tab)
}

scaleVmat <- function(Vmat, FUN = 'pctsum') {
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
    } else if (FUN == 'Zscore') {
        Vmat <- scale(Vmat)
    }
    return(Vmat)
}

normalizeVmat <- function(Vmat1, background = NULL, scale = TRUE, FUN = `/`, normFUN = 'pctsum', roll = 1) {
    if (is.null(background))
        background <- 0
    if (scale == TRUE) {
        Vmat <- scaleVmat(Vmat1 - background, normFUN)
    } else {
        Vmat <- do.call(FUN, Vmat1, background)
    }
    if (roll > 1) {
        Vmat <- zoo::rollmean(Vmat, roll)
    }
    return(Vmat)
}

shuffleVmat <- function(Vmat, SEED = 999) {
    shuffled.Vmat <- c(Vmat)
    set.seed(SEED)
    shuffled.Vmat <- shuffled.Vmat[sample(x = 1:length(shuffled.Vmat), size = length(shuffled.Vmat), replace = F)]
    shuffled.Vmat <- matrix(shuffled.Vmat, nrow = nrow(Vmat), ncol = ncol(Vmat))
    return(shuffled.Vmat)
}

#### ---- PlotVmat function ---- ####

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
) 
{
    # Extract the Vmat to the specified limits
    if (!is.null(xlim) & !is.null(ylim)) {
        if (class(Vmat) == 'matrix' || class(Vmat) == 'table') {
            Vmat <- Vmat[(round(nrow(Vmat)/2) + xlim[1]) : (round(nrow(Vmat)/2) + xlim[2]), (ylim[1]+1) : (ylim[2])]
            row.names(Vmat) <- as.character(xlim[1]:xlim[2])
        } else if (class(Vmat) == 'list') {
            Vmat <- lapply(Vmat, function(V) {
                V <- V[(round(nrow(V)/2) + xlim[1]) : (round(nrow(V)/2) + xlim[2]), (ylim[1]+1) : (ylim[2])]
                row.names(V) <- as.character(xlim[1]:xlim[2])
                return(V)
            })
        }
    }
    # Define breaks and clamp matrix within breaks
    if (class(Vmat) == 'matrix' || class(Vmat) == 'table') {
        if (is.null(breaks)) {
            probs <- quantile(c(Vmat), probs = seq(0, 1, length.out = 101))
            breaks <- seq(probs[(100-HM.COLOR.CUTOFF)/2], probs[HM.COLOR.CUTOFF+(100-HM.COLOR.CUTOFF)/2], length.out = length(colors)+1)
        } 
        Vmat[Vmat < min(breaks)] <- min(breaks)
        Vmat[Vmat > max(breaks)] <- max(breaks)
    } else if (class(Vmat) == 'list') {
        if (is.null(breaks)) {
            probs <- quantile(c(do.call(rbind, Vmat)), probs = seq(0, 1, length.out = 101))
            breaks <- seq(probs[(100-HM.COLOR.CUTOFF)/2], probs[HM.COLOR.CUTOFF+(100-HM.COLOR.CUTOFF)/2], length.out = length(colors)+1)
        }
        Vmat <- lapply(Vmat, function(V) {
            V[V < min(breaks)] <- min(breaks)
            V[V > max(breaks)] <- max(breaks)
            return(V)
        })
    }
    # Plot
    df <- reshape2::melt(Vmat)
    if (class(Vmat) == 'list') {
        colnames(df) <- c('Var1', 'Var2', 'value', 'Cond.')
        df$Cond. <- factor(df$Cond., levels = names(Vmat))
    } else {
        colnames(df) <- c('Var1', 'Var2', 'value')
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = value))
    p <- p + ggplot2::geom_raster()
    p <- p + ggplot2::scale_fill_gradientn(
        breaks = c(min(breaks), max(breaks)), 
        colors = colors
    )
    p <- p + ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(), 
        panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = 1)
    )
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
    p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
    p <- p + ggplot2::labs(
        title = main,
        x = xlab,
        y = ylab, 
        fill = 'Score'
    )
    p <- p + ggplot2::coord_fixed()
    if (!is.null(pdf))
        ggplot2::ggsave(pdf, p)
    return(p)
}

plotVmat.GRanges <- function(
    bam_granges, 
    granges,
    XAXIS.CUTOFF = 1000, 
    YAXIS.CUTOFF = 1000, 
    stranded = TRUE,
    normalize = TRUE,
    Vmat2 = NULL,
    estimate_background = TRUE,
    background.granges = NULL,
    normFun = 'pctsum',
    roll = 3,
    return_Vmat = FALSE,
    genome = 'ce11',
    ...
) 
{
    # Calculate Vmat
    Vmat <- computeVmat(bam_granges, granges, XAXIS.CUTOFF, YAXIS.CUTOFF, stranded)
    # Normalize Vmat using Vmat2
    if (normalize) {
        if (estimate_background == TRUE & is.null(Vmat2)) {
            if (is.null(background.granges)) background.granges <- shuffleGRanges(granges, genome = genome)
            Vmat2 <- computeVmat(bam_granges, background.granges, XAXIS.CUTOFF, YAXIS.CUTOFF)
        }
        Vmat <- normalizeVmat(Vmat, background = Vmat2, scale = TRUE, normFUN = normFun, roll = roll)
    }
    # Replace NA / inf values by 0
    Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
    if (return_Vmat == TRUE) {
        return(Vmat)
    } else {
        plotVmat(Vmat, ...)
    }
}

plotVmat.list <- function(
    list_params, 
    cores = 1,
    XAXIS.CUTOFF = 1000, 
    YAXIS.CUTOFF = 1000, 
    stranded = TRUE,
    normalize = TRUE,
    Vmat2 = NULL,
    estimate_background = TRUE,
    background.granges = NULL,
    normFun = 'pctsum',
    roll = 3,
    return_Vmat = FALSE,
    genome = 'ce11',
    verbose = TRUE,
    ...
) 
{
    # Calculate all the Vmats
    if (is.null(names(list_params))) stop('Please provide a *named* list of parameters. Aborting.')
    if (cores == 1) {
        Vmats_list <- lapply(seq_along(list_params), function(K) {
            if (verbose) message('- Processing sample ', K, '/', length(list_params)) 
            bam_granges <- list_params[[K]][[1]]
            granges <- list_params[[K]][[2]]
            Vmat <- plotVmat(
                bam_granges, 
                granges, 
                XAXIS.CUTOFF, 
                YAXIS.CUTOFF, 
                stranded, 
                normalize, 
                Vmat2, 
                estimate_background,
                background.granges, 
                normFun, 
                roll, 
                return_Vmat = TRUE, 
                genome
            )
            # Replace NA / inf values by 0
            Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
            return(Vmat)
        }) %>% setNames(names(list_params))
    } 
    else if (cores > 1) {
        Vmats_list <- parallel::mclapply(seq_along(list_params), function(K) {
            if (verbose) message('- Processing sample ', K, '/', length(list_params)) 
            bam_granges <- list_params[[K]][[1]]
            granges <- list_params[[K]][[2]]
            Vmat <- plotVmat(
                bam_granges, 
                granges, 
                XAXIS.CUTOFF, 
                YAXIS.CUTOFF, 
                stranded, 
                normalize, 
                Vmat2, 
                estimate_background,
                background.granges, 
                normFun, 
                roll, 
                return_Vmat = TRUE, 
                genome
            )
            # Replace NA / inf values by 0
            Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
            return(Vmat)
        }, mc.cores = cores) %>% setNames(names(list_params))
    }
    # Plot
    if (return_Vmat == TRUE) {
        return(Vmats_list)
    } 
    else {
        plotVmat.default(Vmats_list, ...)
    }
}

#### ---- nucleosomeEnrichment function ---- ####

nucleosomeEnrichment <- function(x, ...) {
    UseMethod("nucleosomeEnrichment")
}

nucleosomeEnrichment.default <- function(Vmat, background, ...) {
    mat <- matrix(getVec(
        Vmat = Vmat, 
        background = background, 
        plus1_nuc_only = plus1_nuc_only
    ), ncol = 2)
    attr(Vmat, "fisher.test") <- fisher.test(mat)
    return(Vmat)
}

nucleosomeEnrichment.GRanges <- function(
    bam_granges, 
    granges, 
    estimate_background = TRUE, 
    plus1_nuc_only = FALSE, 
    verbose = TRUE, 
    ...
) 
{
    if (verbose) message("Computing Vmat...")
    Vmat <- computeVmat(
        bam_granges, 
        granges
    )
    if (estimate_background) {
        if (verbose) message("Computing background...")
        background <- computeVmat(
            bam_granges, 
            shuffleGRanges(granges)
        )
    } 
    else {
        background = NULL
    }
    if (verbose) message("Computing enrichment...")
    res <- computeNucleosomeEnrichmentOverBackground(
        Vmat = Vmat, 
        background = background, 
        plus1_nuc_only = plus1_nuc_only
    )
    return(res)
}

computeNucleosomeEnrichmentOverBackground <- function(
    Vmat, 
    background = NULL, 
    plus1_nuc_only = TRUE, 
    minus1_nuc = list(c(xmin = -150, xmax = -100), c(ymin = 175, ymax = 270)), 
    minus1_nuc_neg = list(c(xmin = -150, xmax = -100), c(ymin = 50, ymax = 145)), 
    plus1_nuc = list(c(xmin = 100, xmax = 150), c(ymin = 175, ymax = 270)), 
    plus1_nuc_neg = list(c(xmin = 100, xmax = 150), c(ymin = 50, ymax = 145))
) 
{
    if (is.null(background)) {
        background <- shuffleVmat(Vmat)
    }
    # Generate plot
    minus1_nuc_df <- data.frame(xmin = minus1_nuc[[1]][[1]], xmax = minus1_nuc[[1]][[2]], ymin = minus1_nuc[[2]][[1]], ymax = minus1_nuc[[2]][[2]]) 
    plus1_nuc_df <- data.frame(xmin = plus1_nuc[[1]][[1]], xmax = plus1_nuc[[1]][[2]], ymin = plus1_nuc[[2]][[1]], ymax = plus1_nuc[[2]][[2]])
    minus1_nuc_neg_df <- data.frame(xmin = minus1_nuc_neg[[1]][[1]], xmax = minus1_nuc_neg[[1]][[2]], ymin = minus1_nuc_neg[[2]][[1]], ymax = minus1_nuc_neg[[2]][[2]])
    plus1_nuc_neg_df <- data.frame(xmin = plus1_nuc_neg[[1]][[1]], xmax = plus1_nuc_neg[[1]][[2]], ymin = plus1_nuc_neg[[2]][[1]], ymax = plus1_nuc_neg[[2]][[2]])
    df_controls <- data.frame(
        rbind(minus1_nuc_df, plus1_nuc_df, minus1_nuc_neg_df, plus1_nuc_neg_df), 
        locus = c('minus1_nuc', 'plus1_nuc', 'minus1_nuc_neg', 'plus1_nuc_neg'), 
        type = factor(c('nuc', 'nuc', 'neg', 'neg'), levels = c('nuc', 'neg'))
    )
    p <- plotVmat(Vmat) + 
        geom_rect(
            data = df_controls, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, col = type), 
            inherit.aes=FALSE, 
            fill = NA
        )
    q <- plotVmat(background) + 
        geom_rect(
            data = df_controls, 
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, col = type), 
            inherit.aes=FALSE, 
            fill = NA
        )
    r <- cowplot::plot_grid(p, q, nrow = 1)
    # Shift ranges to center around the "0"
    minus1_nuc <- list(
        minus1_nuc[[1]] + round(nrow(Vmat)/2), 
        minus1_nuc[[2]] - as.numeric(colnames(Vmat)[1])
    )
    plus1_nuc <- list(
        plus1_nuc[[1]] + round(nrow(Vmat)/2), 
        plus1_nuc[[2]] - as.numeric(colnames(Vmat)[1])
    )
    minus1_nuc_neg <- list(
        minus1_nuc_neg[[1]] + round(nrow(Vmat)/2), 
        minus1_nuc_neg[[2]] - as.numeric(colnames(Vmat)[1])
    )
    plus1_nuc_neg <- list(
        plus1_nuc_neg[[1]] + round(nrow(Vmat)/2), 
        plus1_nuc_neg[[2]] - as.numeric(colnames(Vmat)[1])
    )
    # Count reads in each window
    if (!plus1_nuc_only) {
        vec <- c(
            sum(Vmat[minus1_nuc[[1]][1] : minus1_nuc[[1]][2], minus1_nuc[[2]][1] : minus1_nuc[[2]][2]]) + sum(Vmat[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]), 
            sum(Vmat[minus1_nuc_neg[[1]][1] : minus1_nuc_neg[[1]][2], minus1_nuc_neg[[2]][1] : minus1_nuc_neg[[2]][2]]) + sum(Vmat[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]]),
            sum(background[minus1_nuc[[1]][1] : minus1_nuc[[1]][2], minus1_nuc[[2]][1] : minus1_nuc[[2]][2]]) + sum(background[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]), 
            sum(background[minus1_nuc_neg[[1]][1] : minus1_nuc_neg[[1]][2], minus1_nuc_neg[[2]][1] : minus1_nuc_neg[[2]][2]]) + sum(background[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]])
        )
    } else {
        vec <- c(
            sum(Vmat[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]), 
            sum(Vmat[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]]),
            sum(background[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]), 
            sum(background[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]])
        )
    }
    # Return result
    return(list(
        scores = matrix(vec, ncol = 2),
        fisher_test = fisher.test(matrix(vec, ncol = 2)),
        plot = r
    ))
}
