#' A function to compute Vmat.
#'
#' @param reads.granges GRanges. The paired-end fragments
#' @param target.granges GRanges. The regions to map the fragments onto
#' @param XAXIS.CUTOFF single integer. The x limits (-/+) of the computed Vmat
#' @param YAXIS.CUTOFF single integer. The y upper limit of the computed Vmat
#' @param stranded Boolean. Should the Vmat be stranded? If TRUE, negative-strand 
#' features are flipped.
#' 
#' @return A table object
#' 
#' @import magrittr
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @export

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
    extended.targets <- center.targets + 2*XAXIS.CUTOFF
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
    UP <- -2 * XAXIS.CUTOFF
    DOWN <- 2 * XAXIS.CUTOFF
    dists %<>% factor(levels = c(UP:DOWN))
    widths <- c(IRanges::width(granges)) %>% 
        factor(levels = 0:YAXIS.CUTOFF)
    tab <- table(dists, widths)
    tab <- tab[seq(floor(nrow(tab) / 2) - XAXIS.CUTOFF, floor(nrow(tab) / 2) + XAXIS.CUTOFF, 1)[1:(XAXIS.CUTOFF*2)], 1:YAXIS.CUTOFF]
    return(tab)
}

#' A function to scale a Vmat.
#'
#' @param Vmat 
#' @param FUN string. A Vmat can be scaled relative to its median ('pctmedian'), 
#' to its mean ('pctmean'), to its max ('pctmax'). Otherwise it could be zscore-d
#' entirely ('Zscore') or by rows ('rowZscore') or by columns ('colZscore').
#' 
#' @return A table object
#' 
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @export

scaleVmat <- function(Vmat, FUN = 'pctsum') {
    if (FUN == 'pctmedian') {
        Vmat <- Vmat/median(Vmat)*100
    } else if (FUN == 'pctmean') {
        Vmat <- Vmat/mean(Vmat)*1000000
    } else if (FUN == 'pctsum') {
        Vmat <- Vmat/sum(Vmat)*1000000
    } else if (FUN == 'pctmax') {
        Vmat <- Vmat/max(Vmat)*100
    } else if (FUN == 'colZscore') {
        Vmat <- apply(Vmat, 2, scale)
    } else if (FUN == 'rowZscore') {
        Vmat <- t(apply(Vmat, 1, scale))
    } else if (FUN == 'Zscore') {
        Vmat2 <- matrix(scale(c(Vmat)), byrow = FALSE, nrow = nrow(Vmat))
        colnames(Vmat2) <- colnames(Vmat)
        row.names(Vmat2) <- row.names(Vmat)
        Vmat <- Vmat2
    }
    return(Vmat)
}

#' A function to normalized a Vmat to a given background.
#'
#' @param Vmat1 A Vmat, usually output of computeVmat
#' @param background A Vmat computed from a background. 
#' @param scale Boolean Should the background be removed?
#' @param normFun string A Vmat can be scaled relative to its median ('pctmedian'), 
#' to its mean ('pctmean'), to its max ('pctmax'). Otherwise it could be zscore-d
#' entirely ('Zscore') or by rows ('rowZscore') or by columns ('colZscore').
#' @param roll integer to use as the window to smooth the Vmat rows by rolling mean
#' 
#' @return A normalized Vmat object
#' 
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @importFrom zoo rollmean
#' @export

normalizeVmat <- function(Vmat1, background = NULL, scale = TRUE, FUN = `-`, normFun = 'pctsum', roll = 1) {
    if (is.null(background))
        background <- 0
    if (scale == TRUE) {
        Vmat <- do.call(FUN, list(Vmat1, background))
        Vmat <- scaleVmat(Vmat, normFun)
    }
    if (roll > 1) {
        Vmat <- zoo::rollmean(Vmat, roll)
    }
    return(Vmat)
}

#' A function to shuffle a Vmat
#'
#' @param Vmat A Vmat, usually output of computeVmat
#' @param background A Vmat computed from a background. 
#' @param scale Boolean Should the background be removed?
#' @param normFun string A Vmat can be scaled relative to its median ('pctmedian'), 
#' to its mean ('pctmean'), to its max ('pctmax'). Otherwise it could be zscore-d
#' entirely ('Zscore') or by rows ('rowZscore') or by columns ('colZscore').
#' @param roll integer to use as the window to smooth the Vmat rows by rolling mean
#' 
#' @return A normalized Vmat object
#' 
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @export

shuffleVmat <- function(Vmat, SEED = 999) {
    shuffled.Vmat <- c(Vmat)
    set.seed(SEED)
    shuffled.Vmat <- shuffled.Vmat[sample(x = 1:length(shuffled.Vmat), size = length(shuffled.Vmat), replace = F)]
    shuffled.Vmat <- matrix(shuffled.Vmat, nrow = nrow(Vmat), ncol = ncol(Vmat))
    return(shuffled.Vmat)
}

#' A function to plot a Vmat
#'
#' @return A Vmat ggplot
#' 
#' @export

plotVmat <- function(x, ...) {
    UseMethod("plotVmat")
}

#' A function to plot a computed Vmat
#'
#' @param Vmat A computed Vmat (should be normalized)
#' @param pdf string save the plot in a pdf
#' @param HM.COLOR.CUTOFF Integer should be between 0 and 100. Used to automatically
#' scale the range of colors (ebst to keep between 90 and 100)
#' @param colors a vector of colors
#' @param breaks a vector of breaks. length(breaks) == length(colors) + 1
#' @param xlim vector of two integers  x limits
#' @param ylim vector of two integers y limits
#' @param main string Title of the plot
#' @param xlab string x-axis label
#' @param ylab string y-axis label
#' 
#' @return A Vmat ggplot
#' 
#' @import ggplot2
#' @import RColorBrewer
#' @import reshape2
#' @export


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

#' A function to plot a computed Vmat
#'
#' @param bam_granges GRanges. The paired-end fragments
#' @param granges GRanges. The regions to map the fragments onto
#' @param XAXIS.CUTOFF single integer. The x limits (-/+) of the computed Vmat
#' @param YAXIS.CUTOFF single integer. The y upper limit of the computed Vmat
#' @param stranded Boolean. Should the Vmat be stranded? If TRUE, negative-strand 
#' features are flipped
#' @param normalize Boolean Should the Vmat be normalized?
#' @param Vmat2 table A Vmat pre-computed over a background
#' @param estimate_background Boolean Should the background be estimated automatically?
#' @param background_granges GRanges The genomic loci used to compute the background Vmat
#' @param normFun string A Vmat can be scaled relative to its median ('pctmedian'), 
#' to its mean ('pctmean'), to its max ('pctmax'). Otherwise it could be zscore-d
#' entirely ('Zscore') or by rows ('rowZscore') or by columns ('colZscore').
#' @param roll integer to use as the window to smooth the Vmat rows by rolling mean
#' @param return_Vmat Boolean. Should the function return the computed Vmat 
#' rather than the plot?
#' @param genome a BSgenome object. See getChromSizes for more details
#' 
#' @return A Vmat ggplot
#' 
#' @import ggplot2
#' @import RColorBrewer
#' @import reshape2
#' @export

plotVmat.GRanges <- function(
    bam_granges, 
    granges,
    XAXIS.CUTOFF = 1000, 
    YAXIS.CUTOFF = 1000, 
    stranded = TRUE,
    normalize = TRUE,
    Vmat2 = NULL,
    estimate_background = TRUE,
    background_granges = NULL,
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
            if (is.null(background_granges)) background_granges <- shuffleGRanges(granges, genome = genome)
            Vmat2 <- computeVmat(bam_granges, background_granges, XAXIS.CUTOFF, YAXIS.CUTOFF)
        }
        Vmat <- normalizeVmat(Vmat, background = Vmat2, scale = TRUE, normFun = normFun, roll = roll)
    }
    # Replace NA / inf values by 0
    Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
    if (return_Vmat == TRUE) {
        return(Vmat)
    } else {
        plotVmat(Vmat, ...)
    }
}

#' A function to plot a list of Vmats
#'
#' @param list_params list Each element of the list should be a list, containing
#' the bam_granges as the first sub-element and the granges as the 
#' second sub-element. Can accept the background_granges as a third sub-element.
#' @param XAXIS.CUTOFF single integer. The x limits (-/+) of the computed Vmat
#' @param YAXIS.CUTOFF single integer. The y upper limit of the computed Vmat
#' @param stranded Boolean. Should the Vmat be stranded? If TRUE, negative-strand 
#' features are flipped
#' @param normalize Boolean Should the Vmat be normalized?
#' @param Vmat2 table A Vmat pre-computed over a background
#' @param estimate_background Boolean Should the background be estimated automatically?
#' @param background_granges GRanges The genomic loci used to compute the background Vmat
#' @param normFun string A Vmat can be scaled relative to its median ('pctmedian'), 
#' to its mean ('pctmean'), to its max ('pctmax'). Otherwise it could be zscore-d
#' entirely ('Zscore') or by rows ('rowZscore') or by columns ('colZscore').
#' @param roll integer to use as the window to smooth the Vmat rows by rolling mean
#' @param return_Vmat Boolean. Should the function return the computed Vmat 
#' rather than the plot?
#' @param genome a BSgenome object. See getChromSizes for more details
#' @param verbose Boolean
#' 
#' @return A list of Vmat ggplots
#' 
#' @import parallel
#' @import ggplot2
#' @import RColorBrewer
#' @import reshape2
#' @export


plotVmat.list <- function(
    list_params, 
    cores = 1,
    XAXIS.CUTOFF = 1000, 
    YAXIS.CUTOFF = 1000, 
    stranded = TRUE,
    normalize = TRUE,
    Vmat2 = NULL,
    estimate_background = TRUE,
    background_granges = NULL,
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
            if (length(list_params[[K]]) >= 3) {background_granges <- list_params[[K]][[3]]} else {background_granges}
            Vmat <- plotVmat(
                bam_granges, 
                granges, 
                XAXIS.CUTOFF, 
                YAXIS.CUTOFF, 
                stranded, 
                normalize, 
                Vmat2, 
                estimate_background,
                background_granges, 
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
            if (length(list_params[[K]]) >= 3) {background_granges <- list_params[[K]][[3]]} else {background_granges}
            Vmat <- plotVmat(
                bam_granges, 
                granges, 
                XAXIS.CUTOFF, 
                YAXIS.CUTOFF, 
                stranded, 
                normalize, 
                Vmat2, 
                estimate_background,
                background_granges, 
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

#' A function to compute nucleosome enrichment over a set of GRanges
#'
#' @return list
#' 
#' @export

nucleosomeEnrichment <- function(x, ...) {
    UseMethod("nucleosomeEnrichment")
}

#' A function to compute nucleosome enrichment over a set of GRanges
#'
#' @param Vmat a computed Vmat. Does not need to be normalized. 
#' @param background a background Vmat
#' 
#' @return list
#' 
#' @export

nucleosomeEnrichment.default <- function(Vmat, background, ...) {
    res <- computeNucleosomeEnrichmentOverBackground(
        Vmat = Vmat, 
        background = background, 
        plus1_nuc_only = plus1_nuc_only, 
        ...
    )
    return(res)
}

#' A function to compute nucleosome enrichment over a set of GRanges
#'
#' @param bam_granges GRanges. The paired-end fragments
#' @param granges GRanges. The regions to map the fragments onto
#' @param estimate_background Boolean Should the background be estimated automatically?
#' @param plus1_nuc_only Should compute nucleosome enrichment only for +1 
#' nucleosome?
#' @param verbose Boolean
#' @param genome a BSgenome object. See getChromSizes for more details

#' @return list
#' 
#' @export

nucleosomeEnrichment.GRanges <- function(
    bam_granges, 
    granges, 
    estimate_background = TRUE, 
    plus1_nuc_only = FALSE, 
    verbose = TRUE, 
    genome = 'ce11',
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
            shuffleGRanges(granges, genome = genome)
        )
    } 
    else {
        background = NULL
    }
    if (verbose) message("Computing enrichment...")
    res <- computeNucleosomeEnrichmentOverBackground(
        Vmat = Vmat, 
        background = background, 
        plus1_nuc_only = plus1_nuc_only, 
        ...
    )
    return(res)
}

#' A function to compute nucleosome enrichment of a Vmat. Should only be used 
#' internally.
#'
#' @param Vmat A Vmat computed by nucleosomeEnrichment function
#' @param background a background Vmat
#' @param plus1_nuc_only Should compute nucleosome enrichment only for +1 
#' nucleosome?
#' @param minus1_nuc list where the -1 nucleosome is located
#' @param minus1_nuc_neg where the background of the -1 nucleosome is located
#' @param plus1_nuc where the +1 nucleosome is located
#' @param plus1_nuc where the background of the +1 nucleosome is located

#' @return list
#' 
#' @export

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
    } 
    else {
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
        plot = list(p, q)
    ))
}
