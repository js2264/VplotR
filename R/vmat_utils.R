#' A function to compute Vplot matrix
#' 
#' This function computes the underlying matrix shown as a heatmap
#' in Vplots. For each pair of coordinates (x: distance from fragment
#' midpoint to center of GRanges of interest; y: fragment size), the 
#' function computes how many fragments there are. 
#'
#' @param bam_granges GRanges, paired-end fragments
#' @param granges GRanges, regions to map the fragments onto
#' @param cores Integer, nb of threads to parallelize fragments subsetting
#' @param xlims The x limits of the computed Vmat
#' @param ylims The y limits of the computed Vmat
#' @return A table object
#' 
#' @import magrittr
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @export 
#' 
#' @examples
#' data(bam_test)
#' data(ce11_all_REs)
#' Vmat <- computeVmat(bam_test, ce11_all_REs)
#' dim(Vmat)
#' Vmat[seq(1,5), seq(1,10)]

computeVmat <- function(
    bam_granges, 
    granges, 
    cores = 1,
    xlims = c(-250, 250), 
    ylims = c(50, 300)
) 
{
    # Fragments
    frags <- bam_granges
    
    # Subset frags whose center overlap with extended targets
    targets_extended <- GenomicRanges::resize(
        granges, fix = 'center', width = (xlims[2] - xlims[1])
    )
    breaks <- seq(1, length(frags), length.out = cores + 1)
    frags_l <- parallel::mclapply(
        mc.cores = cores, 
        seq_len(cores), 
        function(K) {
            g <- IRanges::subsetByOverlaps(
                frags[breaks[K]:breaks[K+1]], 
                targets_extended, 
                ignore.strand = TRUE
            )
            return(g)
        }
    )
    frag_ov <- unlist(GRangesList(frags_l))
    
    # Subset frags whose size is within ylims
    frags_subset <- frag_ov[(IRanges::width(frag_ov) >= ylims[1]) & 
        (IRanges::width(frag_ov) <= ylims[2])]
    
    # Fast resize many GRanges
    frags_subset_centers <- GenomicRanges::GRanges(
        GenomicRanges::seqnames(frags_subset),
        IRanges::IRanges(
            GenomicRanges::start(frags_subset) + 
                round(width(frags_subset)/2) - 1,
            width = 1
        )
    )
    
    # Split promoters into forward and reverse
    targets <- deconvolveBidirectionalPromoters(granges)
    targets_plus <- targets[GenomicRanges::strand(targets) == '+']
    targets_minus <- targets[GenomicRanges::strand(targets) == '-']
    
    # Forward promoters
    if (length(targets_plus)) {
        # Resize targets 
        targets_centers <- GenomicRanges::resize(
            targets_plus, fix = 'center', width = 1
        )
        targets_extended <- GenomicRanges::resize(
            targets_centers, fix = 'center', width = (xlims[2] - xlims[1])
        )
        # Get oriented distance between frag center and target center
        filt <- seqnames(
            frags_subset_centers
        ) %in% as.character(unique(seqnames(targets_centers)))
        n <- nearest(
            frags_subset_centers[filt],
            targets_centers
        )
        nearest_targets <- targets_centers[n]
        frags_subset$distance_to_nearest <- NA
        frags_subset[filt]$distance_to_nearest <- GenomicRanges::start(
            frags_subset_centers[filt]
        ) - 
            GenomicRanges::start(nearest_targets)
        frags_subset$strand_nearest_target <- NA
        frags_subset[filt]$strand_nearest_target <- GenomicRanges::strand(
            nearest_targets
        )
        frags_subset$distance_to_nearest <- ifelse(
            frags_subset$strand_nearest_target == '-', 
            -frags_subset$distance_to_nearest, 
            frags_subset$distance_to_nearest
        )
        # dists / widths
        frags_subset_sub <- frags_subset[!is.na(
            frags_subset$distance_to_nearest
        )]
        dists_plus <- frags_subset_sub$distance_to_nearest
        widths_plus <- IRanges::width(frags_subset_sub)
    } 
    else {
        dists_plus <- NULL
        widths_plus <- NULL
    }
    
    # Reverse promoters
    if (length(targets_minus)) {
        # Resize targets 
        targets_centers <- GenomicRanges::resize(
            targets_minus, fix = 'center', width = 1
        )
        targets_extended <- GenomicRanges::resize(
            targets_centers, fix = 'center', width = (xlims[2] - xlims[1])
        )
        # Get oriented distance between frag center and target center
        filt <- seqnames(
            frags_subset_centers
        ) %in% as.character(unique(seqnames(targets_centers)))
        nearest_targets <- targets_centers[
            nearest(frags_subset_centers[filt], targets_centers)
        ]
        frags_subset$distance_to_nearest <- NA
        frags_subset[filt]$distance_to_nearest <- GenomicRanges::start(
            frags_subset_centers[filt]
        ) - 
            GenomicRanges::start(nearest_targets)
        frags_subset$strand_nearest_target <- NA
        frags_subset[filt]$strand_nearest_target <- GenomicRanges::strand(
            nearest_targets
        )
        frags_subset$distance_to_nearest <- ifelse(
            frags_subset$strand_nearest_target == '-', 
            -frags_subset$distance_to_nearest, 
            frags_subset$distance_to_nearest
        )
        # dists / widths
        frags_subset_sub <- frags_subset[!is.na(
            frags_subset$distance_to_nearest
        )]
        dists_minus <- frags_subset_sub$distance_to_nearest
        widths_minus <- IRanges::width(frags_subset_sub)
    }
    else {
        dists_minus <- NULL
        widths_minus <- NULL
    }
    
    # Return table of dists / widths
    tab <- table(
        factor(
            c(dists_plus, dists_minus), 
            levels = xlims[1]:xlims[2]
        ), 
        factor(
            c(widths_plus, widths_minus), 
            levels = c(ylims[1]:ylims[2])
        )
    )
    
    # Replace central Vplot line (dist = 0) by average
    # between two flanking lines
    tab["0", ] = round((tab["-1", ] + tab["1", ])/2)
    
    return(tab)
}

#' A function to shuffle a Vmat
#' 
#' This function works on a Vmat (the output of computeVmat()). It shuffles
#' the matrix to randomize the fragment densities.
#'
#' @param Vmat A Vmat, usually output of computeVmat
#' @return A shuffled Vmat object
#' 
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_all_REs)
#' Vmat <- computeVmat(bam_test, ce11_all_REs)
#' Vmat <- shuffleVmat(Vmat)

shuffleVmat <- function(Vmat) {
    shuffled.Vmat <- c(Vmat)
    shuffled.Vmat <- shuffled.Vmat[sample(
        x = seq_along(shuffled.Vmat), 
        size = length(shuffled.Vmat), 
        replace = FALSE
    )]
    shuffled.Vmat <- matrix(
        shuffled.Vmat, nrow = nrow(Vmat), ncol = ncol(Vmat)
    )
    colnames(shuffled.Vmat) <- colnames(Vmat)
    row.names(shuffled.Vmat) <- row.names(Vmat)
    return(shuffled.Vmat)
}

#' A function to normalized a Vmat
#' 
#' This function normalizes a Vmat. Several different approaches have 
#' been implemented to normalize the Vmats. 
#'
#' @param Vmat A Vmat, usually output of computeVmat
#' @param bam_granges GRanges, the paired-end fragments
#' @param granges GRanges, the regions to map the fragments onto
#' @param normFun character. A Vmat should be scaled either by:
#' \itemize{
#'     \item 'libdepth+nloci', e.g. the library depth and the number of 
#'     loci used to compute the Vmat;
#'     \item zscore, if relative patterns of fragment density 
#'     are more important than density per se;
#'     \item Alternatively, the Vmat can be scaled to % ('pct'), to 
#'     a chosen quantile ('quantile') or to the max Vmat value ('max').
#' }
#' @param s A float indicating which quantile to use if 'quantile'
#' normalization is chosen
#' @param roll integer, to use as the window to smooth the Vmat rows 
#' by rolling mean.
#' @param verbose Boolean 
#' @return A normalized Vmat object
#' 
#' @import S4Vectors
#' @import GenomicRanges
#' @import IRanges
#' @importFrom zoo rollmean
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_all_REs)
#' Vmat <- computeVmat(bam_test, ce11_all_REs)
#' Vmat <- normalizeVmat(
#'     Vmat, 
#'     bam_test, 
#'     ce11_all_REs,
#'     normFun = c('libdepth+nloci')
#' )

normalizeVmat <- function(
    Vmat, 
    bam_granges, 
    granges,
    normFun = c('zscore'),
    s = 0.99, 
    roll = 1,
    verbose = TRUE
)
{
    ## Normalize the matrix 
    if (normFun == 'libdepth+nloci') {
        if (verbose) message('Computing raw library depth')
        Vmat <- Vmat / length(bam_granges) * 1000000
        if (verbose) message('Dividing Vmat by its number of loci')
        Vmat <- Vmat / length(granges)
    } 
    else if (normFun == 'pct') {
        Vmat <- Vmat/sum(Vmat)*100
    }
    else if (normFun == 'quantile') {
        # Cap by quantile
        q <- quantile(Vmat, probs = s)
        Vmat[Vmat >= q] <- q
        # Then scale so that max value is one
        Vmat <- Vmat/q*100
    }
    else if (normFun == 'max') {
        Vmat <- Vmat/max(Vmat)
    }
    else if (normFun == 'zscore') {
        Vmat2 <- matrix(scale(c(Vmat)), byrow = FALSE, nrow = nrow(Vmat))
        colnames(Vmat2) <- colnames(Vmat)
        row.names(Vmat2) <- row.names(Vmat)
        Vmat <- Vmat2
    }
    else if (normFun %in% c('', 'none', 'skip')) {
        if (verbose) message('No normalization applied')
        Vmat <- Vmat
    }
    else {
        if (verbose) message('CAUTION: Normalization function 
        not found. Skipped normalization')
    }
    ## Smooth the heatmap
    if (roll > 1) {
        if (verbose) message('Smoothing the matrix')
        Vmat_rolled <- zoo::rollmean(Vmat, roll)
        row.names(Vmat_rolled) <- zoo::rollmean(
            as.numeric(rownames(Vmat)), 
            roll
        )
        Vmat <- Vmat_rolled
    }
    return(Vmat)
}

