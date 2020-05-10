#' A function to plot footprint of paired-end data at given loci
#' 
#' This function takes paired-end fragments, extract the "cuts" (i.e.
#' extremities) and plot the footprint profile over a set of GRanges.
#'
#' @param frags GRanges, the paired-end fragments
#' @param targets GRanges, the loci to map the fragments onto
#' @param plot_central plot grey rectangle over the loci
#' @param xlim numeric vector of length 2, the x limits of the computed Vmat
#' @param bin Integer, bin used to smooth the gootprint profile
#' @param verbose Integer
#' @return A footprint ggplot
#' 
#' @import ggplot2
#' @import IRanges
#' @import GenomicRanges
#' @importFrom methods as
#' @importFrom zoo rollmean
#' @importFrom stats qt
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' plotFootprint(bam_test, ce11_proms)

plotFootprint <- function(
    frags,
    targets, 
    plot_central = TRUE, 
    xlim = c(-50, 50), 
    bin = 1, 
    verbose = 1
)
{
    # get cut sites (i.e. extremities of fragments)
    if (verbose) {message('- Getting cuts')}
    cuts <- getCuts(frags)
    if (verbose) {message('- Getting cut coverage')}
    cov_cuts <- GenomicRanges::coverage(cuts)
    
    # Get the cuts coverage / target
    if (verbose) {message('- Getting cut coverage / target')}
    w <- xlim[2] - xlim[1]
    targets_resized <- GenomicRanges::resize(
        targets, fix = 'center', width = w
    )
    targets_resized <- IRanges::subsetByOverlaps(
        targets_resized,
        GenomicRanges::reduce(methods::as(cov_cuts, 'GRanges')), 
        type = 'within'
    )
    cov_targets <- cov_cuts[targets_resized]
    
    # Turn it into a rectangular matrix and flip reverse strand scores
    if (verbose) {message('- Reformatting data into matrix')}
    cov_targets <- suppressWarnings(suppressMessages(
        matrix(
            as.vector(unlist(cov_targets)), nrow = length(targets_resized),
            byrow = TRUE
        )
    ))
    cov_targets_flipped <- t(apply(cov_targets, 1, rev))
    cov_targets <- matrix(
        unlist(lapply(
            seq_len(nrow(cov_targets)), 
            function(K) {
                if(
                    (as.vector(
                        GenomicRanges::strand(targets_resized)
                    ) == '-')[K]
                ) {
                    cov_targets_flipped[K,]
                } else {
                    cov_targets[K,]
                }
            }
        )),
        nrow = length(targets_resized), 
        byrow = TRUE
    )
    
    # Wrangle the data
    if (verbose) {message('- Plotting footprint')}
    mat <- cov_targets
    colnames(mat) <- c(
        -unique(GenomicRanges::width(targets_resized))/2, 
        rep('', (ncol(mat) - 3)/2), 
        '', 
        rep('', (ncol(mat) - 2)/2), 
        paste0('+', unique(GenomicRanges::width(targets_resized))/2)
    )
    row.names(mat) <- paste0('locus_', as.character(seq_len(nrow(mat))))
    means <- c(
        rep(NA, bin/2), 
        zoo::rollmean(apply(mat, 2, mean), bin), 
        rep(NA, bin/2)
    )[seq_len(ncol(mat))]
    sums <- c(
        rep(NA, bin/2), 
        zoo::rollmean(colSums(mat), bin), 
        rep(NA, bin/2)
    )[seq_len(ncol(mat))]
    conint <- c(
        rep(NA, bin/2), 
        zoo::rollmean(
            apply(mat, 2, function (n) {
                Q <- stats::qt(0.975,sum(!is.na(n)))
                S <- sd(n,na.rm=TRUE)
                sq <- sqrt(sum(!is.na(n)))
                Q * S / sq
            }), 
            bin
        ), 
        rep(NA, bin/2)
    )[seq_len(ncol(mat))]
    topEE <- means + conint
    bottomEE <- means - conint
    coords <- round(
        seq(
            -unique(GenomicRanges::width(targets_resized))/2, 
            unique(GenomicRanges::width(targets_resized))/2, 
            length.out = unique(GenomicRanges::width(targets_resized))
        )
    )
    EE <- data.frame(
        means = means, 
        meansUp = topEE, 
        meansDown = bottomEE, 
        x = coords,
        y = sums,
        stringsAsFactors = FALSE
    )
    
    # Plot the footprint
    p <- ggplot2::ggplot(col = 'black', fill = 'black')
    p <- p + ggplot2::geom_line(
        data = EE,
        aes(
            x = x, 
            y = sums
        )
    )
    # p <- p + ggplot2::geom_ribbon(alpha = 0.2, col = NA)
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::scale_y_continuous(
        limits = c(min(EE$y), max(EE$y))
    )
    p <- p + ggplot2::labs(
        x = 'Distance from center of motif', 
        y = '# of cuts'
    )
    p <- p + ggplot2::guides(fill = FALSE)
    if (plot_central) {
        p <- p + ggplot2::geom_rect(
            data = data.frame(
                x1 = -width(targets) / 2,
                x2 = width(targets) / 2, 
                y1 = -Inf, 
                y2 = Inf
            )[1,],
            mapping = aes(
                xmin = x1, 
                xmax = x2, 
                ymin = y1, 
                ymax = y2
            ),
            alpha = 0.1, 
            col = NA, 
            fill = 'black'
        )
    }
    
    return(p)
}

#' Internal function
#' 
#' Function to extract cuts (i.e. extremities) of fragments stored
#' as GRanges.
#'
#' @param gr GRanges Paired-end fragments used to extract their extremities
#' @return GRanges
#' 
#' @import GenomicRanges
#' @import IRanges

getCuts <- function(gr) {
    g <- c(
        GenomicRanges::GRanges(
            seqnames = GenomicRanges::seqnames(gr), 
            IRanges::IRanges(start = GenomicRanges::start(gr), width = 1), 
            strand = GenomicRanges::strand(gr)
        ),
        GenomicRanges::GRanges(
            seqnames = GenomicRanges::seqnames(gr), 
            IRanges::IRanges(start = GenomicRanges::end(gr), width = 1), 
            strand = GenomicRanges::strand(gr)
        )
    )
    return(g)
}