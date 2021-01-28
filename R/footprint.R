#' A function to plot footprint of paired-end data at given loci
#' 
#' This function takes paired-end fragments, extract the "cuts" (i.e.
#' extremities) and plot the footprint profile over a set of GRanges.
#'
#' @param frags GRanges, the paired-end fragments
#' @param targets GRanges, the loci to map the fragments onto
#' @param split_strand Boolean, should the + and - strand be splitted?
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
    split_strand = FALSE, 
    plot_central = TRUE, 
    xlim = c(-75, 75), 
    bin = 1,
    verbose = 1
)
{
    if (split_strand) {
        # get cut sites (i.e. extremities of fragments)
        if (verbose) {message('- Getting cuts')}
        cuts <- getCuts(frags)
        s <- GenomicRanges::strand(cuts)
        if (verbose) {message('- Getting cut coverage')}
        cov_cuts_for <- GenomicRanges::coverage(cuts[s == '+'])
        cov_cuts_for_boxed <- GenomicRanges::GRanges(
            names(cov_cuts_for), 
            IRanges(1, lengths(cov_cuts_for))
        )
        cov_cuts_rev <- GenomicRanges::coverage(cuts[s == '-'])
        cov_cuts_rev_boxed <- GenomicRanges::GRanges(
            names(cov_cuts_rev), 
            IRanges(1, lengths(cov_cuts_rev))
        )

        # Get the cuts coverage / target
        if (verbose) {message('- Getting cut coverage / target')}
        w <- xlim[2] - xlim[1] + 1
        targets_resized <- GenomicRanges::resize(
            targets, fix = 'center', width = w
        )
        targets_resized <- IRanges::subsetByOverlaps(
            targets_resized,
            GenomicRanges::intersect(cov_cuts_for_boxed, cov_cuts_rev_boxed), 
            type = 'within'
        )
        cov_targets_for <- cov_cuts_for[
            targets_resized[targets_resized %over% cov_cuts_for_boxed]
        ]
        cov_targets_rev <- cov_cuts_rev[
            targets_resized[targets_resized %over% cov_cuts_rev_boxed]
        ]
        
        # Turn it into a rectangular matrix and flip reverse strand scores
        if (verbose) {message('- Reformatting data into matrix')}
        cov_targets_for <- suppressWarnings(suppressMessages(
            matrix(
                as.vector(unlist(cov_targets_for)), 
                nrow = length(cov_targets_for),
                byrow = TRUE
            )
        ))
        cov_targets_rev <- suppressWarnings(suppressMessages(
            matrix(
                as.vector(unlist(cov_targets_rev)), 
                nrow = length(cov_targets_rev),
                byrow = TRUE
            )
        ))

        # Wrangle the data
        if (verbose) {message('- Plotting footprint')}
        mat <- cov_targets_for
        sums <- zoo::rollmean(colSums(mat), bin, fix = 'center', fill = NA)
        coords <- seq(
            -unique(GenomicRanges::width(targets_resized))/2 + 0.5, 
            unique(GenomicRanges::width(targets_resized))/2 - 0.5, 
            length.out = unique(GenomicRanges::width(targets_resized))
        )
        EE_for <- data.frame(
            x = coords,
            y = sums,
            strand = factor('+', levels = c('+', '-')),
            stringsAsFactors = FALSE
        )
        #
        mat <- cov_targets_rev
        sums <- zoo::rollmean(colSums(mat), bin, fix = 'center', fill = NA)
        coords <- rev(seq(
            -unique(GenomicRanges::width(targets_resized))/2 + 0.5, 
            unique(GenomicRanges::width(targets_resized))/2 - 0.5, 
            length.out = unique(GenomicRanges::width(targets_resized))
        ))
        EE_rev <- data.frame(
            x = coords,
            y = sums,
            strand = factor('-', levels = c('+', '-')),
            stringsAsFactors = FALSE
        )
        #
        EE <- rbind(EE_for, EE_rev)

        # Plot the footprint
        p <- ggplot2::ggplot(col = 'black', fill = 'black')
        p <- p + ggplot2::geom_line(
            data = EE,
            aes(
                x = x, 
                y = y, 
                col = strand
            )
        )
        # p <- p + ggplot2::geom_ribbon(alpha = 0.2, col = NA)
        p <- p + theme_ggplot2()
        p <- p + scale_color_manual(values = c('#f10000', '#0038ff'))
        p <- p + ggplot2::scale_y_continuous(
            limits = c(min(EE$y), max(EE$y))
        )
        p <- p + ggplot2::labs(
            x = 'Distance from center of motif', 
            y = '# of cuts'
        )
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
                    ymax = y2, 
                    fill = ''
                ),
                alpha = 0.1, 
                col = NA
            )
            p <- p + scale_fill_manual('Location of\nthe motif',
                values = '#000000', 
                guide = guide_legend(override.aes = list(alpha = 0.1))
            ) 
        }

        return(p)
    } 
    else {
        # get cut sites (i.e. extremities of fragments)
        if (verbose) {message('- Getting cuts')}
        cuts <- getCuts(frags)
        s <- GenomicRanges::strand(cuts)
        if (verbose) {message('- Getting cut coverage')}
        cov_cuts <- GenomicRanges::coverage(cuts[s == '+'])
        cov_cuts_boxed <- GenomicRanges::GRanges(
            names(cov_cuts), 
            IRanges(1, lengths(cov_cuts))
        )

        # Get the cuts coverage / target
        if (verbose) {message('- Getting cut coverage / target')}
        w <- xlim[2] - xlim[1] + 1
        targets_resized <- GenomicRanges::resize(
            targets, fix = 'center', width = w
        )
        targets_resized <- IRanges::subsetByOverlaps(
            targets_resized,
            cov_cuts_boxed,
            type = 'within'
        )
        cov_targets <- cov_cuts[
            targets_resized[targets_resized %over% cov_cuts_boxed]
        ]
        
        # Turn it into a rectangular matrix and flip reverse strand scores
        if (verbose) {message('- Reformatting data into matrix')}
        cov_targets <- suppressWarnings(suppressMessages(
            matrix(
                as.vector(unlist(cov_targets)), nrow = length(cov_targets),
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
        sums <- zoo::rollmean(colSums(mat), bin, fix = 'center', fill = NA)
        coords <- seq(
            -unique(GenomicRanges::width(targets_resized))/2 + 0.5, 
            unique(GenomicRanges::width(targets_resized))/2 - 0.5, 
            length.out = unique(GenomicRanges::width(targets_resized))
        )
        EE <- data.frame(
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
                y = y
            )
        )
        # p <- p + ggplot2::geom_ribbon(alpha = 0.2, col = NA)
        p <- p + theme_ggplot2()
        p <- p + ggplot2::scale_y_continuous(
            limits = c(min(EE$y), max(EE$y))
        )
        p <- p + ggplot2::labs(
            x = 'Distance from center of motif', 
            y = '# of cuts'
        )
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
                    ymax = y2, 
                    fill = ''
                ),
                alpha = 0.1, 
                col = NA
            )
            p <- p + scale_fill_manual('Location of\nthe motif',
                values = '#000000', 
                guide = guide_legend(override.aes = list(alpha = 0.1))
            ) 
        }
        
        return(p)
        
    }
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
#' @keywords internal

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