#' A function to re-align a GRanges object to TSSs
#' 
#' This function re-aligns ranges (typically regulatory elements)
#' to a set of coordinates, either the TSS column or the
#' TSS.fwd and TSS.rev columns. If none are found, the function
#' assumes the ranges are promoters and that the end or the ranges
#' are the TSSs. 
#'
#' @param granges A stranded GRanges object with a TSS column 
#' or TSS.rev and TSS.fwd columns
#' @param upstream How many bases upstream of the TSS should the GRanges
#' object by extended by? [Default: 0]
#' @param downstream How many bases downstream of the TSS should the GRanges
#' object by extended by? [Default: 1]
#' @return GRanges aligned to the TSS column or to TSS.rev 
#' and TSS.fwd columns, and extended by upstream/downstream bp. 
#' 
#' @import GenomicRanges
#' @import IRanges
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' ce11_proms
#' alignToTSS(ce11_proms)

alignToTSS <- function(granges, upstream = 0, downstream = 1) {
    if (any(GenomicRanges::strand(granges) == '*')) {
        granges <- deconvolveBidirectionalPromoters(granges)
    }
    if (!is.null(granges$TSS)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(
                as.vector(GenomicRanges::strand(granges)) == '+',
                (granges$TSS - upstream), 
                (granges$TSS - downstream + 1)
            ),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else if (!is.null(granges$TSS.fwd) & !is.null(granges$TSS.rev)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(
                as.vector(GenomicRanges::strand(granges)) == '+', 
                (granges$TSS.fwd - upstream), 
                (granges$TSS.rev - downstream + 1)
            ),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else {
        TSSs <- GenomicRanges::start(
            GenomicRanges::resize(granges, fix = 'end', width = 1)
        )
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(
                as.vector(GenomicRanges::strand(granges)) == '+', 
                (TSSs - upstream), 
                (TSSs - downstream + 1)
            ),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    return(granges)
}

#' A function to duplicate bi-directional GRanges
#' 
#' This function splits bi-directional ranges into + and - 
#' stranded ranges. It duplicates the ranges which are '*'.
#'
#' @param granges A stranded GRanges object 
#' @return GRanges with only '+' and '-' strands. GRanges with '*' strand 
#' have been duplicated and split into forward and reverse strands.
#' 
#' @import GenomicRanges
#' @export
#' 
#' @examples
#' data(ce11_all_REs)
#' library(GenomicRanges)
#' proms <- ce11_all_REs[grepl('prom', ce11_all_REs$regulatory_class)]
#' proms
#' table(strand(proms))
#' proms <- deconvolveBidirectionalPromoters(proms)
#' proms
#' table(strand(proms))

deconvolveBidirectionalPromoters <- function(granges) {
    filt <- GenomicRanges::strand(granges) == '+' | 
        GenomicRanges::strand(granges) == '-'
    unid <- granges[filt]
    bid <- granges[GenomicRanges::strand(granges) == '*']
    bid.fwd <- bid
    GenomicRanges::strand(bid.fwd) <- '+'
    bid.rev <- bid
    GenomicRanges::strand(bid.rev) <- '-'
    granges_shifted <- sort(c(unid, bid.fwd, bid.rev), ignore.strand = TRUE)
    return(granges_shifted)
}

#' A function to sample GRanges from GRanges
#'
#' This function takes a given GRanges and returns another GRanges 
#' object. The new GRanges has the same number of ranges and the same
#' chromosome, width and strand distributions than the original 
#' GRanges. 
#'
#' @param x GRanges object
#' @param n Integer, number of sampled GRanges
#' @param width Integer, width of sampled GRanges
#' @param exclude Boolean, should the original GRanges be excluded?
#' @param avoid_overlap Boolean, should the 
#' sampled GRanges not be overlapping?
#' @return A GRanges object of length n
#' 
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' sampleGRanges(ce11_proms, 100)

sampleGRanges <- function(
    x, 
    n = NULL, 
    width = NULL,
    exclude = FALSE, 
    avoid_overlap = FALSE
)
{
    UseMethod("sampleGRanges")
}

#' A function to sample GRanges within GRanges
#'
#' This function takes a given GRanges and returns another GRanges 
#' object. The new GRanges has the same number of ranges and the same
#' chromosome, width and strand distributions than the original 
#' GRanges. 
#'
#' @param x GRanges object
#' @param n Integer, number of sampled GRanges
#' @param width Integer, width of sampled GRanges
#' @param exclude Boolean, should the original GRanges be excluded?
#' @param avoid_overlap Boolean, should the sampled GRanges 
#' not be overlapping?
#' @return A GRanges object of length n
#' 
#' @importFrom methods as
#' @import GenomicRanges
#' @import IRanges
#' @import GenomeInfoDb
#' @export
#' 
#' @examples
#' data(ce11_proms)
#' sampleGRanges(ce11_proms, 100)

sampleGRanges.GRanges <- function(
    x, 
    n = NULL, 
    width = NULL, 
    exclude = FALSE, 
    avoid_overlap = FALSE
)
{
    granges <- x
    maxed_granges <- GenomicRanges::GRanges(
        levels(GenomicRanges::seqnames(granges)), 
        IRanges::IRanges(
            1, 
            width = lengths(GenomicRanges::coverage(granges))
        )
    )
    if (is.null(n)) n <- length(x)
    N <- 2*n 
    g <- GRanges()
    #
    while (length(g) < n) {
        # Sample chrs based on their size
        chrs <- unlist(lapply(
            seq_along(IRanges::width(maxed_granges)), 
            function(K) {
                rep(
                    as.character(GenomicRanges::seqnames(maxed_granges)[K]), 
                    GenomicRanges::end(maxed_granges)[K]
                )
            }
        ))
        chrs <- sample(chrs, N)
        chrs <- as.character(sort(factor(
            chrs, 
            levels = levels(GenomicRanges::seqnames(maxed_granges))
        )))
        # For each chr, get random positions within this chr.
        pos <- unlist(lapply(
            unique(chrs), 
            function(chr) {
                sample(
                    lengths(GenomicRanges::coverage(granges)[chr]), 
                    table(chrs)[chr]
                )
            }
        ))
        # For each chr, get random strands within this chr.
        strands <- unlist(RleList(lapply(
            unique(chrs), 
            function(chr) {
                sample(
                    as.vector(GenomicRanges::strand(
                        granges[GenomicRanges::seqnames(granges) == chr]
                    )), 
                    table(chrs)[chr], 
                    replace = TRUE
                )
            }
        )))
        # For each chr, get widths within this chr.
        if (is.null(width)) {
            widths <- as.vector(unlist(RleList(lapply(
                unique(chrs), 
                function(chr) {
                    sample(
                        IRanges::width(
                            granges[GenomicRanges::seqnames(granges) == chr]
                        ), 
                        table(chrs)[chr], 
                        replace = TRUE
                    )
                }
            ))))
        } 
        else {
            widths <- width
        }
        # Build a GRanges object 
        suppressWarnings({
            g2 <- GenomicRanges::GRanges(
                chrs, 
                IRanges::IRanges(pos, width = widths),
                strand = strands, 
                seqinfo = GenomeInfoDb::seqinfo(granges)
            )
            GenomeInfoDb::seqlengths(g2) <- lengths(maxed_granges)
        })
        if (avoid_overlap) {
            g2 <- reduce(g2)
        }
        # Remove regions overlapping with initial GRanges
        if (exclude) {
            g2 <- g2[!(IRanges::`%over%`(g2, granges))]
        }
        # Remove the extended granges not embedded within the initial granges
        g2 <- IRanges::subsetByOverlaps(
            g2, 
            maxed_granges, 
            type = 'within'
        )
        # g2 <- GenomicRanges::trim(g2)
        if (length(unique(widths)) == 1) {
            g2 <- g2[GenomicRanges::width(g2) == unique(widths)]
        }
        g <- c(g, g2)
    }
    g <- g[sample(seq_len(length(g)), n)]
    g <- sort(g)
    seqlengths(g) <- NA
    return(g)
}

