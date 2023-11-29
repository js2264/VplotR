#' A function to compute sizes distribution for paired-end fragments 
#' 
#' This function takes fragments and compute the distribution of their
#' sizes over a set or multiple sets of GRanges. 
#'
#' @param fragments GRanges object containing paired-end fragments.
#' See importPEBamFiles for more details on how to create such object.
#' @param granges_list GRanges, can be a list of different sets of GRanges.
#' @param extend_granges numeric vector of length 2, how the GRanges 
#' should be extended.
#' @param limits numeric vector of length 2, only consider 
#' fragments within this window of sizes.
#' @param roll Integer, apply a moving average of this size 
#' @param cores Integer, number of threads used to compute 
#' fragment size distribution
#' @return A list of tbl, one for each .bam file.
#' 
#' @import parallel
#' @import IRanges
#' @import GenomicRanges
#' @importFrom magrittr `%>%`
#' @importFrom zoo rollmean
#' @importFrom methods is
#' @importFrom graphics hist
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' df <- getFragmentsDistribution(
#'     bam_test, 
#'     ce11_proms,
#'     extend_granges = c(-500, 500)
#' )
#' head(df)
#' which.max(df$y)

getFragmentsDistribution <- function(
    fragments, 
    granges_list = NULL, 
    extend_granges = c(-500, 500), 
    limits = c(0, 600), 
    roll = 3, 
    cores = 1
) 
{
    `%>%` <- magrittr::`%>%` 
    if (is.null(granges_list)) {
        granges_list <- GenomicRanges::GRanges(
            GenomeInfoDb::seqinfo(fragments)
        )
    }
    if (methods::is(granges_list, 'GRanges')) {
        granges_list <- list(granges_list)
        names(granges_list) <- 'granges'
    }
    df_res <- parallel::mclapply(seq_along(granges_list), function(K) {
        extended_granges <- granges_list[[K]] %>% 
            GenomicRanges::resize(
                width = IRanges::width(.) - extend_granges[1], 
                fix = 'end'
            ) %>%
            GenomicRanges::resize(
                width = IRanges::width(.) - 
                    extend_granges[1] + 
                    extend_granges[2],
                fix = 'start'
            )
        counts <- fragments %>%
            IRanges::subsetByOverlaps(extended_granges) %>% 
            IRanges::width() %>% 
            '['(. >= limits[1] & . <= limits[2]) %>%
            graphics::hist(
                plot = FALSE, 
                breaks = seq(1, limits[2], length.out = limits[2] + 1)
            ) %>% 
            '$'(counts) %>% 
            zoo::rollmean(roll, na.pad = TRUE, 'center')
        res <- data.frame(
            'class' = factor(names(granges_list))[[K]], 
            'x' = seq_len(limits[2]),
            'y' = counts
        )
        return(res)
    }, mc.cores = cores)
    #
    df_res <- do.call(rbind, df_res)
    return(df_res)
}

