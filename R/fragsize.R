#' A function to compute sizes distribution for PE fragments 
#'
#' @param fragments GRanges object containing PE fragments. See importPEBamFiles
#' for more details.
#' @param granges GRanges. cam be a list of different sets of GRanges.
#' should overlap in order to be considered when estimating fragment sizes. 
#' @param extend_granges A numeric vector of length 2. How the GRanges should
#' be extended.
#' @param limits A numeric vector of length 2. Only consider fragments within
#' this window of sizes.
#' 
#' @return A list of tbl, one for each .bam file.
#' 
#' @import parallel
#' @import IRanges
#' @import GenomicRanges
#' @export

getFragmentsDistribution <- function(fragments, granges_list = NULL, extend_granges = c(-500, 500), limits = c(0, 600)) {
    if (any(class(granges_list) == 'GRanges')) {
        granges_list <- list(granges_list)
        names(granges_list) <- 'granges'
    }
    df_res <- parallel::mclapply(seq_along(granges_list), function(K) {
        extended_granges <- granges_list[[K]] %>% 
            GenomicRanges::resize(width = IRanges::width(.) - extend_granges[1], fix = 'end') %>%
            GenomicRanges::resize(width = IRanges::width(.) - extend_granges[1] + extend_granges[2], fix = 'start')
        counts <- fragments %>% 
            IRanges::subsetByOverlaps(extended_granges) %>% 
            IRanges::width() %>% 
            '['(. >= limits[1] & . <= limits[2]) %>%
            hist(plot = F, breaks = seq(1, limits[2], length.out = limits[2] + 1)) %>% 
            '$'(counts)
        res <- data.frame(
            'REs' = factor(names(granges_list))[[K]], 
            'x' = 1:limits[2],
            'y' = counts
        )
        return(res)
    }, mc.cores = length(granges_list)) %>% 
        setNames(names(granges_list)) %>% 
        namedListToLongFormat()
    return(df_res)
}

