#' A fragsize function to compute fragments sizes
#'
#' \code{get.distr(bam.list)} returns the number of fragments for each fragment 
#'   length between \code{limits[1]} and \code{limits[2]} bp. 
#'
#' @param bam.list Output of import.bam function. 
#' @param granges A GRanges object. The genomic ranges to which the fragments 
#'   should overlap in order to be considered when estimating fragment sizes. 
#' @param extend.granges A numeric vector of length 2. How the granges should
#'   be extended.
#' @param limits A numeric vector of length 2. Only consider fragments within
#'   this window of sizes.
#' @return A list of tbl, one for each .bam file.
#' 

getFragmentsDistribution <- function(fragments, list.granges, extend.granges = c(-2000, 2000), limits = c(50, 650)) {
    if(any(class(list.granges) == 'GRanges'))
        list.granges <- list(list.granges)
    list.res <- parallel::mclapply(seq_along(list.granges), function(K) {
        extended.granges <- list.granges[[K]] %>% 
            GenomicRanges::resize(width = IRanges::width(.) - extend.granges[1], fix = 'end') %>%
            GenomicRanges::resize(width = IRanges::width(.) - extend.granges[1] + extend.granges[2], fix = 'start')
        counts <- fragments %>% 
            IRanges::subsetByOverlaps(extended.granges) %>% 
            IRanges::width() %>% 
            '['(. >= limits[1] & . <= limits[2]) %>%
            hist(plot = F, breaks = seq(1, limits[2], length.out = limits[2] + 1)) %>% 
            '$'(counts)
        res <- tibble::as_tibble(data.frame('sample' = attributes(fragments)$names, 'REs' = factor(names(list.granges))[[K]], 'x' = 1:limits[2], 'y' = counts))
        return(res)
    }, mc.cores = length(list.granges)) %>% do.call(rbind, .)
    if (length(list.res) == 1)
        list.res <- list.res[[1]]
    return(list.res)
}

#' A fragsize function to plot fragments sizes
#' 
#' @param distr.list Output of get.distr. A list of tbl, one for each .bam file.
#' @return A list of ggplots

plotFragmentsDistribution <- function(sizes) {
    if(any(class(sizes) == 'data.frame'))
        sizes <- list(sizes)
    list <- lapply(seq_along(sizes), function(K) {
        distr <- sizes[[K]]
        ggplot2::ggplot(distr, ggplot2::aes(x = x, y = y, color = REs)) +
            ggplot2::geom_line() +
            ggplot2::theme_classic() +
            ggplot2::labs(
                title = names(sizes)[[K]], 
                x = 'Fragment size',
                y = '# of fragments',
                color = 'REs classes'
            )
    })
    if (length(list) == 1)
        list <- list[[1]]
    return(list)
}
