setColNames <- function(df, names) {
    colnames(df) <- names
    return(df)
}

#' A function to easily coerce a named list into a long data.frame
#'
#' \code{namedListToLongFormat(x)} returns a data.frame in long format, with 
#' an added 'name' column, containing the names of the input list.  
#'
#' @param x A named list vector.
#' 
#' @export
#' @import magrittr
#' 
#' @return A long data frame

namedListToLongFormat <- function(x) {
    lapply(names(x), function(NAME) {
        L <- x[[NAME]]
        if (is.null(ncol(L))) {
            data.frame(value = L, name = rep(NAME, length(L)))
        } 
        else {
            data.frame(L, name = rep(NAME, nrow(L)))
        }
    }) %>% do.call(rbind, .)
}


#### ---- GRanges utils ---- ####

deconvolveBidirectionalPromoters <- function(granges) {
    unid <- granges[GenomicRanges::strand(granges) == '+' | GenomicRanges::strand(granges) == '-']
    bid <- granges[GenomicRanges::strand(granges) == '*']
    bid.fwd <- bid
    GenomicRanges::strand(bid.fwd) <- '+'
    bid.rev <- bid
    GenomicRanges::strand(bid.rev) <- '-'
    granges_shifted <- sort(c(unid, bid.fwd, bid.rev), ignore.strand = T)
    return(granges_shifted)
}

AlignToTSS <- function(granges, upstream, downstream) {
    if (any(GenomicRanges::strand(granges) == '*')) {
        granges <- deconvolveBidirectionalPromoters(granges)
    }
    if (is.null(granges$TSS.fwd)) {
        granges$TSS.fwd <- IRanges::start(granges)
        granges$TSS.rev <- IRanges::end(granges)
    }
    GenomicRanges::ranges(granges) <- IRanges::IRanges(
        start = ifelse(as.vector(GenomicRanges::strand(granges)) == '+', (granges$TSS.fwd - upstream), (granges$TSS.rev - downstream + 1)),
        width = downstream + upstream,
        names = names(IRanges::ranges(granges))
    )
    return(granges)
}

shuffleGRanges <- function(granges, opt.shuffle = '-chrom -noOverlapping', excl.itself = TRUE, genome = 'ce11', SEED = 222) {
    BEDTOOLS <- Sys.which('bedtools')
    tmp1 <- tempfile()
    tmp2 <- tempfile()
    tmp3 <- tempfile()
    tmp4 <- tempfile()

    write.table(getChromSizes(genome) %>% as.data.frame() %>% '['(c('seqnames', 'end')), tmp1, sep = '\t', col.names = F, row.names = F, quote = F)
    rtracklayer::export.bed(simplifyGRanges(granges), con = tmp2)
    if (excl.itself) {
        excl.granges <- simplifyGRanges(granges)
        rtracklayer::export.bed(excl.granges, con = tmp3)
        system(sprintf("%s shuffle %s -excl %s -seed %i -i %s -g %s | sort -k1,1 -k2.2n > %s", BEDTOOLS, opt.shuffle, tmp3, SEED, tmp2, tmp1, tmp4))
    } else {
        system(sprintf("%s shuffle %s -seed %i  -i %s -g %s | sort -k1,1 -k2.2n > %s", BEDTOOLS, opt.shuffle, SEED, tmp2, tmp1, tmp4))
    }
    t <- rtracklayer::import.bed(tmp4)
    system(sprintf('rm %s %s %s %s', tmp2, tmp4, tmp3, tmp1))
    return(sort(t))
}

simplifyGRanges <- function(granges) {
    GenomicRanges::mcols(granges) <- NULL
    return(granges)
}

#### ---- getChromSizes function ---- ####

getChromSizes <- function(x, ...) {
    UseMethod("getChromSizes")
}

getChromSizes.default <- function(genome = c('ce11', 'dm6', 'mm10', 'hg38')) {
    genome <- switch(genome, 
        'ce11' = BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11, 
        'dm6' = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6,
        'mm10' = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, 
        'hg38' = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    )
    chrom.sizes <- GenomeInfoDb::seqinfo(genome) %>% GenomicRanges::GRanges()
    return(chrom.sizes)
}

getChromSizes.GRanges <- function(granges) {
    message("Approximating chromosome sizes based on the coverage of the input GRanges object...")
    chrs <- levels(GenomicRanges::seqnames(granges))
    starts <- rep(1, length(chrs))
    ends <- sapply(chrs, function(chr) {
        granges[GenomicRanges::seqnames(granges) == chr] %>% 
            sort() %>% 
            tail(1) %>% 
            IRanges::end()
    })
    granges <- GenomicRanges::GRanges(seqnames = chrs, ranges = IRanges::IRanges(starts, width = ends))
    names(granges) <- GenomicRanges::seqnames(granges)
    return(granges)
}


