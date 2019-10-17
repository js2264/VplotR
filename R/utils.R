#' A function to easily set column names for a data.frame
#'
#' \code{setColNames(df, names)} returns a data.frame with new column names.
#' This function can be convinently used in the tidyverse. 
#'
#' @param df A data frame.
#' @param names A character vector, length(names) == ncol(df).

#' @return A data frame with new column names
#' 
#' @export

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
#' @return A long data frame
#' 
#' @import magrittr
#' @export

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

#' A function to duplicate unstranded GRanges into '+' and '-' GRanges 
#'
#'
#' @param granges A GRanges object.
#' @param names A character vector, length(names) == ncol(df).
#' 
#' @return A data frame with new column names
#' 
#' @import GenomicRanges
#' @export

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

#' A function to re-align GRanges to their TSS
#'
#' @param granges A GRanges object with a TSS column or TSS.rev and TSS.fwd columns
#' @param upstream How many bases upstream of the TSS?
#' @param downstream How many bases downstream of the TSS?
#' 
#' @return GRanges aligned to the TSS column or to TSS.rev and TSS.fwd columns
#' 
#' @import GenomicRanges
#' @import IRanges
#' @export

alignToTSS <- function(granges, upstream, downstream) {
    if (any(GenomicRanges::strand(granges) == '*')) {
        granges <- deconvolveBidirectionalPromoters(granges)
    }
    if (!is.null(granges$TSS)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(as.vector(GenomicRanges::strand(granges)) == '+', (granges$TSS - upstream), (granges$TSS - downstream + 1)),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else if (!is.null(granges$TSS.fwd) & !is.null(granges$TSS.rev)) {
        GenomicRanges::ranges(granges) <- IRanges::IRanges(
            start = ifelse(as.vector(GenomicRanges::strand(granges)) == '+', (granges$TSS.fwd - upstream), (granges$TSS.rev - downstream + 1)),
            width = downstream + upstream,
            names = names(IRanges::ranges(granges))
        )
    }
    else {
        stop("No TSS column found. Aborting.")
    }
    return(granges)
}

#' A function to shuffle GRanges along a given genome
#' Internaly, this function relies on bedTools shuffle, so the bedTools
#' suite has to be installed on the workstation. 
#'
#' @param granges A GRanges object to shuffle
#' @param genome a BSgenome object. See getChromSizes for more details
#' @param opt.shuffle string of options (in single quotes) to pass to bedtools shuffle
#' @param exclude_itself Boolean. Should the shuffled GRanges overlap the input GRanges or not? 
#' @param SEED Integer. Make shuffling reproducible
#' 
#' @return a GRanges object with shuffled ranges
#' 
#' @import magrittr
#' @import GenomicRanges
#' @import IRanges
#' @export

shuffleGRanges <- function(granges, genome = NULL, opt.shuffle = '-chrom -noOverlapping', exclude_itself = TRUE, SEED = 222) {
    BEDTOOLS <- Sys.which('bedtools')
    tmp1 <- tempfile()
    tmp2 <- tempfile()
    tmp3 <- tempfile()
    tmp4 <- tempfile()
    if (is.null(genome)) genome <- granges
    getChromSizes(genome) %>% 
        as.data.frame() %>% 
        '['(c('seqnames', 'end')) %>% 
        write.table(tmp1, sep = '\t', col.names = F, row.names = F, quote = F)
    rtracklayer::export.bed(simplifyGRanges(granges), con = tmp2)
    if (exclude_itself) {
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

#' A function to remove metadata columns
#'
#' @param granges A GRanges object 
#' 
#' @return GRanges without any metadata columns
#' 
#' @import GenomicRanges
#' @export

simplifyGRanges <- function(granges) {
    GenomicRanges::mcols(granges) <- NULL
    return(granges)
}

#### ---- getChromSizes function ---- ####

#' A function to get / estimate the size of chromosomes for the input genome
#' 
#' @return GRanges of whole chromosomes for the input genome
#' 
#' @export

getChromSizes <- function(x, ...) {
    UseMethod("getChromSizes")
}

#' A function to get the size of chromosomes for the input genome
#' 
#' @param genome a BSgenome object
#' 
#' @return GRanges of whole chromosomes for the input genome
#' 
#' @import magrittr
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @export

getChromSizes.default <- function(genome = c('ce11', 'dm6', 'mm10', 'hg38', 'sacCer3', 'danRer10')) {
    if (is.character(genome)) {
        genome <- switch(genome, 
            'ce11' = BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11, 
            'dm6' = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6,
            'mm10' = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, 
            'hg38' = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
            'sacCer3' = BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3,
            'danRer10' = BSgenome.Drerio.UCSC.danRer10::BSgenome.Drerio.UCSC.danRer10
        )
    }
    chrom.sizes <- GenomeInfoDb::seqinfo(genome) %>% GenomicRanges::GRanges()
    return(chrom.sizes)
}

#' A function to estimate the size of chromosomes for the input genome
#' 
#' @param genome a GRanges object
#' 
#' @return GRanges of whole chromosomes for the input genome
#' 
#' @import magrittr
#' @import GenomicRanges
#' @import IRanges
#' @export

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


