#' A function to import paired end bam files as GRanges
#' 
#' This function takes bam file paths and read them into GRanges 
#' objects.
#' Note: Can be quite lengthy for .bam files with 5+ millions fragments. 
#'
#' @param files character vector, each element of the vector is the path 
#' of an individual .bam file.
#' @param where GRanges, only import the fragments mapping to the 
#' input GRanges (can fasten the import process a lot). 
#' @param max_insert_size Integer, filter out fragments larger 
#' than this size.
#' @param shift_ATAC_fragments Boolean, if the fragments come 
#' from ATAC-seq, one might want to shift the extremities by +5 / -4 bp. 
#' @param verbose Boolean
#' @return A GRanges object containing fragments from the input .bam file. 
#' 
#' @import parallel
#' @import Rsamtools
#' @import GenomicRanges
#' @import IRanges
#' @import GenomicAlignments
#' @export

importPEBamFiles <- function(
    files, 
    where = NULL, 
    max_insert_size = 1000, 
    shift_ATAC_fragments = FALSE, 
    verbose = TRUE
)
{
    if (any(!file.exists(files))) stop('
        (Some) files are missing. Aborting.'
    )
    if (any(!file.exists(paste0(files, '.bai')))) stop('
        Bam index not found. 
        Please create one with samtools index. Aborting.'
    )
    list.bam <- parallel::mclapply(files, function(FILE) {
        bam <- Rsamtools::BamFile(
            FILE, yieldSize = 500000000, asMates = TRUE
        )
        if (is.null(where)) {
            where <- Rsamtools::scanBamHeader(
                bam, what = c("targets")
            )$targets
            where <- GenomicRanges::GRanges(
                seqnames = names(where), IRanges::IRanges(1, where)
            )
        } else {
            where <- GenomicRanges::reduce(
                IRanges::resize(
                    where, 
                    width = IRanges::width(where) + 1000, 
                    fix = 'center'
                )
            )
        }
        params <- Rsamtools::ScanBamParam(
            which = where, 
            what = c("qname", "mapq", "isize", "flag"), 
            flag = Rsamtools::scanBamFlag(
                isPaired = TRUE, 
                isProperPair = TRUE, 
                isSecondaryAlignment = FALSE,
                isUnmappedQuery = FALSE,
                isNotPassingQualityControls = FALSE,
                isSupplementaryAlignment = FALSE
            )
        )
        if (verbose) message('> Importing ', FILE, ' ...')
        a <- GenomicAlignments::readGAlignmentPairs(bam, param = params)
        if (verbose) message('> Filtering ', FILE, ' ...')
        g <- GenomicAlignments::granges(a)
        # Filter by insert size
        if (max_insert_size > 0 & !is.null(max_insert_size)) 
            g <- g[IRanges::width(g) <= max_insert_size]
        # Shift ATAC fragments
        if (shift_ATAC_fragments) {
            if (verbose) message('> Shifting ', FILE, ' ...')
            g <- shiftATACGranges(g)
        }
        # Return GRanges
        if (verbose) message('> ', FILE, ' import completed.')
        return(g)
    }, mc.cores = length(files))
    if (length(files) == 1) list.bam <- list.bam[[1]]
    return(list.bam)
}

#' A function to shift GRanges fragments by 5/-4. This is useful 
#' when dealing with fragments coming from ATAC-seq. 
#'
#' @param g GRanges of ATAC-seq fragments
#' @param pos_shift Integer. How many bases should fragments on 
#' direct strand be shifted by?
#' @param neg_shift Integer. How many bases should fragments on 
#' negative strand be shifted by?
#' @return A GRanges object containing fragments from the input 
#' .bam file. 
#' 
#' @import parallel
#' @import Rsamtools
#' @import GenomicRanges
#' @import IRanges
#' @import GenomicAlignments
#' @export

shiftATACGranges <- function(g, pos_shift = 5, neg_shift = 4) {
    if (any(GenomicRanges::strand(g) == '*')) stop(
        'Error: the input GRanges has some unstranded fragments. 
        Please read bam files using 
        importPEBamFiles with shift_ATAC_fragments = TRUE.'
    )
    IRanges::shift(
        g, ifelse(GenomicRanges::strand(g) == '+', pos_shift, neg_shift)
    )
}