shuffleGRanges <- function(granges, opt.shuffle = '-chrom -noOverlapping', excl.itself = T, genome = 'ce11', SEED = 222) {
    BEDTOOLS <- Sys.which('bedtools')
    write.table(makeChromSizes(genome) %>% as.data.frame() %>% '['(c('seqnames', 'end')), "chrom.sizes", sep = '\t', col.names = F, row.names = F, quote = F)
    rtracklayer::export.bed(simplifyGRanges(granges), con = 'tmp.bed')
    if (excl.itself) {
        excl.granges <- simplifyGRanges(granges)
        rtracklayer::export.bed(excl.granges, con = 'excl.bed')
        system(sprintf("%s shuffle %s -excl excl.bed -seed %i -i %s -g %s | sort -k1,1 -k2.2n > tmp2.bed", BEDTOOLS, opt.shuffle, SEED, 'tmp.bed', "chrom.sizes"))
    } else {
        system(sprintf("%s shuffle %s -seed %i  -i %s -g %s | sort -k1,1 -k2.2n > tmp2.bed", BEDTOOLS, opt.shuffle, SEED, 'tmp.bed', "chrom.sizes"))
    }
    t <- rtracklayer::import.bed('tmp2.bed')
    system('rm tmp.bed tmp2.bed excl.bed chrom.sizes')
    return(sort(t))
}

simplifyGRanges <- function(granges) {
    GenomicRanges::mcols(granges) <- NULL
    return(granges)
}

makeChromSizes <- function(genome) {
    if (genome == 'ce11')
        genome <- BSgenome.Celegans.UCSC.ce11::BSgenome.Celegans.UCSC.ce11
    chrom.sizes <- GenomeInfoDb::seqinfo(genome) %>% GenomicRanges::GRanges()
    return(chrom.sizes)
}