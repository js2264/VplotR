context("test-import")

test_that("import bam and shift fragments", {
    skip('skip')
    expect_equal({
        PROJECT_PATH <- '~/Rpackages/VplotR/'
        setwd(PROJECT_PATH)
        require(devtools)
        require(tidyverse)
        require(magrittr)
        require(rtracklayer)
        load_all()
        #
        bam_ATAC_hg38 <- readRDS('~/ATAC_hg38_Corces2017.rds')
        bam_ATAC_hg38_merged <- unlist(GRangesList(bam_ATAC_hg38))
        data(CTCF_hg38)
        g <- importPEBamFiles(
            '~/20190718_ATAC-hPGCs/_bam-files/SRR5427806.map_pe_hg38^rm_chrM^rm_blacklist^q10.bam', 
            shift_ATAC_fragments = TRUE
        )
        #
        frags <- bam_ATAC_hg38_merged_shifted
        frags <- g
        data(CTCF_hg38)
        p <- plotFootprint(
            frags,
            CTCF_hg38
        )
        ggsave('footprint_CTCF_sh+1.pdf', width = 7, height = 3)
        #
        #
        bam_MNase_sacCer3 <- importPEBamFiles(
            '~/20200429_yeast_MNase/_bam-files/SRR3193260.map_pe_sacCer3^rm_chrM^rm_blacklist^q10.bam', 
            shift_ATAC_fragments = FALSE
        )
        data(REB1_sacCer3)
        p <- plotFootprint(
            bam_MNase_sacCer3,
            REB1_sacCer3
        )
        ggsave('REB1_footprint.pdf', width = 7, height = 3)
        #
        #
        TRUE
    }, TRUE)
})
