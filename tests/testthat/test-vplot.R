context("test-vplot")

test_that("import bam", {
    expect_equal({
        data(ce11_proms)
        data(bam_test)
        # Test quick Vplot
        p <- plotVmat(
            bam_test,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'][2], 
            normFun = '', 
            colors = colorRampPalette(c('white', 'red'))(10), 
            breaks = c(0, 2)
        )
        ggplot2::ggsave('tmp.pdf')
        # Test shuffling
        p <- plotVmat(
            bam_test,
            sampleGRanges(
                ce11_proms[ce11_proms$which.tissues == 'Germline'][seq(1,2)]
            ), 
            normFun = '', 
            colors = colorRampPalette(c('white', 'red'))(10), 
            breaks = c(0, 1)
        )
        ggplot2::ggsave('tmp.pdf')
        #
        unlink('tmp.pdf')
        TRUE
    }, TRUE)
})

test_that("quick vplot", {
    expect_equal({
        bam_list <- readRDS(
            url('http://ahringerlab.com/VplotR/ATAC_PE_fragments.rds')
        )
        data(ce11_proms)
        # Size distr
        d <- getFragmentsDistribution(
            bam_list[['Muscle']], 
            list(
                'Ub' = ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
                'M' = ce11_proms[ce11_proms$which.tissues == 'Muscle']
            )
        )
        d <- getFragmentsDistribution(
            bam_list[['Muscle']], 
            ce11_proms[ce11_proms$which.tissues == 'Muscle']
        )
        # Test combined Vplot
        V <- plotVmat(
            bam_list[['Germline']][seq_len(1000000)],
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
            normFun = 'libdepth+nloci',
            return_Vmat = TRUE
        )
        V <- shuffleVmat(V)
        p1 <- plotVmat(V)
        ggsave('tmp.pdf')
        p2 <- plotVmat(
            bam_list[['mixed']][seq_len(1000000)],
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
            normFun = 'quantile'
        )
        ggsave('tmp.pdf')
        p2 <- plotVmat(
            bam_list[['mixed']][seq_len(1000000)],
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
            normFun = 'zscore'
        )
        ggsave('tmp.pdf')
        list_params <- list(
            "Germline ATAC-seq over Ubiq. proms." = list(bam_list[['Germline']][seq_len(1000000)], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Germline ATAC-seq over Germline proms." = list(bam_list[['Germline']][seq_len(1000000)], ce11_proms[ce11_proms$which.tissues == 'Germline']),
            "Neurons ATAC-seq over Ubiq. proms." = list(bam_list[['Neurons']][seq_len(1000000)], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Neurons ATAC-seq over Neurons proms." = list(bam_list[['Neurons']][seq_len(1000000)], ce11_proms[ce11_proms$which.tissues == 'Neurons'])
        )
        p3 <- plotVmat(
            list_params, 
            cores = 1,
            return_Vmat = TRUE,
            normFun = 'max'
        )
        p3 <- plotVmat(p3, ncol = 2, nrow = 2)
        ggsave('tmp.pdf')
        p4 <- plotVmat(
            list_params, 
            cores = 1,
            normFun = '',
            roll = 3, 
            nrow = 1, 
            ncol = 4
        )
        ggsave('tmp.pdf')
        #
        ce_TSSs <- alignToTSS(ce11_proms, 200, 200)
        list_params <- list(
            "Germline ATAC-seq over Ubiq. proms." = list(bam_list[['Germline']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Germline ATAC-seq over Germline proms." = list(bam_list[['Germline']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Germline']),
            "Neurons ATAC-seq over Ubiq. proms." = list(bam_list[['Neurons']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Neurons ATAC-seq over Neurons proms." = list(bam_list[['Neurons']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Neurons']),
            "Muscle ATAC-seq over Ubiq. proms." = list(bam_list[['Muscle']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Muscle ATAC-seq over Muscle proms." = list(bam_list[['Muscle']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Muscle']),
            "Hypod. ATAC-seq over Ubiq. proms." = list(bam_list[['Hypod.']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Hypod. ATAC-seq over Hypod. proms." = list(bam_list[['Hypod.']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Hypod.']),
            "Intest. ATAC-seq over Ubiq. proms." = list(bam_list[['Intest.']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Intest. ATAC-seq over Intest. proms." = list(bam_list[['Intest.']][seq_len(1000000)], ce_TSSs[ce_TSSs$which.tissues == 'Intest.'])
        )
        plots <- plotVmat(
            list_params, 
            cores = 1,
            normFun = 'libdepth+nloci', 
            nrow = 2, 
            ncol = 5
        )
        ggsave('tmp.pdf', height = 7.2, width = 18)
        #
        p <- plotFootprint(
            bam_list[['Germline']][seq_len(1000000)],
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
        )
        ggsave('tmp.pdf')
        #
        unlink('tmp.pdf')
        TRUE
    }, TRUE)
})

test_that("REB1 Yeast MNase Vplot", {
    skip('skip')
    expect_equal({
        bam_MNase_sacCer3 <- readRDS(
            url('http://ahringerlab.com/VplotR/MNase_sacCer3_Henikoff2011.rds')
        )
        data(REB1_sacCer3)
        #
        p <- plotVmat(
            bam_MNase_sacCer3,
            REB1_sacCer3, 
            roll = 3,
            ylim = c(25, 200),
            xlim = c(-600, 600)
        )
        ggsave('REB1_sacCer3_Vplot.pdf')
        #
        p <- plotFootprint(
            bam_MNase_sacCer3,
            REB1_sacCer3
        )
        ggsave('REB1_sacCer3_Footprint.pdf')
        #
        TRUE
    }, TRUE)
})

test_that("CTCF Human ATAC Vplot", {
    skip('skip')
    expect_equal({
        bam_ATAC_hg38 <- readRDS(
            url('http://ahringerlab.com/VplotR/ATAC_hg38_Chen2018.rds')
        )
        data(CTCF_hg38)
        #
        p <- plotVmat(
            bam_ATAC_hg38,
            CTCF_hg38, 
            roll = 3,
            ylims = c(80, 700), 
            xlims = c(-600, 600)
        )
        ggsave('CTCF_hg38_Vplot.pdf')
        #
        p <- plotFootprint(
            bam_ATAC_hg38,
            CTCF_hg38
        )
        ggsave('CTCF_hg38_Footprint.pdf')
        #
        TRUE
    }, TRUE)
})

test_that("TSSs Elegans ATAC Vplots", {
    skip('skip')
    expect_equal({
        bam_list <- readRDS(url('http://ahringerlab.com/VplotR/ATAC_PE_fragments.rds'))
        data(ce11_proms)
        ce11_proms <- alignToTSS(
            ce11_proms[strand(ce11_proms) != '*'],
            300, 
            300
        )
        list_params <- list(
            "Germline ATAC-seq over Ubiq. proms." = list(bam_list[['Germline']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Germline ATAC-seq over Germline proms." = list(bam_list[['Germline']], ce11_proms[ce11_proms$which.tissues == 'Germline']),
            "Neurons ATAC-seq over Ubiq. proms." = list(bam_list[['Neurons']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Neurons ATAC-seq over Neurons proms." = list(bam_list[['Neurons']], ce11_proms[ce11_proms$which.tissues == 'Neurons']),
            "Muscle ATAC-seq over Ubiq. proms." = list(bam_list[['Muscle']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Muscle ATAC-seq over Muscle proms." = list(bam_list[['Muscle']], ce11_proms[ce11_proms$which.tissues == 'Muscle']),
            "Hypod. ATAC-seq over Ubiq. proms." = list(bam_list[['Hypod.']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Hypod. ATAC-seq over Hypod. proms." = list(bam_list[['Hypod.']], ce11_proms[ce11_proms$which.tissues == 'Hypod.']),
            "Intest. ATAC-seq over Ubiq. proms." = list(bam_list[['Intest.']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Intest. ATAC-seq over Intest. proms." = list(bam_list[['Intest.']], ce11_proms[ce11_proms$which.tissues == 'Intest.'])
        )
        plots <- plotVmat(
            list_params, 
            cores = 10,
            normFun = 'zscore', 
            nrow = 2, 
            ncol = 5
        )
        ggsave('TSSs_ATAC_ce11_Serizay2020.pdf', height = 7.2, width = 18)
    })
})
