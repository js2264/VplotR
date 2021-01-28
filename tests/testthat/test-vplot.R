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
        data(bam_test)
        data(ce11_proms)
        # Size distr
        d <- getFragmentsDistribution(
            bam_test, 
            list(
                'Ub' = ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
                'M' = ce11_proms[ce11_proms$which.tissues == 'Muscle']
            )
        )
        d <- getFragmentsDistribution(
            bam_test, 
            ce11_proms[ce11_proms$which.tissues == 'Muscle']
        )
        # Test combined Vplot
        V <- plotVmat(
            bam_test,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
            normFun = 'libdepth+nloci',
            return_Vmat = TRUE
        )
        V <- shuffleVmat(V)
        p1 <- plotVmat(V)
        ggsave('tmp.pdf')
        p2 <- plotVmat(
            bam_test,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
            normFun = 'quantile'
        )
        ggsave('tmp.pdf')
        p2 <- plotVmat(
            bam_test,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'],
            normFun = 'zscore'
        )
        ggsave('tmp.pdf')
        list_params <- list(
            "Germline ATAC-seq over Ubiq. proms." = list(bam_test, ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Germline ATAC-seq over Germline proms." = list(bam_test, ce11_proms[ce11_proms$which.tissues == 'Germline']),
            "Neurons ATAC-seq over Ubiq. proms." = list(bam_test, ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "Neurons ATAC-seq over Neurons proms." = list(bam_test, ce11_proms[ce11_proms$which.tissues == 'Neurons'])
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
            "Germline ATAC-seq over Ubiq. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Germline ATAC-seq over Germline proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Germline']),
            "Neurons ATAC-seq over Ubiq. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Neurons ATAC-seq over Neurons proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Neurons']),
            "Muscle ATAC-seq over Ubiq. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Muscle ATAC-seq over Muscle proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Muscle']),
            "Hypod. ATAC-seq over Ubiq. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Hypod. ATAC-seq over Hypod. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Hypod.']),
            "Intest. ATAC-seq over Ubiq. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Ubiq.']),
            "Intest. ATAC-seq over Intest. proms." = list(bam_test, ce_TSSs[ce_TSSs$which.tissues == 'Intest.'])
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
            bam_test,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
        )
        ggsave('tmp.pdf')
        p <- plotFootprint(
            bam_test,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
            split_strand = TRUE
        )
        ggsave('tmp.pdf')
        #
        unlink('tmp.pdf')
        TRUE
    }, TRUE)
})

test_that("profile", {
    expect_equal({
        data(bam_test)
        data(ce11_proms)
        # Size distr
        V <- plotProfile(
            bam_test,
            'chrI:10000-12000',
            loci = ce11_proms,
            annots = ce11_proms,
            min = 80, 
            max = 200
        )
        ggsave('tmp.pdf')
        #
        unlink('tmp.pdf')
        TRUE
    }, TRUE)
})
