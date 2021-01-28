context("test-nucleosome-enrich")

test_that("nucleosome-enrich works", {
    expect_equal({
        data(bam_test)
        data(ce11_proms)
        # 
        frags <- bam_test
        V <- plotVmat(
            frags,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
            normFun = '',
            return_Vmat = TRUE
        )
        V_bg <-  plotVmat(
            frags,
            sampleGRanges(ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], n = 100000), 
            normFun = '',
            return_Vmat = TRUE
        )
        nuc2 <- nucleosomeEnrichment(V, V_bg)
        #
        nuc3 <- nucleosomeEnrichment(
            frags,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.']
        )
        #
        unlink('tmp.pdf')
        TRUE
    }, TRUE)
})

