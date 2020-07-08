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

# Yeast
test_that("Yeast MNase ALL PLOTS", {
    skip('skip')
    expect_equal({
        bam_MNase_sacCer3 <- readRDS('~/MNase_sacCer3_Henikoff2011.rds')
        bam_MNase_sacCer3_merged <- unlist(GRangesList(bam_MNase_sacCer3))
        data(REB1_sacCer3)
        data(ABF1_sacCer3)
        #
        # ABF1 plot
        p <- plotVmat(
            bam_MNase_sacCer3_merged,
            ABF1_sacCer3, 
            ylim = c(25, 200),
            xlim = c(-500, 500)
        )
        ggsave('ABF1_sacCer3_Vplot.pdf')
        #
        # REB1 plot
        p <- plotVmat(
            bam_MNase_sacCer3_merged,
            REB1_sacCer3, 
            ylim = c(25, 200),
            xlim = c(-500, 500), 
        )
        ggsave('REB1_sacCer3_Vplot.pdf')
        #
        # 2 Vplots
        list_params <- list(
            'REB1' = list(bam_MNase_sacCer3_merged, REB1_sacCer3),
            'ABF1' = list(bam_MNase_sacCer3_merged, ABF1_sacCer3)
        )
        p <- plotVmat(
            list_params,
            cores = 2,
            ylim = c(25, 200),
            xlim = c(-500, 500), 
            roll = 3,
            nrow = 2, 
            ncol = 1,
            normFun = 'zscore'
        )
        ggsave('ABF1_REB1_sacCer3_Vplot.pdf')
        # Footprints
        p <- plotFootprint(
            bam_MNase_sacCer3_merged,
            REB1_sacCer3
        )
        ggsave('REB1_sacCer3_Footprint.pdf', width = 7, height = 3)
        #
        p <- plotFootprint(
            bam_MNase_sacCer3_merged,
            ABF1_sacCer3
        )
        ggsave('ABF1_sacCer3_Footprint.pdf', width = 7, height = 3)
        # Genomic profile
        genes_sacCer3 <- GenomicFeatures::genes(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::
            TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
        )
        loc <- "chrXV:186,000-187,800"
        loc <- "chrXV:186,500-187,300"
        p <- plotProfile(
            bam_MNase_sacCer3_merged, 
            window = loc, 
            loci = ABF1_sacCer3, 
            annots = genes_sacCer3,
            min = 20, max = 200, alpha = 0.1, size = 1.5
        )
        ggsave(paste0('sacCer3_', loc, '.pdf'), width = 12, height = 6, limitsize = FALSE)
        # Frag size
        dist <- getFragmentsDistribution(bam_MNase_sacCer3_merged, ABF1_sacCer3)
        ggplot(dist, aes(x = x, y = y)) +
            geom_line() +
            theme_ggplot2() +
            labs(
                title = 'Distribution of fragment sizes', 
                x = 'Fragment size', y = '# of fragments'
            )
        ggsave('ABF1_sacCer3_FragSize.pdf', width = 4, height = 4, limitsize = FALSE)
        #
        TRUE
    }, TRUE)
})

# Human
test_that("CTCF Human ATAC Vplot", {
    skip('skip')
    expect_equal({
        bam_ATAC_hg38 <- readRDS(
            url('http://ahringerlab.com/VplotR/ATAC_hg38_Corces2017.rds')
        )
        bam_ATAC_hg38 <- readRDS('~/ATAC_hg38_Chen2018.rds')
        bam_ATAC_hg38_merged <- unlist(GRangesList(bam_ATAC_hg38))
        data(CTCF_hg38)
        #
        p <- plotVmat(
            bam_ATAC_hg38_merged,
            CTCF_hg38, 
            roll = 3,
            ylims = c(80, 500), 
            xlims = c(-300, 300), 
            normFun = ''
        )
        ggsave('CTCF_hg38_Vplot_ATAC.pdf')
        #
        p2 <- plotFootprint(
            bam_ATAC_hg38_merged,
            CTCF_hg38
        )
        ggsave('CTCF_hg38_Footprint_Chen.pdf', width = 7, height = 3)
        #
        loc <- "chr3:50,117,898-50,126,073"
        p <- plotProfile(bam_ATAC_hg38_merged, loc, min = 20, max = 200, alpha = 0.6)
        ggsave(paste0(loc, '.pdf'), width = 30, height = 3, limitsize = FALSE)
        #
        TRUE
    }, TRUE)
})

test_that("CTCF Human MNase Vplot", {
    skip('skip')
    expect_equal({
        bam_MNase_hg38 <- readRDS(
            url('http://ahringerlab.com/VplotR/MNase_hg38_Mieczkowski2018.rds')
        )
        bam_MNase_hg38_merged <- unlist(GRangesList(bam_MNase_hg38))
        data(CTCF_hg38)
        #
        p <- plotVmat(
            bam_MNase_hg38_merged,
            CTCF_hg38, 
            roll = 3,
            ylims = c(40, 220), 
            xlims = c(-700, 700), 
            normFun = ''
        )
        ggsave('CTCF_hg38_Vplot_MNase.pdf')
        #
        loc <- "chr11:5,269,510-5,303,771"
        p <- plotProfile(bam_MNase_hg38_merged, loc, min = 20, max = 200, alpha = 1)
        ggsave(paste0(loc, '.pdf'), width = 100, height = 3, limitsize = FALSE)
        #
        TRUE
    }, TRUE)
})

test_that("CTCF Human ALL PLOTS", {
    skip('skip')
    expect_equal({
        bam_MNase_hg38 <- readRDS('~/MNase_hg38_Mieczkowski2018.rds')
        bam_DNase_hg38 <- readRDS('~/DNase_hg38_ENCODE.rds')
        bam_ATAC_hg38 <- readRDS('~/ATAC_hg38_Corces2017.rds')
        bam_MNase_hg38_merged <- unlist(GRangesList(bam_MNase_hg38))
        bam_DNase_hg38_merged <- bam_DNase_hg38
        bam_ATAC_hg38_merged <- unlist(GRangesList(bam_ATAC_hg38))
        data(CTCF_hg38)
        CTCF_hg38 <- CTCF_hg38[order(CTCF_hg38$relScore, decreasing = TRUE)][1:10000]
        #
        list_params <- list(
            'MNase' = list(bam_MNase_hg38_merged, CTCF_hg38, c(40, 300), c(-400, 400)), 
            'DNase' = list(bam_DNase_hg38_merged, CTCF_hg38, c(36, 150), c(-200, 200)), 
            'ATAC' = list(bam_ATAC_hg38_merged, CTCF_hg38, c(40, 300), c(-400, 400))
        )
        # 
        # p_vp <- plotVmat(bam_MNase_hg38_merged, CTCF_hg38, ylims = c(40, 220), xlims = c(-600, 600), hm = 90,roll = 3,normFun = 'none', verbose = TRUE)
        # ggsave('CTCF_hg38_MNase_Vplot.pdf', height = 4, width = 8)
        # p_fp <- plotFootprint(bam_DNase_hg38_merged, CTCF_hg38)
        # ggsave('CTCF_hg38_DNase_Footprint.pdf', height = 4, width = 8)
        # p_pr <- plotProfile(bam_MNase_hg38_merged, "chr11:118,926,707-118,929,707", loci = CTCF_hg38, min = 20, max = 260, alpha = 0.3, size = 2)
        # ggsave('CTCF_hg38_MNase_Profile.pdf', height = 4, width = 8)
        # 
        vplots <- parallel::mclapply(mc.cores = 3, list_params, function(params) {
            p <- plotVmat(
                params[[1]], params[[2]], 
                ylims = params[[3]], 
                hm = 90,
                xlims = params[[4]],
                cores = 3,
                roll = 3,
                nrow = 3, ncol = 1, 
                normFun = 'none', 
                verbose = TRUE
            )
        })
        p <- cowplot::plot_grid(plotlist = vplots, nrow = 3, labels = "AUTO")
        ggsave('CTCF_hg38_Vplots.pdf', height = 8, width = 8)
        #
        list_plots <- parallel::mclapply(list_params, function(param) {
            plotFootprint(
                param[[1]],
                param[[2]]
            )
        }, mc.cores = 3)
        p <- cowplot::plot_grid(plotlist = list_plots, nrow = 3, labels = "AUTO")
        ggsave('CTCF_hg38_Footprints.pdf', width = 7, height = 9)
        #
        loc <- "chr11:65,416,250-65,416,750"
        list_plots_profiles <- parallel::mclapply(list_params, function(param) {
            plotProfile(
                param[[1]],
                loc, 
                loci = CTCF_hg38,
                min = 20, max = 250, alpha = 0.3, size = 2
            )
        }, mc.cores = 3)
        p <- cowplot::plot_grid(plotlist = list_plots_profiles, nrow = 3, labels = "AUTO")
        ggsave(paste0('CTCF_hg38_', loc, '.pdf'), width = 8, height = 12, limitsize = FALSE)
        #
        ggMEME <- function(motif) {
            PWM <- motif@profileMatrix
            p <- ggplot() + 
                geom_logo(PWM) + 
                labs(y = '') + 
                theme_ggplot2() + 
                theme(
                    legend.position = 'none', 
                    axis.line.x = element_blank(),
                    axis.line.y = element_line(size = 0.25),
                    axis.text.x = element_text(colour="black"),
                    axis.text.y = element_text(colour="black"),
                    axis.text = element_text(size = 5)
                ) + 
                scale_x_continuous(
                    breaks = seq(1, ncol(PWM), 1),
                    labels = c(c('1', '', '', ''), unlist(lapply(seq(5, ncol(PWM), 5), function(C) c(C, rep('', 4)))[1:ncol(PWM)]))[1:ncol(PWM)],
                    expand = c(0, 0)
                ) + 
                scale_y_continuous(
                    breaks = c(0, 1, 2), 
                    limits = c(0, 23), 
                    expand = c(0, 0)
                )
            return(p)
        }
        CTCF <- toPWM(readJASPARMatrix('http://jaspar.genereg.net/api/v1/matrix/MA0139.1.jaspar'), type='prob')
        p <- ggMEME(CTCF)
        #
        TRUE
    }, TRUE)
})

# Mice
test_that("CTCF Mouse DNAse Vplot", {
    skip('skip')
    expect_equal({
        # bam_DNase_mm10 <- readRDS(url('http://ahringerlab.com/VplotR/DNase_mm10_ENCODE.rds'))
        bam_DNase_mm10 <- readRDS('~/DNase_mm10_ENCODE.rds')
        bam_DNase_mm10_merged <- unlist(GRangesList(bam_DNase_mm10))
        strand(bam_DNase_mm10_merged) <- '*'
        data(CTCF_mm10)
        strand(CTCF_mm10) <- '*'
        #
        p <- plotVmat(
            bam_DNase_mm10_merged,
            CTCF_mm10, 
            roll = 1,
            ylims = c(40, 150), 
            xlims = c(-150, 150)
        )
        ggsave('CTCF_mm10_Vplot.pdf')
        #
        p <- plotFootprint(
            bam_DNase_mm10_merged,
            CTCF_mm10
        )
        ggsave('CTCF_mm10_Footprint.pdf')
        #
        TRUE
    }, TRUE)
})

# Worm
test_that("TSSs Elegans ATAC Vplots", {
    skip('skip')
    expect_equal({
        bam_list <- readRDS('~/ATAC_ce11_Serizay2020.rds')
        #
        data(ce11_proms)
        ce11_proms <- alignToTSS(
            ce11_proms,
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
        ggsave('ce11_TSSs_ATAC_Serizay2020.pdf', height = 7.2, width = 18)
        #
        data(ce11_proms)
        ce11_proms <- resize(
            ce11_proms,
            fix = 'center', 
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
            roll = 3, 
            normFun = 'zscore', 
            nrow = 2, 
            ncol = 5, 
            verbose = 2
        )
        ggsave('ce11_proms_ATAC_Serizay2020.pdf', height = 7.2, width = 18)
        #
        data(ce11_all_REs)
        ce11_all_REs <- resize(
            ce11_all_REs,
            fix = 'center', 
            300
        )
        ce11_enhs <- ce11_all_REs[ce11_all_REs$regulatory_class == 'putative_enhancer']
        list_params <- list(
            "Germline ATAC-seq over Ubiq. proms." = list(bam_list[['Germline']], ce11_enhs[ce11_enhs$which.tissues == 'Ubiq.']),
            "Germline ATAC-seq over Germline proms." = list(bam_list[['Germline']], ce11_enhs[ce11_enhs$which.tissues == 'Germline']),
            "Neurons ATAC-seq over Ubiq. proms." = list(bam_list[['Neurons']], ce11_enhs[ce11_enhs$which.tissues == 'Ubiq.']),
            "Neurons ATAC-seq over Neurons proms." = list(bam_list[['Neurons']], ce11_enhs[ce11_enhs$which.tissues == 'Neurons']),
            "Muscle ATAC-seq over Ubiq. proms." = list(bam_list[['Muscle']], ce11_enhs[ce11_enhs$which.tissues == 'Ubiq.']),
            "Muscle ATAC-seq over Muscle proms." = list(bam_list[['Muscle']], ce11_enhs[ce11_enhs$which.tissues == 'Muscle']),
            "Hypod. ATAC-seq over Ubiq. proms." = list(bam_list[['Hypod.']], ce11_enhs[ce11_enhs$which.tissues == 'Ubiq.']),
            "Hypod. ATAC-seq over Hypod. proms." = list(bam_list[['Hypod.']], ce11_enhs[ce11_enhs$which.tissues == 'Hypod.']),
            "Intest. ATAC-seq over Ubiq. proms." = list(bam_list[['Intest.']], ce11_enhs[ce11_enhs$which.tissues == 'Ubiq.']),
            "Intest. ATAC-seq over Intest. proms." = list(bam_list[['Intest.']], ce11_enhs[ce11_enhs$which.tissues == 'Intest.'])
        )
        plots <- plotVmat(
            list_params, 
            cores = 10,
            roll = 3, 
            normFun = 'libdepth+nloci', 
            nrow = 2, 
            ncol = 5, 
            verbose = 2
        )
        ggsave('ce11_enhs_ATAC_Serizay2020.pdf', height = 7.2, width = 18)
    })
})

test_that("elegans MNase Vplots and profiles", {
    skip('skip')
    expect_equal({
        # bam_MNase_ce11 <- readRDS(url('http://ahringerlab.com/VplotR/MNase_ce11_Janes2018.rds'))
        bam_MNase_ce11 <- readRDS('~/MNase_ce11_Janes2018.rds')
        names(bam_MNase_ce11) <- paste0(names(bam_MNase_ce11) %>% gsub('rep-1_|..map_pe.*', '', .), 'U')
        bam_MNase_ce11_merged <- unlist(GRangesList(bam_MNase_ce11))
        #
        data(ce11_proms)
        ce11_proms <- resize(
            ce11_proms,
            fix = 'center', 
            300
        )
        list_params <- list(
            "MNase-seq 0U over Ubiq. proms." = list(bam_MNase_ce11[['0U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.05U over Ubiq. proms." = list(bam_MNase_ce11[['0.05U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.1U over Ubiq. proms." = list(bam_MNase_ce11[['0.1U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.25U over Ubiq. proms." = list(bam_MNase_ce11[['0.25U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.5U over Ubiq. proms." = list(bam_MNase_ce11[['0.5U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 1U over Ubiq. proms." = list(bam_MNase_ce11[['1U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 8U over Ubiq. proms." = list(bam_MNase_ce11[['8U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 16U over Ubiq. proms." = list(bam_MNase_ce11[['16U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.'])
        )
        plots <- plotVmat(
            list_params, 
            cores = 8,
            ylims = c(50, 200), 
            normFun = 'zscore', 
            nrow = 2, 
            ncol = 4, 
            verbose = 2
        )
        ggsave('ce11_proms_MNase_Serizay2020.pdf', height = 7.5, width = 14)
        #
        loc <- "chrIII:8,083,947-8,087,158"
        p <- plotProfile(bam_MNase_ce11_merged, loc, min = 20, max = 200, alpha = 0.1, size = 2)
        ggsave(paste0(loc, '.pdf'), width = 20, height = 8, limitsize = FALSE)
        #
        loc <- "chrIV:653,277-656,439"
        p <- plotProfile(bam_MNase_ce11_merged, loc, min = 20, max = 200, alpha = 0.1, size = 2)
        ggsave(paste0(loc, '.pdf'), width = 20, height = 8, limitsize = FALSE)
    })
})

test_that("ce11 embryonic promoters", {
    skip('skip')
    expect_equal({
        # bam_MNase_ce11 <- readRDS(url('http://ahringerlab.com/VplotR/MNase_ce11_Janes2018.rds'))
        bam_MNase_ce11 <- readRDS('~/MNase_ce11_Janes2018.rds')
        names(bam_MNase_ce11) <- paste0(names(bam_MNase_ce11) %>% gsub('rep-1_|..map_pe.*', '', .), 'U')
        bam_MNase_ce11_merged <- unlist(GRangesList(bam_MNase_ce11))
        #
        data(ce11_proms)
        ce11_proms <- resize(
            ce11_proms,
            fix = 'center', 
            300
        )
        list_params <- list(
            "MNase-seq 0U over Ubiq. proms." = list(bam_MNase_ce11[['0U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.05U over Ubiq. proms." = list(bam_MNase_ce11[['0.05U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.1U over Ubiq. proms." = list(bam_MNase_ce11[['0.1U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.25U over Ubiq. proms." = list(bam_MNase_ce11[['0.25U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 0.5U over Ubiq. proms." = list(bam_MNase_ce11[['0.5U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 1U over Ubiq. proms." = list(bam_MNase_ce11[['1U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 8U over Ubiq. proms." = list(bam_MNase_ce11[['8U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.']),
            "MNase-seq 16U over Ubiq. proms." = list(bam_MNase_ce11[['16U']], ce11_proms[ce11_proms$which.tissues == 'Ubiq.'])
        )
        plots <- plotVmat(
            list_params, 
            cores = 8,
            ylims = c(50, 200), 
            normFun = 'zscore', 
            nrow = 2, 
            ncol = 4, 
            verbose = 2
        )
        ggsave('ce11_proms_MNase_Serizay2020.pdf', height = 7.5, width = 14)
        #
        loc <- "chrIII:8,083,947-8,087,158"
        p <- plotProfile(bam_MNase_ce11_merged, loc, min = 20, max = 200, alpha = 0.1, size = 2)
        ggsave(paste0(loc, '.pdf'), width = 20, height = 8, limitsize = FALSE)
        #
        loc <- "chrIV:653,277-656,439"
        p <- plotProfile(bam_MNase_ce11_merged, loc, min = 20, max = 200, alpha = 0.1, size = 2)
        ggsave(paste0(loc, '.pdf'), width = 20, height = 8, limitsize = FALSE)
    })
})
