context("test-nucleosome-enrich")

test_that("nucleosome-enrich works", {
    expect_equal({
        bam_list <- readRDS(
            url('http://ahringerlab.com/VplotR/ATAC_PE_fragments.rds')
        )
        data(ce11_proms)
        # 
        frags <- bam_list[['mixed']]
        V <- plotVmat(
            frags,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
            normFun = '',
            return_Vmat = TRUE
        )
        set.seed(42)
        V_bg <-  plotVmat(
            frags,
            sampleGRanges(ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], n = 100000), 
            normFun = '',
            return_Vmat = TRUE
        )
        nuc2 <- nucleosomeEnrichment(V, V_bg)
        #
        set.seed(42)
        nuc3 <- nucleosomeEnrichment(
            frags,
            ce11_proms[ce11_proms$which.tissues == 'Ubiq.']
        )
        #
        unlink('tmp.pdf')
        TRUE
    }, TRUE)
})

test_that("nuc. enrich at ubiq. and tissue-spe. proms", {
    skip('skip')
    expect_equal({
        bam_list <- readRDS(
            url('http://ahringerlab.com/VplotR/ATAC_PE_fragments.rds')
        )
        data(ce11_proms)
        #
        list_scores <- parallel::mclapply(
            c("Germline", "Neurons", "Muscle", "Hypod.", "Intest."), 
            function(TISSUE) {
                message('>> ', TISSUE)
                set.seed(43)
                nucenrich_tissue_spe_proms <- nucleosomeEnrichment(
                    bam_list[[TISSUE]], 
                    ce11_proms[ce11_proms$which.tissues == TISSUE], 
                    verbose = TRUE, 
                    plus1_nuc_only = TRUE
                )
                set.seed(43)
                nucenrich_ubiq_proms <- nucleosomeEnrichment(
                    bam_list[[TISSUE]], 
                    ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
                    verbose = TRUE, 
                    plus1_nuc_only = TRUE
                )
                #
                res <- list(
                    'tissue-spe-proms' = nucenrich_tissue_spe_proms, 
                    'ubiq-proms' = nucenrich_ubiq_proms
                )
                return(res)
            }, 
            mc.cores = 5
        ) %>% setNames(c("Germline", "Neurons", "Muscle", "Hypod.", "Intest."))
        #
        nucenrich_scores <- data.frame(
            tissue = factor(
                rep(names(list_scores), each = 2), levels = names(list_scores)
            ), 
            promoters = c(
                rbind(
                    paste0(names(list_scores), ' promoters'), 
                    rep('Ubiq. promoters', 5)
                )
            ), 
            promoters2 = factor(rep(c(0.5, 0.8), 5)),
            score = unlist(lapply(
                list_scores, 
                function(Vmat) {
                    c(
                        Vmat[[1]]$fisher_test$estimate, 
                        Vmat[[2]]$fisher_test$estimate
                    )
                }
            )), 
            is_ubiq = factor(
                rep(c(
                    'Over tissue-specific promoters', 
                    'Over ubiquitous promoters'
                ), 5), 
                levels = c(
                    'Over ubiquitous promoters', 
                    'Over tissue-specific promoters'
                )
            )
        )
        p <- ggplot(nucenrich_scores, aes(
            x = tissue, 
            y = score, 
            fill = tissue, 
            group = is_ubiq
        )) + 
            geom_col() + 
            facet_wrap(~is_ubiq, nrow = 1) + 
            scale_fill_manual(values = c(
                '#1232D9', 
                '#3B9B46', 
                '#D99B12', 
                '#7e7e7e', 
                '#D912D4'
            )) +
            labs(
                title = 'Flanking nucleosome enrichment score', 
                y = 'Enrichment score', 
                x = 'Tissue-specific data'
            ) + 
            theme_bw() + 
            theme(legend.position = 'none')
        ggsave(
            'Comparison_tissue-specific-nucleosome-enrichment_plus1_seed43.pdf', 
            height = 5, width = 6
        )
        #
        pdf('tmp.pdf')
        lapply(list_scores, function(l) {
            l[['ubiq-proms']]$plot
        })
        lapply(list_scores, function(l) {
            l[['tissue-spe-proms']]$plot
        })
        dev.off()
        TRUE
    }, TRUE)
})

test_that("nuc. enrich at ubiq. and tissue-spe. proms -- with repetitions", {
    skip('skip')
    expect_equal({
        bam_list <- readRDS(
            url('http://ahringerlab.com/VplotR/ATAC_PE_fragments.rds')
        )
        data(ce11_proms)
        #
        l <- mclapply(mc.cores = 20, seq_len(20), function(SEED) {
            list_scores <- parallel::mclapply(
                c("Germline", "Neurons", "Muscle", "Hypod.", "Intest."), 
                function(TISSUE) {
                    message('>> ', TISSUE)
                    set.seed(SEED)
                    nucenrich_tissue_spe_proms <- nucleosomeEnrichment(
                        bam_list[[TISSUE]], 
                        ce11_proms[ce11_proms$which.tissues == TISSUE], 
                        verbose = TRUE, 
                        plus1_nuc_only = TRUE
                    )
                    set.seed(SEED)
                    nucenrich_ubiq_proms <- nucleosomeEnrichment(
                        bam_list[[TISSUE]], 
                        ce11_proms[ce11_proms$which.tissues == 'Ubiq.'], 
                        verbose = TRUE, 
                        plus1_nuc_only = TRUE
                    )
                    #
                    res <- list(
                        'tissue-spe-proms' = nucenrich_tissue_spe_proms, 
                        'ubiq-proms' = nucenrich_ubiq_proms
                    )
                    return(res)
                }, 
                mc.cores = 5
            ) %>% setNames(c("Germline", "Neurons", "Muscle", "Hypod.", "Intest."))
            nucenrich_scores <- data.frame(
                tissue = factor(
                    rep(names(list_scores), each = 2), levels = names(list_scores)
                ), 
                promoters = c(
                    rbind(
                        paste0(names(list_scores), ' promoters'), 
                        rep('Ubiq. promoters', 5)
                    )
                ), 
                promoters2 = factor(rep(c(0.5, 0.8), 5)),
                score = unlist(lapply(
                    list_scores, 
                    function(Vmat) {
                        c(
                            Vmat[[1]]$fisher_test$estimate, 
                            Vmat[[2]]$fisher_test$estimate
                        )
                    }
                )), 
                is_ubiq = factor(
                    rep(c(
                        'Over tissue-specific promoters', 
                        'Over ubiquitous promoters'
                    ), 5), 
                    levels = c(
                        'Over ubiquitous promoters', 
                        'Over tissue-specific promoters'
                    )
                ),
                seed = SEED
            )
        })
        l <- do.call(rbind, l)
        #
        p <- ggplot(l, aes(
            x = tissue, 
            y = score, 
            fill = tissue, 
            group = tissue
        )) + 
            geom_boxplot() + 
            facet_wrap(~is_ubiq, nrow = 1) + 
            scale_fill_manual(values = c(
                '#1232D9', 
                '#3B9B46', 
                '#D99B12', 
                '#7e7e7e', 
                '#D912D4'
            )) +
            labs(
                title = 'Flanking nucleosome enrichment score', 
                y = 'Enrichment score', 
                x = 'Tissue-specific data'
            ) + 
            theme_bw() + 
            theme(legend.position = 'none')
        ggsave(
            'Comparison_tissue-specific-nucleosome-enrichment_plus1_boxplot.pdf', 
            height = 5, width = 6
        )
        #
        pdf('tmp.pdf')
        lapply(list_scores, function(l) {
            l[['ubiq-proms']]$plot
        })
        lapply(list_scores, function(l) {
            l[['tissue-spe-proms']]$plot
        })
        dev.off()
        TRUE
    }, TRUE)
})
