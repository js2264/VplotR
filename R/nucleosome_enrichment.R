#' A function to compute nucleosome enrichment over a set of GRanges
#' 
#' @param x a GRanges or Vmat
#' @param ... additional parameters
#' @return list
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' nucleosomeEnrichment(bam_test, ce11_proms)

nucleosomeEnrichment <- function(x, ...) {
    UseMethod("nucleosomeEnrichment")
}

#' A function to compute nucleosome enrichment over a Vmat
#'
#' @param x a computed Vmat. Should be un-normalized.
#' @param background a background Vmat. Should be un-normalized.
#' @param plus1_nuc_only Boolean, should compute nucleosome 
#' enrichment only for +1 nucleosome?
#' @param ... additional parameters
#' @return list
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' V <- plotVmat(
#'     bam_test,
#'     ce11_proms,
#'     normFun = '',
#'     return_Vmat = TRUE
#' )
#' V_bg <- plotVmat(
#'     bam_test,
#'     sampleGRanges(ce11_proms),
#'     normFun = '',
#'     return_Vmat = TRUE
#' )
#' nucleosomeEnrichment(V, V_bg)

nucleosomeEnrichment.Vmat <- function(
    x, 
    background, 
    plus1_nuc_only = FALSE, 
    ...
) 
{
    computeNucleosomeEnrichmentOverBackground(
        Vmat = x, 
        background = background, 
        plus1_nuc_only = plus1_nuc_only, 
        ...
    )
}

#' A function to compute nucleosome enrichment over a set of GRanges
#'
#' @param x GRanges, paired-end fragments
#' @param granges GRanges, loci to map the fragments onto
#' @param plus1_nuc_only Boolean, should compute nucleosome 
#' enrichment only for +1 nucleosome?
#' @param verbose Boolean
#' @param ... additional parameters
#' @return list
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' nucleosomeEnrichment(bam_test, ce11_proms)

nucleosomeEnrichment.GRanges <- function(
    x, 
    granges, 
    plus1_nuc_only = FALSE, 
    verbose = TRUE, 
    ...
) 
{
    # Get Vmat
    if (verbose) message("Computing Vmat...")
    Vmat <- plotVmat(
        x, 
        granges, 
        normFun = '', 
        roll = 1, 
        return_Vmat = TRUE, 
        verbose = 0
    )
    # Get background 
    if (verbose) message("Computing background...")
    background <- plotVmat(
        x, 
        sampleGRanges(granges), 
        normFun = '', 
        roll = 1, 
        return_Vmat = TRUE, 
        verbose = 0
    )
    # Computing enrichment
    if (verbose) message("Computing enrichment...")
    res <- nucleosomeEnrichment(Vmat, background, plus1_nuc_only, ...)
    return(res)
}

#' Internal function
#' 
#' A function to compute nucleosome enrichment of a Vmat
#'
#' @param Vmat A Vmat computed by nucleosomeEnrichment function
#' @param background a background Vmat
#' @param plus1_nuc_only Boolean Should compute nucleosome enrichment
#' only for +1 nucleosome?
#' @param minus1_nuc list where the -1 nucleosome is located
#' @param minus1_nuc_neg where the background of the -1 nucleosome 
#' is located
#' @param plus1_nuc where the +1 nucleosome is located
#' @param plus1_nuc_neg where the background of the +1 nucleosome
#' is located
#' @param ... additional parameters
#' @return list
#' 
#' @import ggplot2
#' @importFrom cowplot plot_grid

computeNucleosomeEnrichmentOverBackground <- function(
    Vmat, 
    background = NULL, 
    plus1_nuc_only = FALSE, 
    minus1_nuc = list(c(xmin = -150, xmax = -70), c(ymin = 165, ymax = 260)), 
    minus1_nuc_neg = list(
        c(xmin = -150, xmax = -70), c(ymin = 60, ymax = 145)
    ), 
    plus1_nuc = list(c(xmin = 70, xmax = 150), c(ymin = 165, ymax = 260)), 
    plus1_nuc_neg = list(c(xmin = 70, xmax = 150), c(ymin = 50, ymax = 145)),
    ...
) 
{
    if (is.null(background)) {
        background <- shuffleVmat(Vmat)
    }
    # Generate plots
    minus1_nuc_df <- data.frame(
        xmin = minus1_nuc[[1]][[1]], 
        xmax = minus1_nuc[[1]][[2]], 
        ymin = minus1_nuc[[2]][[1]], 
        ymax = minus1_nuc[[2]][[2]]
    ) 
    plus1_nuc_df <- data.frame(
        xmin = plus1_nuc[[1]][[1]], 
        xmax = plus1_nuc[[1]][[2]], 
        ymin = plus1_nuc[[2]][[1]], 
        ymax = plus1_nuc[[2]][[2]]
    )
    minus1_nuc_neg_df <- data.frame(
        xmin = minus1_nuc_neg[[1]][[1]], 
        xmax = minus1_nuc_neg[[1]][[2]], 
        ymin = minus1_nuc_neg[[2]][[1]], 
        ymax = minus1_nuc_neg[[2]][[2]]
    )
    plus1_nuc_neg_df <- data.frame(
        xmin = plus1_nuc_neg[[1]][[1]], 
        xmax = plus1_nuc_neg[[1]][[2]], 
        ymin = plus1_nuc_neg[[2]][[1]], 
        ymax = plus1_nuc_neg[[2]][[2]]
    )
    #
    if (plus1_nuc_only == TRUE) {
        df_controls <- data.frame(
            rbind(plus1_nuc_df, plus1_nuc_neg_df), 
            locus = c('plus1_nuc', 'plus1_nuc_neg'), 
            type = factor(c('nuc', 'neg'), levels = c('nuc', 'neg'))
        )
    }
    else {
        df_controls <- data.frame(
            rbind(
                minus1_nuc_df, 
                plus1_nuc_df,
                minus1_nuc_neg_df, 
                plus1_nuc_neg_df
            ), 
            locus = c(
                'minus1_nuc', 
                'plus1_nuc', 
                'minus1_nuc_neg',
                'plus1_nuc_neg'
             ), 
            type = factor(
                c('nuc', 'nuc', 'neg', 'neg'), levels = c('nuc', 'neg')
            )
        )
    }
    #
    p <- plotVmat(Vmat) + 
        ggplot2::geom_rect(
            data = df_controls, 
            ggplot2::aes(
                xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, col = type
            ), 
            inherit.aes=FALSE, 
            fill = NA, 
            size = 2
        )
    q <- plotVmat(background) + 
        ggplot2::geom_rect(
            data = df_controls, 
            ggplot2::aes(
                xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, col = type
            ), 
            inherit.aes=FALSE, 
            fill = NA, 
            size = 2
        )
    pt <- cowplot::plot_grid(p, q, nrow = 1)
    # Shift ranges to center around the "0"
    minus1_nuc <- list(
        minus1_nuc[[1]] + round(nrow(Vmat)/2), 
        minus1_nuc[[2]] - as.numeric(colnames(Vmat)[1])
    )
    plus1_nuc <- list(
        plus1_nuc[[1]] + round(nrow(Vmat)/2), 
        plus1_nuc[[2]] - as.numeric(colnames(Vmat)[1])
    )
    minus1_nuc_neg <- list(
        minus1_nuc_neg[[1]] + round(nrow(Vmat)/2), 
        minus1_nuc_neg[[2]] - as.numeric(colnames(Vmat)[1])
    )
    plus1_nuc_neg <- list(
        plus1_nuc_neg[[1]] + round(nrow(Vmat)/2), 
        plus1_nuc_neg[[2]] - as.numeric(colnames(Vmat)[1])
    )
    # Count reads in each window
    if (!plus1_nuc_only) {
        vec <- c(
            sum(
                Vmat[minus1_nuc[[1]][1] : minus1_nuc[[1]][2], 
                minus1_nuc[[2]][1] : minus1_nuc[[2]][2]]
            ) + sum(
                Vmat[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], 
                plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]
            ), 
            sum(
                Vmat[minus1_nuc_neg[[1]][1] : minus1_nuc_neg[[1]][2], 
                minus1_nuc_neg[[2]][1] : minus1_nuc_neg[[2]][2]]
            ) + sum(
                Vmat[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], 
                plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]]
            ),
            sum(
                background[minus1_nuc[[1]][1] : minus1_nuc[[1]][2], 
                minus1_nuc[[2]][1] : minus1_nuc[[2]][2]]
            ) + sum(
                background[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], 
                plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]
            ), 
            sum(
                background[minus1_nuc_neg[[1]][1] : minus1_nuc_neg[[1]][2], 
                minus1_nuc_neg[[2]][1] : minus1_nuc_neg[[2]][2]]
            ) + sum(
                background[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], 
                plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]]
            )
        )
    } 
    else {
        vec <- c(
            sum(
                Vmat[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], 
                plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]
            ), 
            sum(
                Vmat[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], 
                plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]]
            ),
            sum(
                background[plus1_nuc[[1]][1] : plus1_nuc[[1]][2], 
                plus1_nuc[[2]][1] : plus1_nuc[[2]][2]]
            ), 
            sum(
                background[plus1_nuc_neg[[1]][1] : plus1_nuc_neg[[1]][2], 
                plus1_nuc_neg[[2]][1] : plus1_nuc_neg[[2]][2]]
            )
        )
    }
    # Test 
    test <- fisher.test(matrix(vec, ncol = 2))
    # Return result
    res <- list(
        Vmat = Vmat, 
        background = background,
        scores = matrix(vec, ncol = 2),
        fisher_test = test,
        plot = pt
    )
    return(res)
}
