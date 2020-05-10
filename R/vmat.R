#' A function to generate a Vplot
#'
#' See individual methods for further detail
#'
#' @param x GRanges or list or Vmat
#' @param ... additional parameters
#' @return A Vmat ggplot
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' V <- plotVmat(
#'     bam_test,
#'     ce11_proms,
#'     normFun = 'libdepth+nloci'
#' )

plotVmat <- function(x, ...) {
    UseMethod("plotVmat")
}

#' A function to plot a computed Vmat
#'
#' The default plotVmat method generates a ggplot representing a 
#' heatmap of fragment density. 
#'
#' @param x A computed Vmat (ideally, should be normalized)
#' @param hm Integer, should be between 0 and 100. 
#' Used to automatically
#' scale the range of colors (best to 
#' keep between 90 and 100)
#' @param colors a vector of colors
#' @param breaks a vector of breaks. 
#' length(breaks) ==  length(colors) + 1
#' @param xlim vector of two integers, x limits
#' @param ylim vector of two integers, y limits
#' @param main character, title of the plot
#' @param xlab character, x-axis label
#' @param ylab character, y-axis label
#' @param ... additional parameters
#' @return A Vmat ggplot
#' 
#' @import ggplot2
#' @import RColorBrewer
#' @import reshape2
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' V <- plotVmat(
#'     bam_test,
#'     ce11_proms,
#'     normFun = 'libdepth+nloci', 
#'     return_Vmat = TRUE
#' )
#' plotVmat(V)

plotVmat.default <- function(
    x, 
    hm = 90, 
    colors = COLORSCALE_VMAT, 
    breaks = NULL, 
    xlim = c(-250, 250),
    ylim = c(50, 300),
    main = '', 
    xlab = 'Distance from center of elements',
    ylab = 'Fragment length',
    ...
) 
{
    # Define breaks and clamp matrix within breaks
    if (is(x, 'Vmat')) {
        if (is.null(breaks)) {
            probs <- quantile(c(x), probs = seq(0, 1, length.out = 101))
            breaks <- seq(
                probs[(100-hm)/2], 
                probs[hm+(100-hm)/2], 
                length.out = length(colors)+1
            )
            if (length(unique(breaks)) == 1) {
                breaks <- seq(
                    min(x), max(x), length.out = length(colors)+1
                )
            }
        }
        x[x < min(breaks)] <- min(breaks)
        x[x > max(breaks)] <- max(breaks)
    } 
    else if (is(x, 'VmatList')) {
        if (is.null(breaks)) {
            probs <- quantile(
                c(do.call(rbind, x)), probs = seq(0, 1, length.out = 101)
            )
            breaks <- seq(
                probs[(100-hm)/2], probs[hm+(100-hm)/2], 
                length.out = length(colors)+1
            )
            if (length(unique(breaks)) == 1) {
                breaks <- seq(
                    do.call(min, x), 
                    do.call(max, x), 
                    length.out = length(colors)+1
                )
            }
        }
        x <- lapply(x, function(V) {
            V[V < min(breaks)] <- min(breaks)
            V[V > max(breaks)] <- max(breaks)
            return(V)
        })
        class(x) <- c('VmatList', class(x))
    }
    
    # Plot
    df <- reshape2::melt(x)
    if (is(x, 'VmatList')) {
        colnames(df) <- c('Var1', 'Var2', 'value', 'Cond.')
        df$Cond. <- factor(df$Cond., levels = names(x))
    } 
    else {
        colnames(df) <- c('Var1', 'Var2', 'value')
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(Var1, Var2, fill = value))
    p <- p + ggplot2::geom_raster()
    p <- p + ggplot2::scale_fill_gradientn(
        breaks = c(min(breaks), max(breaks)), 
        colors = colors
    )
    p <- p + ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        legend.background = ggplot2::element_blank(), 
        panel.border = ggplot2::element_rect(
            colour = "black", fill = NA, size = 1
        )
    )
    p <- p + ggplot2::scale_x_continuous(expand = c(0, 0))
    p <- p + ggplot2::scale_y_continuous(expand = c(0, 0))
    p <- p + ggplot2::labs(
        title = main,
        x = xlab,
        y = ylab, 
        fill = 'Score'
    )
    p <- p + ggplot2::coord_fixed()
    return(p)
}

#' A function to plot a computed Vmat
#'
#' The plotVmat.Vmat() method forwards the Vmat to plotVmat.default(). 
#'
#' @param x A computed Vmat (ideally, should be normalized)
#' @param ... additional parameters
#' @return A Vmat ggplot
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' V <- plotVmat(
#'     bam_test,
#'     ce11_proms,
#'     normFun = 'libdepth+nloci', 
#'     return_Vmat = TRUE
#' )
#' plotVmat(V)

plotVmat.Vmat <- function(x, ...) {
    plotVmat.default(x, ...)
}

#' A function to plot a computed VmatList
#'
#' The plotVmat.VmatList() method forwards the Vmat to plotVmat.default(). 
#'
#' @param x A VmatList (output of plotVmat.list())
#' @param nrow Integer, how many rows in facet?
#' @param ncol Integer, how many cols in facet?
#' @param dir str, direction of facets?
#' @param ... additional parameters
#' @return A Vmat ggplot
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' list_params <- list(
#'     'germline' = list(
#'         bam_test,
#'         ce11_proms[ce11_proms$which.tissues == 'Germline']
#'     ),
#'     'muscle' = list(
#'         bam_test,
#'         ce11_proms[ce11_proms$which.tissues == 'Muscle']
#'     )
#' )
#' V <- plotVmat(
#'     list_params,
#'     normFun = 'libdepth+nloci', 
#'     roll = 5
#' )

plotVmat.VmatList <- function(x, nrow = NULL, ncol = NULL, dir = 'v', ...) {
    p <- plotVmat.default(x, ...)
    if (is.null(nrow) | is.null(ncol)) {
        nrow <- 1
        ncol <- length(x)
    }
    p <- p + 
        facet_wrap(~Cond., nrow = nrow, ncol = ncol, dir = dir) + 
        theme(legend.position = 'bottom') +
        theme(panel.spacing = unit(1, "lines"))
}

#' A function to compute (and plot) a Vmat
#'
#' The plotVmat.GRanges() method computes and normalizes a Vmat 
#' before passing it to plotVmat.Vmat() method.
#'
#' @param x GRanges, paired-end fragments
#' @param granges GRanges, loci to map the fragments onto
#' @param xlims x limits of the computed Vmat
#' @param ylims y limits of the computed Vmat
#' @param normFun character. A Vmat should be scaled either by:
#' \itemize{
#'     \item 'libdepth+nloci', e.g. the library depth and the number of 
#'     loci used to compute the Vmat;
#'     \item zscore, if relative patterns of fragment density 
#'     are more important than density per se;
#'     \item Alternatively, the Vmat can be scaled to % ('pct'), to 
#'     a chosen quantile ('quantile') or to the max Vmat value ('max').
#' }
#' @param s A float indicating which quantile to use if 'quantile'
#' normalization is chosen
#' @param roll integer, to use as the window to smooth the Vmat rows 
#' by rolling mean.
#' @param return_Vmat Boolean, should the function return the computed 
#' Vmat rather than the plot?
#' @param verbose Boolean 
#' @param ... additional parameters
#' @return A Vmat ggplot
#' 
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' V <- plotVmat(
#'     bam_test,
#'     ce11_proms,
#'     normFun = 'libdepth+nloci', 
#'     roll = 5
#' )

plotVmat.GRanges <- function(
    x, 
    granges,
    xlims = c(-250, 250), 
    ylims = c(50, 300), 
    normFun = 'libdepth+nloci',
    s = 0.95, 
    roll = 3,
    return_Vmat = FALSE,
    verbose = 1,
    ...
)
{
    # Calculate Vmat
    if (verbose) message('Computing V-mat') 
    Vmat <- computeVmat(x, granges, xlims, ylims)
    # Normalize Vmat 
    if (verbose) message('Normalizing the matrix') 
    Vmat <- normalizeVmat(
        Vmat, 
        x, 
        granges, 
        normFun, 
        s, 
        roll
    )
    # Replace NA / inf values by 0
    Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
    class(Vmat) <- c("Vmat", class(Vmat))
    #
    if (return_Vmat == TRUE) {
        return(Vmat)
    } else {
        p <- plotVmat(Vmat, ...)
        return(p)
    }
}

#' A function to compute (and plot) several Vmats.
#'
#' The plotVmat.GRanges() method computes and normalizes multiple Vmats
#' before passing them to plotVmat.VmatList() method.
#'
#' @param x list Each element of the list should be a list containing
#' paired-end fragments and GRanges of interest.
#' @param cores Integer, number of cores to parallelize the plots
#' @param nrow Integer, how many rows in facet?
#' @param ncol Integer, how many cols in facet?
#' @param xlims x limits of the computed Vmat
#' @param ylims y limits of the computed Vmat
#' @param normFun character. A Vmat should be scaled either by:
#' \itemize{
#'     \item 'libdepth+nloci', e.g. the library depth and the number of 
#'     loci used to compute the Vmat;
#'     \item zscore, if relative patterns of fragment density 
#'     are more important than density per se;
#'     \item Alternatively, the Vmat can be scaled to % ('pct'), to 
#'     a chosen quantile ('quantile') or to the max Vmat value ('max').
#' }
#' @param s A float indicating which quantile to use if 'quantile'
#' normalization is chosen
#' @param roll integer, to use as the window to smooth the Vmat rows 
#' by rolling mean.
#' @param return_Vmat Boolean, should the function return the computed 
#' Vmat rather than the plot?
#' @param verbose Boolean 
#' @param ... additional parameters
#' @return A list of Vmat ggplots
#' 
#' @import parallel
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' list_params <- list(
#'     'germline' = list(
#'         bam_test,
#'         ce11_proms[ce11_proms$which.tissues == 'Germline']
#'     ),
#'     'muscle' = list(
#'         bam_test,
#'         ce11_proms[ce11_proms$which.tissues == 'Muscle']
#'     )
#' )
#' V <- plotVmat(
#'     list_params,
#'     normFun = 'libdepth+nloci', 
#'     roll = 5
#' )

plotVmat.list <- function(
    x, 
    cores = 1,
    nrow = NULL, 
    ncol = NULL, 
    xlims = c(-250, 250), 
    ylims = c(50, 300), 
    normFun = 'libdepth+nloci',
    s = 0.95, 
    roll = 3,
    return_Vmat = FALSE,
    verbose = 1,
    ...
) 
{
    # Calculate all the Vmats
    if (is.null(names(x))) 
        stop('Please provide a *named* list of parameters. Aborting.')
    Vmats_list <- parallel::mclapply(seq_along(x), function(K) {
        if (verbose) message('- Processing sample ', K, '/', length(x)) 
        bam_granges <- x[[K]][[1]]
        granges <- x[[K]][[2]]
        Vmat <- plotVmat(
            bam_granges, 
            granges, 
            xlims, 
            ylims, 
            normFun, 
            s,
            roll, 
            return_Vmat = TRUE
        )
        # Replace NA / inf values by 0
        Vmat[is.infinite(Vmat) | is.na(Vmat)] <- 0
        return(Vmat)
    }, mc.cores = cores)
    names(Vmats_list) <- names(x)
    class(Vmats_list) <- "VmatList"
    # Plot
    if (return_Vmat == TRUE) {
        return(Vmats_list)
    } 
    else {
        plotVmat(Vmats_list, nrow, ncol, ...)
    }
}
