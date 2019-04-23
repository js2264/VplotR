# ------------ getPlottingMetrics function ------------

getPlottingMetrics <- function(x, ...) {
    UseMethod("getPlottingMetrics")
}

getPlottingMetrics.default <- function(x, ...) {
    return()
}

getPlottingMetrics.CompressedGRangesList <- function(granges, list.bw.files, k, norm = 'none', flank = c(500, 500), stranded = T, BIN = 10) {
    names_granges <- names(granges)
    granges <- granges[[k]]
    scores <- list.bw.files
    # Center and resize GRanges
    granges <- GenomicRanges::resize(granges, width = 1, fix = 'center')
    granges <- GenomicRanges::resize(granges, fix = 'end', width = flank[1])
    granges <- GenomicRanges::resize(granges, fix = 'start', width = flank[1] + flank[2])
    # Subset scores over GRanges
    scores.subset <- scores[granges]
    # Turn it into a rectangular matrix and correct for strandness
    scores.subset <- suppressWarnings(matrix(as.vector(unlist(scores.subset)), nrow = length(granges), byrow = T))
    scores.subset.flipped <- t(apply(scores.subset, 1, rev))
    if (stranded) {
        scores.subset <- matrix(sapply(1:nrow(scores.subset), function(K) {
            if((as.vector(strand(granges)) == '-')[K]) {
                scores.subset.flipped[K,]
            } else {
                scores.subset[K,]
            }
        }), nrow = length(granges), byrow = T)
    }
    # Normalize matrix
    if (norm == 'zscore') {
        scores.subset <- apply(scores.subset, 1, scale)
    } else if (norm == 'log2') {
        scores.subset <- log2(scores.subset)
    }
    # Compute specific metrics to plot
    medians <- apply(scores.subset, 2, median)
    means <- apply(scores.subset, 2, mean)
    smoothed_line <- if (use.mean) {zoo::rollmean(means, BIN)} else {zoo::rollmean(medians, BIN)}
    stderror <- zoo::rollmean(apply(scores.subset, 2, function(n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n)) ) }), BIN)
    conint  <- zoo::rollmean(apply(scores.subset, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), BIN)
    topEE <- smoothed_line + conint
    bottomEE <- smoothed_line - conint
    # Return matrix
    df <- data.frame(
        "sample" = ifelse(!is.null(names_granges), names_granges[k], k),
        "x" = c(-flank[1]:flank[2])[1:length(smoothed_line)],
        "smoothed_line" = smoothed_line, 
        "topEE" = topEE, 
        "bottomEE" = bottomEE
    )
    return(df)
}

getPlottingMetrics.GRanges <- function(granges, list.bw.files, k, norm = 'none', flank = c(500, 500), stranded = T, BIN = 10) {
    scores <- list.bw.files[[k]]
    # Center and resize GRanges
    granges <- GenomicRanges::resize(granges, width = 1, fix = 'center')
    granges <- GenomicRanges::resize(granges, fix = 'end', width = flank[1])
    granges <- GenomicRanges::resize(granges, fix = 'start', width = flank[1] + flank[2])
    # Subset scores over GRanges
    scores.subset <- scores[granges]
    # Turn it into a rectangular matrix and correct for strandness
    scores.subset <- suppressWarnings(matrix(as.vector(unlist(scores.subset)), nrow = length(granges), byrow = T))
    scores.subset.flipped <- t(apply(scores.subset, 1, rev))
    if (stranded) {
        scores.subset <- matrix(sapply(1:nrow(scores.subset), function(K) {
            if((as.vector(strand(granges)) == '-')[K]) {
                scores.subset.flipped[K,]
            } else {
                scores.subset[K,]
            }
        }), nrow = length(granges), byrow = T)
    }
    # Normalize matrix
    if (norm == 'zscore') {
        scores.subset <- apply(scores.subset, 1, scale)
    } else if (norm == 'log2') {
        scores.subset <- log2(scores.subset)
    }
    # Compute specific metrics to plot
    medians <- apply(scores.subset, 2, median)
    means <- apply(scores.subset, 2, mean)
    smoothed_line <- if (use.mean) {zoo::rollmean(means, BIN)} else {zoo::rollmean(medians, BIN)}
    stderror <- zoo::rollmean(apply(scores.subset, 2, function(n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n)) ) }), BIN)
    conint  <- zoo::rollmean(apply(scores.subset, 2, function (n) { qt(0.975,sum(!is.na(n)))*sd(n,na.rm=TRUE)/sqrt(sum(!is.na(n))) }), BIN)
    topEE <- smoothed_line + conint
    bottomEE <- smoothed_line - conint
    # Return matrix
    df <- data.frame(
        "sample" = ifelse(!is.null(names(list.bw.files)), names(list.bw.files)[k], k),
        "x" = c(-flank[1]:flank[2])[1:length(smoothed_line)],
        "smoothed_line" = smoothed_line, 
        "topEE" = topEE, 
        "bottomEE" = bottomEE
    )
    return(df)
}

# ------------ plotMetaCoverage function ------------

coverageMetaplot <- function(x, ...) {
    UseMethod("coverageMetaplot")
}

coverageMetaplot.default <- function(x, ...) {
    stop("Provide a bw.list in first argument")
}

coverageMetaplot.character <- function(file, name = "Sample", ...) {
    list.bw <- importBigWig(file, verbose = F) %>% stats::setNames(name)
    coverageMetaplot(list.bw, ...)
}

coverageMetaplot.SimpleRleList <- function(x, name = "Sample", ...) {
    x <- list(x) %>% stats::setNames(name)
    class(x) <- c("listBigWig", class(x))
    coverageMetaplot(x, ...)
}

coverageMetaplot.listBigWig <- function(
    list.bw.files, # Can be either a bw.file or an alreday imported Rle
    list.granges, # The regions of interest
    colors = NULL,
    stranded = T,
    use.mean = T,
    BIN = 10,
    ylim = NULL, 
    xlim = NULL,
    xlab = NULL,
    ylab = NULL,
    by.granges = T,
    ...
    ) {
    # Check that there is 1 unique set of GRanges in list.granges
    if (!any(grepl("GRangesList", class(list.granges)[[1]]))) {
        if (any(class(list.granges) == 'GRanges')) {
            list.granges <- GenomicRanges::GRangesList(list(list.granges)) %>% setNames("Segments")
        } else {
            stop('Please provide a GRangesList object in second argument')
        }
    }
    # Define plotting variables
    xlab <- if (is.null(xlab)) {'Segments coordinates'} else {xlab}
    ylab <- ifelse(is.null(ylab), 'Score', ylab)
    if (!is.null(colors)) {
        colors <- rep(colors, 5)
    } else {
        colors <- c( '#1232D9', '#3B9B46', '#D99B12', '#9e9e9e', '#D912D4', '#9E7FCC', '#B0E797', '#D1B3B3', '#991919', '#23A4A3', '#000000', '#dbdbdb')
    }
    # Loop through Granges
    if (by.granges) {
        list.plots <- parallel::mclapply(1:length(list.granges), function(k.g) {
            # Get the metrics
            df <- lapply(seq_along(list.bw.files), function(k.b) {
                getPlottingMetrics(list.granges[[k.g]], list.bw.files, k.b, flank = flank, stranded = stranded, BIN = BIN)
            }) %>% do.call(rbind, .)
            if (is.null(xlim)) {
                xlim2 <- c(-flank[1], flank[2])
            } else {
                xlim2 <- xlim
            }
            if (is.null(ylim)) {
                ylim2 <- c(min(df$smoothed_line) * 0.9, max(df$smoothed_line) * 1.1)
            } else {
                ylim2 <- ylim
            }
            # Plot
            ggplot2::ggplot(df, ggplot2::aes(x = x, y = smoothed_line, ymin = bottomEE, ymax = topEE, group = sample, col = sample, fill = sample)) +
                ggplot2::geom_line() + 
                ggplot2::geom_ribbon(alpha = 0.2, col = NA) + 
                ggplot2::geom_vline(xintercept = 0, lty = 3, alpha = 0.4) +
                ggplot2::theme_classic() + 
                ggplot2::labs(
                    title = names(list.granges)[k.g],
                    x = xlab, 
                    y = ylab
                ) + 
                ggplot2::xlim(xlim2) + 
                ggplot2::ylim(ylim2) + 
                ggplot2::scale_color_manual(values = colors) +
                ggplot2::scale_fill_manual(values = colors) + 
                ggplot2::theme(legend.position = "none")
        }, mc.cores = min(4, length(list.granges)))
        df <- data.frame("Group" = {if (!is.null(names(list.bw.files))) {factor(names(list.bw.files), levels = names(list.bw.files))} else {factor(1:length(list.bw.files))}})
        for_legend <- ggplot2::ggplot(df, ggplot2::aes(x = NA, y = NA, ymin = NA, ymax = NA, group = Group, col = Group, fill = Group)) +
            ggplot2::geom_line() + 
            ggplot2::geom_ribbon(alpha = 0.2, col = NA) + 
            ggplot2::scale_color_manual(values = colors) +
            ggplot2::scale_fill_manual(values = colors)
        list.plots[[length(list.plots) + 1]] <- suppressMessages(cowplot::get_legend(for_legend))
    } else {
        list.plots <- parallel::mclapply(1:length(list.bw.files), function(k.b) {
            # Get the metrics
            df <- lapply(seq_along(list.granges), function(k.g) {
                getPlottingMetrics(list.granges, list.bw.files[[k.b]], k.g, flank = flank, stranded = stranded, BIN = BIN)
            }) %>% do.call(rbind, .)
            if (is.null(xlim)) {
                xlim2 <- c(-flank[1], flank[2])
            } else {
                xlim2 <- xlim
            }
            if (is.null(ylim)) {
                ylim2 <- c(min(df$smoothed_line) * 0.9, max(df$smoothed_line) * 1.1)
            } else {
                ylim2 <- ylim
            }
            # Plot
            ggplot2::ggplot(df, ggplot2::aes(x = x, y = smoothed_line, ymin = bottomEE, ymax = topEE, group = sample, col = sample, fill = sample)) +
                ggplot2::geom_line() + 
                ggplot2::geom_ribbon(alpha = 0.2, col = NA) + 
                ggplot2::geom_vline(xintercept = 0, lty = 3, alpha = 0.4) +
                ggplot2::theme_classic() + 
                ggplot2::labs(
                    title = names(list.bw.files)[k.b],
                    x = xlab, 
                    y = ylab
                ) + 
                ggplot2::xlim(xlim2) + 
                ggplot2::ylim(ylim2) + 
                ggplot2::scale_color_manual(values = colors) +
                ggplot2::scale_fill_manual(values = colors) + 
                ggplot2::theme(legend.position = "none")
        }, mc.cores = min(4, length(list.bw.files)))
        df <- data.frame("Group" = {if (!is.null(names(list.granges))) {factor(names(list.granges), levels = names(list.granges))} else {factor(1:length(list.granges))}})
        for_legend <- ggplot2::ggplot(df, ggplot2::aes(x = NA, y = NA, ymin = NA, ymax = NA, group = Group, col = Group, fill = Group)) +
            ggplot2::geom_line() + 
            ggplot2::geom_ribbon(alpha = 0.2, col = NA) + 
            ggplot2::scale_color_manual(values = colors) +
            ggplot2::scale_fill_manual(values = colors)
        list.plots[[length(list.plots) + 1]] <- suppressMessages(cowplot::get_legend(for_legend))
    }
    if (length(list.plots) == 2)
        list.plots <- list.plots[[1]] + ggplot2::theme(legend.position = "bottom")
    return(list.plots)
}