#' A function to generate a Vplot along chromosome coordinates
#'
#' @param fragments GRanges
#' @param window character, chromosome location
#' @param loci GRanges, optional genomic locus. Fragments overlapping 
#' this locus will be in red.
#' @param annots GRanges, optional gene annotations
#' @param min integer, minimum fragment size
#' @param max integer, maximum fragment size
#' @param alpha float, transparency value
#' @param size float, dot size
#' @param with_densities Boolean, should the densities be plotted?
#' @param verbose Boolean
#' @return A ggplot
#' 
#' @import GenomicRanges
#' @import IRanges
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @importFrom zoo rollmean
#' @export
#' 
#' @examples
#' data(bam_test)
#' data(ce11_proms)
#' V <- plotProfile(
#'     bam_test,
#'     'chrI:10000-12000',
#'     loci = ce11_proms,
#'     min = 80, 
#'     max = 200
#' )

plotProfile <- function(
    fragments, 
    window = loc, 
    loci = NULL,
    annots = NULL,
    min = 50,
    max = 200, 
    alpha = 0.5,
    size = 1,
    with_densities = TRUE,
    verbose = TRUE
) 
{
    if (verbose) message("Filtering and resizing fragments")
    loc <- GenomicRanges::GRanges(gsub(',', '', window))
    fragments <- IRanges::subsetByOverlaps(fragments, loc)
    if (verbose) message(
        length(fragments), 
        " fragments mapped over ", 
        GenomicRanges::width(loc), 
        " bases"
    )
    
    # Map fragments
    if (verbose) message("Filtering and resizing fragments")
    fragments_1bp <- GenomicRanges::resize(
        fragments, fix = 'center', width = 1
    )
    fragments_1bp$original_width <- GenomicRanges::width(fragments)
    df <- data.frame(
        pos = GenomicRanges::start(fragments_1bp), 
        width = GenomicRanges::width(fragments)
    )
    if (!is.null(loci)) {
        sub <- subsetByOverlaps(fragments, loci, ignore.strand = TRUE)
        df$col <- factor(
            ifelse(
                fragments %in% sub, 
                'red', 'black'
            ), 
            levels = c('black', 'red')
        )
    } 
    else {
        df$col <- factor('black', levels = c('black', 'red'))
    }
    
    # Dot plot
    if (verbose) message("Generating plot")
    p <- ggplot2::ggplot(
        df, ggplot2::aes(pos, width, col = col)
    )
    p <- p + ggplot2::geom_point(shape = 16, size = size, alpha = alpha)
    p <- p + ggplot2::scale_color_manual(values = c('#000000', '#ff0000'))
    p <- p + ggplot2::scale_x_continuous(
        limits = c(GenomicRanges::start(loc), GenomicRanges::end(loc)),
        expand = c(0, 0)
    )
    p <- p + ggplot2::scale_y_continuous(
        limits = c(min, max),
        expand = c(0, 0)
    )
    p <- p + theme_ggplot2()
    p <- p + ggplot2::theme(
        plot.background = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
            colour = "black", fill = NA, size = 1
        ), 
        legend.position = 'none'
    )
    p <- p + ggplot2::labs(
        x = paste0(
            "Coordinates (", as.character(GenomeInfoDb::seqnames(loc)), ")"
        ), 
        y = "Fragment sizes", 
        title = loc
    )
    # p <- p + ggplot2::coord_fixed()
    if (!is.null(loci)) {
        g <- loci[loci %over% fragments]
        p <- p + ggplot2::geom_rect(
            data = data.frame(
                x1 = GenomicRanges::start(g),
                x2 = GenomicRanges::end(g), 
                y1 = min, 
                y2 = max
            ),
            mapping = aes(
                xmin = x1, 
                xmax = x2, 
                ymin = y1, 
                ymax = y2
            ),
            inherit.aes = FALSE,
            alpha = 0.25, 
            col = NA, 
            fill = '#ff0000'
        )
    }
    
    # Gene tracks
    if (!is.null(annots)) {
        g <- annots[annots %over% fragments]
        d <- data.frame(
            x1 = unlist(lapply(
                GenomicRanges::start(g), 
                function(s) max(s, GenomicRanges::start(loc))
            )),
            x2 = unlist(lapply(
                GenomicRanges::end(g), 
                function(s) min(s, GenomicRanges::end(loc))
            )),
            y1 = 0, 
            y2 = 1
        )
        gene_track <- ggplot2::ggplot(d) + 
            ggplot2::geom_rect(aes(
                xmin = x1, 
                xmax = x2, 
                ymin = y1, 
                ymax = y2
            ),
            alpha = 0.25, 
            col = NA, 
            fill = '#000000'
        )
        gene_track <- gene_track + ggplot2::scale_x_continuous(
            limits = c(GenomicRanges::start(loc), GenomicRanges::end(loc)),
            expand = c(0, 0)
        )
        gene_track <- gene_track + theme_ggplot2()
        gene_track <- gene_track + ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(), 
            panel.grid = ggplot2::element_blank(), 
            panel.grid.major.y = ggplot2::element_blank(), 
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(), 
            panel.grid.minor.x = ggplot2::element_blank(), 
            plot.margin = ggplot2::margin(t=0, b=0, l=0, r=0, unit='pt')
        )
        gene_track <- gene_track + ggplot2::labs(
            x = "", 
            y = "", 
            title = ''
        )
    } 
    else {
        gene_track <- ggplot2::ggplot() + ggplot2::theme_void()
    }
    
    # Densities
    marg_y <- data.frame(
        pos = IRanges::start(loc):IRanges::end(loc), 
        dens = as.vector(table(factor(
            df$pos, 
            levels = IRanges::start(loc):IRanges::end(loc)
        )))
    )
    marg_y$dens <- zoo::rollmean(marg_y$dens, 50, fill = NA, fix = 'center')
    marg_y$dens <- marg_y$dens/sum(marg_y$dens, na.rm = TRUE)
    q1 <- ggplot2::ggplot(marg_y, ggplot2::aes(x = pos, y = dens))
    q1 <- q1 + ggplot2::geom_line()
    q1 <- q1 + ggplot2::scale_x_continuous(
        limits = c(GenomicRanges::start(loc), GenomicRanges::end(loc)),
        expand = c(0, 0)
    )
    q1 <- q1 + theme_ggplot2()
    q1 <- q1 + ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(), 
        panel.grid.major.y = ggplot2::element_blank(), 
        panel.grid.minor.y = ggplot2::element_blank(), 
        plot.margin = ggplot2::margin(t=0, b=0, l=0, r=0, unit='pt')
    )
    q1 <- q1 + ggplot2::labs(
        x = "", 
        y = "", 
        title = ''
    )
    
    marg_x <- data.frame(
        widths = min:max, 
        dens = as.vector(table(factor(
            df$width, 
            levels = min:max
        )))
    )
    marg_x$dens <- zoo::rollmean(marg_x$dens, 20, fill = NA, fix = 'center')
    marg_x$dens <- marg_x$dens/sum(marg_x$dens, na.rm = TRUE)
    q2 <- ggplot2::ggplot(marg_x, ggplot2::aes(x = widths, y = dens))
    q2 <- q2 + ggplot2::geom_line()
    q2 <- q2 + ggplot2::scale_x_continuous(
        limits = c(min, max),
        expand = c(0, 0)
    )
    q2 <- q2 + ggplot2::coord_flip()
    q2 <- q2 + theme_ggplot2()
    q2 <- q2 + ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(), 
        panel.grid = ggplot2::element_blank(), 
        panel.grid.major.x = ggplot2::element_blank(), 
        panel.grid.minor.x = ggplot2::element_blank(), 
        plot.margin = margin(t=12, b=20, r=30, l=3, unit='pt')
    )
    q2 <- q2 + ggplot2::labs(
        x = "", 
        y = "", 
        title = ''
    )
    
    if (with_densities) {
        # Assemble everything
        if (!is.null(annots)) {
            p_left <- cowplot::plot_grid(
                gene_track, q1, p,
                nrow = 3, 
                align = "v", 
                rel_heights = c(1, 3, 6)
            )
            p_right <- cowplot::plot_grid(
                ggplot2::ggplot() + ggplot2::theme_void(),
                ggplot2::ggplot() + ggplot2::theme_void(),
                q2, 
                nrow = 3, 
                align = "v", 
                rel_heights = c(1, 3, 6)
            )
            p <- cowplot::plot_grid(
                p_left, p_right,
                align = "h", 
                nrow = 1, 
                rel_widths = c(8, 2)
            )
        }
        else {
            p_left <- cowplot::plot_grid(
                q1, p,
                nrow = 2, 
                align = "v", 
                rel_heights = c(3, 7)
            )
            p_right <- cowplot::plot_grid(
                ggplot2::ggplot() + ggplot2::theme_void(), q2, 
                nrow = 2, 
                align = "v", 
                rel_heights = c(3, 7)
            )
            p <- cowplot::plot_grid(
                p_left, p_right,
                align = "h", 
                nrow = 1, 
                rel_widths = c(8, 2)
            )
        }
    }
    
    return(p)
}

#' Personal ggplot2 theming function, adapted from roboto-condensed 
#' at https://github.com/hrbrmstr/hrbrthemes/
#'
#' @param base_family,base_size base font family and size
#' @param plot_title_family,plot_title_face, plot title family, face
#' @param plot_title_size,plot_title_margin, plot title size and margin
#' @param subtitle_face,subtitle_size plot subtitle family, 
#' face and size
#' @param subtitle_margin plot subtitle margin bottom (single numeric value)
#' @param strip_text_family,strip_text_face,strip_text_size facet label font 
#' family, face and size
#' @param caption_face,caption_size,caption_margin plot caption
#' family, face, size and margin
#' @param axis_title_family,axis_title_face,axis_title_size axis title font 
#' family, face and size
#' @param axis_title_just axis title font justificationk one of `[blmcrt]`
#' @param axis_text_size font size of axis text
#' @param plot_margin plot margin (specify with [ggplot2::margin])
#' @param panel_spacing panel spacing (use `unit()`)
#' @param grid_col grid color
#' @param grid panel grid (`TRUE`, `FALSE`, or a combination of 
#' `X`, `x`, `Y`, `y`)
#' @param axis_col axis color
#' @param axis add x or y axes? `TRUE`, `FALSE`, "`xy`"
#' @param ticks ticks if `TRUE` add ticks
#' @param border border if `TRUE` add border
#' @return theme A ggplot theme
#' 
#' @import ggplot2
#' @export
#' 
#' @examples
#' library(ggplot2)
#'
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   labs(x="Fuel effiency (mpg)", y="Weight (tons)",
#'        title="Seminal ggplot2 scatterplot example") +
#'   theme_ggplot2()

theme_ggplot2 <- function(
    grid = TRUE,
    border = TRUE, 
    base_family = NULL, base_size = 8,
    plot_title_family = base_family, plot_title_size = 12,
    plot_title_face = "plain", plot_title_margin = 5,
    subtitle_size = 11,
    subtitle_face = "plain", subtitle_margin = 5,
    strip_text_family = base_family, strip_text_size = 10,
    strip_text_face = "bold",
    caption_size = 9,
    caption_face = "plain", caption_margin = 3,
    axis_text_size = base_size,
    axis_title_family = base_family,
    axis_title_size = 9,
    axis_title_face = "plain",
    axis_title_just = "rt",
    panel_spacing = grid::unit(2, "lines"),
    grid_col = "#cccccc", 
    plot_margin = margin(12, 12, 12, 12),
    axis_col = "#cccccc", 
    axis = FALSE, 
    ticks = FALSE
) 
{
    
    ret <- ggplot2::theme_minimal(
        base_family = base_family, base_size = base_size
    )
    
    ret <- ret + theme(legend.background = element_blank())
    ret <- ret + theme(legend.key = element_blank())
    ret <- ret + theme(legend.position = 'bottom')
    ret <- ret + theme(legend.title = element_text(size = 10, face="bold")) 
    ret <- ret + theme(legend.text = element_text(size = 9))
    
    ret <- ret + theme(plot.margin = plot_margin)
    
    ret <- ret + theme(panel.spacing = panel_spacing)
    
    if (inherits(grid, "character") | grid == TRUE) {
        ret <- ret + theme(
            panel.grid = element_line(color = grid_col, size = 0.2), 
            panel.grid.major = element_line(color = grid_col, size = 0.2), 
            panel.grid.minor = element_line(color = grid_col, size = 0.15)
        )
        
        if (inherits(grid, "character")) {
            if (regexpr("X", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.major.x = element_blank())
            if (regexpr("Y", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.major.y = element_blank())
            if (regexpr("x", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.minor.x = element_blank())
            if (regexpr("y", grid)[1] < 0) 
                ret <- ret + theme(panel.grid.minor.y = element_blank())
        }
    } else {
        ret <- ret + theme(panel.grid = element_blank())
        ret <- ret + theme(panel.grid.major  = element_blank())
        ret <- ret + theme(panel.grid.major.x  = element_blank())
        ret <- ret + theme(panel.grid.major.y  = element_blank())
        ret <- ret + theme(panel.grid.minor  = element_blank())
        ret <- ret + theme(panel.grid.minor.x  = element_blank())
        ret <- ret + theme(panel.grid.minor.y  = element_blank())
    }
    
    if (border == TRUE) {
        ret <- ret + theme(
            panel.border = element_rect(
                colour = "black", fill = NA, size = 0.5
            )
        )
    }
    
    if (inherits(axis, "character") | axis == TRUE) {
        ret <- ret + theme(
            axis.line = element_line(color = axis_col, size = 0.15)
        )
        if (inherits(axis, "character")) {
            axis <- tolower(axis)
            if (regexpr("x", axis)[1] < 0) {
                ret <- ret + theme(axis.line.x = element_blank())
            } else {
                ret <- ret + theme(
                    axis.line.x = element_line(color = axis_col, size = 0.15)
                )
            }
            if (regexpr("y", axis)[1] < 0) {
                ret <- ret + theme(axis.line.y = element_blank())
            } else {
                ret <- ret + theme(
                    axis.line.y = element_line(color = axis_col, size = 0.15)
                )
            }
        } else {
            ret <- ret + theme(
                axis.line.x = element_line(color = axis_col, size = 0.15), 
                axis.line.y = element_line(color = axis_col, size = 0.15)
            )
        }
    } else {
        ret <- ret + theme(axis.line = element_blank())
    }
    
    if (!ticks) {
        ret <- ret + theme(axis.ticks = element_blank())
        ret <- ret + theme(axis.ticks.x = element_blank())
        ret <- ret + theme(axis.ticks.y = element_blank())
    } else {
        ret <- ret + theme(axis.ticks = element_line(size = 0.15))
        ret <- ret + theme(axis.ticks.x = element_line(size = 0.15))
        ret <- ret + theme(axis.ticks.y = element_line(size = 0.15))
        ret <- ret + theme(axis.ticks.length = grid::unit(5, "pt"))
    }
    
    xj <- switch(
        tolower(substr(axis_title_just, 1, 1)), 
        b = 0, l = 0, m = 0.5, c = 0.5, r = 1, t = 1
    )
    yj <- switch(
        tolower(substr(axis_title_just, 2, 2)), 
        b = 0, l = 0, m = 0.5, c = 0.5, r = 1, t = 1
    )
    
    ret <- ret + theme(axis.text = element_text(
        size = axis_text_size, margin = margin(t = 0, r = 0)
    ))
    ret <- ret + theme(axis.text.x = element_text(
        size = axis_text_size, margin = margin(t = 0)
    ))
    ret <- ret + theme(axis.text.y = element_text(
        size = axis_text_size, margin = margin(r = 0)
    ))
    
    ret <- ret + theme(axis.title = element_text(
        size = axis_title_size, family = axis_title_family)
    )
    ret <- ret + theme(axis.title.x = element_text(
        # hjust = xj, 
        size = axis_title_size,
        family = axis_title_family, face = axis_title_face
    ))
    ret <- ret + theme(axis.title.y = element_text(
        # hjust = yj, 
        size = axis_title_size,
        family = axis_title_family, face = axis_title_face
    ))
    ret <- ret + theme(axis.title.y.right = element_text(
        # hjust = yj, 
        size = axis_title_size, angle = 90,
        family = axis_title_family, face = axis_title_face
    ))
    
    ret <- ret + theme(strip.text = element_text(
        hjust = 0, size = strip_text_size,
        face = strip_text_face, family = strip_text_family
    ))
    
    ret <- ret + theme(plot.title = element_text(
        hjust = 0, size = plot_title_size,
        margin = margin(b = plot_title_margin),
        family = plot_title_family, face = plot_title_face
    ))
    ret <- ret + theme(plot.subtitle = element_text(
        hjust = 0, size = subtitle_size,
        margin = margin(b = subtitle_margin),
        face = subtitle_face
    ))
    ret <- ret + theme(plot.caption = element_text(
        hjust = 1, size = caption_size,
        margin = margin(t = caption_margin),
        face = caption_face
    ))
    
    ret
    
}