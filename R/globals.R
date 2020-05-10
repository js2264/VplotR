utils::globalVariables(c(
    ".",
    "x",
    "x1",
    "x2",
    "y1", 
    "y2",
    "xmin", 
    "xmax", 
    "ymin",
    "ymax",
    "value",
    "type",
    "fisher.test",
    "colorRampPalette",
    "Var1",
    "Var2"
))

COLORSCALE_VMAT <- c(
    colorRampPalette(
        rev(RColorBrewer::brewer.pal('Spectral', n = 10))[seq(1, 5)]
    )(30), 
    colorRampPalette(
        rev(RColorBrewer::brewer.pal('Spectral', n = 10))[seq(6, 10)]
    )(30)
)