image.scale <- function(breaks, col = heat.colors(12), z= NULL, zlim = NULL, axis.pos=1, add.axis=TRUE, ...) {
     if(!missing(breaks)){
      if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
     }
     if(missing(breaks) & !is.null(zlim)){
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
     if(missing(breaks) & is.null(zlim)){
      zlim <- range(z, na.rm=TRUE)
      zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
      zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
      breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
     }
     poly <- vector(mode="list", length(col))
     for(i in seq(poly)){
      poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
     }
     if(axis.pos %in% c(1,3)){ylim<-c(0,1); xlim<-range(breaks)}
     if(axis.pos %in% c(2,4)){ylim<-range(breaks); xlim<-c(0,1)}
     plot(1,1,t="n",ylim=ylim, xlim=xlim, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i", ...)
     for(i in seq(poly)){
      if(axis.pos %in% c(1,3)){
       polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
      }
      if(axis.pos %in% c(2,4)){
       polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
      }
     }
     box()
     if(add.axis) {axis(axis.pos)}
}
