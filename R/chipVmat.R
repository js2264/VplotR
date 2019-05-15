
estimateFragmentsLength <- function(reads, samples=1000, window=1000, min.shift=1, max.shift=100, mc.cores = 1, as.shift = FALSE) {
    
    # Inspired from fragmentLenDetect function in the nucleR package
    
    # Get =/- coverage of each chromosome
    chromSizes <- getChromSizes(reads)
    .doCover <- function(x) {
        lapply(
            list(pos="+", neg="-"),
            function (s)
            coverage(x[strand(x) == s, ])[[1]]
        )
    }
    chrs <- levels(seqnames(chromSizes))
    splt <- split(reads, chrs)
    cover <- lapply(splt, .doCover)
    
    # Randomly select regions in the available chromosome bounds
    chrLength <- IRanges::end(chromSizes) %>% stats::setNames(seqnames(chromSizes))
    probs <- chrLength / sum(chrLength)
    chrSample <- sample(chrs, samples, replace=TRUE, prob=probs)
    position <- round(runif(samples, max = chrLength[chrSample] - window))
    dd <- data.frame(chrSample, position)
    
    # For each sampled region, look for the shift with higher correlation
    # between strands
    shiftPos <- function(i) {
        chr <- as.character(dd[i, "chrSample"])
        sta <- as.numeric(dd[i, "position"])
        end <- sta + window
        cpos <- try(as.vector(cover[[chr]][["pos"]][sta:end]), silent=TRUE)
        cneg <- try(as.vector(cover[[chr]][["neg"]][sta:end]), silent=TRUE)
        if (class(cpos) == "try-error" | class(cneg) == "try-error") {
            return(NA)
        }
        if (sum(cpos) == 0 | sum(cneg) == 0) {
            return(NA)
        }
        cpos[is.na(cpos)] <- 0
        cneg[is.na(cneg)] <- 0
        x <- sapply(
            min.shift:max.shift,
            function (s)
            cor(cpos[1:(window + 1 - s)], cneg[s:window])
        )
        res <- which(x == max(x, na.rm=TRUE))[1]
        # We only shifted the negative strand, but both strands will be
        # shifted the half of this amount
        res <- res / 2
        # Discard NAs
        if (is.na(res) | !is.numeric(res)) {
            return(numeric(0))
        }
        return(res + min.shift - 1)
    }
    
    shift <- sapply(1:nrow(dd), shiftPos) %>% 
        median(na.rm=TRUE) %>%
        round()
    
    # Fragment length is the shift * 2 + the length of the read
    fragLen <- shift * 2 + width(reads)[1]
    
    if (as.shift) {
        return (shift)
    } else {
        return (fragLen)
    }
}
