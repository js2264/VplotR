importBigWig <- function(files, verbose = T) {
    # Check that all the files exist
    stopifnot(all(file.exists(files)))
    list.bw <- parallel::mclapply(files, function(FILE) {
        if (verbose) message('Importing ', FILE, '...')
        rtracklayer::import(FILE, as = 'Rle')
    }, mc.cores = min(10, length(files)))
    class(list.bw) <- c("listBigWig", class(list.bw))
    return(list.bw)
}
