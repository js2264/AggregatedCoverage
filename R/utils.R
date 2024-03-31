.set_seqinfo <- function(tracks) {
    n <- names(tracks)
    sls <- lapply(tracks, lengths)
    has_seqlengths <- any(!is.na(lapply(sls, sum)))
    if (!has_seqlengths) {
        sis <- lapply(tracks, function(rle) {
            Seqinfo(
                names(lengths(rle)), 
                lengths(rle)
            )
        }) |> unique()
        if (length(unique(sis)) > 1) 
            warning("EXPERIMENTAL: More than 1 seqinfo inferred from the tracks. Using `seqinfo` from the first track.")
        si <- sis[[1]]
    } 
    else {
        sl <- sls[[which(!is.na(lapply(sls, sum)))[[1]]]]
        si <- list(Seqinfo(
            names(sl), 
            sl
        ))
    }
    tracks <- lapply(tracks, function(rle) {
        GenomeInfoDb::seqinfo(rle) <- si
        return(rle)
    })
    names(tracks) <- n
    return(tracks)
}

.set_seqinfo_bwfl <- function(tracks) {
    n <- names(tracks)
    sis <- lapply(tracks, seqinfo) |> unique()
    if (length(sis) > 1) 
        warning("EXPERIMENTAL: More than 1 seqinfo inferred from the tracks. Using `seqinfo` from the first track.")
    si <- sis[[1]]
    return(tracks)
} 

.resize_granges <- function(gr, width, seqinfo) {
    GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(seqinfo)
    GenomeInfoDb::seqinfo(gr) <- seqinfo
    gr <- suppressWarnings(resize(gr, fix = 'center', width = width))
    w <- width
    gr <- trim(gr)
    gr <- gr[width(gr) == w]
    sort(gr)
}

#' @importFrom IRanges NumericList
#' @importFrom stats qt

.compute_cov <- function(scores, gr, center, scale, ignore.strand = TRUE) {
    m <- as.matrix(scores)
    if (!ignore.strand) {
        m <- data.frame(m)
        which.flip <- which(as.vector(strand(gr)) == '-')
        m[which.flip, ] <- rev(m[which.flip, ])
        m <- as.matrix(m)
    }
    m <- t(scale(t(m), center = center, scale = scale))
    m
}

.summarize_cov <- function(scores) {
    mean <- colMeans(scores, na.rm = TRUE)
    median <- apply(scores, 2, median, na.rm = TRUE)
    min <- apply(scores, 2, min, na.rm = TRUE)
    max <- apply(scores, 2, max, na.rm = TRUE)
    sd <- apply(scores, 2, sd, na.rm = TRUE)
    se <- sd/sqrt(nrow(scores))
    ci_low <- mean - stats::qt(1 - (0.05 / 2), nrow(scores) - 1) * se
    ci_high <- mean + stats::qt(1 - (0.05 / 2), nrow(scores) - 1) * se
    data.frame(
        coord = seq(-ncol(scores)/2, ncol(scores)/2-1, 1), 
        mean = mean, median = median, min = min, max = max, 
        sd = sd, se = se, ci_low = ci_low, ci_high = ci_high
    )
}

.coarsen_vec <- function(x, bin, FUN, ...) {
    if ({length(x) %% bin} != 0) stop(
        "The length of the provided vector should be divided by the bin size without remainder.
Please adjust `bin` argument."
    )
    stats::aggregate(
        x, 
        by = list(rep(seq(1, length(x)/bin), each = bin)), 
        FUN = FUN, ...
    )$x
}

.coarsen_mat <- function(x, bin, FUN = mean, ...) {
    if ({ncol(x) %% bin} != 0) stop(
        "The column number of the provided matrix should be divided by the window size without remainder.
Please adjust `window` argument."
    )
    apply(x, 1, function(vec) {
        stats::aggregate(
            vec, 
            by = list(rep(seq(1, length(vec)/bin), each = bin)), 
            FUN = FUN, ...
        )$x
    }) |> t()
}
