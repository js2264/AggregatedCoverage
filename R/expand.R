#' Expand a CoverageExperiment object
#' 
#' @description
#' A `CoverageExperiment` object can be coerced into a `tibble` using the 
#' `tidySummarizedExperiment` package, but this will not turn 
#' each coverage matrix into a "long" format. The `expand` function 
#' provided here allows one to coerce a `CoverageExperiment`
#' object into a long data frame, and adds the `ranges` 
#' and `seqnames` to the resulting `tibble`. 
#' 
#' @name expand,CoverageExperiment
#' @aliases expand,CoverageExperiment-method
#' @rdname expand
#' 
#' @param data a `CoverageExperiment` object
#' @param ...,.name_repair ignored
#' @return a `tibble` object
#' 
#' @importFrom tidyr expand
#' @export
#' @examples 
#' data(ce)
#' ce
#' 
#' expand(ce)

expand.CoverageExperiment <- function(data, ..., .name_repair = NULL) {
    tracks <- colData(data)$track
    features <- rowData(data)$features
    w <- width(rowRanges(data)[[1]])[[1]]
    bin <- w / ncol(assay(data, "coverage")[1, 1][[1]])
    df <- lapply(features, function(f) {
        rr <- rowRanges(data)[[f]]
        rrdf <- as.data.frame(rr)
        coord <- lapply(seq_len(nrow(rrdf)), function(K) seq(rrdf[K, 'start'], rrdf[K, 'end'], by = bin)) |> unlist()
        coord.scaled <- lapply(seq_len(nrow(rrdf)), function(K) seq(-w/2, w/2-1, by = bin)) |> unlist()
        lapply(tracks, function(t) {
            m <- assay(data, "coverage")[f, t][[1]] |> 
                as.data.frame() |>
                dplyr::mutate(
                    track = factor(t, levels = tracks), 
                    features = factor(f, levels = features), 
                    chr = as.vector(seqnames(rr)), 
                    ranges = as.character(rr), 
                    strand = as.vector(strand(rr))
                ) |> 
                dplyr::relocate(tidyr::all_of(c("track", "features", "chr", "ranges", "strand")))
            d <- tidyr::pivot_longer(
                m, 
                !tidyr::any_of(c("track", "features", "chr", "ranges", "strand", "coord.scaled")), 
                names_to = "coord", values_to = "coverage"
            )
            d$coord <- coord
            d$coord.scaled <- coord.scaled
            d
        }) |> dplyr::bind_rows()
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::left_join(colData(data) |> as.data.frame(), by = 'track') |> 
        dplyr::group_by(track, features, ranges)
    return(df)
}
