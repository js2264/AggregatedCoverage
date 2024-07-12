#' as_tibble
#' 
#' Coerce an `CoverageExperiment` or `AggregatedCoverage` object into a `tibble`
#' 
#' @name as_tibble-methods
#' @rdname as_tibble-methods
#' @param x an `CoverageExperiment` object
#' @param ... ignored
#' @return `tibble`
#' 
#' @importFrom tidyr any_of
#' @importFrom tidyr all_of
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr as_tibble
#' @examples 
#' data(ac)
#' as_tibble(ac)
NULL

#' @rdname as_tibble-methods
#' @export

as_tibble.AggregatedCoverage <- function(x, ...) {
    tracks <- colData(x)$track
    features <- rowData(x)$features
    assays <- names(assays(x))
    w <- width(rowRanges(x)[[1]])[[1]]
    bin <- w / length(assay(x, "mean")[1, 1][[1]])
    lapply(tracks, function(t) {
        lapply(features, function(f) {
            lapply(assays(x)[assays], `[[`, f, t) |> 
                stats::setNames(assays) |> 
                dplyr::bind_cols() |> 
                dplyr::mutate(
                    coord = seq(-w/2, w/2-1, by = bin), 
                    track = factor(t, levels = tracks), 
                    features = factor(f, levels = features)
                ) |> 
                dplyr::relocate(coord) |> 
                dplyr::relocate(features) |> 
                dplyr::relocate(track)
        }) |> dplyr::bind_rows()
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::left_join(colData(x) |> as.data.frame(), by = 'track') |> 
        dplyr::mutate(.feature = features, .sample = track) |> 
        dplyr::relocate(.feature) |> 
        dplyr::relocate(.sample)
}
