#' geom_aggrcoverage
#'
#' #' @description
#' 
#' `geom_aggrcoverage()`
#' 
#' @name geom_aggrcoverage
#' @rdname geom_aggrcoverage
#' 
#' @param mapping mapping
#' @param data data
#' @param ... ...
#' @param ci ci
#' @param na.rm na.rm
#' @param show.legend show.legend
#' @param inherit.aes inherit.aes
#' @return A `ggplot` object`
#'
#' @import ggplot2
#'
#' @examples
#' library(rtracklayer)
#' library(plyranges)
#' library(ggplot2)
#' library(purrr)
#' TSSs_bed <- system.file("extdata", "TSSs.bed", package = "tidyCoverage")
#' features <- list(
#'     TSS_fwd = import(TSSs_bed) |> filter(strand == '+'), 
#'     TSS_rev = import(TSSs_bed) |> filter(strand == '-')
#' )
#' tracks <- list(
#'     RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
#'     RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
#' ) |> map(import, as = 'Rle')
#' df <- CoverageExperiment(tracks, features, width = 5000, ignore.strand = FALSE) |> 
#'  aggregate() |> 
#'  as_tibble()
NULL 

GeomAggrCoverage <- ggplot2::ggproto("GeomAggrCoverage", ggplot2::Geom,
    setup_params = function(data, params) {
        params$ci <- params$ci
        params
    },
    extra_params = c("na.rm"),
    required_aes = c("x", "y", "ymin", "ymax"), 
    default_aes = ggplot2::aes(
        colour = "black", 
        linewidth = 1, 
        linetype = 1, 
        alpha = 0.4
    ),
    draw_group = function(data, params, coord, ci = TRUE, ...) {
        forLine <- transform(data, alpha = 1)
        forRibbon <- transform(
            data, 
            alpha = data$alpha, 
            fill = data$colour, 
            colour = NA
        )
        grid::gList(
            if (ci) ggplot2::GeomRibbon$draw_panel(forRibbon, params, coord, ...),
            ggplot2::GeomLine$draw_panel(forLine, params, coord, ...)
        )
    }, 
    draw_key = ggplot2::draw_key_smooth
)

#' @rdname geom_aggrcoverage
#' @export

geom_aggrcoverage <- function(
    mapping = NULL, 
    data = NULL, 
    ..., 
    ci = TRUE, 
    na.rm = FALSE, 
    show.legend = NA, 
    inherit.aes = TRUE
) {
    m <- ggplot2::aes(x = coord, y = mean, ymin = ci_low, ymax = ci_high, group = interaction(.sample, .feature))
    if (!is.null(mapping)) m <- utils::modifyList(m, mapping)

    ggplot2::layer(
        data = data, 
        mapping = m,  
        stat = "identity", 
        geom = GeomAggrCoverage, 
        position = "identity", 
        show.legend = show.legend, 
        inherit.aes = inherit.aes,
        params = list(na.rm = na.rm, ci = ci, ...)
    )
}

