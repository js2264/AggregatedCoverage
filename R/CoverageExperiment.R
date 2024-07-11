#' CoverageExperiment
#'
#' #' @description
#' 
#' `CoverageExperiment` objects store coverages for individual
#' tracks over different sets of features. The `coverage` 
#' assay contains a separate matrix for each combination of 
#' track x features. `CoverageExperiment` objects are instantiated 
#' using the `CoverageExperiment()` #' function, and can be 
#' coarsened using the `coarsen()` function.
#'
#' @name CoverageExperiment
#' @rdname CoverageExperiment
#' 
#' @param tracks A genomic track imported as a `RleList` or a *named* list of 
#' genomic tracks.
#' @param features A set of features imported as `GRanges` or a *named* 
#' `GRangesList`. 
#' @param width Width to resize each set of genomic features
#' @param ignore.strand Logical, whether to not take the features strand
#' information
#' @param scale,center Logical, whether to scale and/or center tracks prior to 
#' summarization
#' @param BPPARAM Passed to BiocParallel
#' @param ... Passed to the relevant method
#' @return A `CoverageExperiment` object
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import SummarizedExperiment
#' @import methods
#' @importFrom rtracklayer BigWigFile
#' @importFrom rtracklayer BigWigFileList
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom BiocParallel bpparam
#' @importFrom BiocParallel bplapply
#'
#' @examples
#' library(rtracklayer)
#' library(purrr)
#' library(plyranges)
#' TSSs_bed <- system.file("extdata", "TSSs.bed", package = "tidyCoverage")
#' features <- import(TSSs_bed) |> filter(strand == '+')
#' 
#' #############################################################################
#' ## 1. Creating a `CoverageExperiment` object from a single BigWigFile
#' #############################################################################
#' 
#' RNA_fwd <- system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage")
#' tracks <- BigWigFile(RNA_fwd)
#' CoverageExperiment(tracks, features, width = 5000)
#' 
#' #############################################################################
#' ## 2. Creating a `CoverageExperiment` object from a BigWigFileList
#' #############################################################################
#' 
#' RNA_rev <- system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
#' tracks <- BigWigFileList(list(RNA_fwd = RNA_fwd, RNA_rev = RNA_rev))
#' CoverageExperiment(tracks, features, width = 5000)
#' 
#' #############################################################################
#' ## 3. Creating a `CoverageExperiment` object from imported bigwig files
#' #############################################################################
#' 
#' tracks <- list(
#'     RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
#'     RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
#' ) |> map(import, as = 'Rle')
#' CoverageExperiment(tracks, features, width = 5000)
#' 
#' #############################################################################
#' ## 4. Correct for strandness when recovering coverage
#' #############################################################################
#' 
#' TSSs_bed <- system.file("extdata", "TSSs.bed", package = "tidyCoverage")
#' features <- list(
#'     TSS_fwd = import(TSSs_bed) |> filter(strand == '+'), 
#'     TSS_rev = import(TSSs_bed) |> filter(strand == '-')
#' )
#' tracks <- list(
#'     RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
#'     RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
#' ) |> map(import, as = 'Rle')
#' CoverageExperiment(tracks, features, width = 5000, ignore.strand = FALSE)
#' 
NULL

################################################################################
################################################################################
########################  MAIN METHOD FOR BIGWIGS  #############################
################################################################################
################################################################################

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "BigWigFileList", features = "GRangesList"), 
    function(
        tracks, features, width = NULL, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        window = 1, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        ## Check that input args are valid
        stopifnot(length(unique(names(tracks))) == length(tracks)) # named tracks list
        stopifnot(length(unique(names(features))) == length(features)) # named features GRangesList

        ## Extend and filter features
        tracks <- .set_seqinfo_bwfl(tracks)
        si <- seqinfo(tracks[[1]])
        if (is.null(width)) .check_granges_widths(features)
        features <- lapply(features, .resize_granges, width, si) 

        ## Prepare cData and rData
        cData <- data.frame(
            track = names(tracks)
        )
        rData <- GRangesList(features)
        mcols(rData) <- data.frame(
            features = names(features), 
            n = lengths(features)
        )

        ## Extract coverage scores
        combs <- expand.grid(names(tracks), names(features)) |> 
            setNames(c("tracks", "features"))
        l <- BiocParallel::bplapply(seq_len(nrow(combs)), function(K) {
            t <- combs[K, "tracks"]
            f <- combs[K, "features"]
            scores <- rtracklayer::import(
                tracks[[t]], 
                selection = rtracklayer::BigWigSelection(ranges = features[[f]]), 
                as = "NumericList", 
                format = "bigWig"
            )
            m <- .compute_cov(
                scores, features[[f]], 
                scale = scale, center = center, ignore.strand = ignore.strand
            )
            if (window > 1) {
                m <- .coarsen_mat(m, window)
            }
            return(m)
        }, BPPARAM = BPPARAM)
        names(l) <- paste(combs$tracks, combs$features, sep = '^')

        ## Fill out different assay matrices
        m <- matrix(
            list(), 
            nrow = length(features), ncol = length(tracks)
        )
        colnames(m) <- names(tracks)
        rownames(m) <- names(features)
        for (t in names(tracks)) {
            for (f in names(features)) {
                name <- paste(t, f, sep = '^')
                m[f, t][[1]] <- l[[name]]
            }
        }
        l_assays <- list(coverage = m)

        ## Instantiate and return the CoverageExperiment final object
        CE <- methods::new(
            "CoverageExperiment",
            SummarizedExperiment::SummarizedExperiment(
                rowRanges = rData,
                colData = cData,
                assays = l_assays
            )
        )
        return(CE)
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "BigWigFileList", features = "GRanges"), 
    function(
        tracks, features, ...
    ) {
        features <- GRangesList(features = features)
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "BigWigFileList", features = "list"), 
    function(
        tracks, features, ...
    ) {
        features <- as(features, 'GRangesList')
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @importFrom BiocIO resource
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "BigWigFile", features = "GRangesList"), 
    function(
        tracks, features, ...
    ) {
        tracks <- BigWigFileList(list(track = BiocIO::resource(tracks)))
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "BigWigFile", features = "GRanges"), 
    function(
        tracks, features, ...
    ) {
        tracks <- BigWigFileList(list(track = BiocIO::resource(tracks)))
        features <- GRangesList(features = features)
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "BigWigFile", features = "list"), 
    function(
        tracks, features, ...
    ) {
        tracks <- BigWigFileList(list(track = BiocIO::resource(tracks)))
        features <- as(features, 'GRangesList')
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

################################################################################
################################################################################
########################  MAIN METHOD FOR RleLists  ############################
################################################################################
################################################################################

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "list", features = "GRangesList"), 
    function(
        tracks, features, width = NULL, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        window = 1, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        ## Check that input args are valid
        stopifnot(length(unique(names(tracks))) == length(tracks)) # named tracks list
        stopifnot(length(unique(names(features))) == length(features)) # named features GRangesList

        ## Extend and filter features
        tracks <- .set_seqinfo(tracks)
        si <- seqinfo(tracks[[1]])[[1]]
        if (is.null(width)) .check_granges_widths(features)
        features <- lapply(features, .resize_granges, width, si) 

        ## Prepare cData and rData
        cData <- data.frame(
            track = names(tracks)
        )
        rData <- GRangesList(features)
        mcols(rData) <- data.frame(
            features = names(features),
            n = lengths(features)
        )

        ## Extract coverage scores
        combs <- expand.grid(names(tracks), names(features)) |> 
            setNames(c("tracks", "features"))
        l <- BiocParallel::bplapply(seq_len(nrow(combs)), function(K) {
            t <- combs[K, "tracks"]
            f <- combs[K, "features"]
            scores <- IRanges::NumericList(tracks[[t]][features[[f]]])
            m <- .compute_cov(
                scores, features[[f]], 
                scale = scale, center = center, ignore.strand = ignore.strand
            )
            if (window > 1) {
                m <- .coarsen_mat(m, window)
            }
            return(m)
        }, BPPARAM = BPPARAM)
        names(l) <- paste(combs$tracks, combs$features, sep = '^')

        ## Fill out different assay matrices
        m <- matrix(
            list(), 
            nrow = length(features), ncol = length(tracks)
        )
        colnames(m) <- names(tracks)
        rownames(m) <- names(features)
        for (t in names(tracks)) {
            for (f in names(features)) {
                name <- paste(t, f, sep = '^')
                m[f, t][[1]] <- l[[name]]
            }
        }
        l_assays <- list(coverage = m)

        ## Instantiate and return the CoverageExperiment final object
        CE <- methods::new(
            "CoverageExperiment",
            SummarizedExperiment::SummarizedExperiment(
                rowRanges = rData,
                colData = cData,
                assays = l_assays
            )
        )
        return(CE)
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "list", features = "GRanges"), 
    function(
        tracks, features, ...
    ) {
        features <- GRangesList(features = features)
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "list", features = "list"), 
    function(
        tracks, features, ...
    ) {
        features <- as(features, 'GRangesList')
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "RleList", features = "GRangesList"), 
    function(
        tracks, features, ...
    ) {
        tracks <- list(track = tracks)
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "RleList", features = "GRanges"), 
    function(
        tracks, features, ...
    ) {
        tracks <- list(track = tracks)
        features <- GRangesList(features = features)
        CoverageExperiment(
            tracks, features, ...
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "RleList", features = "list"), 
    function(
        tracks, features, ...
    ) {
        tracks <- list(track = tracks)
        features <- as(features, 'GRangesList')
        CoverageExperiment(
            tracks, features, ...
        )
    }
)
