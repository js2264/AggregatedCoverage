test_that("CoverageExperiment works", {

    library(rtracklayer)
    library(plyranges)
    library(purrr)

    # ~~~~~~~~~~~~~~~ Import genomic features into a named list ~~~~~~~~~~~~~~~ #
    features <- list(
        TSSs = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
        `Convergent transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage")
    ) |> map(import) |> map(filter, strand == '+') 

    # ~~~~~~~~~~~~~~~ Import coverage tracks into a named list ~~~~~~~~~~~~~~~ #
    tracks <- list(
        Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage"), 
        RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
        RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
    ) |> map(import, as = 'Rle')

    ## ~~~~~~~~~~~~~~~ TEST GETTERS ~~~~~~~~~~~~~~~ ##
    expect_s4_class(
        CE <- CoverageExperiment(
            tracks, features, width = 100, scale = TRUE, center = TRUE
        ), 
        "CoverageExperiment"
    )
    expect_true(all(dim(CE) == c(2, 3)))
    expect_true(names(assays(CE)) == 'coverage')
    expect_true(length(assay(CE, 'coverage')) == 6L)
    expect_true(length(assay(CE, 'coverage')[[1]]) == 87300L)
    expect_true(length(assay(CE, 'coverage')[[2]]) == 46900L)
    expect_true(length(assay(CE, 'coverage')[[3]]) == 87300L)
    expect_true(length(assay(CE, 'coverage')[[4]]) == 46900L)
    expect_true(length(assay(CE, 'coverage')[[5]]) == 87300L)
    expect_true(length(assay(CE, 'coverage')[[6]]) == 46900L)
    expect_true(all(dim(assay(CE, 'coverage')[[1]]) == c(873, 100)))
    expect_true(all(dim(assay(CE, 'coverage')[[2]]) == c(469, 100)))
    expect_true(all(rowData(CE)$n == c(873, 469)))
    expect_true(all(colnames(CE) == c('Scc1', 'RNA_fwd', 'RNA_rev')))

})

test_that("other CoverageExperiment methods work", {
    # ~~~~~~~~ BigWigFile + GRanges
    x <- BigWigFile(system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"))
    y <- import(system.file("extdata", "TSSs.bed", package = "tidyCoverage"))[1:200]
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")
    expect_s4_class(CoverageExperiment(x, y), "CoverageExperiment")

    # ~~~~~~~~ BigWigFile + GRangesList
    y <- GRangesList(tss = y, tss2 = y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ BigWigFile + list
    y <- as.list(y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ BigWigFileList + GRanges
    x <- BigWigFileList(list(
        track = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"), 
        track2 = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage")
    ))
    y <- import(system.file("extdata", "TSSs.bed", package = "tidyCoverage"))
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ BigWigFileList + GRangesList
    y <- GRangesList(tss = y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ BigWigFileList + list
    y <- as.list(y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ RleList + GRanges
    x <- import(system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"), as = 'Rle')
    y <- import(system.file("extdata", "TSSs.bed", package = "tidyCoverage"))
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")
    expect_s4_class(CoverageExperiment(x, y), "CoverageExperiment")

    # ~~~~~~~~ RleList + GRangesList
    y <- GRangesList(tss = y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ RleList + list
    y <- as.list(y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ list + GRanges
    x <- list(
        track = import(system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"), as = 'Rle'), 
        track2 = import(system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"), as = 'Rle')
    )
    y <- import(system.file("extdata", "TSSs.bed", package = "tidyCoverage"))[1:200]
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ list + GRangesList
    y <- GRangesList(tss = y, tss2 = y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")

    # ~~~~~~~~ list + list
    y <- as.list(y)
    expect_s4_class(CoverageExperiment(x, y, width = 100), "CoverageExperiment")
})

test_that("aggregate works", {

    features <- list(
        TSSs = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
        `Convergent transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage")
    ) |> map(import) |> map(filter, strand == '+') 
    tracks <- list(
        Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage"), 
        RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
        RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
    ) |> map(import, as = 'Rle')
    CE <- CoverageExperiment(
        tracks, features, width = 100, scale = TRUE, center = TRUE
    )

    expect_s4_class(
        AC <- aggregate(CE), 
        "AggregatedCoverage"
    )
    expect_true(all(dim(AC) == dim(CE)))
    expect_true(all(names(assays(AC)) == c("mean", "median", "min", "max", "sd", "se", "ci_low", "ci_high")))
    expect_true(length(assay(AC, 'mean')) == 6L)
    expect_true(length(assay(AC, 'mean')[[1]]) == 100L)
    expect_true(length(assay(AC, 'mean')[[2]]) == 100L)
    expect_true(all(rowData(AC)$n == c(873, 469)))
    expect_true(all(colnames(AC) == c('Scc1', 'RNA_fwd', 'RNA_rev')))
})

test_that("print/show work", {
    features <- list(
        TSSs = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
        `Convergent transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage")
    ) |> map(import) |> map(filter, strand == '+') 
    tracks <- list(
        Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage"), 
        RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
        RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
    ) |> map(import, as = 'Rle')
    CE <- CoverageExperiment(
        tracks, features, width = 100, scale = TRUE, center = TRUE
    )
    AC <- aggregate(CE)

    options(restore_SummarizedExperiment_show = TRUE)
    expect_no_error(show(CE))
    expect_no_error(show(AC))
    expect_no_error(print(CE))
    expect_no_error(print(AC))
    options(restore_SummarizedExperiment_show = FALSE)
    expect_no_error(show(CE))
    expect_no_error(show(AC))
    expect_no_error(print(CE))
    expect_no_error(print(AC))
})

test_that("as_tibble methods work", {
    features <- list(
        TSSs = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
        `Convergent transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage")
    ) |> map(import) |> map(filter, strand == '+') 
    tracks <- list(
        Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage"), 
        RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
        RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
    ) |> map(import, as = 'Rle')
    CE <- CoverageExperiment(
        tracks, features, width = 100, scale = TRUE, center = TRUE
    )
    AC <- aggregate(CE)

    expect_true(all(dim(as_tibble(AC)) == c(600, 13)))
    expect_true(all(colnames(as_tibble(AC)) == c(".sample", ".feature", "track", "features", "coord", "mean", "median", "min", "max", "sd", "se", "ci_low", "ci_high")))
})

test_that("expand method works", {
    features <- list(
        TSSs = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
        `Convergent transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage")
    ) |> map(import) |> map(filter, strand == '+') 
    tracks <- list(
        Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage"), 
        RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
        RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
    ) |> map(import, as = 'Rle')
    CE <- CoverageExperiment(
        tracks, features, width = 100, scale = TRUE, center = TRUE
    )
    expect_true(all(dim(expand(CE)) == c(402600, 8)))
    expect_true(all(colnames(expand(CE)) == c("track", "features", "chr", "ranges", "strand", "coord", "coverage", "coord.scaled")))
})

test_that("coarsen method works", {
    features <- list(
        TSSs = system.file("extdata", "TSSs.bed", package = "tidyCoverage"),
        `Convergent transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "tidyCoverage")
    ) |> map(import) |> map(filter, strand == '+') 
    tracks <- list(
        Scc1 = system.file("extdata", "Scc1.bw", package = "tidyCoverage"), 
        RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
        RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
    ) |> map(import, as = 'Rle')
    CE <- CoverageExperiment(
        tracks, features, width = 100, scale = TRUE, center = TRUE
    )
    AC <- aggregate(CE)

    expect_s4_class(
        CCE <- coarsen(CE, 10), 
        "CoverageExperiment"
    )
    expect_true(all(dim(assay(CCE, 'coverage')[[1]]) == c(873, 10)))
    expect_true(all(dim(assay(CCE, 'coverage')[[2]]) == c(469, 10)))
    expect_true(all(rowData(CCE)$n == c(873, 469)))
})
