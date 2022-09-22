# This tests the stageMultiAssayExperiment function.
# library(testthat); library(alabaster.mae); source("test-stage-mae.R")

set.seed(100)
mat <- matrix(rpois(2000, 10), ncol=10)
colnames(mat) <- paste0("SAMPLE_", seq_len(ncol(mat)))

library(SummarizedExperiment)
se <- SummarizedExperiment(list(counts=mat, cpm=mat/10))
se$stuff <- LETTERS[1:10]
se$blah <- runif(10)
rowData(se)$whee <- runif(nrow(se))

library(MultiAssayExperiment)
mae <- MultiAssayExperiment(list(rnaseq=se))

test_that("stageMultiAssayExperiment works as expected", {
    temp <- tempfile()
    dir.create(temp)
    info <- stageObject(mae, temp, path="thing")

    map.df <- read.csv(file.path(temp, info$dataset$sample_mapping$resource$path), stringsAsFactors=FALSE)
    transformed <- data.frame(assay = map.df$experiment, primary = map.df$sample, colname = map.df$column)
    samp <- as.data.frame(sampleMap(mae))
    samp$assay <- as.character(samp$assay)
    expect_identical(transformed, samp)

    cdf <- read.csv(file.path(temp, info$dataset$sample_data$resource$path), stringsAsFactors=FALSE)
    expect_identical(cdf[,1], rownames(colData(mae)))

    # Assays are saved correctly.
    expect_true(file.exists(file.path(temp, "thing", "experiment-1", "assay-1", "array.h5")))
    expect_true(file.exists(file.path(temp, "thing", "experiment-1", "assay-2", "array.h5")))

    # Round trip works correctly.
    out <- loadMultiAssayExperiment(info, temp)
    expect_identical(sampleMap(mae), sampleMap(out))
    expect_identical(colData(mae), colData(out))
    expect_identical(names(experiments(out)), names(experiments(out)))
})

test_that("stageMultiAssayExperiment works as expected with non-empty colData", {
    colData(mae)$foo <- runif(10)

    temp <- tempfile()
    dir.create(temp)
    info <- stageObject(mae, temp, path="thing")

    cdf <- read.csv(file.path(temp, info$dataset$sample_data$resource$path), row.names=1)
    expect_equal(cdf, as.data.frame(colData(mae)))

    # Round trip works correctly.
    out <- loadMultiAssayExperiment(info, temp)
    expect_equal(colData(mae), colData(out))
})

test_that("stageMultiAssayExperiment works as expected with non-trivial sample mappings", {
    temp <- tempfile()
    dir.create(temp)
    mae <- MultiAssayExperiment(list(rnaseq=se[,1:5], chipseq=se[1:20,3:10]))
    info <- stageObject(mae, temp, "whee")

    map.df <- read.csv(file.path(temp, info$dataset$sample_mapping$resource$path))
    transformed <- data.frame(assay = map.df$experiment, primary = map.df$sample, colname = map.df$column)
    samp <- as.data.frame(sampleMap(mae))
    samp$assay <- as.character(samp$assay)
    expect_equal(transformed, samp)

    expect_identical(info$dataset$experiments[[1]]$name, "rnaseq")
    expect_identical(info$dataset$experiments[[2]]$name, "chipseq")

    # Round trip works correctly.
    out <- loadMultiAssayExperiment(info, temp)

    s1 <- sampleMap(mae)
    s1$assay <- as.character(s1$assay)
    s2 <- sampleMap(out)
    s2$assay <- as.character(s2$assay)
    expect_identical(s1, s2)

    expect_identical(names(experiments(out)), c("rnaseq", "chipseq"))
})

test_that("stageMultiAssayExperiment respects entries in the metadata", {
    temp <- tempfile()
    dir.create(temp)
    metadata(mae)$WHEE <- "foo"
    info <- stageObject(mae, temp, path="thing")

    # Something was saved here.
    lines <- readLines(file.path(temp, info$dataset$other_data$resource$path))
    expect_true(any(grepl("WHEE", lines)))
    expect_true(any(grepl("foo", lines)))

    # Loading recovers it.
    mae2 <- loadMultiAssayExperiment(info, temp)
    expect_identical(metadata(mae2)$WHEE, "foo")
})

test_that("loadMultiAssayExperiment works as expected with restricted experiments", {
    temp <- tempfile()
    dir.create(temp)
    mae <- MultiAssayExperiment(list(rnaseq=se[,1:5], chipseq=se[1:20,3:10]))
    info <- stageObject(mae, temp, "whee")

    # Single loads work as expected.
    out <- loadMultiAssayExperiment(info, temp, experiments="rnaseq")
    expect_identical(names(experiments(out)), "rnaseq")
    expect_equal(colData(out[[1]]), colData(se)[1:5,])
    expect_equal(rowData(out[[1]]), rowData(se))

    out <- loadMultiAssayExperiment(info, temp, experiments=2)
    expect_identical(names(experiments(out)), "chipseq")
    expect_equal(colData(out[[1]]), colData(se)[3:10,])
    expect_equal(rowData(out[[1]]), rowData(se)[1:20,,drop=FALSE])

    # Multiple loads work as expected.
    out <- loadMultiAssayExperiment(info, temp, experiments=1:2)
    expect_identical(names(experiments(out)), c("rnaseq", "chipseq"))

    out <- loadMultiAssayExperiment(info, temp, experiments=c("chipseq", "rnaseq"))
    expect_identical(names(experiments(out)), c("chipseq", "rnaseq"))
})

test_that("loadMultiAssayExperiment works as expected with parallelization", {
    temp <- tempfile()
    dir.create(temp)
    mae <- MultiAssayExperiment(list(rnaseq=se[,1:5], chipseq=se[1:20,3:10]))
    info <- stageObject(mae, temp, "whee")

    serial <- loadMultiAssayExperiment(info, temp)
    parallel <- loadMultiAssayExperiment(info, temp, BPPARAM=BiocParallel::SnowParam())
    expect_identical(serial, parallel)
})

test_that("loadMultiAssayExperiment allows us to request no-nesting", {
    mae$yyy <- "here"
    colData(mae)[["xxx"]] <- DataFrame(stuff=1:10)
                  
    temp <- tempfile()
    dir.create(temp)
    info <- stageObject(mae, temp, "whee")

    regular <- loadMultiAssayExperiment(info, temp)
    expect_identical(regular$yyy, mae$yyy)
    expect_s4_class(regular$xxx, "DataFrame")
    expect_identical(regular$xxx$stuff, 1:10)

    nonest <- loadMultiAssayExperiment(info, temp, include.nested=FALSE)
    expect_identical(nonest$yyy, mae$yyy)
    expect_null(nonest$xxx)
})
