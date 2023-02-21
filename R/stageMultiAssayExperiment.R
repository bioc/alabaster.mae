#' Stage a dataset
#'
#' Save the metadata and annotations of a \linkS4class{MultiAssayExperiment} in a staging directory.
#' 
#' @param x A \linkS4class{MultiAssayExperiment} object or one of its subclasses.
#' @inheritParams alabaster.base::stageObject
#' @param meta.name String containing the name of the metadata file.
#' @param sm.name String containing the prefix of the sample mapping file.
#' @param sd.name String containing the prefix of the sample data file.
#'
#' @return A named list containing the metadata for this dataset.
#' The contents of \code{x} are saved into a \code{path} subdirectory inside \code{dir}.
#'
#' @details
#' \code{meta.name} is only needed to set up the output \code{path}, for consistency with the \code{\link{stageObject}} contract.
#' Callers should make sure to write the metadata to the same document by using \code{\link{.writeMetadata}} with \code{meta.only=TRUE}.
#'
#' @author Aaron Lun
#' 
#' @examples
#' # Mocking up an MAE
#' mat <- matrix(rnorm(1000), ncol=10)
#' colnames(mat) <- letters[1:10]
#' rownames(mat) <- sprintf("GENE_%i", seq_len(nrow(mat)))
#' 
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(list(counts=mat))
#' 
#' library(MultiAssayExperiment)
#' mae <- MultiAssayExperiment(list(gene=se))
#'
#' # Staging it:
#' tmp <- tempfile()
#' dir.create(tmp)
#' stageObject(mae, tmp, "dataset")
#' 
#' @export
#' @rdname stageMultiAssayExperiment
#' @aliases stageObject,SampleMapFrame-method
#' @import methods alabaster.base alabaster.se
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom MultiAssayExperiment colData sampleMap
#' @importFrom S4Vectors DataFrame
setMethod("stageObject", "MultiAssayExperiment", function(x, dir, path, child=FALSE, sm.name="sample_mapping", sd.name="sample_data", meta.name="dataset.json") {
    dir.create(file.path(dir, path), showWarnings=FALSE)

    # Saving the sample mappings.
    sm <- sampleMap(x)
    df <- DataFrame(
        sample = as.character(sm$primary),
        experiment = as.character(sm$assay),
        column = as.character(sm$colname)
    )

    map.info <- tryCatch({
        meta <- .stageObject(df, dir, file.path(path, sm.name), child=TRUE)
        .writeMetadata(dir=dir, meta)
    }, error=function(e) {
        stop("failed to stage 'sampleMap(<", class(x)[1], ">)'\n  - ", e$message)
    })

    # Saving the sample data.
    cd <- colData(x)
    if (is.null(rownames(cd)) || anyDuplicated(rownames(cd))) {
        stop("row names of 'colData(<", class(x)[1], ">)' should be non-NULL and unique")
    }

    sdpath <- file.path(path, sd.name)
    sd.info <- tryCatch({
        meta <- .stageObject(cd, dir, sdpath, child=TRUE)
        .writeMetadata(dir=dir, meta)
    }, error=function(e) {
        stop("failed to stage 'colData(<", class(x)[1], ">)'\n  - ", e$message)
    })

    # Saving the experiments.
    exp.info <- .stage_experiments(x, dir, path)

    # Saving other metadata.
    meta.info <- .processMetadata(x, dir, path, "metadata")

    list(
        `$schema`="dataset/v1.json",
        path=file.path(path, meta.name),
        dataset=list(
            experiments=exp.info,
            sample_mapping=list(resource=map.info),
            sample_data=list(resource=sd.info),
            other_data=meta.info
        ),
        is_child=child
    )
})

#' @import alabaster.base 
#' @importFrom MultiAssayExperiment experiments 
#' @importMethodsFrom alabaster.se stageObject
.stage_experiments <- function(x, dir, path) {
    exp.names <- names(experiments(x))
    if (anyDuplicated(exp.names)) {
        stop("detected duplicated experiment names in a ", class(x)[1], " object")
    }
    if (any(exp.names == "")) {
        stop("detected empty experiment name in a ", class(x)[1], " object")
    }

    collected <- list()
    for (i in seq_along(exp.names)) {
        newname <- file.path(path, paste0("experiment-", i))
        exp.name <- exp.names[i]
        se <- x[[exp.name]]
        if (is.null(colnames(se)) || anyDuplicated(colnames(se))) {
            stop("column names of '<", class(x)[1], ">[[", i, "]]' should be non-NULL and unique")
        }

        processed <- tryCatch({
            meta <- .stageObject(se, dir=dir, path=newname, child=TRUE)
            meta$path <- file.path(newname, "experiment.json")
            .writeMetadata(meta=meta, dir=dir)
        }, error=function(e) {
            stop("failed to stage experiment '", exp.name, "'\n  - ", e$message)
        })

        collected[[i]] <- list(name=exp.name, resource=processed)
    }

    collected
}
