.onLoad <- function(libname, pkgname) {
    registerReadObjectFunction("multi_sample_dataset", readMultiAssayExperiment)
}

.onUnload <- function(libname, pkgname) {
    registerReadObjectFunction("multi_sample_dataset", NULL)
}
