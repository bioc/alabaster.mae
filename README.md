# Save `MultiAssayExperiment`s to file

The **alabaster.mae** package implements methods for saving and loading `MultiAssayExperiment` objects under the **alabaster** framework.
It provides a language-agnostic method for serializing the sample mappings, sample data, and the various `SummarizedExperiment` instances nested within the MAE.
To get started, install the package and its dependencies from GitHub:

```r
devtools::install_github("ArtifactDB/alabaster.schemas")
devtools::install_github("ArtifactDB/alabaster.base")
devtools::install_github("ArtifactDB/alabaster.ranges")
devtools::install_github("ArtifactDB/alabaster.matrix")
devtools::install_github("ArtifactDB/alabaster.se")
devtools::install_github("ArtifactDB/alabaster.mae")
```

To demonstrate, let's create a mildly complicated MAE containing RNA-seq and ChIP-seq data with partial overlaps:

```r
library(SummarizedExperiment)
rna.counts <- matrix(rpois(60, 10), ncol=6)
colnames(rna.counts) <- c("disease1", "disease2", "disease3", "control1", "control2", "control3")
rownames(rna.counts) <- c("ENSMUSG00000000001", "ENSMUSG00000000003", "ENSMUSG00000000028", 
    "ENSMUSG00000000031", "ENSMUSG00000000037", "ENSMUSG00000000049",  "ENSMUSG00000000056", 
    "ENSMUSG00000000058", "ENSMUSG00000000078",  "ENSMUSG00000000085")
rna.se <- SummarizedExperiment(list(counts=rna.counts))
colData(rna.se)$disease <- rep(c("disease", "control"), each=3)

chip.counts <- matrix(rpois(100, 10), ncol=4)
colnames(chip.counts) <- c("disease1", "disease2", "control1", "control3")
chip.peaks <- GRanges("chr1", IRanges(1:25*100+1, 1:25*100+100))
chip.se <- SummarizedExperiment(list(counts=chip.counts), rowRanges=chip.peaks)

library(MultiAssayExperiment)
mapping <- DataFrame(
    primary = c(colnames(rna.se), colnames(chip.se)), # sample identifiers
    assay = rep(c("rnaseq", "chipseq"), c(ncol(rna.se), ncol(chip.se))), # experiment name
    colname = c(colnames(rna.se), colnames(chip.se)) # column names inside each experiment
)
mae <- MultiAssayExperiment(list(rnaseq=rna.se, chipseq=chip.se), sampleMap=mapping)
```

Now we can just save it to file:

```r
library(alabaster.mae)
tmp <- tempfile()
dir.create(tmp)
meta <- stageObject(mae, tmp, "my_mae")
meta[["$schema"]]
## [1] "dataset/v1.json"
```

And easily load it back:

```r
roundtrip <- loadObject(meta, tmp)
roundtrip
## A MultiAssayExperiment object of 2 listed
##  experiments with user-defined names and respective classes.
##  Containing an ExperimentList class object of length 2:
##  [1] rnaseq: SummarizedExperiment with 10 rows and 6 columns
##  [2] chipseq: RangedSummarizedExperiment with 25 rows and 4 columns
## Functionality:
##  experiments() - obtain the ExperimentList instance
##  colData() - the primary/phenotype DataFrame
##  sampleMap() - the sample coordination DataFrame
##  `$`, `[`, `[[` - extract colData columns, subset, or experiment
##  *Format() - convert into a long or wide DataFrame
##  assays() - convert ExperimentList to a SimpleList of matrices
##  exportClass() - save data to flat files
```
