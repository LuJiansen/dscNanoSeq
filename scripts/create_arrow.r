library(parallel)
library(ArchR)

args <- commandArgs(TRUE)
frag <- args[1]
ref <- args[2]
prefix <- args[3]

set.seed(1)

print(ref)
if (!ref %in% c("hg38", "hg19", "mm10", "mm9")) {
  if (ref == "GRCh38") {
    ref <- "hg38"
  } else if (ref == "GRCh37") {
    ref <- "hg19"
  } else if (ref == "GRCm37") {
    ref <- "mm9"
  } else if (ref == "GRCm38") {
    ref <- "mm10"
  } else {
    stop("Invalid reference genome provided.")
  }
}

addArchRGenome(ref)

ArrowFiles <- createArrowFiles(
  inputFiles = frag,
  sampleNames = prefix,
  minTSS = 0, minFrags = 0,
  maxFrags = 1e+08,
  addTileMat = TRUE,
  TileMatParams = list(blacklist = NULL),
  addGeneScoreMat = TRUE,
  excludeChr = NULL,
  GeneScoreMatParams = list(blacklist = NULL),
  force = TRUE,
  cleanTmp = FALSE # DO NOT remove tmp, which may lead to disruption of other arrow processes!
)
