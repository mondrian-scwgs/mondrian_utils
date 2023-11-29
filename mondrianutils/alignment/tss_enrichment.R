#!/usr/bin/env Rscript

library(getopt)
suppressPackageStartupMessages(library("ATACseqQC"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg18.knownGene"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg38.knownGene"))
options(error=traceback)


spec <- matrix(
  c(
    "bamfile","bamfile", 1, "character", "input bam",
    "chromosomes","chromosomes", 1, "character", "comma separated chromosome list",
    "output","output", 1, "character", "output file",
    "tempdir","tempdir", 1, "character", "temppdir",
    "genome_version","genome_version", 1, "character", "genome identifier"
  ),
  ncol = 5,
  byrow = TRUE,
)
opt <- getopt(spec)


if (is.null(opt$bamfile)) {
  stop("--bamfile is required.")
}
if (is.null(opt$chromosomes)) {
  stop("--chromosomes is required.")
}
if (is.null(opt$output)) {
  stop("--output is required.")
}
if (is.null(opt$tempdir)) {
  stop("--output is required.")
}
if (is.null(opt$genome_version)) {
  stop("--genome version is required.")
}

if (!(opt$genome_version %in% c("grch37", "grch38"))) {
  stop("Invalid genome version. Please use either 'grch37' or 'grch38'.\n")
}

if (!dir.exists(opt$tempdir)) {
  dir.create(opt$tempdir, recursive = TRUE)
}

shiftedBamfile <- file.path(opt$tempdir, "shifted.bam")
if (file.exists(shiftedBamfile)) {
  file.remove(shiftedBamfile)
}


chromosomes <-  strsplit(opt$chromosomes, ",")[[1]]


if (opt$genome_version == "grch37"){
  seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg18.knownGene)
  which <- as(seqinformation, "GRanges")
  seqlevelsStyle(which) <- "NCBI"
  seqlevels(which, pruning.mode="coarse") <- chromosomes

  txs <- transcripts(TxDb.Hsapiens.UCSC.hg18.knownGene)
  seqlevelsStyle(txs) <- "NCBI"
  seqlevels(txs, pruning.mode="coarse") <- chromosomes
} else {
  seqinformation <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
  which <- as(seqinformation, "GRanges")
  seqlevels(which, pruning.mode="coarse") <- chromosomes

  txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
  seqlevels(txs, pruning.mode="coarse") <- chromosomes
}

gal <- readBamFile(opt$bamfile, which=which, asMates=TRUE, bigFile=TRUE)
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)

tsse <- TSSEscore(gal1, txs)
score <- tsse$TSSEscore


writeLines(as.character(score), opt$output)