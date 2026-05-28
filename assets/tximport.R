#!/usr/bin/env Rscript

# Summarise Salmon/Kallisto transcript-level quantification to gene level
# with tximport, ready to feed DESeq2 / edgeR.
#
# Usage:
#   Rscript tximport.R <samplesheet.csv> <gtf> <salmon|kallisto> <quant_parent_dir>
#
# <quant_parent_dir> contains one sub-directory per sample (named after the
# sample id), each holding the quantification output.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: tximport.R <samplesheet.csv> <gtf> <salmon|kallisto> <quant_parent_dir>")
}

samplesheet <- args[1]
gtf_file    <- args[2]
aligner     <- tolower(args[3])
quant_dir   <- args[4]

library(tximport)

# ---- 1. Samples -----------------------------------------------------------
samples <- read.csv(samplesheet, stringsAsFactors = FALSE)
if (!all(c("sample", "condition") %in% colnames(samples))) {
  stop("Samplesheet must contain 'sample' and 'condition' columns")
}

# ---- 2. Locate the per-sample quantification files ------------------------
quant_name <- if (aligner == "salmon") "quant.sf" else "abundance.tsv"
files <- file.path(quant_dir, samples$sample, quant_name)
names(files) <- samples$sample

missing <- !file.exists(files)
if (any(missing)) {
  stop(paste0("Quantification file (", quant_name, ") not found for sample(s): ",
              paste(samples$sample[missing], collapse = ", ")))
}

# ---- 3. Build the transcript-to-gene map from the GTF ---------------------
message("Building transcript-to-gene map from: ", gtf_file)
con   <- if (grepl("\\.gz$", gtf_file)) gzfile(gtf_file) else file(gtf_file)
lines <- readLines(con)
close(con)

# Keep only annotation lines that carry both a transcript_id and a gene_id.
lines <- lines[!grepl("^#", lines)]
lines <- lines[grepl('transcript_id "', lines) & grepl('gene_id "', lines)]
if (length(lines) == 0) {
  stop("No transcript_id/gene_id attributes found in the GTF - cannot build tx2gene")
}

txid  <- sub('.*transcript_id "([^"]+)".*', "\\1", lines)
gene  <- sub('.*gene_id "([^"]+)".*',       "\\1", lines)
tx2gene <- unique(data.frame(TXNAME = txid, GENEID = gene, stringsAsFactors = FALSE))
message("tx2gene: ", nrow(tx2gene), " transcript-to-gene pairs")

# ---- 4. Import and summarise to gene level --------------------------------
txi <- tximport(
  files,
  type            = if (aligner == "salmon") "salmon" else "kallisto",
  tx2gene         = tx2gene,
  ignoreTxVersion = TRUE   # match ENST00000456328.2 <-> ENST00000456328
)

# ---- 5. Save outputs ------------------------------------------------------
saveRDS(txi, "txi.rds")
write.csv(round(txi$counts), "tximport_gene_counts.csv")

message("tximport complete: ", nrow(txi$counts), " genes x ",
        ncol(txi$counts), " samples")
