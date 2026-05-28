#!/usr/bin/env Rscript

# Differential transcript usage (DTU) testing with DEXSeq.
#
# Each transcript is treated as an "exon" feature of its gene, so DEXSeq's
# differential-exon-usage test becomes a test of whether transcript usage
# (the proportions of isoforms within a gene) changes between conditions.
# This complements gene-level differential expression: a gene can show a
# significant isoform-usage switch even when its total expression is unchanged.
#
# Usage: dtu_dexseq.R <samplesheet.csv> <gtf> <salmon|kallisto> <quant_dir> [cpus]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: dtu_dexseq.R <samplesheet.csv> <gtf> <salmon|kallisto> <quant_dir> [cpus]")
}
samplesheet <- args[1]
gtf_file    <- args[2]
aligner     <- tolower(args[3])
quant_dir   <- args[4]
cpus        <- if (length(args) >= 5) suppressWarnings(as.integer(args[5])) else 1L
if (is.na(cpus) || cpus < 1L) cpus <- 1L

suppressPackageStartupMessages({
  library(tximport)
  library(DEXSeq)
  library(BiocParallel)
})

dir.create("dtu_output", showWarnings = FALSE)
bpp <- MulticoreParam(workers = cpus)

emit_empty <- function(msg) {
  message(msg, " Writing empty DTU result.")
  write.csv(data.frame(), file.path("dtu_output", "dtu_transcript_results.csv"))
  write.csv(data.frame(), file.path("dtu_output", "dtu_gene_qvalues.csv"))
  quit(save = "no", status = 0)
}

# ---- samples & quantification files ---------------------------------------
samples <- read.csv(samplesheet, stringsAsFactors = FALSE)
quant_name <- if (aligner == "salmon") "quant.sf" else "abundance.tsv"
files <- file.path(quant_dir, samples$sample, quant_name)
names(files) <- samples$sample
if (!all(file.exists(files))) {
  stop("Quantification files not found for: ",
       paste(samples$sample[!file.exists(files)], collapse = ", "))
}

# ---- transcript-to-gene map from the GTF ----------------------------------
strip <- function(x) sub("\\.[0-9]+$", "", x)   # drop version suffixes
con   <- if (grepl("\\.gz$", gtf_file)) gzfile(gtf_file) else file(gtf_file)
lines <- readLines(con)
close(con)
lines <- lines[!grepl("^#", lines)]
lines <- lines[grepl('transcript_id "', lines) & grepl('gene_id "', lines)]
if (length(lines) == 0) stop("No transcript_id/gene_id found in the GTF.")
tx2gene <- unique(data.frame(
  TXNAME = strip(sub('.*transcript_id "([^"]+)".*', "\\1", lines)),
  GENEID = sub('.*gene_id "([^"]+)".*',       "\\1", lines),
  stringsAsFactors = FALSE))

# gene-symbol / biotype map, used to annotate the result tables
g_id <- sub('.*gene_id "([^"]+)".*', "\\1", lines)
g_nm <- ifelse(grepl('gene_name "', lines),
               sub('.*gene_name "([^"]+)".*', "\\1", lines), g_id)
g_bt <- ifelse(grepl('gene_biotype "', lines),
               sub('.*gene_biotype "([^"]+)".*', "\\1", lines),
        ifelse(grepl('gene_type "', lines),
               sub('.*gene_type "([^"]+)".*', "\\1", lines), "NA"))
geneinfo <- data.frame(gene_id = g_id, gene_name = g_nm, gene_biotype = g_bt,
                       stringsAsFactors = FALSE)
geneinfo <- geneinfo[!duplicated(geneinfo$gene_id), ]

# ---- import transcript-level counts (DTU-scaled) --------------------------
txi <- tximport(files, type = aligner, txOut = TRUE,
                countsFromAbundance = "dtuScaledTPM", ignoreTxVersion = TRUE)
cts <- txi$counts
rownames(cts) <- strip(rownames(cts))
cts <- cts[rowSums(cts) > 0, , drop = FALSE]

gene_of <- tx2gene$GENEID[match(rownames(cts), tx2gene$TXNAME)]
ok <- !is.na(gene_of)
cts <- cts[ok, , drop = FALSE]
gene_of <- gene_of[ok]

# DTU is only defined for genes with >= 2 transcripts.
keep_multi <- function(g) g %in% names(which(table(g) > 1L))
m <- keep_multi(gene_of); cts <- cts[m, , drop = FALSE]; gene_of <- gene_of[m]

# Light expression filter: transcript seen (>= 5 counts) in >= 2 samples,
# then re-restrict to genes that still have >= 2 transcripts.
expr <- rowSums(cts >= 5) >= 2
cts <- cts[expr, , drop = FALSE]; gene_of <- gene_of[expr]
m <- keep_multi(gene_of); cts <- cts[m, , drop = FALSE]; gene_of <- gene_of[m]

if (nrow(cts) < 2 || length(unique(gene_of)) < 1) {
  emit_empty("No multi-transcript genes passed filtering.")
}
message(sprintf("DTU input: %d transcripts across %d multi-transcript genes.",
                nrow(cts), length(unique(gene_of))))

# ---- DEXSeq DTU test ------------------------------------------------------
condition <- factor(samples$condition)
if ("REF" %in% levels(condition)) condition <- relevel(condition, ref = "REF")
sampleData <- data.frame(condition = condition, row.names = samples$sample)

dxd <- DEXSeqDataSet(
  countData  = round(cts),
  sampleData = sampleData,
  design     = ~ sample + exon + condition:exon,
  featureID  = rownames(cts),
  groupID    = gene_of)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM = bpp)
dxd <- testForDEU(dxd, BPPARAM = bpp)
dxr <- DEXSeqResults(dxd)

# ---- write results --------------------------------------------------------
cols <- intersect(c("groupID", "featureID", "exonBaseMean",
                    "dispersion", "stat", "pvalue", "padj"), colnames(dxr))
res <- as.data.frame(dxr[, cols])
res <- res[order(res$padj), ]
gm  <- match(res$groupID, geneinfo$gene_id)
res <- data.frame(groupID      = res$groupID,
                  gene_name    = geneinfo$gene_name[gm],
                  gene_biotype = geneinfo$gene_biotype[gm],
                  res[, setdiff(colnames(res), "groupID"), drop = FALSE],
                  check.names = FALSE, stringsAsFactors = FALSE)
write.csv(res, file.path("dtu_output", "dtu_transcript_results.csv"),
          row.names = FALSE)

gq   <- perGeneQValue(dxr)
gqdf <- data.frame(gene = names(gq), gene_qvalue = as.numeric(gq),
                   stringsAsFactors = FALSE)
gqdf <- gqdf[order(gqdf$gene_qvalue), ]
gm   <- match(gqdf$gene, geneinfo$gene_id)
gqdf <- data.frame(gene         = gqdf$gene,
                   gene_name    = geneinfo$gene_name[gm],
                   gene_biotype = geneinfo$gene_biotype[gm],
                   gene_qvalue  = gqdf$gene_qvalue,
                   stringsAsFactors = FALSE)
write.csv(gqdf, file.path("dtu_output", "dtu_gene_qvalues.csv"),
          row.names = FALSE)

message(sprintf("DTU complete: %d transcripts and %d genes with adjusted p/q < 0.05.",
                sum(res$padj < 0.05, na.rm = TRUE),
                sum(gqdf$gene_qvalue < 0.05, na.rm = TRUE)))
