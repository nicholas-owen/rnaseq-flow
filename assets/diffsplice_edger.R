#!/usr/bin/env Rscript

# Differential splicing with edgeR's diffSpliceDGE.
#
# diffSpliceDGE tests, for each feature (an exon or a transcript), whether its
# log-fold-change between conditions differs from the overall log-fold-change
# of the gene it belongs to. A feature that shifts more (or less) than its gene
# indicates differential usage:
#   * exon level       -> differential exon usage       (STAR / HISAT2)
#   * transcript level -> differential transcript usage  (Salmon / Kallisto)
# It complements gene-level DE (DESeq2/edgeR) and rMATS' event-based test.
#
# Usage:
#   diffsplice_edger.R <samplesheet.csv> <level> <gtf> <aligner> <data_dir> [cpus]
#     level    : 'exon' or 'transcript'
#     gtf      : annotation GTF (gene symbols/biotype; transcript mode also
#                builds the transcript-to-gene map)
#     aligner  : star | hisat2 | salmon | kallisto
#     data_dir : directory holding the staged inputs --
#                  exon mode       -> <sample>.exon.featureCounts.txt files
#                  transcript mode -> per-sample Salmon/Kallisto quant directories

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: diffsplice_edger.R <samplesheet.csv> <level> <gtf> <aligner> <data_dir> [cpus]")
}
samplesheet <- args[1]
level       <- tolower(args[2])
gtf_file    <- args[3]
aligner     <- tolower(args[4])
data_dir    <- args[5]
if (!level %in% c("exon", "transcript")) {
  stop("level must be 'exon' or 'transcript'")
}

suppressPackageStartupMessages(library(edgeR))

outdir <- "diffsplice_output"
dir.create(outdir, showWarnings = FALSE)
feature_csv <- file.path(outdir, paste0("diffsplice_", level, "_results.csv"))
gene_csv    <- file.path(outdir, "diffsplice_gene_results.csv")

emit_empty <- function(msg) {
  message(msg, " Writing empty diffSplice result.")
  write.csv(data.frame(), feature_csv, row.names = FALSE)
  write.csv(data.frame(), gene_csv,    row.names = FALSE)
  quit(save = "no", status = 0)
}

strip <- function(x) sub("\\.[0-9]+$", "", x)   # drop version suffixes

# ---- samples --------------------------------------------------------------
samples   <- read.csv(samplesheet, stringsAsFactors = FALSE)
condition <- factor(samples$condition)
if ("REF" %in% levels(condition)) condition <- relevel(condition, ref = "REF")
if (nlevels(condition) < 2) emit_empty("Fewer than two conditions.")

# ---- GTF: gene-symbol / biotype map (transcript mode also reuses `lines`) --
con   <- if (grepl("\\.gz$", gtf_file)) gzfile(gtf_file) else file(gtf_file)
lines <- readLines(con); close(con)
lines <- lines[!grepl("^#", lines) & grepl('gene_id "', lines)]
if (length(lines) == 0) stop("No gene_id attributes found in the GTF.")
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
annotate_genes <- function(df) {
  # Prepend gene_name + gene_biotype, matched on the GeneID column.
  if (nrow(df) == 0) return(df)
  m <- match(df$GeneID, geneinfo$gene_id)
  data.frame(GeneID       = df$GeneID,
             gene_name    = geneinfo$gene_name[m],
             gene_biotype = geneinfo$gene_biotype[m],
             df[, setdiff(colnames(df), "GeneID"), drop = FALSE],
             check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- assemble the per-feature count matrix --------------------------------
if (level == "exon") {
  # featureCounts -f output: one row per GTF exon, columns
  # Geneid, Chr, Start, End, Strand, Length, <count>.
  files <- file.path(data_dir, paste0(samples$sample, ".exon.featureCounts.txt"))
  names(files) <- samples$sample
  if (!all(file.exists(files))) {
    stop("Exon featureCounts files missing for: ",
         paste(samples$sample[!file.exists(files)], collapse = ", "))
  }
  geneid <- NULL; exonid <- NULL; mats <- list()
  for (s in samples$sample) {
    df  <- read.table(files[s], header = TRUE, comment.char = "#",
                       sep = "\t", stringsAsFactors = FALSE)
    key <- paste(df$Geneid, df$Chr, df$Start, df$End, sep = ":")
    # An exon shared by several transcripts is listed once per transcript;
    # collapse to one row per genomic exon region (counts are identical).
    uniq <- !duplicated(key)
    if (is.null(geneid)) { geneid <- df$Geneid[uniq]; exonid <- key[uniq] }
    v <- df[[ncol(df)]][uniq]; names(v) <- key[uniq]
    mats[[s]] <- v[exonid]
  }
  counts <- do.call(cbind, mats)
} else {
  suppressPackageStartupMessages(library(tximport))
  # reuse the gene_id-bearing GTF lines read above; keep the transcript lines
  txlines <- lines[grepl('transcript_id "', lines)]
  if (length(txlines) == 0) stop("No transcript_id/gene_id pairs found in the GTF.")
  tx2gene <- unique(data.frame(
    TXNAME = strip(sub('.*transcript_id "([^"]+)".*', "\\1", txlines)),
    GENEID = sub('.*gene_id "([^"]+)".*', "\\1", txlines),
    stringsAsFactors = FALSE))

  quant_name <- if (aligner == "salmon") "quant.sf" else "abundance.tsv"
  files <- file.path(data_dir, samples$sample, quant_name)
  names(files) <- samples$sample
  if (!all(file.exists(files))) {
    stop("Quantification files missing for: ",
         paste(samples$sample[!file.exists(files)], collapse = ", "))
  }
  txi <- tximport(files, type = aligner, txOut = TRUE,
                  countsFromAbundance = "scaledTPM", ignoreTxVersion = TRUE)
  counts <- txi$counts
  rownames(counts) <- strip(rownames(counts))
  geneid <- tx2gene$GENEID[match(rownames(counts), tx2gene$TXNAME)]
  ok     <- !is.na(geneid)
  counts <- counts[ok, , drop = FALSE]
  geneid <- geneid[ok]
  exonid <- rownames(counts)
}
colnames(counts) <- samples$sample
storage.mode(counts) <- "double"

# ---- restrict to multi-feature genes (splicing is undefined otherwise) -----
multi  <- function(g) g %in% names(which(table(g) > 1L))
keep   <- multi(geneid)
counts <- counts[keep, , drop = FALSE]; geneid <- geneid[keep]; exonid <- exonid[keep]
if (nrow(counts) < 2) emit_empty("No multi-feature genes to test.")
message(sprintf("diffSplice input: %d %s features across %d genes.",
                nrow(counts), level, length(unique(geneid))))

# ---- edgeR fit ------------------------------------------------------------
y <- DGEList(counts = round(counts),
             genes  = data.frame(GeneID = geneid, FeatureID = exonid,
                                 stringsAsFactors = FALSE))
keep <- filterByExpr(y, group = condition)
y    <- y[keep, , keep.lib.sizes = FALSE]
if (sum(table(y$genes$GeneID) > 1L) < 1L) {
  emit_empty("No multi-feature genes passed expression filtering.")
}
y <- calcNormFactors(y)

design <- model.matrix(~ condition)
colnames(design) <- make.names(colnames(design))
y   <- estimateDisp(y, design)
fit <- glmFit(y, design)

# ---- diffSpliceDGE, once per condition coefficient (vs the reference) ------
feature_tables <- list(); gene_tables <- list()
for (cf in colnames(design)[-1]) {
  cmp <- sub("^condition", "", cf)
  ds  <- diffSpliceDGE(fit, geneid = "GeneID", exonid = "FeatureID", coef = cf)
  ft  <- topSpliceDGE(ds, test = "exon",  number = Inf)
  gt  <- topSpliceDGE(ds, test = "Simes", number = Inf)
  if (nrow(ft) > 0) ft$comparison <- cmp
  if (nrow(gt) > 0) gt$comparison <- cmp
  feature_tables[[cf]] <- ft
  gene_tables[[cf]]    <- gt
  message(sprintf("  %s: %d significant %s features, %d significant genes (FDR < 0.05).",
                  cmp, sum(ft$FDR < 0.05, na.rm = TRUE), level,
                  sum(gt$FDR < 0.05, na.rm = TRUE)))
}

write.csv(annotate_genes(do.call(rbind, feature_tables)), feature_csv, row.names = FALSE)
write.csv(annotate_genes(do.call(rbind, gene_tables)),    gene_csv,    row.names = FALSE)
message("diffSplice complete.")
