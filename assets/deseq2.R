#!/usr/bin/env Rscript

# Differential expression with DESeq2.
#
# Usage:
#   Rscript deseq2.R <samplesheet.csv> <design> <gene_info.tsv> <featureCounts files...>
#   Rscript deseq2.R <samplesheet.csv> <design> <gene_info.tsv> <txi.rds>
#   Rscript deseq2.R <samplesheet.csv> <design> <gene_info.tsv> --matrix <count_matrix>
#
# The --matrix form reads a single pre-computed gene x sample count matrix
# (featureCounts wide output, or a plain gene-id + samples table) and enters the
# analysis directly at differential expression - see read_count_matrix() below.
#
# <design> is the model formula (e.g. "~ condition" or "~ batch + condition");
# an empty value falls back to "~ condition". <gene_info.tsv> is the
# gene_id/gene_name/gene_biotype table (from gtf2geneinfo.py) used to annotate
# the result tables with gene symbols. The single-.rds form is a tximport
# object built from Salmon/Kallisto transcript quantification; it is imported
# with DESeqDataSetFromTximport, which uses the average transcript length as a
# normalisation offset. Each contrast's log2 fold changes are shrunk with the
# apeglm estimator (lfcShrink) for better gene ranking and cleaner MA / volcano
# plots; the Wald-test stat / p-value / FDR are kept from the unshrunken fit.

args <- commandArgs(trailingOnly = TRUE)

# The optional --matrix flag selects single-count-matrix mode. It is stripped
# out here so the remaining positional arguments are unchanged from the other
# two input modes (featureCounts files / txi.rds).
matrix_mode <- "--matrix" %in% args
args        <- args[args != "--matrix"]

if (length(args) < 4) {
  stop("Usage: deseq2.R <samplesheet.csv> <design> <gene_info.tsv> [--matrix] <featureCounts files... | txi.rds | count_matrix>")
}

samplesheet_path <- args[1]
design_str       <- args[2]
gene_info_path   <- args[3]
inputs           <- args[4:length(args)]
rds_input        <- inputs[grepl("\\.rds$", inputs)]

# Load libraries. DESeq2 and ggplot2 are always present in the DESeq2
# container; pheatmap is optional and the heatmap is skipped if it is missing.
library(DESeq2)
library(ggplot2)
have_pheatmap <- requireNamespace("pheatmap", quietly = TRUE)

# Gene-symbol / biotype lookup (from the GTF-derived gene_info table). Used to
# prepend readable gene_name + gene_biotype columns to the result tables.
geneinfo <- NULL
if (!is.na(gene_info_path) && nzchar(gene_info_path) &&
    file.exists(gene_info_path) && file.info(gene_info_path)$size > 0) {
    geneinfo <- tryCatch(read.delim(gene_info_path, stringsAsFactors = FALSE),
                         error = function(e) NULL)
}
annotate_genes <- function(df) {
    # Prepend gene_id + gene_name + gene_biotype, matched on the row names (gene
    # IDs); falls back to a version-insensitive match for any unmatched IDs, then
    # to a gene_name match (for a symbol-keyed --matrix input, whose row names are
    # already symbols rather than gene IDs -- in that mode gene_id carries those
    # symbols, since they are the key the matrix is actually indexed by).
    #
    # gene_id is emitted as a named column and the callers write with
    # row.names = FALSE. Previously the IDs were left as row names, which
    # write.csv() emits as a leading column with an empty header: the
    # authoritative identifier was the one field without a label, read.csv()
    # renamed it to 'X', pandas to 'Unnamed: 0', and a reader skimming the sheet
    # saw gene_name as the first real column. Gene symbols are neither unique nor
    # stable across annotation releases (two Ensembl IDs can share a symbol, and
    # symbols get renamed), so the stable ID has to be plainly labelled.
    ids <- rownames(df)
    if (is.null(geneinfo)) return(cbind(gene_id = ids, df))
    strip <- function(x) sub("\\.[0-9]+$", "", x)
    m  <- match(ids, geneinfo$gene_id)
    na <- is.na(m)
    if (any(na)) m[na] <- match(strip(ids[na]), strip(geneinfo$gene_id))
    na <- is.na(m)
    if (any(na)) m[na] <- match(ids[na], geneinfo$gene_name)
    cbind(gene_id      = ids,
          gene_name    = geneinfo$gene_name[m],
          gene_biotype = geneinfo$gene_biotype[m],
          df)
}

# Volcano plot: log2 fold change vs -log10 p-value, coloured by significance
# (padj < 0.05 and |log2FC| > 1; up in red, down in blue, the rest grey).
draw_volcano <- function(lfc, pval, padj, title) {
    ok   <- !is.na(lfc) & !is.na(pval)
    lfc  <- lfc[ok]; pval <- pval[ok]; padj <- padj[ok]
    sig  <- !is.na(padj) & padj < 0.05 & abs(lfc) > 1
    cols <- ifelse(sig, ifelse(lfc > 0, "#c2255c", "#1c7ed6"), "#ced4da")
    plot(lfc, -log10(pval), pch = 20, cex = 0.55, col = cols,
         xlab = "log2 fold change", ylab = expression(-log[10] ~ italic(p)),
         main = title)
    abline(v = c(-1, 1), lty = 2, col = "#868e96")
    if (any(sig)) abline(h = -log10(max(pval[sig])), lty = 2, col = "#868e96")
}

# ---------------------------------------------------------------------------
# Write one figure as both PNG and SVG.
#
# `draw` is a function of no arguments that renders the plot. It is called once
# per device, because a plot can only be written while its device is open: a
# ggplot object must be print()ed inside it, base and grid (pheatmap) plots are
# simply drawn.
#
# SVG is written with grDevices::svg() rather than ggsave(..., ".svg"):
# ggsave() delegates SVG to the svglite package, which none of the pipeline
# containers carry (verified on the DESeq2, edgeR and gProfiler images), while
# svg() needs only cairo, which all of them have. Both devices are given the
# same inch dimensions so the raster and vector copies are identical in layout.
#
# The SVG is what the Quarto report embeds: it is vector, so it stays sharp at
# any zoom. Note that cairo renders text as glyph outlines, not <text> elements
# -- so the labels are not selectable, searchable, or editable in Illustrator or
# Inkscape. That makes the file font-independent (it renders identically
# everywhere, with no substitution risk) at the cost of post-editing. Producing
# editable text would need the svglite package, which none of the pipeline
# containers carry; a user who wants it can swap the device in the reproduce
# script. The PNG is kept for pasting into email, slides and documents, and for
# anything that cannot consume SVG.
#
# NOTE: kept identical in deseq2.R and edger.R (as annotate_genes and
# draw_volcano already are). Edit both together.
# ---------------------------------------------------------------------------
FIG_RES <- 150   # px per inch for the raster copy

save_figure <- function(dir, name, draw, width = 7, height = 7) {
    png(file.path(dir, paste0(name, ".png")),
        width = width, height = height, units = "in", res = FIG_RES)
    draw()
    dev.off()
    svg(file.path(dir, paste0(name, ".svg")), width = width, height = height)
    draw()
    dev.off()
}

# ---------------------------------------------------------------------------
# Count-matrix reader (--matrix mode).
#
# Reads a single pre-computed gene x sample count matrix and returns an integer
# matrix whose columns are ordered to match `sample_ids` (the samplesheet
# order). Two layouts are auto-detected:
#   * featureCounts wide output - carries the annotation columns Geneid, Chr,
#     Start, End, Strand, Length (these are dropped); count columns are named
#     after BAM paths.
#   * plain matrix - first column gene IDs, every other column a sample.
# Count columns are matched to the samplesheet by exact (cleaned) name or exact
# path token - never a raw substring, so a numeric id '100' cannot match '1100'.
# Non-integer values are rejected: they indicate pseudo-aligner estimates, which
# must be imported via tximport (a Salmon/Kallisto run), not as a flat matrix.
#
# NOTE: kept identical in deseq2.R and edger.R (as annotate_genes / draw_volcano
# already are). Edit both together.
# ---------------------------------------------------------------------------
FEATURECOUNTS_ANNO_COLS <- c("Chr", "Start", "End", "Strand", "Length")

path_basename <- function(x) sub("^.*[/\\\\]", "", x)         # portable basename (/ or \)
path_tokens   <- function(x) unlist(strsplit(x, "[/\\\\.]"))  # split on / \ .

clean_sample_name <- function(h, clean_exts) {
    b <- path_basename(h)
    b <- sub("\\.bam$", "", b, ignore.case = TRUE)
    for (ext in clean_exts) b <- sub(ext, "", b, fixed = TRUE)
    sub("[._-]+$", "", b)   # trim separators left behind by suffix removal
}

read_count_matrix <- function(path, sample_ids,
                              clean_exts = c(".fastp", ".featureCounts",
                                             "Aligned.sortedByCoord.out",
                                             "Aligned.out")) {
    # Sniff the delimiter from the first non-comment line (tab or comma).
    # readLines / read.table transparently decompress a gzipped (.gz) path, so
    # GEO-style gzipped matrices are read without a separate decompression step.
    hdr_lines <- readLines(path, n = 50)
    hdr_lines <- hdr_lines[!grepl("^#", hdr_lines)]
    if (!length(hdr_lines)) stop("count matrix '", path, "' has no data lines")
    first <- hdr_lines[1]
    n_tab <- lengths(regmatches(first, gregexpr("\t", first)))
    n_com <- lengths(regmatches(first, gregexpr(",",  first)))
    sep   <- if (n_tab >= n_com) "\t" else ","

    # A header row with one fewer field than the data rows makes read.table
    # promote the first data column to row names -- common in GEO matrices whose
    # gene-ID column has no header. Detect it by comparing header and first-data
    # field counts, so the gene IDs are taken from the row names rather than
    # mistaken for a sample column (which would otherwise fail later with a
    # confusing "sample not found" error).
    n_hdr  <- length(strsplit(first, sep, fixed = TRUE)[[1]])
    n_data <- if (length(hdr_lines) >= 2)
                  length(strsplit(hdr_lines[2], sep, fixed = TRUE)[[1]]) else n_hdr
    promoted <- n_data == n_hdr + 1

    df <- read.table(path, header = TRUE, sep = sep, comment.char = "#",
                     check.names = FALSE, stringsAsFactors = FALSE, quote = "")
    if (ncol(df) < 1 || (!promoted && ncol(df) < 2))
        stop("count matrix '", path, "' has too few columns for gene IDs plus samples")

    # Gene IDs come from the row names when the gene column was header-less
    # (promoted), otherwise from the first column. Then drop the featureCounts
    # annotation columns, but only when all five are present (a plain matrix has
    # none of them).
    if (promoted) {
        gene_ids   <- rownames(df)
        count_cols <- df
    } else {
        gene_ids   <- as.character(df[[1]])
        count_cols <- df[, -1, drop = FALSE]
    }
    anno <- intersect(FEATURECOUNTS_ANNO_COLS, colnames(count_cols))
    if (length(anno) == length(FEATURECOUNTS_ANNO_COLS)) {
        count_cols <- count_cols[, !(colnames(count_cols) %in% anno), drop = FALSE]
    }
    headers <- colnames(count_cols)

    # Resolve each samplesheet sample to exactly one count column.
    cleaned  <- vapply(headers, clean_sample_name, character(1), clean_exts)
    col_idx  <- integer(length(sample_ids))
    problems <- character(0)
    for (i in seq_along(sample_ids)) {
        s   <- sample_ids[i]
        hit <- which(cleaned == s)                                   # exact cleaned name
        if (!length(hit))                                            # else exact path token
            hit <- which(vapply(headers,
                                function(h) s %in% path_tokens(h), logical(1)))
        if (length(hit) == 0) {
            problems <- c(problems,
                          sprintf("  - sample '%s': no matching count column", s))
        } else if (length(hit) > 1) {
            problems <- c(problems,
                          sprintf("  - sample '%s': ambiguous, matches %d columns (%s)",
                                  s, length(hit), paste(headers[hit], collapse = ", ")))
        } else {
            col_idx[i] <- hit
        }
    }
    if (length(problems))
        stop("could not map samplesheet samples to count columns:\n",
             paste(problems, collapse = "\n"))

    unused <- setdiff(seq_along(headers), col_idx)
    if (length(unused))
        message("Note: ", length(unused),
                " count column(s) not listed in the samplesheet were ignored: ",
                paste(headers[unused], collapse = ", "))

    # Order to samplesheet, coerce to integer, reject non-integer (pseudo-aligner) input.
    m <- as.matrix(count_cols[, col_idx, drop = FALSE])
    storage.mode(m) <- "double"
    if (any(is.na(m)))
        stop("count matrix contains missing or non-numeric values")
    if (any(abs(m - round(m)) > 1e-6))
        stop("count matrix has non-integer values. DESeq2 and edgeR need raw ",
             "integer counts, so normalised values (TPM, FPKM / RPKM, or ",
             "DESeq2- / CPM-normalised counts) and Salmon/Kallisto estimated ",
             "counts are not accepted here. For pseudo-aligner output, re-run via ",
             "--aligner salmon|kallisto (tximport imports it correctly); for a ",
             "normalised matrix, obtain the raw integer counts instead.")
    m <- round(m)
    storage.mode(m) <- "integer"
    rownames(m) <- gene_ids
    colnames(m) <- sample_ids
    m
}

# Report how well the matrix gene IDs match the GTF-derived annotation, and stop
# on an obvious mismatch (wrong genome build, or an unexpected key). Only
# meaningful for --matrix input, where the counts and the GTF are supplied
# independently; on an aligned run they always match by construction.
check_annotation_match <- function(ids, geneinfo, floor = 0.5) {
    if (is.null(geneinfo)) return(invisible(NULL))
    strip    <- function(x) sub("\\.[0-9]+$", "", x)
    id_rate  <- mean(strip(ids) %in% strip(geneinfo$gene_id))
    sym_rate <- mean(ids %in% geneinfo$gene_name)
    message(sprintf("Annotation match: %.1f%% of %d matrix genes match gene_id, %.1f%% match gene_name",
                    100 * id_rate, length(ids), 100 * sym_rate))
    if (max(id_rate, sym_rate) < floor) {
        stop(sprintf(paste0(
            "only %.1f%% of matrix gene IDs match the --gtf annotation (and %.1f%% match ",
            "gene symbols). The count matrix and the GTF look mismatched - check they are ",
            "the same genome / annotation build, and that the matrix's first column is the ",
            "gene identifier."), 100 * id_rate, 100 * sym_rate))
    }
    if (id_rate < floor && sym_rate >= floor)
        message("Note: matrix appears to be gene-symbol-keyed; results are annotated via gene_name.")
    invisible(NULL)
}

# 1. Read Samplesheet
samples <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
# Force sample IDs to character so numeric-looking IDs (e.g. '100') are matched
# as strings, never coerced to integers.
samples$sample <- as.character(samples$sample)
rownames(samples) <- samples$sample

# Treat a condition named 'REF' as the baseline (denominator) of contrasts.
if ("REF" %in% samples$condition) {
    samples$condition <- relevel(factor(samples$condition), ref = "REF")
} else {
    warning("Condition 'REF' not found. Using default (alphabetical) level ordering.")
    samples$condition <- factor(samples$condition)
}

# Resolve the model formula. An empty <design> falls back to '~ condition';
# covariate columns (everything except 'condition') are modelled as factors.
if (is.na(design_str) || !nzchar(trimws(design_str))) design_str <- "~ condition"
if (!grepl("^\\s*~", design_str)) design_str <- paste("~", design_str)
design_formula <- as.formula(design_str)
for (v in all.vars(design_formula)) {
    if (!v %in% colnames(samples)) {
        stop(sprintf("design variable '%s' is not a column of the samplesheet", v))
    }
    if (v != "condition") samples[[v]] <- factor(samples[[v]])
}
message("DESeq2 design: ", design_str)

# 2. Build the DESeqDataSet - from a count matrix, a tximport object, or
#    per-sample featureCounts files.
if (matrix_mode) {
    message("Input mode: pre-computed count matrix")
    count_matrix <- read_count_matrix(inputs[1], samples$sample)
    check_annotation_match(rownames(count_matrix), geneinfo)

    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData   = samples,
                                  design    = design_formula)
} else if (length(rds_input) >= 1) {
    message("Input mode: tximport (transcript-level quantification)")
    txi <- readRDS(rds_input[1])

    ord <- samples$sample
    if (!all(ord %in% colnames(txi$counts))) {
        stop("tximport object is missing samples listed in the samplesheet")
    }
    # Reorder tximport matrices to match the samplesheet/colData order.
    txi$counts    <- txi$counts[,    ord, drop = FALSE]
    txi$abundance <- txi$abundance[, ord, drop = FALSE]
    txi$length    <- txi$length[,    ord, drop = FALSE]

    dds <- DESeqDataSetFromTximport(txi, colData = samples, design = design_formula)
} else {
    message("Input mode: featureCounts gene counts")
    count_files <- inputs
    counts_list <- list()
    for (s in samples$sample) {
        # featureCounts files are named "<sample>.featureCounts.txt".
        match_file <- grep(paste0(s, ".featureCounts.txt"), count_files, value = TRUE)
        if (length(match_file) == 0) {
            stop(paste("Count file for sample", s, "not found in inputs"))
        }
        # featureCounts format: Geneid, Chr, Start, End, Strand, Length, <count>
        df <- read.table(match_file[1], header = TRUE, comment.char = "#",
                         stringsAsFactors = FALSE)
        counts <- df[, ncol(df)]          # count is the last column
        names(counts) <- df$Geneid
        counts_list[[s]] <- counts
    }
    count_matrix <- do.call(cbind, counts_list)
    rownames(count_matrix) <- names(counts_list[[1]])

    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData   = samples,
                                  design    = design_formula)
}

# 3. Filter low-count genes (applied identically to both input modes)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]

# 4. Run DESeq2
dds <- DESeq(dds)

# 5. Results & Visualisations
dir.create("deseq2_output", showWarnings = FALSE)

# Variance Stabilizing Transformation for PCA/Heatmap
vsd <- vst(dds, blind = FALSE)

# PCA
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = 1.5, hjust = 1.5) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

save_figure("deseq2_output", "pca_plot", function() print(p))

# The plotted coordinates, so the PCA can be redrawn (or restyled) without
# rerunning DESeq2. percentVar is carried as columns because it belongs to the
# axis labels, not to any single sample.
pca_out <- pcaData
pca_out$percentVar_PC1 <- percentVar[1]
pca_out$percentVar_PC2 <- percentVar[2]
write.csv(pca_out, file.path("deseq2_output", "pca_data.csv"), row.names = FALSE)

# Heatmap of the 20 most variable genes (variance computed with base R so no
# extra package dependency is needed).
gene_var    <- apply(assay(vsd), 1, var)
topVarGenes <- head(order(gene_var, decreasing = TRUE), 20)
mat <- assay(vsd)[topVarGenes, ]
df <- as.data.frame(colData(dds)[, c("condition")])
rownames(df) <- colnames(mat)
colnames(df) <- "condition"

# The variance-stabilised values behind the heatmap, annotated like the results
# tables. check.names = FALSE keeps sample names verbatim. Written whether or
# not pheatmap is available, so the data exists even if the figure does not.
write.csv(annotate_genes(data.frame(mat, check.names = FALSE)),
          file.path("deseq2_output", "heatmap_top_var.csv"), row.names = FALSE)
# The column annotation (sample -> condition), so the heatmap is reproducible
# from these two files alone.
write.csv(data.frame(sample = rownames(df), condition = df$condition),
          file.path("deseq2_output", "heatmap_top_var_annotation.csv"),
          row.names = FALSE)

if (have_pheatmap) {
    save_figure("deseq2_output", "heatmap_top_var",
                function() pheatmap::pheatmap(mat, annotation_col = df,
                                              show_rownames = TRUE),
                width = 8, height = 7)
} else {
    message("pheatmap not available - skipping heatmap_top_var.png/.svg")
}

# 6. Contrasts
#
# Each contrast's log2 fold changes are shrunk with the apeglm estimator
# (lfcShrink), which pulls low-count / high-variance estimates toward zero for
# better gene ranking and cleaner MA / volcano plots. apeglm shrinks a model
# *coefficient*, so 'condA vs condB' must be a single coefficient: when condB
# is not already the model's reference level, the condition factor is releveled
# and the GLM refitted (nbinomWaldTest; releveling does not change the
# dispersions, so they are reused). The Wald-test stat / p-value / FDR come
# from the unshrunken fit -- shrinkage changes the effect-size estimate, not
# the test -- so the 'stat' column is preserved for downstream GSEA ranking.
cond_levels <- levels(samples$condition)

run_contrast <- function(condA, condB) {
    coef_name <- paste0("condition_", condA, "_vs_", condB)

    if (condB == levels(dds$condition)[1]) {
        # condB is already the model's reference level.
        dds_use <- dds
    } else {
        # Relevel so condB is the reference, then refit the GLM so that
        # 'condA vs condB' becomes a single coefficient apeglm can shrink.
        dds_use <- dds
        dds_use$condition <- relevel(dds_use$condition, ref = condB)
        dds_use <- nbinomWaldTest(dds_use)
    }

    # Unshrunken Wald results (supply stat / pvalue / padj) and the
    # apeglm-shrunken log2 fold changes (used for ranking and the plots).
    res_mle <- results(dds_use, name = coef_name)
    res     <- lfcShrink(dds_use, coef = coef_name, type = "apeglm",
                         res = res_mle)

    # CSV: apeglm-shrunken log2FoldChange + lfcSE, with the Wald 'stat' carried
    # over from the unshrunken fit so GSEA can still rank genes by it.
    res_df <- as.data.frame(res)
    res_df$stat <- res_mle$stat
    res_df <- res_df[, c("baseMean", "log2FoldChange", "lfcSE",
                         "stat", "pvalue", "padj")]
    res_df <- res_df[order(res_df$pvalue), ]

    filename <- paste0("deseq2_results_", condA, "_vs_", condB, ".csv")
    # row.names = FALSE: annotate_genes() has already promoted the row names to
    # an explicit gene_id column, so writing them again would duplicate the IDs
    # under a blank header.
    write.csv(annotate_genes(res_df),
              file = file.path("deseq2_output", filename), row.names = FALSE)

    # No separate CSV for the MA and volcano plots: both are drawn entirely from
    # baseMean / log2FoldChange / pvalue / padj, every one of which is already a
    # column of the deseq2_results_*.csv written just above. Duplicating a
    # whole-transcriptome table twice more would add nothing.
    save_figure("deseq2_output", paste0("maplot_", condA, "_vs_", condB),
                function() plotMA(res, main = paste(condA, "vs", condB,
                                                    "(apeglm-shrunk LFC)"),
                                  ylim = c(-2, 2)))

    save_figure("deseq2_output", paste0("volcano_", condA, "_vs_", condB),
                function() draw_volcano(res$log2FoldChange, res$pvalue, res$padj,
                                        paste(condA, "vs", condB)),
                width = 8, height = 6.5)
}

pairs <- combn(cond_levels, 2)

for (i in 1:ncol(pairs)) {
    c1 <- pairs[1, i]
    c2 <- pairs[2, i]

    # Orient each contrast so 'REF' is the denominator where present.
    if (c1 == "REF") {
        run_contrast(c2, c1)
    } else if (c2 == "REF") {
        run_contrast(c1, c2)
    } else {
        run_contrast(c1, c2)
    }
}
