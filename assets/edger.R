#!/usr/bin/env Rscript

# Differential expression with edgeR.
#
# Usage:
#   Rscript edger.R <samplesheet.csv> <design> <gene_info.tsv> <featureCounts files...>
#   Rscript edger.R <samplesheet.csv> <design> <gene_info.tsv> <txi.rds>
#   Rscript edger.R <samplesheet.csv> <design> <gene_info.tsv> --matrix <count_matrix>
#
# The --matrix form reads a single pre-computed gene x sample count matrix
# (featureCounts wide output, or a plain gene-id + samples table) and enters the
# analysis directly at differential expression - see read_count_matrix() below.
#
# <design> is the model formula (e.g. "~ condition" or "~ batch + condition");
# an empty value falls back to "~ condition". <gene_info.tsv> is the
# gene_id/gene_name/gene_biotype table (from gtf2geneinfo.py) used to annotate
# the result tables with gene symbols. Every contrast is tested with the
# quasi-likelihood GLM (glmQLFTest), so any covariates in the design are
# properly adjusted for. The single-.rds form is a tximport object built from
# Salmon/Kallisto transcript quantification, imported via the official
# tximport -> edgeR recipe (transcript length becomes a model offset).

args <- commandArgs(trailingOnly = TRUE)

# The optional --matrix flag selects single-count-matrix mode. It is stripped
# out here so the remaining positional arguments are unchanged from the other
# two input modes (featureCounts files / txi.rds).
matrix_mode <- "--matrix" %in% args
args        <- args[args != "--matrix"]

if (length(args) < 4) {
  stop("Usage: edger.R <samplesheet.csv> <design> <gene_info.tsv> [--matrix] <featureCounts files... | txi.rds | count_matrix>")
}

samplesheet_path <- args[1]
design_str       <- args[2]
gene_info_path   <- args[3]
inputs           <- args[4:length(args)]
rds_input        <- inputs[grepl("\\.rds$", inputs)]

# Load libraries (edgeR pulls in limma; all plots here use base graphics).
library(edgeR)

# Gene-symbol / biotype lookup (from the GTF-derived gene_info table). Used to
# prepend readable gene_name + gene_biotype columns to the result tables.
geneinfo <- NULL
if (!is.na(gene_info_path) && nzchar(gene_info_path) &&
    file.exists(gene_info_path) && file.info(gene_info_path)$size > 0) {
    geneinfo <- tryCatch(read.delim(gene_info_path, stringsAsFactors = FALSE),
                         error = function(e) NULL)
}
annotate_genes <- function(df) {
    # Prepend gene_name + gene_biotype, matched on the row names (gene IDs);
    # falls back to a version-insensitive match for any unmatched IDs, then to a
    # gene_name match (for a symbol-keyed --matrix input, whose row names are
    # already symbols rather than gene IDs).
    if (is.null(geneinfo)) return(df)
    strip <- function(x) sub("\\.[0-9]+$", "", x)
    m  <- match(rownames(df), geneinfo$gene_id)
    na <- is.na(m)
    if (any(na)) m[na] <- match(strip(rownames(df)[na]), strip(geneinfo$gene_id))
    na <- is.na(m)
    if (any(na)) m[na] <- match(rownames(df)[na], geneinfo$gene_name)
    cbind(gene_name    = geneinfo$gene_name[m],
          gene_biotype = geneinfo$gene_biotype[m],
          df)
}

# Volcano plot: log2 fold change vs -log10 p-value, coloured by significance
# (FDR < 0.05 and |logFC| > 1; up in red, down in blue, the rest grey).
draw_volcano <- function(lfc, pval, fdr, title) {
    ok   <- !is.na(lfc) & !is.na(pval)
    lfc  <- lfc[ok]; pval <- pval[ok]; fdr <- fdr[ok]
    sig  <- !is.na(fdr) & fdr < 0.05 & abs(lfc) > 1
    cols <- ifelse(sig, ifelse(lfc > 0, "#c2255c", "#1c7ed6"), "#ced4da")
    plot(lfc, -log10(pval), pch = 20, cex = 0.55, col = cols,
         xlab = "log2 fold change", ylab = expression(-log[10] ~ italic(p)),
         main = title)
    abline(v = c(-1, 1), lty = 2, col = "#868e96")
    if (any(sig)) abline(h = -log10(max(pval[sig])), lty = 2, col = "#868e96")
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

# 1. Read samplesheet
samples <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
# Force sample IDs to character so numeric-looking IDs (e.g. '100') are matched
# as strings, never coerced to integers.
samples$sample <- as.character(samples$sample)
rownames(samples) <- samples$sample

# 'condition' is the variable of interest; 'REF' is its baseline level.
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
message("edgeR design: ", design_str)

# 2. Build the DGEList - from a count matrix, a tximport object, or per-sample
#    featureCounts files.
if (matrix_mode) {
    message("Input mode: pre-computed count matrix")
    count_matrix <- read_count_matrix(inputs[1], samples$sample)
    check_annotation_match(rownames(count_matrix), geneinfo)

    y <- DGEList(counts = count_matrix, group = samples$condition)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
} else if (length(rds_input) >= 1) {
    message("Input mode: tximport (transcript-level quantification)")
    txi <- readRDS(rds_input[1])

    ord <- samples$sample
    if (!all(ord %in% colnames(txi$counts))) {
        stop("tximport object is missing samples listed in the samplesheet")
    }
    cts     <- txi$counts[, ord, drop = FALSE]
    normMat <- txi$length[, ord, drop = FALSE]

    # Official tximport -> edgeR recipe: turn transcript lengths into an
    # offset that corrects for length and library composition.
    normMat <- normMat / exp(rowMeans(log(normMat)))
    normCts <- cts / normMat
    eff.lib <- calcNormFactors(normCts) * colSums(normCts)
    normMat <- sweep(normMat, 2, eff.lib, "*")
    normMat <- log(normMat)

    y <- DGEList(cts, group = samples$condition)
    y <- scaleOffset(y, normMat)
    keep <- filterByExpr(y)
    y <- y[keep, ]
} else {
    message("Input mode: featureCounts gene counts")
    count_files <- inputs
    counts_list <- list()
    for (s in samples$sample) {
        match_file <- grep(paste0(s, ".featureCounts.txt"), count_files, value = TRUE)
        if (length(match_file) == 0) stop(paste("File not found for", s))

        df <- read.table(match_file[1], header = TRUE, comment.char = "#",
                         stringsAsFactors = FALSE)
        counts <- df[, ncol(df)]
        names(counts) <- df$Geneid
        counts_list[[s]] <- counts
    }
    count_matrix <- do.call(cbind, counts_list)
    rownames(count_matrix) <- names(counts_list[[1]])

    y <- DGEList(counts = count_matrix, group = samples$condition)
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
}

# 3. Design matrix and dispersion
design <- model.matrix(design_formula, data = samples)
y <- estimateDisp(y, design)

# 4. Quasi-likelihood fit
fit <- glmQLFit(y, design)

dir.create("edger_output", showWarnings = FALSE)

# MDS Plot
png(file.path("edger_output", "mds_plot.png"))
plotMDS(y)
dev.off()

# 5. Contrasts.
# Each contrast is built over the 'condition' coefficients of the design and
# tested with glmQLFTest, so any covariates in the design are adjusted for.
cond_levels <- levels(samples$condition)
pairs <- combn(cond_levels, 2)

for (i in seq_len(ncol(pairs))) {
    c1 <- pairs[1, i]
    c2 <- pairs[2, i]

    # Orient so 'REF' is the denominator where present.
    numerator   <- c1
    denominator <- c2
    if (c1 == "REF") { numerator <- c2; denominator <- c1 }

    # Contrast over the design columns: +1 on the numerator's condition
    # coefficient, -1 on the denominator's. A level that is the model's
    # reference (intercept) level has no column and is simply left at 0.
    contrast <- setNames(rep(0, ncol(design)), colnames(design))
    num_col  <- paste0("condition", numerator)
    den_col  <- paste0("condition", denominator)
    if (num_col %in% names(contrast)) contrast[num_col] <-  1
    if (den_col %in% names(contrast)) contrast[den_col] <- -1
    if (all(contrast == 0)) next

    test_res  <- glmQLFTest(fit, contrast = contrast)
    res_table <- annotate_genes(topTags(test_res, n = Inf)$table)

    filename <- paste0("edger_results_", numerator, "_vs_", denominator, ".csv")
    write.csv(res_table, file = file.path("edger_output", filename))

    png(file.path("edger_output", paste0("smear_", numerator, "_vs_", denominator, ".png")))
    plotSmear(test_res, de.tags = rownames(res_table)[res_table$FDR < 0.05])
    dev.off()

    png(file.path("edger_output", paste0("volcano_", numerator, "_vs_", denominator, ".png")),
        width = 760, height = 640)
    draw_volcano(res_table$logFC, res_table$PValue, res_table$FDR,
                 paste(numerator, "vs", denominator))
    dev.off()
}
