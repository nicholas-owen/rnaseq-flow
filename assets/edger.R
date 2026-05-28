#!/usr/bin/env Rscript

# Differential expression with edgeR.
#
# Usage:
#   Rscript edger.R <samplesheet.csv> <design> <gene_info.tsv> <featureCounts files...>
#   Rscript edger.R <samplesheet.csv> <design> <gene_info.tsv> <txi.rds>
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

if (length(args) < 4) {
  stop("Usage: edger.R <samplesheet.csv> <design> <gene_info.tsv> <featureCounts files... | txi.rds>")
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
    # falls back to a version-insensitive match for any unmatched IDs.
    if (is.null(geneinfo)) return(df)
    strip <- function(x) sub("\\.[0-9]+$", "", x)
    m  <- match(rownames(df), geneinfo$gene_id)
    na <- is.na(m)
    if (any(na)) m[na] <- match(strip(rownames(df)[na]), strip(geneinfo$gene_id))
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

# 1. Read samplesheet
samples <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
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

# 2. Build the DGEList - from a tximport object or from featureCounts files
if (length(rds_input) >= 1) {
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
