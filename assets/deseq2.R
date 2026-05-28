#!/usr/bin/env Rscript

# Differential expression with DESeq2.
#
# Usage:
#   Rscript deseq2.R <samplesheet.csv> <design> <gene_info.tsv> <featureCounts files...>
#   Rscript deseq2.R <samplesheet.csv> <design> <gene_info.tsv> <txi.rds>
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

if (length(args) < 4) {
  stop("Usage: deseq2.R <samplesheet.csv> <design> <gene_info.tsv> <featureCounts files... | txi.rds>")
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

# 1. Read Samplesheet
samples <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
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

# 2. Build the DESeqDataSet - from a tximport object or from featureCounts files
if (length(rds_input) >= 1) {
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

ggsave(file.path("deseq2_output", "pca_plot.png"), plot = p)

# Heatmap of the 20 most variable genes (variance computed with base R so no
# extra package dependency is needed).
gene_var    <- apply(assay(vsd), 1, var)
topVarGenes <- head(order(gene_var, decreasing = TRUE), 20)
mat <- assay(vsd)[topVarGenes, ]
df <- as.data.frame(colData(dds)[, c("condition")])
rownames(df) <- colnames(mat)
colnames(df) <- "condition"

if (have_pheatmap) {
    png(file.path("deseq2_output", "heatmap_top_var.png"))
    pheatmap::pheatmap(mat, annotation_col = df, show_rownames = TRUE)
    dev.off()
} else {
    message("pheatmap not available - skipping heatmap_top_var.png")
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
    write.csv(annotate_genes(res_df),
              file = file.path("deseq2_output", filename))

    png(file.path("deseq2_output", paste0("maplot_", condA, "_vs_", condB, ".png")))
    plotMA(res, main = paste(condA, "vs", condB, "(apeglm-shrunk LFC)"),
           ylim = c(-2, 2))
    dev.off()

    png(file.path("deseq2_output", paste0("volcano_", condA, "_vs_", condB, ".png")),
        width = 760, height = 640)
    draw_volcano(res$log2FoldChange, res$pvalue, res$padj, paste(condA, "vs", condB))
    dev.off()
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
