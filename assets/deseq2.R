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
# padj is the adjusted p-value (DESeq2 padj, edgeR FDR). The cutoffs are
# arguments rather than literals so that the generated reproduce/ scripts can
# expose them as editable settings -- inlined with them hardcoded, a user could
# change PADJ_CUT / LFC_CUT and see the figure not move.
draw_volcano <- function(lfc, pval, padj, title, padj_cut = 0.05, lfc_cut = 1) {
    ok   <- !is.na(lfc) & !is.na(pval)
    lfc  <- lfc[ok]; pval <- pval[ok]; padj <- padj[ok]
    sig  <- !is.na(padj) & padj < padj_cut & abs(lfc) > lfc_cut
    cols <- ifelse(sig, ifelse(lfc > 0, "#c2255c", "#1c7ed6"), "#ced4da")
    plot(lfc, -log10(pval), pch = 20, cex = 0.55, col = cols,
         xlab = "log2 fold change", ylab = expression(-log[10] ~ italic(p)),
         main = title)
    abline(v = c(-lfc_cut, lfc_cut), lty = 2, col = "#868e96")
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
FIG_RES <- 300   # px per inch for the raster copy. 300 is the usual journal minimum for halftone / combination figures. Affects the PNG only -- the SVG is vector, so resolution is meaningless there.

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
# Shrinkage diagnostic: the same contrast before and after apeglm.
#
# DESeq2-specific, so deliberately NOT part of the block kept identical with
# edger.R -- edgeR's quasi-likelihood test moderates the dispersion rather than
# the effect size, so it has no equivalent.
#
# `df` is the annotated results table, which carries both log2FoldChange (the
# apeglm-shrunken estimate the rest of the pipeline uses) and
# log2FoldChange_MLE (the raw maximum-likelihood one). Significance comes from
# the unshrunken Wald test in both panels -- shrinkage changes the estimate, not
# the test -- which is why a gene can move between panels without changing
# colour.
#
# Every gene is drawn, with no thinning: a published figure has to be the same
# every time it is made, and sampling the non-significant cloud would make it
# depend on the RNG. (The Quarto report thins its interactive version, where
# widget size matters and exact reproducibility does not.)
# ---------------------------------------------------------------------------
draw_shrinkage_ma <- function(df, padj_cut = 0.05, lfc_cut = 1) {
    if (!("log2FoldChange_MLE" %in% names(df))) return(NULL)
    d <- df[!is.na(df$baseMean) & df$baseMean > 0 &
            !is.na(df$log2FoldChange) & !is.na(df$log2FoldChange_MLE), , drop = FALSE]
    if (!nrow(d)) return(NULL)

    sig <- !is.na(d$padj) & d$padj < padj_cut & abs(d$log2FoldChange) > lfc_cut
    dir <- ifelse(sig, ifelse(d$log2FoldChange > 0, "Up", "Down"), "Not significant")

    # rbind rather than a pivot: tidyr is not in this container.
    long <- rbind(
        data.frame(baseMean = d$baseMean, lfc = d$log2FoldChange_MLE,
                   direction = dir, estimate = "Unshrunken (MLE)",
                   stringsAsFactors = FALSE),
        data.frame(baseMean = d$baseMean, lfc = d$log2FoldChange,
                   direction = dir, estimate = "apeglm-shrunken",
                   stringsAsFactors = FALSE))
    long$estimate <- factor(long$estimate,
                            levels = c("Unshrunken (MLE)", "apeglm-shrunken"))

    ggplot(long, aes(baseMean, lfc, colour = direction)) +
        geom_point(size = 0.5, alpha = 0.5) +
        geom_hline(yintercept = 0, colour = "#868e96", linewidth = 0.3) +
        scale_x_log10() +
        facet_wrap(~ estimate) +
        scale_colour_manual(values = c("Up" = "#c2255c", "Down" = "#1c7ed6",
                                       "Not significant" = "#ced4da")) +
        labs(x = "mean of normalised counts", y = "log2 fold change",
             colour = NULL) +
        theme_bw(base_size = 11)
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

# ---------------------------------------------------------------------------
# reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R objects the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits deseq2_output as a whole directory.
#
# The guiding split is that the parent folder is for *reading* and reproduce/ is
# for *running*. Parent CSVs therefore stay plain, so they open by double-click
# in Excel; the reproduce/ copies are gzipped, which read.csv() handles
# transparently and which matters once a results table is tens of MB. Anyone
# wanting to eyeball that data uses the plain copy one directory up.
# ---------------------------------------------------------------------------
REPRO_DIR <- file.path("deseq2_output", "reproduce")
dir.create(REPRO_DIR, showWarnings = FALSE, recursive = TRUE)

# Write a data frame into reproduce/ as a gzipped CSV.
write_repro_csv <- function(df, name) {
    con <- gzfile(file.path(REPRO_DIR, paste0(name, ".csv.gz")), "w")
    on.exit(close(con))
    write.csv(df, con, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Emitter for the reproduce/ scripts.
#
# Each generated script redraws exactly one figure from the artefacts beside
# it. Two rules keep generated code correct:
#
#   * Helper functions are emitted by deparse()-ing the live function object,
#     never hand-copied into a string. The copy in the script therefore cannot
#     drift from the one the pipeline just used. deparse() drops comments --
#     they are not part of a function object -- so the script carries its own
#     header instead.
#   * Every injected value goes through deparse() as well, so quoting and
#     escaping are R's problem rather than ours. Nothing is pasted into code as
#     a raw string.
#
# Generated scripts assume the working directory is their own folder. They
# check their inputs up front and stop with the available contrasts listed,
# rather than letting R emit a bare "cannot open file".
# ---------------------------------------------------------------------------

# Deparse a value into source text: handles quoting, escaping and vectors.
emit_value <- function(x) paste(deparse(x), collapse = "\n")

# Deparse a live function into `name <- function(...) {...}`.
emit_function <- function(name, fn) {
    src    <- deparse(fn)
    src[1] <- paste0(name, " <- ", src[1])
    src
}

# file        path of the script to write
# description one-line summary for the header
# libraries   packages to attach (also reported as versions in the header)
# selectors   named list of the values that pick which figure to draw, e.g.
#             list(CONTRAST = c("A_vs_REF", "B_vs_REF")) or, where a figure is
#             indexed by more than one thing, list(CONTRAST = ..., DIRECTION =
#             c("UP", "DOWN")). Each becomes a constant defaulting to its first
#             value, overridable by positional command-line arguments in the
#             order given. NULL for a figure that is produced once per run.
# constants   named list emitted as the editable settings block, in order
# functions   named list of live functions to inline
# inputs      R expressions (as strings) evaluating to the files it reads
# body        the plotting code
write_repro_script <- function(file, description, libraries = character(0),
                               selectors = NULL, constants = list(),
                               functions = list(), inputs = character(0),
                               body = character(0)) {
    rule <- paste0("# ", strrep("-", 74))
    # Selector constants lead the settings block: they are what a reader changes
    # most often, so they should be the first thing under the heading.
    #
    # Every generated figure is prefixed, so a regenerated one is never mistaken
    # for the pipeline's own output. This matters because the settings above it
    # are editable: change a cutoff and the result is a different figure that
    # would otherwise carry the published figure's exact filename. Exposed as a
    # constant so anyone deliberately replacing a published figure can clear it.
    constants <- c(lapply(selectors, function(v) v[1]), constants,
                   list(OUTPUT_PREFIX = "repro_"))
    # Tolerant version lookup: an optional package (pheatmap) may be absent, and
    # failing to record a version must never abort the analysis task.
    vers <- c(paste0("R ", getRversion()),
              vapply(libraries, function(p)
                  paste0(p, " ", tryCatch(as.character(packageVersion(p)),
                                          error = function(e) "not installed")),
                  character(1)))

    out <- c(
        "#!/usr/bin/env Rscript",
        "#",
        paste0("# ", description),
        "#",
        "# Redraws this figure from the files in this folder. It does not rerun",
        "# the pipeline -- everything needed is here, so the folder can be copied",
        "# anywhere and still work.",
        "#",
        "# It reads its inputs from the working directory, so that must be the",
        "# folder this file is in.",
        "#",
        "#   From a terminal:",
        "#",
        "#     cd <the folder containing this script>",
        paste0("#     Rscript ", basename(file)),
        "#",
        "#   In RStudio: open this file, then",
        "#",
        "#     Session > Set Working Directory > To Source File Location",
        "#",
        "#   and click Source.")

    if (length(selectors)) {
        out <- c(out,
            "#",
            paste0("# To draw a different ",
                   paste(tolower(names(selectors)), collapse = " or "),
                   ", edit the settings below. From a terminal you can instead"),
            "# pass them as arguments, in this order (this does not apply when",
            "# sourcing the file in RStudio):",
            "#",
            paste0("#     Rscript ", basename(file), " ",
                   paste(vapply(selectors, function(v) v[1], character(1)),
                         collapse = " ")),
            "#",
            "# Values available in this folder. Copy a line over the matching",
            "# setting below:",
            "#")
        for (nm in names(selectors))
            out <- c(out, paste0("#     ", nm, " <- ",
                                 vapply(selectors[[nm]], emit_value, character(1))))
    }

    out <- c(out,
        "#",
        paste0("# Generated ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               " by the rnaseq-flow pipeline."),
        paste0("# The published figure was drawn with: ", paste(vers, collapse = "; ")),
        "",
        rule,
        "# Settings -- edit these to adapt the figure.",
        rule)

    if (length(constants)) {
        pad <- max(nchar(names(constants)))
        out <- c(out, vapply(names(constants), function(n)
            paste0(formatC(n, width = -pad), " <- ", emit_value(constants[[n]])),
            character(1)))
    }

    out <- c(out, "", rule, "# No need to edit below here.", rule)

    if (length(libraries))
        out <- c(out, paste0("suppressMessages(library(", libraries, "))"))

    if (length(selectors)) {
        out <- c(out, "",
            "# Optional arguments override the settings above, in this order, so",
            "# the pipeline can drive this same script without editing it.",
            ".args <- commandArgs(trailingOnly = TRUE)")
        for (i in seq_along(selectors))
            out <- c(out, sprintf(
                "if (length(.args) >= %dL && nzchar(.args[%d])) %s <- .args[%d]",
                i, i, names(selectors)[i], i))
        out <- c(out, paste0(".SELECTORS <- ", emit_value(selectors)))
    }

    if (length(inputs)) {
        out <- c(out, "",
            "# Fail with something actionable rather than R's \"cannot open file\".",
            paste0(".inputs  <- c(", paste(inputs, collapse = ",\n              "), ")"),
            ".missing <- .inputs[!file.exists(.inputs)]",
            "if (length(.missing)) {",
            "    .msg <- paste0(\"cannot find: \", paste(.missing, collapse = \", \"),",
            "                   \"\\n  Run this script with the working directory set\",",
            "                   \" to its own folder.\")")
        if (length(selectors))
            out <- c(out,
            "    .msg <- paste0(.msg, \"\\n  Values available here:\")",
            "    for (.n in names(.SELECTORS))",
            "        .msg <- paste0(.msg, \"\\n    \", .n, \" : \",",
            "                       paste(.SELECTORS[[.n]], collapse = \", \"))")
        out <- c(out,
            "    stop(.msg, call. = FALSE)",
            "}")
    }

    # Provenance, fixed at generation time, plus the settings to record.
    out <- c(out, "",
        paste0(".SETTINGS  <- ", emit_value(names(constants))),
        paste0(".LIBS      <- ", emit_value(unname(libraries))),
        paste0(".PIPELINE  <- ", emit_value("rnaseq-flow")),
        paste0(".URL       <- ",
               emit_value("https://github.com/nicholas-owen/rnaseq-flow")),
        paste0(".GENERATED <- ",
               emit_value(format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
        paste0(".ORIGINAL  <- ", emit_value(paste(vers, collapse = "; "))))

    out <- c(out, "",
        "# Alongside each figure, record what was written, with which settings,",
        "# and where it came from -- so a figure that leaves this folder can",
        "# still be traced back and its settings checked.",
        "write_info <- function(base) {",
        "    .pad <- max(nchar(.SETTINGS))",
        "    vals <- vapply(.SETTINGS, function(n)",
        "        paste0(\"  \", formatC(n, width = -.pad), \" = \",",
        "               paste(deparse(get(n)), collapse = \"\")), character(1))",
        "    libs <- if (length(.LIBS))",
        "        vapply(.LIBS, function(p) paste0(\"  \", p, \" \",",
        "            tryCatch(as.character(packageVersion(p)),",
        "                     error = function(e) \"not installed\")), character(1))",
        "        else character(0)",
        "    writeLines(c(",
        "        paste0(\"Figure written : \", base, \".png\"),",
        "        paste0(\"                 \", base, \".svg\"),",
        "        paste0(\"Date           : \", format(Sys.time(), \"%Y-%m-%d %H:%M:%S\")),",
        "        \"\",",
        "        \"Settings used\",",
        "        vals,",
        "        \"\",",
        "        \"Regenerated with\",",
        "        paste0(\"  R \", getRversion()),",
        "        libs,",
        "        \"\",",
        "        \"Provenance\",",
        "        paste0(\"  Pipeline         : \", .PIPELINE),",
        "        paste0(\"  Repository       : \", .URL),",
        "        paste0(\"  Script generated : \", .GENERATED),",
        "        paste0(\"  Published figure drawn with : \", .ORIGINAL),",
        "        \"\",",
        "        \"This figure was regenerated from the published data by a\",",
        "        \"reproduce/ script. It is not the pipeline's own output: if the\",",
        "        \"settings above have been edited, it will differ from the\",",
        "        \"published figure of the same name in the parent folder.\"),",
        "        paste0(base, \".info.txt\"))",
        "}",
        "",
        "# Prefix the name, draw both formats, then write the .info.txt beside it.",
        "save_and_document <- function(name, draw, width = 7, height = 7) {",
        "    base <- paste0(OUTPUT_PREFIX, name)",
        "    save_figure(\".\", base, draw, width = width, height = height)",
        "    write_info(base)",
        "    message(\"wrote \", base, \".png, \", base, \".svg and \", base, \".info.txt\")",
        "    invisible(base)",
        "}")

    for (nm in names(functions))
        out <- c(out, "", emit_function(nm, functions[[nm]]))

    out <- c(out, "", body, "")
    writeLines(out, file)
    invisible(file)
}

# Variance Stabilizing Transformation for PCA/Heatmap
vsd <- vst(dds, blind = FALSE)

# PCA. The plot is built by a function rather than inline so that the same code
# can be inlined into reproduce/pca_plot.R -- the published figure and the
# regenerated one are then the same drawing, not two that happen to look alike.
# `pca` needs PC1, PC2, condition and name columns; the variance percentages are
# passed separately because they belong to the axis labels, not to any sample.
draw_pca <- function(pca, pct1, pct2, point_size = 3) {
  ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
    geom_point(size = point_size) +
    geom_text(vjust = 1.5, hjust = 1.5) +
    xlab(paste0("PC1: ", pct1, "% variance")) +
    ylab(paste0("PC2: ", pct2, "% variance")) +
    coord_fixed() +
    theme_bw()
}

pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

save_figure("deseq2_output", "pca_plot",
            function() print(draw_pca(pcaData, percentVar[1], percentVar[2])))

# The plotted coordinates, so the PCA can be redrawn (or restyled) without
# rerunning DESeq2. percentVar is carried as columns because a flat CSV has
# nowhere else to put it.
pca_out <- pcaData
pca_out$percentVar_PC1 <- percentVar[1]
pca_out$percentVar_PC2 <- percentVar[2]
write_repro_csv(pca_out, "pca_data")

# Heatmap of the 20 most variable genes (variance computed with base R so no
# extra package dependency is needed).
gene_var    <- apply(assay(vsd), 1, var)
topVarGenes <- head(order(gene_var, decreasing = TRUE), 20)
mat <- assay(vsd)[topVarGenes, ]
df <- as.data.frame(colData(dds)[, c("condition")])
rownames(df) <- colnames(mat)
colnames(df) <- "condition"

# As with the PCA, the drawing is a function so the identical code can be
# inlined into reproduce/heatmap_top_var.R. The clustering arguments are named
# explicitly rather than left to pheatmap's defaults: hclust is deterministic,
# so stating the distance and linkage is what makes the row/column ordering
# reproducible on purpose instead of by accident.
draw_heatmap <- function(mat, ann, dist_method = "euclidean",
                         hclust_method = "complete", scale = "none") {
    pheatmap::pheatmap(mat, annotation_col = ann, show_rownames = TRUE,
                       clustering_distance_rows = dist_method,
                       clustering_distance_cols = dist_method,
                       clustering_method        = hclust_method,
                       scale                    = scale)
}

# The variance-stabilised values behind the heatmap, annotated like the results
# tables. check.names = FALSE keeps sample names verbatim. Written whether or
# not pheatmap is available, so the data exists even if the figure does not.
write_repro_csv(annotate_genes(data.frame(mat, check.names = FALSE)),
                "heatmap_top_var")
# The column annotation (sample -> condition), so the heatmap is reproducible
# from these two files alone.
write_repro_csv(data.frame(sample = rownames(df), condition = df$condition),
                "heatmap_top_var_annotation")

if (have_pheatmap) {
    save_figure("deseq2_output", "heatmap_top_var",
                function() draw_heatmap(mat, df),
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

# Contrasts actually produced, appended by run_contrast() as it goes. The
# reproduce/ scripts list these in their headers as copy-pasteable CONTRAST
# lines, so the list has to reflect what was really written to disk.
repro_contrasts <- character(0)

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
    # The unshrunken (maximum-likelihood) fold change, kept beside the shrunken
    # one so the effect of shrinkage can actually be seen: the report draws the
    # two as a before/after MA pair, which is the standard check on how far
    # apeglm pulled in the low-count, high-variance genes. Without this column
    # the published table records only the outcome of shrinkage, never its size.
    # res and res_mle come from the same fit and share a row order.
    res_df$log2FoldChange_MLE <- res_mle$log2FoldChange
    res_df <- res_df[, c("baseMean", "log2FoldChange", "log2FoldChange_MLE",
                         "lfcSE", "stat", "pvalue", "padj")]
    res_df <- res_df[order(res_df$pvalue), ]

    stem <- paste0("deseq2_results_", condA, "_vs_", condB)
    # row.names = FALSE: annotate_genes() has already promoted the row names to
    # an explicit gene_id column, so writing them again would duplicate the IDs
    # under a blank header.
    res_out <- annotate_genes(res_df)
    write.csv(res_out, file = file.path("deseq2_output", paste0(stem, ".csv")),
              row.names = FALSE)

    # reproduce/: both result objects, plus a gzipped copy of the table. `res`
    # is the apeglm-shrunken object the MA plot is drawn from; `res_mle` is the
    # unshrunken Wald fit that supplies the `stat` column, kept so a reworked
    # figure has everything without rerunning DESeq2. The object also carries
    # provenance the CSV cannot: the mcols() column descriptions, the alpha and
    # independent-filtering threshold, and the apeglm prior. Named after the
    # data, not the figure, since the MA and volcano scripts share them.
    saveRDS(list(res = res, res_mle = res_mle),
            file.path(REPRO_DIR, paste0(stem, ".rds")))
    write_repro_csv(res_out, stem)
    repro_contrasts <<- c(repro_contrasts, paste0(condA, "_vs_", condB))

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

    # Before/after apeglm, from the table just written -- so the figure, the
    # reproduce script and the report's interactive version all draw the same
    # thing from the same source.
    shrink <- draw_shrinkage_ma(res_out)
    if (!is.null(shrink)) {
        save_figure("deseq2_output", paste0("shrinkage_ma_", condA, "_vs_", condB),
                    function() print(shrink), width = 9, height = 5)
    }
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

# ---------------------------------------------------------------------------
# 7. reproduce/ scripts -- one per figure type, with the contrast in a single
# constant at the top rather than a loop, so the code stays flat and literal
# for readers who are not comfortable in R. The width/height constants match
# what the pipeline used above, so a regenerated figure matches the published
# one.
# ---------------------------------------------------------------------------

write_repro_script(
    file        = file.path(REPRO_DIR, "pca_plot.R"),
    description = "DESeq2 PCA of variance-stabilised counts.",
    libraries   = "ggplot2",
    constants   = list(POINT_SIZE = 3, WIDTH = 7, HEIGHT = 7, FIG_RES = FIG_RES),
    functions   = list(draw_pca = draw_pca, save_figure = save_figure),
    inputs      = '"pca_data.csv.gz"',
    body        = c(
        '# The coordinates were computed by DESeq2::plotPCA and saved; this',
        '# redraws them, so no variance-stabilising transform is repeated.',
        'pca <- read.csv("pca_data.csv.gz", stringsAsFactors = FALSE)',
        '',
        'save_and_document("pca_plot",',
        '                  function() print(draw_pca(pca, pca$percentVar_PC1[1],',
        '                                            pca$percentVar_PC2[1],',
        '                                            point_size = POINT_SIZE)),',
        '                  width = WIDTH, height = HEIGHT)'))

write_repro_script(
    file        = file.path(REPRO_DIR, "heatmap_top_var.R"),
    description = "Heatmap of the 20 most variable genes (variance-stabilised counts).",
    libraries   = "pheatmap",
    constants   = list(DIST_METHOD = "euclidean", HCLUST_METHOD = "complete",
                       SCALE = "none", WIDTH = 8, HEIGHT = 7, FIG_RES = FIG_RES),
    functions   = list(draw_heatmap = draw_heatmap, save_figure = save_figure),
    inputs      = c('"heatmap_top_var.csv.gz"',
                    '"heatmap_top_var_annotation.csv.gz"'),
    body        = c(
        '# Rebuild the matrix pheatmap expects: gene ids as row names, samples as',
        '# columns. check.names = FALSE keeps sample names exactly as written.',
        'tbl <- read.csv("heatmap_top_var.csv.gz", check.names = FALSE,',
        '                stringsAsFactors = FALSE)',
        'ann <- read.csv("heatmap_top_var_annotation.csv.gz", stringsAsFactors = FALSE)',
        '',
        '.meta <- c("gene_id", "gene_name", "gene_biotype")',
        'mat   <- as.matrix(tbl[, setdiff(names(tbl), .meta), drop = FALSE])',
        '# Row labels are the gene ids, matching the published figure. To label',
        '# by gene symbol instead, use tbl$gene_name here.',
        'rownames(mat) <- tbl$gene_id',
        '',
        'ann_df <- data.frame(condition = ann$condition, row.names = ann$sample)',
        '',
        'save_and_document("heatmap_top_var",',
        '                  function() draw_heatmap(mat, ann_df,',
        '                                          dist_method   = DIST_METHOD,',
        '                                          hclust_method = HCLUST_METHOD,',
        '                                          scale         = SCALE),',
        '                  width = WIDTH, height = HEIGHT)'))

if (length(repro_contrasts)) {
    write_repro_script(
        file        = file.path(REPRO_DIR, "maplot.R"),
        description = "DESeq2 MA plot: apeglm-shrunken log2 fold change vs mean expression.",
        libraries   = "DESeq2",
        selectors   = list(CONTRAST = repro_contrasts),
        constants   = list(YLIM = c(-2, 2), WIDTH = 7, HEIGHT = 7,
                           FIG_RES = FIG_RES),
        functions   = list(save_figure = save_figure),
        inputs      = 'paste0("deseq2_results_", CONTRAST, ".rds")',
        body        = c(
            '# The saved object holds both fits: $res is apeglm-shrunken (what the',
            '# figure shows) and $res_mle is the unshrunken Wald fit.',
            'fits <- readRDS(paste0("deseq2_results_", CONTRAST, ".rds"))',
            '',
            'save_and_document(paste0("maplot_", CONTRAST),',
            '                  function() plotMA(fits$res, ylim = YLIM,',
            '                                    main = paste(sub("_vs_", " vs ", CONTRAST),',
            '                                                 "(apeglm-shrunk LFC)")),',
            '                  width = WIDTH, height = HEIGHT)'))

    write_repro_script(
        file        = file.path(REPRO_DIR, "shrinkage_ma.R"),
        description = "DESeq2 shrinkage diagnostic: the same contrast before and after apeglm.",
        libraries   = "ggplot2",
        selectors   = list(CONTRAST = repro_contrasts),
        constants   = list(PADJ_CUT = 0.05, LFC_CUT = 1, WIDTH = 9, HEIGHT = 5,
                           FIG_RES = FIG_RES),
        functions   = list(draw_shrinkage_ma = draw_shrinkage_ma,
                           save_figure       = save_figure),
        inputs      = 'paste0("deseq2_results_", CONTRAST, ".csv.gz")',
        body        = c(
            '# Both estimates are columns of the results table, so this needs no',
            '# saved object: log2FoldChange is the apeglm-shrunken fold change and',
            '# log2FoldChange_MLE the raw maximum-likelihood one.',
            'res <- read.csv(paste0("deseq2_results_", CONTRAST, ".csv.gz"),',
            '                stringsAsFactors = FALSE)',
            '',
            'save_and_document(paste0("shrinkage_ma_", CONTRAST),',
            '                  function() print(draw_shrinkage_ma(res,',
            '                                                     padj_cut = PADJ_CUT,',
            '                                                     lfc_cut  = LFC_CUT)),',
            '                  width = WIDTH, height = HEIGHT)'))

    write_repro_script(
        file        = file.path(REPRO_DIR, "volcano.R"),
        description = "DESeq2 volcano plot: log2 fold change vs -log10 p-value.",
        selectors   = list(CONTRAST = repro_contrasts),
        constants   = list(PADJ_CUT = 0.05, LFC_CUT = 1,
                           WIDTH = 8, HEIGHT = 6.5, FIG_RES = FIG_RES),
        functions   = list(draw_volcano = draw_volcano, save_figure = save_figure),
        inputs      = 'paste0("deseq2_results_", CONTRAST, ".csv.gz")',
        body        = c(
            '# Drawn entirely from the results table -- no R object needed.',
            'res <- read.csv(paste0("deseq2_results_", CONTRAST, ".csv.gz"),',
            '                stringsAsFactors = FALSE)',
            '',
            'save_and_document(paste0("volcano_", CONTRAST),',
            '                  function() draw_volcano(res$log2FoldChange,',
            '                                          res$pvalue, res$padj,',
            '                                          sub("_vs_", " vs ", CONTRAST),',
            '                                          padj_cut = PADJ_CUT,',
            '                                          lfc_cut  = LFC_CUT),',
            '                  width = WIDTH, height = HEIGHT)'))
}

message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
