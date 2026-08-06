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
# (FDR < 0.05 and |logFC| > 1; up in red, down in blue, the rest grey).
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

# ---------------------------------------------------------------------------
# reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R objects the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits edger_output as a whole directory.
#
# The guiding split is that the parent folder is for *reading* and reproduce/ is
# for *running*. Parent CSVs therefore stay plain, so they open by double-click
# in Excel; the reproduce/ copies are gzipped, which read.csv() handles
# transparently and which matters once a results table is tens of MB. Anyone
# wanting to eyeball that data uses the plain copy one directory up.
# ---------------------------------------------------------------------------
REPRO_DIR <- file.path("edger_output", "reproduce")
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
# contrasts   character vector, or NULL for a figure that is not per-contrast
# constants   named list emitted as the editable settings block, in order
# functions   named list of live functions to inline
# inputs      R expressions (as strings) evaluating to the files it reads
# body        the plotting code
write_repro_script <- function(file, description, libraries = character(0),
                               contrasts = NULL, constants = list(),
                               functions = list(), inputs = character(0),
                               body = character(0)) {
    rule <- paste0("# ", strrep("-", 74))
    # Every generated figure is prefixed, so a regenerated one is never mistaken
    # for the pipeline's own output. This matters because the settings above it
    # are editable: change a cutoff and the result is a different figure that
    # would otherwise carry the published figure's exact filename. Exposed as a
    # constant so anyone deliberately replacing a published figure can clear it.
    constants <- c(constants, list(OUTPUT_PREFIX = "repro_"))
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

    if (!is.null(contrasts)) {
        out <- c(out,
            "#",
            "# To draw a different contrast, edit the CONTRAST setting below.",
            "# From a terminal you can instead pass it as the first argument",
            "# (this does not apply when sourcing the file in RStudio):",
            "#",
            paste0("#     Rscript ", basename(file), " ", contrasts[1]),
            "#",
            "# Contrasts available in this folder. Copy one of these lines over",
            "# the CONTRAST setting below:",
            "#",
            paste0("#     CONTRAST <- ", vapply(contrasts, emit_value, character(1))))
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

    if (!is.null(contrasts)) {
        out <- c(out, "",
            "# An optional first argument overrides CONTRAST, so the pipeline can",
            "# drive this same script without editing it.",
            ".args <- commandArgs(trailingOnly = TRUE)",
            "if (length(.args) >= 1L && nzchar(.args[1])) CONTRAST <- .args[1]",
            paste0(".CONTRASTS <- ", emit_value(contrasts)))
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
        if (!is.null(contrasts))
            out <- c(out,
            "    .msg <- paste0(.msg, \"\\n  Contrasts available here: \",",
            "                   paste(.CONTRASTS, collapse = \", \"))")
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

# MDS Plot. plotMDS() returns the computed coordinates invisibly, so run it once
# with plot = FALSE to capture them for the CSV, then let save_figure() redraw
# it per device.
mds <- plotMDS(y, plot = FALSE)
save_figure("edger_output", "mds_plot", function() plotMDS(y))

# The MDS object itself goes to reproduce/, not the DGEList: limma has a plot
# method for the MDS class, so re-plotting this object reproduces the figure
# exactly without recomputing distances -- verified as byte-identical PNG
# output. A DGEList would carry the entire count matrix (~28x larger at 2000
# genes, and the gap widens with gene count, since an MDS object scales only
# with the number of samples).
saveRDS(mds, file.path(REPRO_DIR, "mds_plot.rds"))

# The plotted coordinates, so the MDS can be redrawn without rerunning edgeR.
# var.explained is only present in newer edgeR, so it is added when available.
# Matched on sample name rather than position: every input path orders the
# columns by samples$sample, but a mismatch here would mislabel the plot
# silently rather than fail.
mds_out <- data.frame(sample    = colnames(y),
                      condition = as.character(
                          samples$condition[match(colnames(y), samples$sample)]),
                      x         = as.numeric(mds$x),
                      y         = as.numeric(mds$y),
                      stringsAsFactors = FALSE)
if (!is.null(mds$var.explained)) {
    mds_out$var_explained_dim1 <- round(100 * mds$var.explained[1], 1)
    mds_out$var_explained_dim2 <- round(100 * mds$var.explained[2], 1)
}
# Figure data, so it lives in reproduce/ rather than the parent folder: the
# parent carries the analysis tables, reproduce/ carries what redraws a figure.
write_repro_csv(mds_out, "mds_data")

# 5. Contrasts.
# Each contrast is built over the 'condition' coefficients of the design and
# tested with glmQLFTest, so any covariates in the design are adjusted for.
cond_levels <- levels(samples$condition)
pairs <- combn(cond_levels, 2)

# Contrasts actually produced, collected as the loop runs. The reproduce/
# scripts list these in their headers as copy-pasteable CONTRAST lines, so the
# list has to reflect what exists rather than what was theoretically possible
# (a contrast whose coefficients are all zero is skipped below).
repro_contrasts <- character(0)

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

    stem <- paste0("edger_results_", numerator, "_vs_", denominator)
    # row.names = FALSE: annotate_genes() has already promoted the row names to
    # an explicit gene_id column, so writing them again would duplicate the IDs
    # under a blank header.
    write.csv(res_table, file = file.path("edger_output", paste0(stem, ".csv")),
              row.names = FALSE)

    # reproduce/: the DGELRT that plotSmear consumes (and that a rewritten
    # volcano would need), plus a gzipped copy of the results table. Both are
    # named after the *data*, not the figure, because the smear and volcano
    # scripts share them -- naming per figure would store the same object twice
    # for every contrast.
    saveRDS(test_res, file.path(REPRO_DIR, paste0(stem, ".rds")))
    write_repro_csv(res_table, stem)
    repro_contrasts <- c(repro_contrasts,
                         paste0(numerator, "_vs_", denominator))

    # No separate CSV for the smear and volcano plots: both are drawn entirely
    # from logCPM / logFC / PValue / FDR, every one of which is already a column
    # of the edger_results_*.csv written just above.
    de_tags <- rownames(res_table)[res_table$FDR < 0.05]
    save_figure("edger_output", paste0("smear_", numerator, "_vs_", denominator),
                function() plotSmear(test_res, de.tags = de_tags))

    save_figure("edger_output", paste0("volcano_", numerator, "_vs_", denominator),
                function() draw_volcano(res_table$logFC, res_table$PValue,
                                        res_table$FDR,
                                        paste(numerator, "vs", denominator)),
                width = 8, height = 6.5)
}

# ---------------------------------------------------------------------------
# 6. reproduce/ scripts -- one per figure type, with the contrast in a single
# constant at the top rather than a loop, so the code stays flat and literal
# for readers who are not comfortable in R. The width/height constants match
# what the pipeline used above, so a regenerated figure matches the published
# one.
# ---------------------------------------------------------------------------

write_repro_script(
    file        = file.path(REPRO_DIR, "mds_plot.R"),
    description = "edgeR MDS plot: samples in two dimensions of leading log2 fold change.",
    libraries   = "limma",
    constants   = list(WIDTH = 7, HEIGHT = 7, FIG_RES = 150),
    functions   = list(save_figure = save_figure),
    inputs      = '"mds_plot.rds"',
    body        = c(
        '# The saved MDS object already holds the computed coordinates, so this',
        '# redraws the figure without recomputing distances from the counts.',
        'mds <- readRDS("mds_plot.rds")',
        '',
        'save_and_document("mds_plot", function() plotMDS(mds),',
        '                  width = WIDTH, height = HEIGHT)'))

if (length(repro_contrasts)) {
    write_repro_script(
        file        = file.path(REPRO_DIR, "smear.R"),
        description = "edgeR smear plot: log fold change vs average abundance.",
        libraries   = "edgeR",
        contrasts   = repro_contrasts,
        constants   = list(CONTRAST = repro_contrasts[1], FDR_CUT = 0.05,
                           WIDTH = 7, HEIGHT = 7, FIG_RES = 150),
        functions   = list(save_figure = save_figure),
        inputs      = c('paste0("edger_results_", CONTRAST, ".rds")',
                        'paste0("edger_results_", CONTRAST, ".csv.gz")'),
        body        = c(
            '# The test object supplies the plot; the results table supplies which',
            '# genes to highlight, so the FDR cutoff above is a real setting.',
            'test_res <- readRDS(paste0("edger_results_", CONTRAST, ".rds"))',
            'res      <- read.csv(paste0("edger_results_", CONTRAST, ".csv.gz"),',
            '                     stringsAsFactors = FALSE)',
            'de_tags  <- res$gene_id[!is.na(res$FDR) & res$FDR < FDR_CUT]',
            '',
            'save_and_document(paste0("smear_", CONTRAST),',
            '                  function() plotSmear(test_res, de.tags = de_tags),',
            '                  width = WIDTH, height = HEIGHT)'))

    write_repro_script(
        file        = file.path(REPRO_DIR, "volcano.R"),
        description = "edgeR volcano plot: log2 fold change vs -log10 p-value.",
        contrasts   = repro_contrasts,
        constants   = list(CONTRAST = repro_contrasts[1], FDR_CUT = 0.05,
                           LFC_CUT = 1, WIDTH = 8, HEIGHT = 6.5, FIG_RES = 150),
        functions   = list(draw_volcano = draw_volcano, save_figure = save_figure),
        inputs      = 'paste0("edger_results_", CONTRAST, ".csv.gz")',
        body        = c(
            '# Drawn entirely from the results table -- no R object needed.',
            'res <- read.csv(paste0("edger_results_", CONTRAST, ".csv.gz"),',
            '                stringsAsFactors = FALSE)',
            '',
            'save_and_document(paste0("volcano_", CONTRAST),',
            '                  function() draw_volcano(res$logFC, res$PValue,',
            '                                          res$FDR,',
            '                                          sub("_vs_", " vs ", CONTRAST),',
            '                                          padj_cut = FDR_CUT,',
            '                                          lfc_cut  = LFC_CUT),',
            '                  width = WIDTH, height = HEIGHT)'))
}

message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
