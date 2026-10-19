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

# ---------------------------------------------------------------------------
# Write one figure as both PNG and SVG.
#
# `draw` is a function of no arguments that renders the plot. It is called once
# per device, because a plot can only be written while its device is open.
# plotSpliceDGE() is base graphics and draws straight to the current device,
# so wrapping the call is all that is needed.
#
# SVG is written with grDevices::svg() rather than ggsave(..., ".svg"):
# ggsave() delegates SVG to the svglite package, which none of the pipeline
# containers carry, while svg() needs only cairo, which all of them have. Both
# devices are given the same inch dimensions so the raster and vector copies
# are identical in layout.
#
# Mirrors save_figure() in deseq2.R, edger.R, gsea.R, gprofiler.R and
# isoform_switch.R. The copies are character-identical and are edited together.
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
# Baseline level from --reference_level (default "REF"); see deseq2.R for why it
# arrives as an environment variable rather than an argument.
REF_LEVEL <- Sys.getenv("RNASEQ_FLOW_REFERENCE_LEVEL", "REF")
if (REF_LEVEL %in% levels(condition)) condition <- relevel(condition, ref = REF_LEVEL)
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
  # dropInfReps: only txi$counts is used below, and importing the inferential
  # replicates needs jsonlite for Salmon or rhdf5 for kallisto. See the note in
  # dtu_dexseq.R.
  txi <- tximport(files, type = aligner, txOut = TRUE,
                  countsFromAbundance = "scaledTPM",
                  dropInfReps = TRUE, ignoreTxVersion = TRUE)
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

# ---------------------------------------------------------------------------
# reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R objects the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits diffsplice_output as a whole directory.
#
# The guiding split is that the parent folder is for *reading* and reproduce/ is
# for *running*. Parent CSVs therefore stay plain, so they open by double-click
# in Excel; the reproduce/ copies are gzipped, which read.csv() handles
# transparently and which matters once a results table is tens of MB. Anyone
# wanting to eyeball that data uses the plain copy one directory up.
# ---------------------------------------------------------------------------
REPRO_DIR <- file.path(outdir, "reproduce")
# Created only once a figure exists, unlike the four scripts this is copied
# from. An analysis that plots nothing has nothing to reproduce.

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

# ---------------------------------------------------------------------------
# Figure settings.
#
# plotSpliceDGE is edgeR's own plot for this test and needs nothing the
# diffSpliceDGE object does not already carry. It draws one point per feature
# of a gene, at the feature's log-fold-change against the gene's overall
# log-fold-change, with the features that pass FDR marked in red. That is
# exactly what diffSpliceDGE tests, which is why this is the figure rather
# than a volcano: a volcano of feature FDR would hide the within-gene pattern
# the test is about.
#
# Wide and short, because the x-axis labels are feature IDs written vertically
# at cex 0.6, and a gene with forty exons needs the width more than the height.
# ---------------------------------------------------------------------------
PLOT_TOP_N  <- 10     # genes per comparison, by Simes FDR
PLOT_WIDTH  <- 9
PLOT_HEIGHT <- 5.5
PLOT_FDR    <- 0.05   # features at or below this are marked on the plot

# The baseline every coefficient is measured against. relevel() above moved
# --reference_level to the front when the samplesheet has it; otherwise it is
# whichever level sorted first. Either way it is what model.matrix dropped, so
# naming the figure folder <condition>_vs_<baseline> is accurate, and matches
# rmats_output/ and isoform_switch_output/.
BASELINE <- levels(condition)[1]

# '/' and ':' turn up in gene symbols and are not legal in a path everywhere.
path_safe <- function(x) gsub("[/:]", "-", x)

# ---- diffSpliceDGE, once per condition coefficient (vs the reference) ------
feature_tables <- list(); gene_tables <- list(); plotted <- list()
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

  # ---- figures for this comparison ----------------------------------------
  # topSpliceDGE already returns the gene table sorted by p-value, so the head
  # is the top of the ranking.
  if (nrow(gt) == 0) next
  contrast <- paste0(cmp, "_vs_", BASELINE)
  top      <- head(gt, PLOT_TOP_N)

  pad   <- formatC(seq_len(nrow(top)), width = nchar(nrow(top)), flag = "0")
  gname <- geneinfo$gene_name[match(top$GeneID, geneinfo$gene_id)]
  figure <- paste0(pad, "_splice_plot_", path_safe(top$GeneID))
  # geneinfo falls back to the gene ID when the GTF carries no gene_name, so
  # the suffix is added only where it says something the ID does not.
  named <- !is.na(gname) & nzchar(gname) & gname != top$GeneID
  figure[named] <- paste0(figure[named], "_aka_", path_safe(gname[named]))

  cmp_dir <- file.path(outdir, contrast)
  dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)

  for (i in seq_len(nrow(top))) {
    gid <- top$GeneID[i]; stem <- figure[i]
    tryCatch(
      save_figure(cmp_dir, stem,
                  function() plotSpliceDGE(ds, geneid = gid,
                                           genecolname = "GeneID",
                                           FDR = PLOT_FDR),
                  width = PLOT_WIDTH, height = PLOT_HEIGHT),
      error = function(e) {
        # The tables are the result; one unplottable gene is not worth failing
        # the task for. save_figure() opens a device before it draws, so an
        # error part way through leaves that device open and the next gene
        # would be drawn into it. Rscript starts with only the null device,
        # device 1, so anything above that is this loop's leak. png() also
        # creates its file before anything is drawn, so remove the blank
        # image rather than publish it.
        while (dev.cur() > 1) dev.off()
        unlink(file.path(cmp_dir, paste0(stem, c(".png", ".svg"))))
        message("Could not plot ", gid, ": ", conditionMessage(e))
      })
  }

  # Both formats have to be there. The SVG is written second, so a figure with
  # only a PNG is one that failed between the two devices.
  stems <- file.path(cmp_dir, figure)
  ok    <- file.exists(paste0(stems, ".png")) & file.exists(paste0(stems, ".svg"))
  if (!any(ok)) next
  message(sprintf("  %s: wrote %d splice %s.", contrast, sum(ok),
                  if (sum(ok) == 1) "figure" else "figures"))

  # The diffSpliceDGE object is what plotSpliceDGE redraws from, and no table
  # substitutes for it: the plot needs every feature's coefficient in gene
  # order, which the published tables do not preserve. It is affordable whole,
  # unlike edger.R's DGEList, because diffSpliceDGE returns coefficients and
  # p-values and keeps no count matrix.
  dir.create(REPRO_DIR, showWarnings = FALSE, recursive = TRUE)
  saveRDS(ds, file.path(REPRO_DIR, paste0("diffsplice_", contrast, ".rds")))

  plotted[[contrast]] <- data.frame(
    comparison = contrast,
    rank       = seq_len(nrow(top))[ok],
    GeneID     = top$GeneID[ok],
    gene_name  = gname[ok],
    NExons     = top$NExons[ok],
    FDR        = top$FDR[ok],
    figure     = figure[ok],
    stringsAsFactors = FALSE)
}

feature_results <- annotate_genes(do.call(rbind, feature_tables))
gene_results    <- annotate_genes(do.call(rbind, gene_tables))
write.csv(feature_results, feature_csv, row.names = FALSE)
write.csv(gene_results,    gene_csv,    row.names = FALSE)

# ---- reproduce/ -----------------------------------------------------------
if (length(plotted)) {
  index <- do.call(rbind, plotted)
  rownames(index) <- NULL

  # One row per published figure. splice_plot.R reads this rather than parsing
  # file names, and it doubles as the index of what is in the comparison
  # folders.
  write_repro_csv(index, "splice_plot_genes")

  # The two published tables, gzipped.
  write_repro_csv(feature_results, paste0("diffsplice_", level, "_results"))
  write_repro_csv(gene_results,    "diffsplice_gene_results")

  write_repro_script(
    file        = file.path(REPRO_DIR, "splice_plot.R"),
    description = paste("edgeR diffSplice plot: each feature's log-fold-change",
                        "against its gene's, for one gene."),
    libraries   = "edgeR",
    # A GENE selector rather than the CONTRAST the DE and enrichment scripts
    # use, because these figures are per gene. COMPARISON picks the condition
    # pair and only does any work on a design with more than two conditions.
    selectors   = list(GENE       = unique(index$GeneID),
                       COMPARISON = unique(index$comparison)),
    constants   = list(FDR_CUT = PLOT_FDR, WIDTH = PLOT_WIDTH,
                       HEIGHT = PLOT_HEIGHT, FIG_RES = FIG_RES),
    functions   = list(save_figure = save_figure),
    inputs      = c('"splice_plot_genes.csv.gz"',
                    'paste0("diffsplice_", COMPARISON, ".rds")'),
    body        = c(
      '# splice_plot_genes.csv.gz lists every figure this folder can redraw:',
      '# one row per gene per comparison, with the name the pipeline published',
      '# it under, the number of features tested and the gene-level FDR.',
      'genes <- read.csv("splice_plot_genes.csv.gz", stringsAsFactors = FALSE)',
      '',
      '# GENE takes the gene ID or the symbol, because anyone arriving from a',
      '# figure filename has the symbol to hand. The version suffix is',
      '# optional on either side, since some annotations carry one.',
      '.bare <- function(x) sub("\\\\.[0-9]+$", "", x)',
      'fig <- genes[genes$comparison == COMPARISON &',
      '             (genes$GeneID == GENE |',
      '              .bare(genes$GeneID) == .bare(GENE) |',
      '              (!is.na(genes$gene_name) & genes$gene_name == GENE)), ,',
      '             drop = FALSE]',
      '',
      'if (nrow(fig) != 1) {',
      '    stop("no single figure for GENE = ", GENE,',
      '         " and COMPARISON = ", COMPARISON,',
      '         "\\n  Figures in this folder:\\n",',
      '         paste0("    ", genes$comparison, "  ", genes$GeneID,',
      '                ifelse(is.na(genes$gene_name), "",',
      '                       paste0(" (", genes$gene_name, ")")),',
      '                collapse = "\\n"), call. = FALSE)',
      '}',
      '',
      '# The diffSpliceDGE result for this comparison. It carries every',
      '# feature coefficient in gene order, which is what plotSpliceDGE walks',
      '# and what no published table preserves.',
      'ds <- readRDS(paste0("diffsplice_", COMPARISON, ".rds"))',
      '',
      '# Features at or below FDR_CUT are marked in red, sized by how far',
      '# below. Raising it marks more.',
      'save_and_document(paste0(fig$comparison, "_", fig$figure),',
      '                  function() plotSpliceDGE(ds, geneid = fig$GeneID,',
      '                                           genecolname = "GeneID",',
      '                                           FDR = FDR_CUT),',
      '                  width = WIDTH, height = HEIGHT)'))

  message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
}

message("diffSplice complete.")
