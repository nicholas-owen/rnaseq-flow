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

outdir <- "dtu_output"
dir.create(outdir, showWarnings = FALSE)
bpp <- MulticoreParam(workers = cpus)

# ---------------------------------------------------------------------------
# Write one figure as both PNG and SVG.
#
# `draw` is a function of no arguments that renders the plot. It is called once
# per device, because a plot can only be written while its device is open.
#
# SVG is written with grDevices::svg() rather than ggsave(..., ".svg"):
# ggsave() delegates SVG to the svglite package, which none of the pipeline
# containers carry, while svg() needs only cairo, which all of them have. Both
# devices are given the same inch dimensions so the raster and vector copies
# are identical in layout.
#
# Mirrors save_figure() in deseq2.R, edger.R, gsea.R, gprofiler.R,
# isoform_switch.R and diffsplice_edger.R. The copies are character-identical
# and are edited together.
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

# Isoform fractions: each transcript's share of its gene, within one sample.
# Any per-sample scalar cancels, so it does not matter that the counts went
# through dtuScaledTPM, and size factors are not needed either. A sample with
# no reads on the gene gives NA rather than 0/0, and drops out of the figure.
dtu_fractions <- function(counts) {
    total <- colSums(counts)
    sweep(counts, 2, ifelse(total > 0, total, NA_real_), "/")
}

# Isoform usage for one gene: a bar per condition per transcript, at that
# condition's mean isoform fraction, with every sample drawn as a point on top.
#
# This is the quantity DEXSeq tests here. The gene-level q-value says a gene's
# transcript proportions moved; this says which transcripts moved and by how
# much, on the 0 to 1 scale the proportions live on.
#
# DEXSeq's own plotDEXSeq is not used. It is written for exon bins laid out on
# genomic coordinates, and draws a gene model keyed to them; with whole
# transcripts as features those spans overlap into a single block, and the
# connector lines between the panels cross at random because transcripts have
# no left-to-right order. It also plots fitted counts on a log axis rather than
# proportions, so a switch has to be read off two lines instead of an axis.
#
# counts  features x samples, for one gene
# conds   condition of each column, in column order. A factor fixes the order
#         the bars and the legend appear in, which is how the baseline from
#         --reference_level ends up on the left of every group.
# padj    per-feature adjusted p-value, in row order; NA allowed
draw_dtu_usage <- function(counts, conds, padj, title,
                           max_isoforms = 12, alpha = 0.05, cols = NULL) {
    frac <- dtu_fractions(counts)
    levs  <- if (is.factor(conds)) levels(droplevels(conds)) else unique(conds)
    conds <- as.character(conds)
    cmean <- vapply(levs, function(l)
        rowMeans(frac[, conds == l, drop = FALSE], na.rm = TRUE),
        numeric(nrow(frac)))
    cmean <- matrix(cmean, nrow = nrow(frac),
                    dimnames = list(rownames(frac), levs))

    # A gene can carry fifty annotated transcripts and the panel cannot. Show
    # every significant one, fill the rest of the panel with the most used, and
    # sum whatever is left into a single "other" bar so the bars still add to
    # one. Drawn in order of usage, so the dominant isoform is on the left.
    ranked <- order(rowMeans(cmean, na.rm = TRUE), decreasing = TRUE)
    sig    <- ranked[ranked %in% which(!is.na(padj) & padj < alpha)]
    keep   <- head(union(sig, ranked), max_isoforms)
    keep   <- keep[order(match(keep, ranked))]
    rest   <- setdiff(seq_len(nrow(frac)), keep)

    bars <- cmean[keep, , drop = FALSE]
    pts  <- frac[keep, , drop = FALSE]
    labs <- rownames(frac)[keep]
    star <- !is.na(padj[keep]) & padj[keep] < alpha
    if (length(rest)) {
        bars <- rbind(bars, colSums(cmean[rest, , drop = FALSE], na.rm = TRUE))
        pts  <- rbind(pts,  colSums(frac[rest, , drop = FALSE], na.rm = TRUE))
        labs <- c(labs, sprintf("other (%d)", length(rest)))
        star <- c(star, FALSE)
    }

    if (is.null(cols))
        cols <- c("#4C72B0", "#DD8452", "#55A868", "#C44E52",
                  "#8172B3", "#937860", "#DA8BC3", "#8C8C8C")
    if (length(levs) > length(cols)) cols <- grDevices::rainbow(length(levs))
    cols <- cols[seq_along(levs)]

    n <- nrow(bars); k <- length(levs)
    at <- matrix(0, nrow = k, ncol = n)
    for (i in seq_len(n)) at[, i] <- (i - 1) * (k + 1) + seq_len(k)

    # Feature IDs are written vertically, so the bottom margin has to follow
    # the longest one: Ensembl transcript IDs need about nine lines, a short
    # exon label needs four, and a fixed margin wastes a third of the panel on
    # one or clips the other.
    cex_lab <- 0.68
    op <- par(mar = c(min(14, max(4, 0.5 * max(nchar(labs)) * cex_lab + 2.5)),
                      4.5, 4, 2))
    on.exit(par(op))
    plot.new()
    plot.window(xlim = c(0.3, max(at) + 1),
                ylim = c(0, max(bars, pts, na.rm = TRUE) * 1.15))
    abline(h = seq(0, 1, 0.1), col = "#EEEEEE")
    rect(at - 0.42, 0, at + 0.42, t(bars), col = cols, border = NA)

    for (i in seq_len(n)) for (j in seq_len(k)) {
        v <- pts[i, conds == levs[j]]
        v <- v[!is.na(v)]
        if (!length(v)) next
        jit <- if (length(v) == 1) 0 else seq(-0.12, 0.12, length.out = length(v))
        points(at[j, i] + jit, v, pch = 16, cex = 0.6, col = "#00000099")
    }

    axis(2, at = seq(0, 1, 0.1), las = 2, cex.axis = 0.8)
    mtext("isoform fraction", side = 2, line = 3, cex = 0.95)
    mid  <- colMeans(at)
    lab2 <- ifelse(star, paste0(labs, " *"), labs)
    # Two passes, because axis() takes one colour for the whole call.
    if (any(star))
        axis(1, at = mid[star], labels = lab2[star], las = 2, cex.axis = cex_lab,
             col.axis = "#C0392B", tick = FALSE, line = -0.4)
    if (any(!star))
        axis(1, at = mid[!star], labels = lab2[!star], las = 2, cex.axis = cex_lab,
             col.axis = "black", tick = FALSE, line = -0.4)
    box(col = "#999999")
    title(main = title, adj = 0, cex.main = 1.1)
    mtext(sprintf("* DEXSeq padj < %s;  points are individual samples", alpha),
          side = 3, adj = 0, line = 0.2, cex = 0.72, col = "#555555")
    legend("topright", legend = levs, fill = cols, border = NA,
           bty = "n", cex = 0.85)
    invisible(NULL)
}

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

# ---- transcripts quantified but absent from the GTF -----------------------
# The transcript FASTA the index was built from and the GTF need not cover the
# same transcripts, and Ensembl's own defaults do not: cdna.all.fa.gz carries
# transcripts on patch and haplotype scaffolds (CHR_HSCHR*), while
# Homo_sapiens.GRCh38.<release>.gtf.gz covers the primary assembly only. Those
# are the two files --download_refs fetches, so the mismatch is the default.
#
# With countsFromAbundance = "dtuScaledTPM" tximport asserts that every
# quantified transcript appears in tx2gene, and the extras abort the import:
#   all(txId %in% tx2gene$tx) is not TRUE
# Map each to itself. It becomes a single-isoform gene, DTU is undefined for
# those, and keep_multi() below drops it before DEXSeq sees it. No real gene's
# scaling changes, because dtuScaledTPM scales within a gene.
#
# tximport strips versions with sub("\\..*", ...), not strip()'s trailing
# .<digits>; match tximport, since it is what performs the check.
quantified <- sub("\\..*", "", read.delim(files[1], stringsAsFactors = FALSE)[[1]])
orphans    <- setdiff(quantified, tx2gene$TXNAME)
if (length(orphans) > 0) {
  message(sprintf(paste("tx2gene: %d of %d quantified transcripts are not in",
                        "the GTF (patch and haplotype scaffolds, typically);",
                        "each is treated as a single-isoform gene and dropped",
                        "before testing."),
                  length(orphans), length(quantified)))
  tx2gene <- rbind(tx2gene, data.frame(TXNAME = orphans, GENEID = orphans,
                                       stringsAsFactors = FALSE))
}

# ---- import transcript-level counts (DTU-scaled) --------------------------
# tx2gene is required here, not optional: dtuScaledTPM scales each transcript
# against the others in its own gene, so tximport has to know the grouping.
# Without it the call fails with "'dtuScaledTPM' requires 'tx2gene' input".
# scaledTPM, which diffsplice_edger.R uses, scales by library size and needs no
# such map, which is why only this script is affected.
# dropInfReps: neither this script nor diffsplice_edger.R reads the inferential
# replicates, only txi$counts, and importing them drags in a format-specific
# dependency that differs per quantifier. Salmon stores them as JSON and needs
# jsonlite; kallisto stores them in abundance.h5 and needs rhdf5. Skipping the
# import avoids both and is faster.
txi <- tximport(files, type = aligner, txOut = TRUE,
                countsFromAbundance = "dtuScaledTPM", tx2gene = tx2gene,
                dropInfReps = TRUE, ignoreTxVersion = TRUE)
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

# ---------------------------------------------------------------------------
# reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R objects the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits dtu_output as a whole directory.
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

# ---- DEXSeq DTU test ------------------------------------------------------
condition <- factor(samples$condition)
# Baseline level from --reference_level (default "REF"); see deseq2.R for why it
# arrives as an environment variable rather than an argument.
REF_LEVEL <- Sys.getenv("RNASEQ_FLOW_REFERENCE_LEVEL", "REF")
if (REF_LEVEL %in% levels(condition)) condition <- relevel(condition, ref = REF_LEVEL)
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

# ---- figures --------------------------------------------------------------
# The genes whose transcript proportions moved most, one figure each. There is
# no contrast dimension here the way there is in rmats_output/ or
# diffsplice_output/: the design is ~ sample + exon + condition:exon, so DEXSeq
# tests the condition:exon interaction once across every condition at once and
# returns one q-value per gene. The figure therefore shows all conditions side
# by side and the files sit flat in dtu_output/.
PLOT_TOP_N       <- 10
PLOT_WIDTH       <- 9
PLOT_HEIGHT      <- 6
PLOT_ALPHA       <- 0.05   # padj below which a transcript is starred
PLOT_MAX_ISOFORM <- 12     # per panel, before the rest are summed into "other"

path_safe <- function(x) gsub("[/:]", "-", x)

top <- head(gqdf[order(gqdf$gene_qvalue), , drop = FALSE], PLOT_TOP_N)
pad <- formatC(seq_len(nrow(top)), width = nchar(nrow(top)), flag = "0")
top$figure <- paste0(pad, "_dtu_plot_", path_safe(top$gene))
# geneinfo falls back to the gene ID when the GTF carries no gene_name, so the
# suffix is added only where it says something the ID does not.
named <- !is.na(top$gene_name) & nzchar(top$gene_name) & top$gene_name != top$gene
top$figure[named] <- paste0(top$figure[named], "_aka_", path_safe(top$gene_name[named]))

# Everything the figures are drawn from, in long form, for the plotted genes
# only. Counts rather than fractions, because the fraction is what the inlined
# helper computes and a reproduce/ script that recomputed nothing would hide
# the one step worth seeing.
plot_rows <- list()

for (i in seq_len(nrow(top))) {
  g        <- top$gene[i]
  features <- rownames(cts)[gene_of == g]
  if (length(features) < 2) next
  # Columns ordered by condition level, so the baseline is the left bar of
  # every group. It also means the reproduce/ table can recover that order
  # from the row order alone, with no extra column to carry it.
  o    <- order(condition)
  m    <- cts[features, o, drop = FALSE]
  # Matched on the gene as well as the feature: a transcript ID is unique, but
  # matching on it alone would silently depend on that staying true.
  padj <- res$padj[match(paste(g, features), paste(res$groupID, res$featureID))]
  lab  <- if (named[i]) paste0(top$gene_name[i], "  (", g, ")") else g

  tryCatch(
    save_figure(outdir, top$figure[i],
                function() draw_dtu_usage(m, condition[o], padj,
                                          title = paste(lab, " isoform usage"),
                                          max_isoforms = PLOT_MAX_ISOFORM,
                                          alpha = PLOT_ALPHA),
                width = PLOT_WIDTH, height = PLOT_HEIGHT),
    error = function(e) {
      # The tables are the result; one unplottable gene is not worth failing
      # the task for. save_figure() opens a device before it draws, so an error
      # part way through leaves that device open and the next gene would be
      # drawn into it. Rscript starts with only the null device, device 1, so
      # anything above that is this loop's leak. png() also creates its file
      # before anything is drawn, so remove the blank image rather than ship it.
      while (dev.cur() > 1) dev.off()
      unlink(file.path(outdir, paste0(top$figure[i], c(".png", ".svg"))))
      message("Could not plot ", g, ": ", conditionMessage(e))
    })

  plot_rows[[g]] <- data.frame(
    gene_id    = g,
    gene_name  = top$gene_name[i],
    feature_id = rep(features, times = ncol(m)),
    sample     = rep(colnames(m), each = length(features)),
    condition  = rep(as.character(condition)[o], each = length(features)),
    count      = as.vector(m),
    padj       = rep(padj, times = ncol(m)),
    stringsAsFactors = FALSE)
}

# Both formats have to be there. The SVG is written second, so a figure with
# only a PNG is one that failed between the two devices.
stems <- file.path(outdir, top$figure)
ok    <- file.exists(paste0(stems, ".png")) & file.exists(paste0(stems, ".svg"))
top   <- top[ok, , drop = FALSE]
plot_rows <- plot_rows[names(plot_rows) %in% top$gene]
message(sprintf("Wrote %d DTU %s.", nrow(top),
                if (nrow(top) == 1) "figure" else "figures"))

# ---- reproduce/ -----------------------------------------------------------
if (nrow(top) > 0 && length(plot_rows)) {
  dir.create(REPRO_DIR, showWarnings = FALSE, recursive = TRUE)

  # No .rds anywhere in this folder. Unlike the switch and splice plots, every
  # number behind a DTU figure is a count and a p-value, so a gzipped table
  # carries the lot: ten genes of a few transcripts across a handful of samples
  # is a few hundred rows.
  write_repro_csv(do.call(rbind, plot_rows), "dtu_plot_data")

  # One row per published figure.
  write_repro_csv(top[, c("gene", "gene_name", "gene_biotype", "gene_qvalue",
                          "figure")],
                  "dtu_plot_genes")

  # The two published tables, gzipped.
  write_repro_csv(res,  "dtu_transcript_results")
  write_repro_csv(gqdf, "dtu_gene_qvalues")

  write_repro_script(
    file        = file.path(REPRO_DIR, "dtu_plot.R"),
    description = paste("DEXSeq DTU plot: each transcript's share of its gene,",
                        "per condition, for one gene."),
    libraries   = character(0),   # base graphics only; DEXSeq is not needed to redraw
    # A GENE selector rather than the CONTRAST the DE and enrichment scripts
    # use, because these figures are per gene. There is no COMPARISON here:
    # DEXSeq tests one interaction across every condition at once.
    selectors   = list(GENE = unique(top$gene)),
    constants   = list(ALPHA = PLOT_ALPHA, MAX_ISOFORMS = PLOT_MAX_ISOFORM,
                       WIDTH = PLOT_WIDTH, HEIGHT = PLOT_HEIGHT,
                       FIG_RES = FIG_RES),
    functions   = list(dtu_fractions  = dtu_fractions,
                       draw_dtu_usage = draw_dtu_usage,
                       save_figure    = save_figure),
    inputs      = c('"dtu_plot_data.csv.gz"', '"dtu_plot_genes.csv.gz"'),
    body        = c(
      '# dtu_plot_genes.csv.gz lists every figure this folder can redraw, with',
      '# the name the pipeline published it under and the gene-level q-value.',
      '# dtu_plot_data.csv.gz holds the counts and adjusted p-values behind',
      '# them, one row per transcript per sample.',
      'genes <- read.csv("dtu_plot_genes.csv.gz", stringsAsFactors = FALSE)',
      'rows  <- read.csv("dtu_plot_data.csv.gz",  stringsAsFactors = FALSE)',
      '',
      '# GENE takes the gene ID or the symbol, because anyone arriving from a',
      '# figure filename has the symbol to hand. The version suffix is',
      '# optional on either side, since some annotations carry one.',
      '.bare <- function(x) sub("\\\\.[0-9]+$", "", x)',
      'fig <- genes[genes$gene == GENE |',
      '             .bare(genes$gene) == .bare(GENE) |',
      '             (!is.na(genes$gene_name) & genes$gene_name == GENE), ,',
      '             drop = FALSE]',
      '',
      'if (nrow(fig) != 1) {',
      '    stop("no single figure for GENE = ", GENE,',
      '         "\\n  Genes in this folder:\\n",',
      '         paste0("    ", genes$gene,',
      '                ifelse(is.na(genes$gene_name), "",',
      '                       paste0(" (", genes$gene_name, ")")),',
      '                collapse = "\\n"), call. = FALSE)',
      '}',
      '',
      '# Back to the features x samples matrix draw_dtu_usage expects. The',
      '# sample order has to match the condition vector, so both come off the',
      '# same rows in the same order.',
      'd    <- rows[rows$gene_id == fig$gene, , drop = FALSE]',
      'feat <- unique(d$feature_id)',
      'samp <- unique(d$sample)',
      'm    <- matrix(d$count[match(paste(rep(feat, times = length(samp)),',
      '                                   rep(samp, each = length(feat))),',
      '                             paste(d$feature_id, d$sample))],',
      '               nrow = length(feat), dimnames = list(feat, samp))',
      '# The rows were written with the baseline condition first, so taking the',
      '# levels in order of appearance puts it on the left of every bar group,',
      '# exactly where the published figure has it.',
      'conds <- d$condition[match(samp, d$sample)]',
      'conds <- factor(conds, levels = unique(conds))',
      'padj  <- d$padj[match(feat, d$feature_id)]',
      '',
      'lab <- if (!is.na(fig$gene_name) && fig$gene_name != fig$gene)',
      '           paste0(fig$gene_name, "  (", fig$gene, ")") else fig$gene',
      '',
      'save_and_document(fig$figure,',
      '                  function() draw_dtu_usage(m, conds, padj,',
      '                                            title = paste(lab, " isoform usage"),',
      '                                            max_isoforms = MAX_ISOFORMS,',
      '                                            alpha = ALPHA),',
      '                  width = WIDTH, height = HEIGHT)'))

  message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
}
