#!/usr/bin/env Rscript

# Usage: Rscript isoform_switch.R <samplesheet.csv> <transcript_fasta> <gtf> <salmon_parent_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: isoform_switch.R <samplesheet.csv> <transcript_fasta> <gtf> <salmon_parent_dir>")
}

samplesheet_file <- args[1]
transcript_fasta <- args[2]
gtf_file <- args[3]
salmon_dir <- args[4]

library(IsoformSwitchAnalyzeR)
# dplyr was attached here and never used. Its forty-odd masking messages were
# the first thing in every failure log, which is what made the real errors in
# this script hard to find.

# theme_bw(), which section 6 hands to switchPlot(). IsoformSwitchAnalyzeR
# Depends on ggplot2, so this attaches nothing that was not already on the
# search path; it names the dependency the code below actually uses.
library(ggplot2)

# ---------------------------------------------------------------------------
# Write one figure as both PNG and SVG.
#
# `draw` is a function of no arguments that renders the plot. It is called once
# per device, because a plot can only be written while its device is open.
# switchPlot() suits this exactly: it takes no filename, pathToOutput or
# fileType argument and instead composes its four panels onto the current
# device with grid viewports, so wrapping the call is all that is needed.
#
# SVG is written with grDevices::svg() rather than ggsave(..., ".svg"):
# ggsave() delegates SVG to the svglite package, which none of the pipeline
# containers carry, while svg() needs only cairo, which all of them have. Both
# devices are given the same inch dimensions so the raster and vector copies
# are identical in layout.
#
# Mirrors save_figure() in deseq2.R, edger.R, gsea.R and gprofiler.R. The four
# copies are character-identical and are edited together.
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

# isoform_switch_output is a declared Nextflow output, so it has to exist
# whatever happens below, and a dataset with no significant switches is a
# legitimate result rather than a failed task. Same pattern as dtu_dexseq.R.
dir.create("isoform_switch_output", showWarnings = FALSE)

emit_empty <- function(msg) {
    message(msg, " Writing empty isoform switch result.")
    write.csv(data.frame(),
              file.path("isoform_switch_output", "isoform_switches.csv"))
    quit(save = "no", status = 0)
}

# preFilter raises an error rather than returning an empty object when nothing
# survives the multi-isoform filter, and importRdata calls it internally at its
# last step to estimate the differential range. So the condition has to be
# caught at both call sites, not checked after them.
#
# An annotation with one transcript per gene reaches this legitimately. Ensembl
# R64-1-1 has 7127 genes and 7127 transcripts, so the yeast example has no gene
# that isoform analysis is defined for, and the correct outcome is an empty
# result rather than a failed task.
#
# Only this one condition is treated as empty; every other error is re-raised,
# because reporting a real failure as "nothing found" would be a wrong answer
# rather than a missing one. The match is deliberately loose: the package's own
# wording is "No genes were left after filtering for mutlipe transcrip genes",
# and those typos may be corrected upstream. If the wording changes beyond
# recognition this stops catching and the task fails loudly, which is the safe
# direction.
empty_if_no_multi_isoform <- function(expr, where) {
    tryCatch(expr, error = function(e) {
        if (grepl("no genes were left", conditionMessage(e), ignore.case = TRUE)) {
            emit_empty(paste0(where, ": no gene has two or more isoforms. ",
                              "On an annotation with one transcript per gene ",
                              "this is expected."))
        }
        stop(e)
    })
}

# 1. Prepare Design
# We need to map samples to their quant files.
# The 'salmon_dir' is expected to contain subdirectories named after samples
# OR we pass the list of files. 
# For simplicity with Nextflow, we might have staged all sample dirs into 'salmon_dir'.

samples <- read.csv(samplesheet_file, stringsAsFactors = FALSE)
# Expect columns: sample, condition
# Construct path to quant.sf
# We assume the directory name matches the sample ID from the samplesheet.
# Update: In Nextflow we will likely stage: salmon_dir/sample_id/quant.sf

samples$path <- file.path(salmon_dir, samples$sample, "quant.sf")
samples$sampleID <- samples$sample
samples$condition <- samples$condition

# Check existence
missing <- !file.exists(samples$path)
if (any(missing)) {
    warning(paste("Missing quant files for:", paste(samples$sample[missing], collapse=", ")))
    samples <- samples[!missing, ]
}

if (nrow(samples) < 2) {
    stop("Not enough samples found.")
}

# The isoform_id column that importRdata needs on both matrices is added by
# addIsofomIdAsColumn (the typo is the package's own), which defaults to TRUE,
# so nothing has to be named here.
#
# This call used to pass addIsolevelInfo = TRUE. No release of
# IsoformSwitchAnalyzeR has ever had an argument by that name, so R raised
# "unused argument" and --isoform_switch failed at its first call on every
# dataset and every aligner. It had been there since the first commit.
#
# ignoreAfterPeriod trims the version suffix off the quantifier's transcript
# IDs. Ensembl cDNA FASTA headers are versioned (ENST00000641871.1) but Ensembl
# GTFs carry the version as a separate transcript_version attribute, leaving
# transcript_id bare, so without this nothing in the count matrix matches the
# annotation importRdata reads below. It defaults to FALSE.
#
# sampleVector is named, because its names become the column names of the
# returned matrices and importRdata matches those against designMatrix$sampleID.
# An unnamed vector falls back to names derived from the paths, which do not
# match the samplesheet, and importRdata stops with "Each sample stored in
# 'designMatrix$sampleID' must have a corresponding expression column". Naming
# the vector also rules out renaming by position, which would quietly swap
# conditions if the function ever returned the columns in another order.
quant <- importIsoformExpression(
    sampleVector      = setNames(samples$path, samples$sampleID),
    ignoreAfterPeriod = TRUE
)

missing_cols <- setdiff(samples$sampleID, colnames(quant$counts))
if (length(missing_cols) > 0) {
    stop("No expression column for: ", paste(missing_cols, collapse = ", "),
         ". importIsoformExpression returned: ",
         paste(colnames(quant$counts), collapse = ", "))
}

# importRdata refuses to run when the quantification and the annotation
# describe different isoform sets, reporting a Jaccard similarity below 0.925.
# Ensembl's defaults guarantee that: cdna.all.fa.gz carries transcripts on
# patch and haplotype scaffolds and Homo_sapiens.GRCh38.<release>.gtf.gz does
# not. On GRCh38.116 that was 14212 of 453553 quantified isoforms.
#
# IsoformSwitchAnalyzeR suggests supplying the chr_patch_hapl_scaff GTF
# instead. That is not available here, because the same --gtf feeds
# featureCounts, DESeq2 annotation and rMATS, where patch and haplotype
# duplicates of real genes would corrupt gene-level counts. Drop the
# unannotated isoforms instead: without exon structure there is nothing to
# analyse, and diffsplice_edger.R and dtu_dexseq.R handle them the same way.
#
# Only the feature-type "transcript" rows are read. Matching transcript_id
# across every exon line too would mean a regex over millions of lines for the
# same answer.
gtf_lines <- readLines(gtf_file)
tx_lines  <- grep("\ttranscript\t", gtf_lines, value = TRUE, fixed = TRUE)
if (length(tx_lines) == 0) {
    # Not every GTF carries transcript rows; fall back to any line with the tag.
    tx_lines <- grep('transcript_id "', gtf_lines, value = TRUE)
}
annotated <- unique(sub('.*transcript_id "([^"]+)".*', "\\1", tx_lines))
rm(gtf_lines, tx_lines)

n_quant <- nrow(quant$counts)
quant$counts    <- quant$counts[quant$counts$isoform_id %in% annotated, , drop = FALSE]
quant$abundance <- quant$abundance[quant$abundance$isoform_id %in% annotated, , drop = FALSE]
if (nrow(quant$counts) < n_quant) {
    message(sprintf(paste("Dropped %d of %d quantified isoforms with no entry",
                          "in the GTF (patch and haplotype scaffolds)."),
                    n_quant - nrow(quant$counts), n_quant))
}
if (nrow(quant$counts) == 0) {
    emit_empty("No quantified isoform appears in the GTF.")
}

# 2. Create SwitchAnalyzeRlist
# Modify column names of design to match required: sampleID, condition
design <- samples[, c("sampleID", "condition")]

# Which condition becomes condition_1 sets the sign of every result, because
# dIF is IF2 - IF1 (verified against a GSE151251 run). Left to itself,
# importRdata builds all unique non-directional pairwise comparisons and
# documents nothing about how it picks the order, so the direction is whatever
# the package happens to do.
#
# On the examples to hand it lands on the pipeline's reference level first,
# since REF precedes TGFb both alphabetically and in samplesheet order. That is
# luck. Rename the conditions, or point --reference_level at the other one, and
# this analysis would invert while DESeq2, edgeR, rMATS, DTU and diffSplice all
# kept the requested direction: those five read RNASEQ_FLOW_REFERENCE_LEVEL and
# this script did not.
#
# State the comparisons instead, with the reference as condition_1 wherever it
# appears. dIF then reads as change relative to the baseline, matching the sign
# of DESeq2's log2 fold change. All pairwise comparisons are still produced, so
# this fixes the direction without narrowing the analysis.
REF_LEVEL <- Sys.getenv("RNASEQ_FLOW_REFERENCE_LEVEL", "REF")
conds <- unique(as.character(design$condition))
comparisons <- NULL
if (!REF_LEVEL %in% conds) {
    message("Reference level '", REF_LEVEL, "' is not a condition in the ",
            "samplesheet. Leaving the comparison direction to ",
            "IsoformSwitchAnalyzeR.")
} else {
    others <- sort(setdiff(conds, REF_LEVEL))
    comparisons <- data.frame(condition_1 = REF_LEVEL, condition_2 = others,
                              stringsAsFactors = FALSE)
    # Pairs not involving the reference still have to be listed or they are
    # dropped. Neither side is a baseline, so the sign carries no special
    # meaning, but the pair still has to be oriented one way and section 6
    # names the figure folder after it. run_rmats.py sorts the conditions and
    # calls the pair <first>_vs_<second>, with the first as numerator. dIF is
    # IF2 - IF1, so putting the second alphabetically in condition_1 makes the
    # first the numerator here too, and the folder name the same string in both
    # directories. A three-condition run then has rmats_output/A_vs_B beside
    # isoform_switch_output/A_vs_B, measuring the same direction.
    if (length(others) > 1) {
        rest <- t(utils::combn(others, 2))
        comparisons <- rbind(comparisons,
                             data.frame(condition_1 = rest[, 2],
                                        condition_2 = rest[, 1],
                                        stringsAsFactors = FALSE))
    }
}

aSwitchList <- empty_if_no_multi_isoform(importRdata(
    isoformCountMatrix   = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix         = design,
    comparisonsToMake    = comparisons,
    isoformExonAnnoation = gtf_file,
    isoformNtFasta       = transcript_fasta,
    showProgress = FALSE,
    # Needed again here, and for a different input: the count matrix arrives
    # already trimmed from importIsoformExpression, and the GTF's transcript_id
    # carries no version, but isoformNtFasta is the cDNA FASTA, whose headers
    # are versioned. Without this the sequences match no isoform.
    ignoreAfterPeriod = TRUE,
    fixStringTieAnnotationProblem = FALSE # Ensembl GTF, not a StringTie assembly
), "importRdata")

# 3. Filtering
aSwitchList <- empty_if_no_multi_isoform(preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 1, # CPM 1
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE,
    keepIsoformInAllConditions = FALSE
), "preFilter")

# removeSingleIsoformGenes removes everything on an annotation where no gene has
# two isoforms, which is exactly the yeast example: Ensembl R64-1-1 has 7127
# genes and 7127 transcripts. That is a correct outcome, not a failure, and the
# guard has to sit here rather than only after the test, because nothing
# downstream of an empty list would get that far.
if (is.null(aSwitchList) || nrow(aSwitchList$isoformFeatures) == 0) {
    emit_empty(paste("No multi-isoform genes survived filtering. On an",
                     "annotation with one transcript per gene this is",
                     "expected."))
}

# 4. Test for Switches
# Pairwise comparisons
aSwitchList <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchList,
    reduceToSwitchingGenes = TRUE
)

# reduceToSwitchingGenes = TRUE discards everything that did not switch, so no
# surviving features is a possible outcome. Errors are deliberately not caught
# here: reporting a genuine failure as "no switches found" would be a wrong
# answer rather than a missing one.
if (is.null(aSwitchList) || nrow(aSwitchList$isoformFeatures) == 0) {
    emit_empty("No isoform switches passed testing.")
}

# 5. Extract Results
# Extract Top Switching Genes
extractSwitchSummary(aSwitchList)

# extractTopSwitches' own defaults, named here because section 6 selects the
# genes to plot with a second call and has to use the same two numbers. A
# figure for a gene absent from isoform_switches.csv would have no explanation.
SWITCH_ALPHA <- 0.05
SWITCH_DIF   <- 0.1

# extractTopSwitches() stops rather than returning an empty frame when no gene
# clears the cutoffs. isoformSwitchTestDEXSeq above filters on the same two
# numbers, so reaching here with nothing left should not happen, but "no
# switches" is a legitimate result and must not surface as a failed task.
switches <- tryCatch(
    extractTopSwitches(
        aSwitchList,
        filterForConsequences = FALSE,
        n = 50,
        alpha = SWITCH_ALPHA,
        dIFcutoff = SWITCH_DIF,
        sortByQvals = TRUE
    ),
    error = function(e) NULL
)

if (is.null(switches) || nrow(switches) == 0) {
    emit_empty(sprintf("No gene cleared q < %s and |dIF| > %s.",
                       SWITCH_ALPHA, SWITCH_DIF))
}

# row.names = FALSE: extractTopSwitches() already returns gene_id / isoform_id
# as named columns, so the row names are just 1..N and write.csv() would emit
# them as a leading column with a blank header.
write.csv(switches, file.path("isoform_switch_output", "isoform_switches.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Figures
#
# switchPlotTopSwitches() used to write these. switchPlot() takes no filename,
# pathToOutput or fileType argument in any release, and this script called it
# in a loop with filename = <path>, so every call raised "unused argument" and
# no plot was ever written, in the same way addIsolevelInfo stopped the import.
# switchPlotTopSwitches() does write files, and fixed that, but it writes PNG
# only, at its own size, into a folder it names condition_1 first. With the
# reference as condition_1 that reads REF_vs_TGFb: the convention this pipeline
# dropped after v1.5.1 and renamed out of rmats_output/ (OUTPUTS.md section 11).
#
# The sign was never wrong. dIF is IF2 - IF1, so it already measured change in
# the test condition, the same direction as DESeq2's log2FoldChange. That is
# the dangerous case, because the folder name said the opposite of what the
# numbers in it did, and nothing in the files resolved the disagreement.
#
# So the loop below places the files itself. switchPlot() draws to the current
# device, so save_figure() wraps it: the figures gain an SVG at the shared
# FIG_RES like every other figure the pipeline publishes, the folder is named
# <condition>_vs_REF like rmats_output/ and the DE tables, and section 8 can
# ship a reproduce/ script that redraws any one of them.
#
# File names are the ones switchPlotTopSwitches used, so a result from before
# this change differs only in which folder it sits in. The rank prefix is
# zero-padded, which makes the folder sort by q-value.
#
# 9 x 6 inches is the size the package chose for this plot, and it needs the
# room: four panels sharing one transcript axis.
# ---------------------------------------------------------------------------
PLOT_TOP_N     <- 10    # genes per comparison
PLOT_WIDTH     <- 9
PLOT_HEIGHT    <- 6
PLOT_IF_CUTOFF <- 0.05  # switchPlot's own default: isoform fraction below which an isoform is drawn as not expressed
PLOT_BASE_SIZE <- 8     # ggplot2 base font size, the one switchPlotTopSwitches passed

top <- tryCatch(
    extractTopSwitches(
        aSwitchList,
        filterForConsequences = FALSE,
        n = PLOT_TOP_N,
        alpha = SWITCH_ALPHA,
        dIFcutoff = SWITCH_DIF,
        inEachComparison = TRUE,
        sortByQvals = TRUE
    ),
    error = function(e) NULL
)

if (is.null(top) || nrow(top) == 0) {
    message("No gene passed the plotting cutoffs. No switch figures written.")
    top <- NULL
} else {
    # dIF is IF2 - IF1, so condition_2 is the numerator and goes first.
    top$comparison <- paste0(top$condition_2, "_vs_", top$condition_1)

    # extractTopSwitches numbers its Rank column 1..N across every comparison
    # at once, so on a three-condition design the second comparison would start
    # at 11. Rank within the comparison instead, which is what a file name
    # sitting in one comparison's folder means.
    top <- top[order(top$comparison, top$gene_switch_q_value), , drop = FALSE]
    top$rank <- as.integer(stats::ave(seq_len(nrow(top)), top$comparison,
                                      FUN = seq_along))

    # Zero-padded only in the file name. Keeping rank itself an integer means
    # the reproduce/ index reads back as a number rather than as "01", which
    # read.csv would turn into 1 anyway and so disagree with the file beside it.
    pad <- formatC(top$rank, width = max(nchar(top$rank)), flag = "0")

    # IsoformSwitchAnalyzeR keeps the Ensembl gene version on gene_id, so a
    # gene arrives here as ENSG00000168710.20. diffsplice_edger.R and
    # dtu_dexseq.R read gene_id straight out of the GTF, where it is bare, so
    # without this the same gene is named two ways across three result
    # directories and a cross-reference by ID quietly matches nothing. AHCYL1,
    # MYLK and NR1H2 all appeared in both folders on the GSE151251 run.
    #
    # Only the published name changes. gene_id itself stays as the object has
    # it, because that is what switchPlot() matches on and what the subset
    # below selects by.
    bare <- function(x) sub("\\.[0-9]+$", "", x)
    top$gene_id_bare <- bare(top$gene_id)

    # Gene symbols carry '/' and ':' often enough to matter, and neither is
    # legal in a path everywhere. The package stripped them the same way.
    path_safe <- function(x) gsub("[/:]", "-", x)
    top$figure <- paste0(pad, "_switch_plot_", path_safe(top$gene_id_bare))
    named <- !is.na(top$gene_name) & nzchar(top$gene_name)
    top$figure[named] <- paste0(top$figure[named], "_aka_",
                                path_safe(top$gene_name[named]))

    for (i in seq_len(nrow(top))) {
        fig     <- top[i, ]
        cmp_dir <- file.path("isoform_switch_output", fig$comparison)
        dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)
        tryCatch(
            save_figure(
                cmp_dir, fig$figure,
                function() switchPlot(
                    aSwitchList,
                    gene       = fig$gene_id,
                    condition1 = fig$condition_1,
                    condition2 = fig$condition_2,
                    IFcutoff   = PLOT_IF_CUTOFF,
                    dIFcutoff  = SWITCH_DIF,
                    # Topology needs analyzeDeepTMHMM(), which this script does
                    # not run. Left at its TRUE default switchPlot prints
                    # "Omitting toplogy visualization" (the typo is the
                    # package's) once per device, so twice per figure and
                    # twenty times over a ten-gene run, for a panel that was
                    # never going to be drawn.
                    plotTopology = FALSE,
                    localTheme = theme_bw(base_size = PLOT_BASE_SIZE)
                ),
                width = PLOT_WIDTH, height = PLOT_HEIGHT),
            error = function(e) {
                # The tables above are the result; one unplottable gene is not
                # worth failing the task for. save_figure() opens a device
                # before it draws, though, so an error inside switchPlot()
                # leaves that device open and the next gene would be drawn into
                # it. Rscript starts with only the null device, which is
                # device 1, so anything above that is this loop's leak.
                while (dev.cur() > 1) dev.off()
                # png() creates its file before draw() is called, so a failure
                # part way through leaves a blank image on disk. Remove it. A
                # figure that is present but empty is worse than one that is
                # missing, and it would otherwise be counted below as written.
                unlink(file.path(cmp_dir,
                                 paste0(fig$figure, c(".png", ".svg"))))
                message("Could not plot ", fig$gene_id, ": ",
                        conditionMessage(e))
            }
        )
    }

    # What reached disk, which is what section 8 may offer to redraw. Both
    # formats have to be there: the SVG is written second, so a figure with
    # only a PNG is one that failed between the two devices.
    drawn <- file.path("isoform_switch_output", top$comparison, top$figure)
    top <- top[file.exists(paste0(drawn, ".png")) &
               file.exists(paste0(drawn, ".svg")), , drop = FALSE]
    message("Wrote ", nrow(top), " switch ",
            if (nrow(top) == 1) "figure" else "figures", ".")
}

# ---------------------------------------------------------------------------
# 7. The analysis object
#
# The whole switchAnalyzeRlist, for anyone taking this further: consequence
# prediction, ORF analysis, or just a different set of genes plotted. Section 8
# saves a trimmed copy for the figures; this is the complete one.
# ---------------------------------------------------------------------------
saveRDS(aSwitchList, file.path("isoform_switch_output", "switchList_analyzed.rds"))

# ---------------------------------------------------------------------------
# 8. reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R object the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits isoform_switch_output as a whole
# directory.
#
# The guiding split is that the parent folder is for *reading* and reproduce/ is
# for *running*. Parent CSVs therefore stay plain, so they open by double-click
# in Excel; the reproduce/ copies are gzipped, which read.csv() handles
# transparently and which matters once a results table is tens of MB. Anyone
# wanting to eyeball that data uses the plain copy one directory up.
#
# Unlike deseq2_output/ and the rest, the folder is created only once a figure
# exists. A run that finds no switch has nothing to redraw, and an empty
# reproduce/ next to an empty result would only raise the question.
# ---------------------------------------------------------------------------
REPRO_DIR <- file.path("isoform_switch_output", "reproduce")

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

# switchPlot() wants the whole switchAnalyzeRlist, not a table: exon structures,
# per-sample expression and the test result all feed one figure. The object
# saved in section 7 carries the entire analysis, which on a human annotation
# runs to hundreds of MB, so trim it to the plotted genes. Same reasoning as
# edger.R keeping its MDS object rather than the DGEList, and gsea.R keeping the
# plotted pathways rather than the whole GMT: the smallest thing that redraws.
plot_list <- if (is.null(top) || nrow(top) == 0) NULL else tryCatch(
    subsetSwitchAnalyzeRlist(
        aSwitchList,
        aSwitchList$isoformFeatures$gene_id %in% unique(top$gene_id)),
    error = function(e) {
        message("Could not trim the switchAnalyzeRlist, so no reproduce/ ",
                "folder was written: ", conditionMessage(e))
        NULL
    })

if (!is.null(plot_list)) {
    # Trimmed to the plotted genes is not the same as small. On a GRCh38 run
    # with ten genes the subset still came to 79 MB, about sixty times every
    # figure in the directory combined, and it is published and hashed into the
    # RO-Crate on every run.
    #
    # So drop whatever the figure cannot read. switchPlot() and the panels it
    # composes reach for nine slots and no others: isoformFeatures, exons,
    # conditions, sourceId, and the five optional annotation slots this script
    # never computes. Checked by taking every switchAnalyzeRlist$<slot> access
    # in isoform_plots.R of 2.10.0, the pinned version: ntSequence, aaSequence,
    # isoformCountMatrix, isoformRepExpression, isoformRepIF,
    # isoformSwitchAnalysis and switchConsequence appear nowhere in that file.
    #
    # designMatrix and runInfo are a few hundred bytes and are kept anyway:
    # they say which samples and which cutoffs produced the figure, and
    # summary() reads them.
    #
    # Nothing is lost. The whole object is still published next door as
    # switchList_analyzed.rds; this copy exists only to redraw these figures.
    PLOT_SLOTS <- c("isoformFeatures", "exons", "conditions", "sourceId",
                    "designMatrix", "runInfo",
                    "orfAnalysis", "domainAnalysis", "idrAnalysis",
                    "signalPeptideAnalysis", "topologyAnalysis",
                    "AlternativeSplicingAnalysis", "subCellLocationAnalysis")
    # Computed before the loop, because the names change as slots are removed.
    # [[<- NULL drops a slot and keeps the S4 wrapper; subsetting with [ would
    # return a plain list and switchPlot() rejects anything whose class is not
    # exactly 'switchAnalyzeRlist'.
    for (s in setdiff(names(plot_list), PLOT_SLOTS)) plot_list[[s]] <- NULL

    dir.create(REPRO_DIR, showWarnings = FALSE, recursive = TRUE)
    rds <- file.path(REPRO_DIR, "switch_plots.rds")
    saveRDS(plot_list, rds)
    # Reported so a run says for itself whether this stayed small, rather than
    # anyone having to go and measure it.
    message(sprintf("reproduce/switch_plots.rds: %.1f MB, %d slots (%s)",
                    file.size(rds) / 1024^2, length(names(plot_list)),
                    paste(names(plot_list), collapse = ", ")))

    # One row per published figure: the gene, the comparison, and the name it
    # was written under. switch_plot.R reads this rather than parsing file
    # names, and it doubles as the index of what is in the comparison folders.
    # gene_id is the bare form, so this table joins to dtu_output/ and
    # diffsplice_output/ on it. gene_id_versioned is what the saved object and
    # isoform_switches.csv carry, and what switchPlot() has to be given.
    write_repro_csv(
        data.frame(comparison          = top$comparison,
                   rank                = top$rank,
                   gene_id             = top$gene_id_bare,
                   gene_id_versioned   = top$gene_id,
                   gene_name           = top$gene_name,
                   condition_1         = top$condition_1,
                   condition_2         = top$condition_2,
                   gene_switch_q_value = top$gene_switch_q_value,
                   figure              = top$figure,
                   stringsAsFactors    = FALSE),
        "switch_plot_genes")

    # The published table, gzipped.
    write_repro_csv(switches, "isoform_switches")

    write_repro_script(
        file        = file.path(REPRO_DIR, "switch_plot.R"),
        description = paste("Isoform switch plot: transcript structures, gene",
                            "expression, isoform expression and isoform usage",
                            "for one gene."),
        libraries   = c("IsoformSwitchAnalyzeR", "ggplot2"),
        # A GENE selector rather than the CONTRAST one the other reproduce/
        # scripts carry, because these figures are per gene. COMPARISON comes
        # second and only does any work on a design with more than two
        # conditions, where one gene can switch in several of them.
        selectors   = list(GENE       = unique(top$gene_id_bare),
                           COMPARISON = unique(top$comparison)),
        constants   = list(IF_CUTOFF = PLOT_IF_CUTOFF, DIF_CUTOFF = SWITCH_DIF,
                           BASE_SIZE = PLOT_BASE_SIZE, WIDTH = PLOT_WIDTH,
                           HEIGHT = PLOT_HEIGHT, FIG_RES = FIG_RES),
        functions   = list(save_figure = save_figure),
        inputs      = c('"switch_plots.rds"', '"switch_plot_genes.csv.gz"'),
        body        = c(
            '# switch_plot_genes.csv.gz lists every figure this folder can',
            '# redraw: one row per gene per comparison, with the name the',
            '# pipeline published it under and the two conditions switchPlot',
            '# needs. gene_switch_q_value is the value the ranking used.',
            'genes <- read.csv("switch_plot_genes.csv.gz", stringsAsFactors = FALSE)',
            '',
            '# GENE takes the gene ID or the gene symbol, because anyone',
            '# arriving from a figure filename has the symbol to hand. The',
            '# version suffix is optional on either side: the figures and the',
            '# gene_id column are bare, so they line up with dtu_output/ and',
            '# diffsplice_output/, while the object carries the versioned form.',
            '.bare <- function(x) sub("\\\\.[0-9]+$", "", x)',
            'fig <- genes[genes$comparison == COMPARISON &',
            '             (genes$gene_id == GENE |',
            '              genes$gene_id_versioned == GENE |',
            '              .bare(genes$gene_id) == .bare(GENE) |',
            '              (!is.na(genes$gene_name) & genes$gene_name == GENE)), ,',
            '             drop = FALSE]',
            '',
            'if (nrow(fig) != 1) {',
            '    stop("no single figure for GENE = ", GENE,',
            '         " and COMPARISON = ", COMPARISON,',
            '         "\\n  Figures in this folder:\\n",',
            '         paste0("    ", genes$comparison, "  ", genes$gene_id,',
            '                ifelse(is.na(genes$gene_name), "",',
            '                       paste0(" (", genes$gene_name, ")")),',
            '                collapse = "\\n"), call. = FALSE)',
            '}',
            '',
            '# The plotted genes only, cut down to the slots switchPlot reads:',
            '# the feature table, the exon structures and the conditions. The',
            '# sequence slots, the replicate matrices and the per-isoform test',
            '# detail are not in here, because the figure never looks at them.',
            '# For any of that, or to run consequence or ORF analysis, use',
            '# switchList_analyzed.rds one directory up, which is the whole',
            '# object.',
            'sal <- readRDS("switch_plots.rds")',
            '',
            '# dIF is IF2 - IF1, so condition2 is the one the figure measures',
            '# the change in. The published figure sits in',
            '# <condition2>_vs_<condition1>/ for the same reason.',
            '# gene_id_versioned, not gene_id: the object keeps the Ensembl',
            '# version and switchPlot matches on it exactly.',
            'save_and_document(paste0(fig$comparison, "_", fig$figure),',
            '                  function() switchPlot(sal,',
            '                                        gene       = fig$gene_id_versioned,',
            '                                        condition1 = fig$condition_1,',
            '                                        condition2 = fig$condition_2,',
            '                                        IFcutoff   = IF_CUTOFF,',
            '                                        dIFcutoff  = DIF_CUTOFF,',
            '                                        # needs analyzeDeepTMHMM(), which the',
            '                                        # pipeline does not run',
            '                                        plotTopology = FALSE,',
            '                                        localTheme = theme_bw(base_size = BASE_SIZE)),',
            '                  width = WIDTH, height = HEIGHT)'))

    message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
}
