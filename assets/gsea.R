#!/usr/bin/env Rscript

# Usage: Rscript gsea.R <gmt_file> <deseq2_output_dir>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: gsea.R <gmt_file> <deseq2_output_dir>")
}

gmt_file <- args[1]
de_dir <- args[2]

library(fgsea)
library(ggplot2)
library(dplyr)
library(tibble)

# Build a named, de-duplicated, descending ranking vector from a value vector and
# a matching identifier vector. Genes with a missing/empty id or value are
# dropped; when several genes share an id (e.g. two Ensembl IDs collapsing to one
# gene symbol) the one with the largest absolute statistic is kept.
build_ranks <- function(values, ids) {
    ids  <- as.character(ids)
    keep <- !is.na(values) & !is.na(ids) & nzchar(ids)
    values <- values[keep]; ids <- ids[keep]
    if (!length(values)) return(numeric(0))
    ord    <- order(abs(values), decreasing = TRUE)   # keep the most extreme per id
    values <- values[ord]; ids <- ids[ord]
    keep2  <- !duplicated(ids)
    values <- values[keep2]
    names(values) <- ids[keep2]
    sort(values, decreasing = TRUE)
}

# Write one figure as both PNG and SVG. `draw` renders the plot and is called
# once per device. SVG uses grDevices::svg() rather than ggsave(..., ".svg"),
# because ggsave() delegates SVG to the svglite package, which the pipeline
# containers do not carry, whereas svg() needs only cairo, which they do. The
# report embeds the SVG; the PNG is kept for slides, email and documents.
# Mirrors save_figure() in deseq2.R / edger.R.
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

# 1. Load Pathways
pathways <- fgsea::gmtPathways(gmt_file)
# Union of all genes across the gene sets, used to pick the ranking identifier
# (gene symbol vs Ensembl ID) that actually overlaps the pathways.
pathway_genes <- unique(unlist(pathways, use.names = FALSE))

# 2. Find DE Result files
res_files <- list.files(de_dir, pattern = "deseq2_results_.*\\.csv", full.names = TRUE)

dir.create("gsea_output", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R objects the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits gsea_output as a whole directory.
#
# The guiding split is that the parent folder is for *reading* and reproduce/ is
# for *running*. Parent CSVs therefore stay plain, so they open by double-click
# in Excel; the reproduce/ copies are gzipped, which read.csv() handles
# transparently and which matters once a results table is tens of MB. Anyone
# wanting to eyeball that data uses the plain copy one directory up.
# ---------------------------------------------------------------------------
REPRO_DIR <- file.path("gsea_output", "reproduce")
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

# Contrasts that actually produced a figure, collected as the loop runs. A
# contrast can be skipped entirely (no gene overlap with the gene sets), so this
# has to reflect what was written rather than every contrast attempted.
repro_contrasts <- character(0)

for (f in res_files) {
    # Extract contrast name from filename
    # filename format: deseq2_results_COND_A_vs_COND_B.csv
    basename_f <- basename(f)
    contrast_name <- sub("deseq2_results_", "", sub(".csv", "", basename_f))
    
    message(paste("Processing contrast:", contrast_name))
    
    # Read DE results
    res <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
    
    # Per-gene ranking statistic: DESeq2's Wald stat, or a signed -log10 p
    # fallback when a shrunken table lacks it.
    if ("stat" %in% colnames(res)) {
        stat <- res$stat
    } else {
        res$pvalue[is.na(res$pvalue)] <- 1
        stat <- sign(res$log2FoldChange) * -log10(res$pvalue + 1e-300)
    }

    # GMT gene sets (MSigDB, and the pipeline's own download) are keyed by gene
    # SYMBOL, but the DESeq2 rows are keyed by Ensembl gene ID. Ranking by the row
    # names would overlap the pathways at essentially zero genes and silently
    # empty the GSEA result. Rank by whichever identifier actually overlaps the
    # gene sets: prefer the annotated gene_name (symbol) column, fall back to the
    # Ensembl IDs -- so both a symbol GMT and an Ensembl-keyed GMT work.
    cand <- list("gene ID" = build_ranks(stat, rownames(res)))
    if ("gene_name" %in% colnames(res))
        cand[["gene symbol"]] <- build_ranks(stat, res$gene_name)
    overlaps <- vapply(cand, function(r) length(intersect(names(r), pathway_genes)),
                       integer(1))
    best  <- names(cand)[which.max(overlaps)]
    ranks <- cand[[best]]
    message(sprintf("  ranking by %s - %d of %d genes overlap the gene sets",
                    best, max(overlaps), length(ranks)))

    # The ranked list fgsea was run on. gsea_stats_*.csv carries the per-pathway
    # results, but the enrichment figure also needs the ranking itself -- and
    # which identifier it is keyed by (symbol or gene ID) is decided above at
    # run time, so it cannot be reconstructed from the DE table alone.
    ranks_out <- data.frame(id = names(ranks), rank_stat = as.numeric(ranks),
                            ranked_by = best, stringsAsFactors = FALSE)
    write.csv(ranks_out,
              file.path("gsea_output", paste0("gsea_ranks_", contrast_name, ".csv")),
              row.names = FALSE)

    # Zero overlap would make fgsea error (or return nothing); skip with a
    # diagnostic instead of failing the run or emitting a silent empty table.
    if (length(ranks) < 1 || max(overlaps) < 1) {
        warning(sprintf(paste0("contrast '%s': no ranked genes overlap the gene sets - skipping ",
            "GSEA. The --gmt is expected to use gene symbols; ranking by symbol needs a --gtf so ",
            "the DE tables carry a gene_name column, and the --gmt organism must match the data."),
            contrast_name))
        next
    }
    if (max(overlaps) < 15)
        warning(sprintf(paste0("contrast '%s': only %d genes overlap the gene sets (< minSize 15); ",
            "GSEA results will be sparse or empty. Check that the --gmt identifiers match the data."),
            contrast_name, max(overlaps)))

    # Run FGSEA.
    # fgsea() dispatches to the multilevel algorithm, which estimates its
    # p-values by Monte Carlo sampling. Without a fixed seed the same ranking
    # yields different p-values -- and so a different pathway order -- on every
    # run, including on -resume. Seeding per contrast (rather than once before
    # the loop) keeps each contrast reproducible on its own, whatever else the
    # run contains and in whatever order the contrasts happen to be processed.
    set.seed(42)
    fgseaRes <- fgsea(pathways = pathways,
                      stats    = ranks,
                      minSize  = 15,
                      maxSize  = 500)
    
    # Filter and Sort
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(padj)
    
    # Save Results.
    # leadingEdge is a list-column of the genes driving each enrichment -- the
    # part of a GSEA result that is actually followed up -- so it is collapsed
    # to a '/'-delimited string rather than dropped. (The previous
    # -which(... %in% ...) form also had a trap: were the column ever absent,
    # -integer(0) selects *zero* columns and writes an empty table.)
    fgseaResTidy$leadingEdge <- vapply(fgseaResTidy$leadingEdge, paste,
                                       character(1), collapse = "/")
    # row.names = FALSE: the row names here are just 1..N, which write.csv()
    # would emit as a leading column with a blank header. The pathway identifier
    # is already a named column.
    write.csv(fgseaResTidy,
              file = file.path("gsea_output", paste0("gsea_stats_", contrast_name, ".csv")),
              row.names = FALSE)
    
    # Plot top pathways
    topPathwaysUp <- fgseaResTidy %>% filter(ES > 0) %>% head(10) %>% pull(pathway)
    topPathwaysDown <- fgseaResTidy %>% filter(ES < 0) %>% head(10) %>% pull(pathway)
    topPathways <- c(topPathwaysUp, topPathwaysDown)
    
    if (length(topPathways) > 0) {
        save_figure("gsea_output", paste0("gsea_plot_", contrast_name),
                    function() print(plotGseaTable(pathways[topPathways], ranks,
                                                   fgseaRes, gseaParam = 0.5)),
                    width = 8.3, height = 6.25)

        # reproduce/: this is the one figure in the pipeline whose inputs cannot
        # be rebuilt from the published tables at all. plotGseaTable needs the
        # gene-set *membership* as well as the ranking, and membership lives in
        # the GMT -- an input, never published with the results. So the bundle is
        # saved: the plotted gene sets, the ranking they were scored against, and
        # the fgsea result. Only the plotted subset of pathways is kept, which is
        # why this is named after the figure rather than after the data (the
        # DESeq2 / edgeR objects are whole-analysis objects; this one is not).
        saveRDS(list(pathways = pathways[topPathways],
                     ranks    = ranks,
                     fgseaRes = fgseaRes),
                file.path(REPRO_DIR, paste0("gsea_plot_", contrast_name, ".rds")))

        # The two tables behind the result, as the inspectable copies. Both stay
        # in the parent folder as well: gsea_stats is the headline result, and
        # gsea_ranks is a genuine analysis output (the ranked gene list), not
        # merely figure scaffolding.
        write_repro_csv(fgseaResTidy, paste0("gsea_stats_", contrast_name))
        write_repro_csv(ranks_out,    paste0("gsea_ranks_", contrast_name))
        repro_contrasts <- c(repro_contrasts, contrast_name)
    }
}

# ---------------------------------------------------------------------------
# reproduce/ script, one per figure type, with the contrast in a single constant
# at the top rather than a loop so the code stays flat and literal. The
# width/height match what the pipeline used above, so a regenerated figure
# matches the published one.
# ---------------------------------------------------------------------------
if (length(repro_contrasts)) {
    write_repro_script(
        file        = file.path(REPRO_DIR, "gsea_plot.R"),
        description = "fgsea summary table: enrichment of the top pathways.",
        libraries   = "fgsea",
        selectors   = list(CONTRAST = repro_contrasts),
        constants   = list(GSEA_PARAM = 0.5, WIDTH = 8.3, HEIGHT = 6.25,
                           FIG_RES = FIG_RES),
        functions   = list(save_figure = save_figure),
        inputs      = 'paste0("gsea_plot_", CONTRAST, ".rds")',
        body        = c(
            '# The bundle holds the three things plotGseaTable needs: the plotted',
            '# gene sets, the ranking they were scored against, and the fgsea',
            '# result. Gene-set membership comes from the GMT, which is a pipeline',
            '# input and is not published with the results, so this figure cannot',
            '# be rebuilt from the CSVs alone -- hence the object.',
            'gsea <- readRDS(paste0("gsea_plot_", CONTRAST, ".rds"))',
            '',
            'save_and_document(paste0("gsea_plot_", CONTRAST),',
            '                  function() print(plotGseaTable(gsea$pathways,',
            '                                                 gsea$ranks,',
            '                                                 gsea$fgseaRes,',
            '                                                 gseaParam = GSEA_PARAM)),',
            '                  width = WIDTH, height = HEIGHT)'))

    message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
}
