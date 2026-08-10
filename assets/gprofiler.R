#!/usr/bin/env Rscript

# Usage: Rscript gprofiler.R <organism> <deseq2_output_dir>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: gprofiler.R <organism> <deseq2_output_dir>")
}

organism <- args[1]
de_dir <- args[2]

library(gprofiler2)
library(ggplot2)

# Find DE Result files
res_files <- list.files(de_dir, pattern = "deseq2_results_.*\\.csv", full.names = TRUE)

dir.create("gprofiler_output", showWarnings = FALSE)

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
# reproduce/ : everything an end user needs to redraw any figure in this
# directory without rerunning the pipeline -- the R objects the plot calls
# consume, the tables behind them, and one script per figure type. Published
# automatically, because the module emits gprofiler_output as a whole directory.
#
# The guiding split is that the parent folder is for *reading* and reproduce/ is
# for *running*. Parent CSVs therefore stay plain, so they open by double-click
# in Excel; the reproduce/ copies are gzipped, which read.csv() handles
# transparently and which matters once a results table is tens of MB. Anyone
# wanting to eyeball that data uses the plain copy one directory up.
# ---------------------------------------------------------------------------
REPRO_DIR <- file.path("gprofiler_output", "reproduce")
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

# Contrast / direction combinations that actually produced a figure, collected
# as the loop runs. A contrast can yield an UP result and no DOWN one (or
# neither), so these have to reflect what was really written rather than every
# combination that was attempted.
repro_contrasts  <- character(0)
repro_directions <- character(0)

for (f in res_files) {
    # Extract contrast name
    basename_f <- basename(f)
    contrast_name <- sub("deseq2_results_", "", sub(".csv", "", basename_f))
    
    message(paste("Processing contrast:", contrast_name))
    
    # Read DE results
    res <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
    
    # Identify significant genes
    # Filter: padj < 0.05 and |log2FC| > 1 (optional, stricter is usually better for GO)
    if ("padj" %in% colnames(res) & "log2FoldChange" %in% colnames(res)) {
       sig_genes_up <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 1)]
       sig_genes_down <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < -1)]
    } else {
       next
    }
    
    # Statistical background for the over-representation test: the genes that
    # were actually tested for this contrast (everything DESeq2 kept after
    # filtering), not every gene g:Profiler has an annotation for. Testing a
    # tissue's expressed genes against the whole annotated genome inflates
    # enrichment p-values -- most severely for tissue-specific categories,
    # which are exactly the ones being looked for.
    bg <- rownames(res)

    # Run GOST for Up and Down separately or together
    run_gost <- function(genes, label) {
        if (length(genes) > 0) {
            gostres <- gost(query = genes, 
                            organism = organism, 
                            ordered_query = FALSE, 
                            multi_query = FALSE, 
                            significant = TRUE, 
                            exclude_iea = FALSE, 
                            measure_underrepresentation = FALSE, 
                            evcodes = FALSE, 
                            user_threshold = 0.05, 
                            correction_method = "g_SCS", 
                            domain_scope = "custom_annotated",
                            custom_bg = bg,
                            numeric_ns = "", 
                            sources = NULL, 
                            as_short_link = FALSE)
            
            if (!is.null(gostres) && !is.null(gostres$result)) {
                stem <- paste0("gprofiler_", label, "_", contrast_name)

                # Save Data.
                #
                # gost() returns `parents` as a list-column, which write.csv()
                # cannot represent, so it has to be flattened. Do NOT use
                # apply(result, 2, as.character) for this: apply() returns a
                # matrix for a multi-row result but drops to a *named vector*
                # when there is exactly one row, and write.csv() then writes a
                # single column named "x" holding one value per field. A
                # contrast with a single enriched term therefore produced a
                # silently unreadable table -- real data, structurally corrupt.
                #
                # Collapsing only the list-columns keeps this a data frame
                # whatever the row count, and preserves the numeric columns as
                # numbers rather than round-tripping everything through
                # character. '/' matches the delimiter gsea.R uses for its
                # leadingEdge column.
                res_flat  <- gostres$result
                list_cols <- vapply(res_flat, is.list, logical(1))
                res_flat[list_cols] <- lapply(res_flat[list_cols], function(col)
                    vapply(col, function(x) paste(as.character(x), collapse = "/"),
                           character(1)))

                # row.names = FALSE: the row names here are just 1..N, which
                # write.csv() would emit as a leading column with a blank
                # header. The term identifier is already a named column.
                write.csv(res_flat,
                          file = file.path("gprofiler_output", paste0(stem, ".csv")),
                          row.names = FALSE)

                # Plot, via the shared save_figure so the PNG and SVG match in
                # layout and this module behaves like the others. (It previously
                # used ggsave for the raster copy, which defaults to 300 dpi;
                # the figure is now written at FIG_RES like every other one.)
                save_figure("gprofiler_output",
                            paste0("gostplot_", label, "_", contrast_name),
                            function() print(gostplot(gostres, capped = TRUE,
                                                      interactive = FALSE)),
                            width = 10, height = 6)

                # reproduce/: gostplot() needs the gost result object -- it reads
                # $result and $meta, and cannot be driven from a flat CSV -- so
                # the object is the reproduction input here, with the table
                # alongside it as the durable, inspectable copy.
                saveRDS(gostres, file.path(REPRO_DIR, paste0(stem, ".rds")))
                write_repro_csv(res_flat, stem)
                repro_contrasts  <<- union(repro_contrasts, contrast_name)
                repro_directions <<- union(repro_directions, label)
            }
        }
    }
    
    run_gost(sig_genes_up, "UP")
    run_gost(sig_genes_down, "DOWN")
}

# ---------------------------------------------------------------------------
# reproduce/ script. This figure is indexed by two things -- the contrast and
# the direction -- so it takes two selector constants rather than one. Not every
# combination exists (a contrast may have significant genes in one direction
# only), which is why the values come from what was actually written above and
# why the input check names the file it could not find.
# ---------------------------------------------------------------------------
if (length(repro_contrasts) && length(repro_directions)) {
    write_repro_script(
        file        = file.path(REPRO_DIR, "gostplot.R"),
        description = "gProfiler Manhattan plot of enriched terms.",
        libraries   = "gprofiler2",
        selectors   = list(CONTRAST  = repro_contrasts,
                           DIRECTION = repro_directions),
        constants   = list(CAPPED = TRUE, WIDTH = 10, HEIGHT = 6, FIG_RES = FIG_RES),
        functions   = list(save_figure = save_figure),
        inputs      = 'paste0("gprofiler_", DIRECTION, "_", CONTRAST, ".rds")',
        body        = c(
            '# gostplot() draws from the gost result object, not from the flat',
            '# table -- the matching .csv.gz beside this file is the same data in',
            '# readable form, for anything other than redrawing this figure.',
            'gostres <- readRDS(paste0("gprofiler_", DIRECTION, "_", CONTRAST, ".rds"))',
            '',
            'save_and_document(paste0("gostplot_", DIRECTION, "_", CONTRAST),',
            '                  function() print(gostplot(gostres, capped = CAPPED,',
            '                                            interactive = FALSE)),',
            '                  width = WIDTH, height = HEIGHT)'))

    message("reproduce/ written: ", paste(list.files(REPRO_DIR), collapse = ", "))
}
