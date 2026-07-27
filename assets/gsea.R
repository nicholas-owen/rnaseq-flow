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

# 1. Load Pathways
pathways <- fgsea::gmtPathways(gmt_file)
# Union of all genes across the gene sets, used to pick the ranking identifier
# (gene symbol vs Ensembl ID) that actually overlaps the pathways.
pathway_genes <- unique(unlist(pathways, use.names = FALSE))

# 2. Find DE Result files
res_files <- list.files(de_dir, pattern = "deseq2_results_.*\\.csv", full.names = TRUE)

dir.create("gsea_output", showWarnings = FALSE)

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

    # Run FGSEA
    fgseaRes <- fgsea(pathways = pathways, 
                      stats    = ranks,
                      minSize  = 15,
                      maxSize  = 500)
    
    # Filter and Sort
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(padj)
    
    # Save Results
    write.csv(fgseaResTidy[, -which(names(fgseaResTidy) %in% c("leadingEdge"))], 
              file = file.path("gsea_output", paste0("gsea_stats_", contrast_name, ".csv")))
    
    # Plot top pathways
    topPathwaysUp <- fgseaResTidy %>% filter(ES > 0) %>% head(10) %>% pull(pathway)
    topPathwaysDown <- fgseaResTidy %>% filter(ES < 0) %>% head(10) %>% pull(pathway)
    topPathways <- c(topPathwaysUp, topPathwaysDown)
    
    if (length(topPathways) > 0) {
        png(file.path("gsea_output", paste0("gsea_plot_", contrast_name, ".png")), width=800, height=600)
        p <- plotGseaTable(pathways[topPathways], ranks, fgseaRes, gseaParam = 0.5)
        print(p)
        dev.off()
    }
}
