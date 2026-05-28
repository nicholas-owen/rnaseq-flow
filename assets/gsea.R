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

# 1. Load Pathways
pathways <- fgsea::gmtPathways(gmt_file)

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
    
    # Rank genes
    # We use stat column from DESeq2
    # If stat is missing (e.g. shrinking was used without recalculating), use signed -log10 pvalue
    if ("stat" %in% colnames(res)) {
        ranks <- res$stat
    } else {
        # Fallback: sign(log2FoldChange) * -log10(pvalue)
        # Avoid infinite pvalues
        res$pvalue[is.na(res$pvalue)] <- 1
        ranks <- sign(res$log2FoldChange) * -log10(res$pvalue + 1e-300)
    }
    
    names(ranks) <- rownames(res)
    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    
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
