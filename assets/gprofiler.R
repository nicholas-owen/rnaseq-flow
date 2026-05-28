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
                            domain_scope = "annotated", 
                            custom_bg = NULL, 
                            numeric_ns = "", 
                            sources = NULL, 
                            as_short_link = FALSE)
            
            if (!is.null(gostres) && !is.null(gostres$result)) {
                # Save Data
                write.csv(apply(gostres$result,2,as.character), file = file.path("gprofiler_output", paste0("gprofiler_", label, "_", contrast_name, ".csv")))
                
                # Plot
                p <- gostplot(gostres, capped = TRUE, interactive = FALSE)
                ggsave(file.path("gprofiler_output", paste0("gostplot_", label, "_", contrast_name, ".png")), plot = p, width = 10, height = 6)
            }
        }
    }
    
    run_gost(sig_genes_up, "UP")
    run_gost(sig_genes_down, "DOWN")
}
