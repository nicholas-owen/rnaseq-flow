#!/usr/bin/env Rscript

# Script to download gene sets (GMT) for a specified organism using msigdbr.
# Usage: Rscript download_gmt.R <organism> <outdir>

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: download_gmt.R <organism> <outdir>")
}

organism     <- args[1] # e.g. "Homo sapiens", "Mus musculus"
out_dir      <- args[2]

# Ensure output directory exists
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Install/Load msigdbr
if (!require("msigdbr", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos="http://cran.us.r-project.org")
    BiocManager::install("msigdbr", update = FALSE, ask = FALSE)
    library(msigdbr)
}
library(dplyr)

message(paste("Querying msigdbr for organism:", organism))

# Simple mapping for common gProfiler codes to msigdbr scientific names
get_scientific_name <- function(code) {
  code <- tolower(code)
  if (code == "hsapiens" || code == "human") return("Homo sapiens")
  if (code == "mmusculus" || code == "mouse") return("Mus musculus")
  if (code == "rnorvegicus" || code == "rat") return("Rattus norvegicus")
  if (code == "drerio" || code == "zebrafish") return("Danio rerio")
  if (code == "celegans") return("Caenorhabditis elegans")
  if (code == "dmelanogaster") return("Drosophila melanogaster")
  if (code == "scerevisiae") return("Saccharomyces cerevisiae")
  return(code) # Return valid scientific name if passed directly
}

organism_sci <- get_scientific_name(organism)
message(paste("Normalized organism name:", organism_sci))

# Check if organism is available
all_species <- msigdbr::msigdbr_species()
if (!(organism_sci %in% all_species$species_name)) {
    stop(paste("Organism '", organism_sci, "' not found in msigdbr database."))
}

# Update organism variable for subsequent calls
organism <- organism_sci

# Function to write GMT
write_gmt <- function(gene_sets, filename) {
    # gene_sets is a data frame with gs_name, gene_symbol
    # GMT format: Name\tDescription\tGene1\tGene2...
    
    # We split by gene set name
    gs_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
    
    conn <- file(filename, "w")
    for (gs_name in names(gs_list)) {
        genes <- gs_list[[gs_name]]
        # Description can be NA or URL
        cat(paste(c(gs_name, "na", genes), collapse = "\t"), file = conn)
        cat("\n", file = conn)
    }
    close(conn)
    message(paste("Saved:", filename))
}

# 1. Download Hallmark (H)
tryCatch({
    h_df <- msigdbr(species = organism, category = "H")
    if (nrow(h_df) > 0) {
        write_gmt(h_df, file.path(out_dir, "hallmark.gmt"))
    } else {
        message("No Hallmark gene sets found.")
    }
}, error = function(e) { message("Error fetching Hallmark sets: ", e$message) })

# 2. Download C2 (Curated - KEGG, Reactome, etc.)
tryCatch({
    c2_df <- msigdbr(species = organism, category = "C2")
    # Sub-split into KEGG, REACTOME if desired, but C2 all is common. 
    # Let's save C2 complete.
    if (nrow(c2_df) > 0) {
        write_gmt(c2_df, file.path(out_dir, "c2_curated.gmt"))
        
        # Also try to extract just KEGG and Reactome for convenience
        kegg_df <- c2_df[grep("^KEGG_", c2_df$gs_name), ]
        if (nrow(kegg_df) > 0) write_gmt(kegg_df, file.path(out_dir, "c2_kegg.gmt"))
        
        reactome_df <- c2_df[grep("^REACTOME_", c2_df$gs_name), ]
        if (nrow(reactome_df) > 0) write_gmt(reactome_df, file.path(out_dir, "c2_reactome.gmt"))
    }
}, error = function(e) { message("Error fetching C2 sets: ", e$message) })

# 3. Download C5 (Ontology - GO)
tryCatch({
    c5_df <- msigdbr(species = organism, category = "C5")
    if (nrow(c5_df) > 0) {
        write_gmt(c5_df, file.path(out_dir, "c5_go.gmt"))
        
        # Split GO BP/MF/CC
        bp_df <- c5_df[c5_df$gs_subcat == "GO:BP", ]
        if (nrow(bp_df) > 0) write_gmt(bp_df, file.path(out_dir, "c5_go_bp.gmt"))
    }
}, error = function(e) { message("Error fetching C5 sets: ", e$message) })

message("Download complete.")
