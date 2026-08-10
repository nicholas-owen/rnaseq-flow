#!/usr/bin/env Rscript

# Download MSigDB gene sets as GMT files for a given organism, using msigdbr.
#
# Usage: Rscript download_gmt.R <organism> <outdir>
#   <organism> may be a gProfiler-style code (hsapiens, mmusculus, ...), a
#   common name (human, mouse), or a scientific name ("Homo sapiens").
#
# Requires msigdbr >= 10. The API changed at 10.0.0:
#   msigdbr(category = "H")   ->  msigdbr(collection = "H")
#   df$gs_subcat              ->  df$gs_subcollection
# (the subcollection *values* are unchanged, e.g. "GO:BP"). The `category` /
# `subcategory` arguments still exist but are deprecated and warn.
#
# NOTE ON NETWORK ACCESS. Up to msigdbr 7.5.1 the gene sets were bundled in the
# package. From msigdbr 24.1.0 they are retrieved over the network on first use
# instead, so this step needs outbound HTTPS from wherever the process runs.
# That is consistent with the rest of the --download_refs workflow (which
# fetches from Ensembl), but on an HPC cluster whose compute nodes have no
# internet route it will fail -- run the download workflow somewhere that does,
# and pass the resulting GMT with --gmt on later runs.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: download_gmt.R <organism> <outdir>")
}

organism <- args[1]
out_dir  <- args[2]

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# msigdbr is a CRAN package and is provided by the process environment (see
# modules/local/download_gmt.nf). It is deliberately NOT installed at runtime:
# a container that installs its own dependencies is neither reproducible nor
# guaranteed to have network or a writable library path. (The previous version
# tried BiocManager::install("msigdbr") here, which could not have worked in any
# case -- msigdbr is not a Bioconductor package.)
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  stop("The 'msigdbr' package is not available in this environment. It is ",
       "provided by the DOWNLOAD_GMT process definition; if you are running ",
       "this script by hand, install it with install.packages('msigdbr').")
}
library(msigdbr)

msigdbr_version <- as.character(utils::packageVersion("msigdbr"))
message("msigdbr version: ", msigdbr_version)
if (utils::compareVersion(msigdbr_version, "10.0.0") < 0) {
  stop("msigdbr >= 10 is required (this script uses the 'collection' argument ",
       "and the 'gs_subcollection' column, introduced at 10.0.0). Found ",
       msigdbr_version, ".")
}

# ---- resolve the organism to a scientific name ----------------------------
# msigdbr_species() lists species_name (scientific) and species_common_name.
all_species <- msigdbr::msigdbr_species()

# Resolve whatever the user passed to msigdbr's canonical scientific name.
# Matching is case-insensitive and covers three input styles: a gProfiler-style
# code (hsapiens), a common name (human / zebrafish), or the scientific name
# itself ("Danio rerio"). Comparisons are done on a lower-cased key, but the
# value returned is always msigdbr's own capitalisation -- returning the
# lower-cased input instead is a bug that made every directly-supplied
# scientific name fail validation ("Danio rerio" -> "danio rerio" -> not found).
resolve_species <- function(input, species_tbl) {
  key <- tolower(trimws(input))

  aliases <- c(
    hsapiens = "Homo sapiens",         human     = "Homo sapiens",
    mmusculus = "Mus musculus",        mouse     = "Mus musculus",
    rnorvegicus = "Rattus norvegicus", rat       = "Rattus norvegicus",
    drerio = "Danio rerio",            zebrafish = "Danio rerio",
    celegans = "Caenorhabditis elegans",
    dmelanogaster = "Drosophila melanogaster",
    scerevisiae = "Saccharomyces cerevisiae"
  )
  if (key %in% names(aliases)) return(unname(aliases[[key]]))

  # Otherwise match against msigdbr's own scientific / common name columns.
  i <- match(key, tolower(species_tbl$species_name))
  if (!is.na(i)) return(species_tbl$species_name[i])
  j <- match(key, tolower(species_tbl$species_common_name))
  if (!is.na(j)) return(species_tbl$species_name[j])

  trimws(input)   # unrecognised; rejected by the check below
}

organism_sci <- resolve_species(organism, all_species)
message("Querying msigdbr for: ", organism, " -> ", organism_sci)

if (!(organism_sci %in% all_species$species_name)) {
  stop("Organism '", organism, "' is not available in msigdbr. Available: ",
       paste(all_species$species_name, collapse = ", "))
}

# Gene sets are taken from the human MSigDB and mapped to the requested species
# by orthology (msigdbr's default, db_species = "HS"). This preserves the
# behaviour of earlier versions of this pipeline. MSigDB also publishes native
# mouse collections (db_species = "MM", collections MH/M1/M2/...), which would
# be a different -- and for mouse, arguably better -- set of gene sets; that is
# left as a deliberate future choice rather than a silent change.

# ---- GMT writer -----------------------------------------------------------
write_gmt <- function(gene_sets, filename) {
  # GMT: <name> <tab> <description> <tab> <gene> <tab> <gene> ...
  gs_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
  conn <- file(filename, "w")
  on.exit(close(conn), add = TRUE)
  for (gs_name in names(gs_list)) {
    genes <- unique(gs_list[[gs_name]])
    cat(paste(c(gs_name, "na", genes), collapse = "\t"), "\n",
        file = conn, sep = "")
  }
  message("Saved: ", filename, " (", length(gs_list), " gene sets)")
  invisible(length(gs_list))
}

# Track what actually got written: a run that fetches nothing must not report
# success. Each collection is still attempted independently so that one
# unavailable collection does not lose the others.
written <- character(0)
failures <- character(0)

fetch_collection <- function(label, collection) {
  tryCatch(
    {
      df <- msigdbr(species = organism_sci, collection = collection)
      if (nrow(df) == 0) {
        message("No ", label, " gene sets returned for ", organism_sci, ".")
        return(NULL)
      }
      df
    },
    error = function(e) {
      failures <<- c(failures, sprintf("%s (%s)", label, conditionMessage(e)))
      message("Error fetching ", label, " sets: ", conditionMessage(e))
      NULL
    }
  )
}

# 1. Hallmark (H)
h_df <- fetch_collection("Hallmark", "H")
if (!is.null(h_df)) {
  write_gmt(h_df, file.path(out_dir, "hallmark.gmt"))
  written <- c(written, "hallmark.gmt")
}

# 2. C2 curated (KEGG, Reactome, ...), plus KEGG/Reactome convenience subsets
c2_df <- fetch_collection("C2", "C2")
if (!is.null(c2_df)) {
  write_gmt(c2_df, file.path(out_dir, "c2_curated.gmt"))
  written <- c(written, "c2_curated.gmt")

  kegg_df <- c2_df[grep("^KEGG", c2_df$gs_name), ]
  if (nrow(kegg_df) > 0) {
    write_gmt(kegg_df, file.path(out_dir, "c2_kegg.gmt"))
    written <- c(written, "c2_kegg.gmt")
  }

  reactome_df <- c2_df[grep("^REACTOME_", c2_df$gs_name), ]
  if (nrow(reactome_df) > 0) {
    write_gmt(reactome_df, file.path(out_dir, "c2_reactome.gmt"))
    written <- c(written, "c2_reactome.gmt")
  }
}

# 3. C5 ontology (GO), plus the GO:BP subset most often used for enrichment
c5_df <- fetch_collection("C5", "C5")
if (!is.null(c5_df)) {
  write_gmt(c5_df, file.path(out_dir, "c5_go.gmt"))
  written <- c(written, "c5_go.gmt")

  # msigdbr >= 10 renamed the column gs_subcat -> gs_subcollection; the values
  # ("GO:BP", "GO:CC", "GO:MF") are unchanged.
  bp_df <- c5_df[c5_df$gs_subcollection == "GO:BP", ]
  if (nrow(bp_df) > 0) {
    write_gmt(bp_df, file.path(out_dir, "c5_go_bp.gmt"))
    written <- c(written, "c5_go_bp.gmt")
  }
}

# ---- report ---------------------------------------------------------------
if (length(written) == 0) {
  stop("No gene sets could be downloaded for '", organism_sci, "'. ",
       if (length(failures)) paste0("Failures: ", paste(failures, collapse = "; "), ". ") else "",
       "msigdbr >= 24 retrieves gene sets over the network, so check that this ",
       "machine has outbound HTTPS access.")
}

if (length(failures) > 0) {
  message("Completed with ", length(failures), " collection(s) unavailable: ",
          paste(failures, collapse = "; "))
}
message("Download complete: ", length(written), " GMT file(s) -> ",
        paste(written, collapse = ", "))
