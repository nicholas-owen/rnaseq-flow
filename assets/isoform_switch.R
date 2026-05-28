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
library(dplyr)

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

# Simplify design for IsoformSwitchAnalyzeR (needs: sampleID, condition, path? No, path is used in import)
# importIsoformExpression wrapper
quant <- importIsoformExpression(
    sampleVector = samples$path,
    addIsolevelInfo = TRUE
)

# 2. Create SwitchAnalyzeRlist
# Modify column names of design to match required: sampleID, condition
design <- samples[, c("sampleID", "condition")]

aSwitchList <- importRdata(
    isoformCountMatrix   = quant$counts,
    isoformRepExpression = quant$abundance,
    designMatrix         = design,
    isoformExonAnnoation = gtf_file,
    isoformNtFasta       = transcript_fasta,
    showProgress = FALSE,
    fixStringTieAnnotationProblem = FALSE # Assuming standard GTF
)

# 3. Filtering
aSwitchList <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = 1, # CPM 1
    isoformExpressionCutoff = 0,
    removeSingleIsoformGenes = TRUE,
    keepIsoformInAllConditions = FALSE
)

# 4. Test for Switches
# Pairwise comparisons
aSwitchList <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchList,
    reduceToSwitchingGenes = TRUE
)

# 5. Extract Results
# Extract Top Switching Genes
extractSwitchSummary(aSwitchList)

switches <- extractTopSwitches(
    aSwitchList, 
    filterForConsequences = FALSE, 
    n = 50, 
    sortByQvals = TRUE
)

dir.create("isoform_switch_output", showWarnings = FALSE)
write.csv(switches, file.path("isoform_switch_output", "isoform_switches.csv"))

# 6. Plots
# Generate plots for top 10 genes
# SwitchPlot requires open reading frames or external analysis for full effect (consequences),
# but usage and splicing plots work without.
top_genes <- head(switches$gene_id, 10)

for (gene in top_genes) {
    tryCatch({
        switchPlot(
            aSwitchList,
            gene = gene,
            condition1 = as.character(design$condition[1]), # Use first two conditions implicitly? 
            # Ideally loop through contrasts. DEXSeq output has condition_1 and condition_2
            # Let's just create generic plots for the top hits which usually picks the relevant contrast.
            filename = file.path("isoform_switch_output", paste0("switch_plot_", gene))
        )
    }, error = function(e) { message(paste("Error plotting gene:", gene)) })
}

# Save Object
saveRDS(aSwitchList, file.path("isoform_switch_output", "switchList_analyzed.rds"))
