process EDGER {
    label 'process_medium'
    // Biocontainer with EdgeR
    // 3.42.0 was never published as a biocontainer; 3.42.4 is the nearest real
    // build of the same minor series. Deliberately not 4.x: edgeR 4 changed the
    // glmQLFit default (legacy=FALSE), which would silently shift every result.
    container 'quay.io/biocontainers/bioconductor-edger:4.8.2--r45h262fe30_1'

    input:
    path samplesheet
    val  design
    path gene_info
    path counts_files
    path r_script
    val  ref_level
    val  padj_cutoff   // what counts as significant; see nextflow.config
    val  lfc_cutoff

    output:
    path "edger_output", emit: results
    path "versions.yml", emit: versions

    script:
    // task.ext.args carries '--matrix' when the run started from a pre-computed
    // count matrix (--counts), and is empty otherwise. It is appended with its
    // own trailing space (not a separate ${} slot) so that when the flag is
    // absent the rendered command is byte-identical to the pre-feature version.
    // That keeps the task hash stable, so a normal run still resumes (-resume)
    // across this change. edger.R strips the flag out wherever it appears.
    def args = task.ext.args ? "${task.ext.args} " : ''
    // Baseline exported from a process input so it reaches the task hash. See
    // modules/local/deseq2.nf and the note in nextflow.config.
    """
    export RNASEQ_FLOW_REFERENCE_LEVEL='${ref_level}'
    export RNASEQ_FLOW_PADJ_CUTOFF='${padj_cutoff}'
    export RNASEQ_FLOW_LFC_CUTOFF='${lfc_cutoff}'
    Rscript ${r_script} ${samplesheet} '${design}' ${gene_info} ${args}${counts_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        edger: \$(Rscript -e "library(edgeR); cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir edger_output
    touch versions.yml
    """
}
