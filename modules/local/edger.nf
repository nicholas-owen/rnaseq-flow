process EDGER {
    label 'process_medium'
    // Biocontainer with EdgeR
    // 3.42.0 was never published as a biocontainer; 3.42.4 is the nearest real
    // build of the same minor series. Deliberately not 4.x: edgeR 4 changed the
    // glmQLFit default (legacy=FALSE), which would silently shift every result.
    container 'quay.io/biocontainers/bioconductor-edger:3.42.4--r43hf17093f_0'

    input:
    path samplesheet
    val  design
    path gene_info
    path counts_files
    path r_script

    output:
    path "edger_output", emit: results
    path "versions.yml", emit: versions

    script:
    // task.ext.args carries '--matrix' when the run started from a pre-computed
    // count matrix (--counts), and is empty otherwise. It is appended with its
    // own trailing space (not a separate ${} slot) so that when the flag is
    // absent the rendered command is byte-identical to the pre-feature version
    // — that keeps the task hash stable, so a normal run still resumes (-resume)
    // across this change. edger.R strips the flag out wherever it appears.
    def args = task.ext.args ? "${task.ext.args} " : ''
    """
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
