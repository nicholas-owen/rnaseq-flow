process EDGER {
    label 'process_medium'
    // Biocontainer with EdgeR
    container 'quay.io/biocontainers/bioconductor-edger:3.42.0--r43hdfd78af_0'

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
    """
    Rscript ${r_script} ${samplesheet} '${design}' ${gene_info} ${counts_files}

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
