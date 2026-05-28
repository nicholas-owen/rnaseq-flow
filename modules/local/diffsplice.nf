process DIFFSPLICE {
    label 'process_high'
    // edgeR (+ tximport for the transcript level), provisioned on the fly by
    // Nextflow Wave from the Conda packages below - no fixed image required.
    conda 'bioconda::bioconductor-edger bioconda::bioconductor-tximport'

    input:
    path samplesheet
    path count_input   // exon featureCounts files OR per-sample quant directories
    path gtf           // transcript-to-gene map (used by the transcript level)
    val  level         // 'exon' or 'transcript'
    val  aligner       // star | hisat2 | salmon | kallisto
    path r_script      // assets/diffsplice_edger.R

    output:
    path "diffsplice_output", emit: results
    path "versions.yml"     , emit: versions

    script:
    """
    Rscript ${r_script} ${samplesheet} ${level} ${gtf} ${aligner} . ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-edger: \$(Rscript -e "cat(as.character(packageVersion('edgeR')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir diffsplice_output
    touch diffsplice_output/diffsplice_${level}_results.csv
    touch diffsplice_output/diffsplice_gene_results.csv
    touch versions.yml
    """
}
