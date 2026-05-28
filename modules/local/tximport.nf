process TXIMPORT {
    label 'process_medium'
    // tximport summarises Salmon/Kallisto transcript quantification to gene level.
    container 'quay.io/biocontainers/bioconductor-tximport:1.28.0--r43hdfd78af_0'

    input:
    path samplesheet
    path quant_dirs   // per-sample Salmon/Kallisto output directories
    path gtf          // used to build the transcript-to-gene map
    val  aligner      // 'salmon' or 'kallisto'
    path r_script     // assets/tximport.R

    output:
    path "txi.rds"                 , emit: txi
    path "tximport_gene_counts.csv", emit: counts
    path "versions.yml"            , emit: versions

    script:
    """
    Rscript ${r_script} ${samplesheet} ${gtf} ${aligner} .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximport: \$(Rscript -e "library(tximport); cat(as.character(packageVersion('tximport')))")
    END_VERSIONS
    """

    stub:
    """
    touch txi.rds
    touch tximport_gene_counts.csv
    touch versions.yml
    """
}
