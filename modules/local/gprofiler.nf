process GPROFILER {
    label 'process_medium'
    // Container with gprofiler2. rocker/verse usually does NOT have it by default, 
    // but we can install it or use a bioconductor container if available. 
    // 'gprofiler2' is on CRAN. 
    // Let's use a container that likely has it or install it on fly (not recommended for NF)
    // or use a safer bet: custom docker or a container from biocontainers that includes it.
    // 'quay.io/biocontainers/r-gprofiler2:0.2.1--r41h7d875b9_0' exists.
    container 'quay.io/biocontainers/r-gprofiler2:0.2.3--r43hc247a5b_0'

    input:
    val organism
    path deseq2_dir
    path r_script

    output:
    path "gprofiler_output", emit: results
    path "versions.yml"    , emit: versions

    script:
    """
    Rscript ${r_script} ${organism} ${deseq2_dir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gprofiler2: \$(Rscript -e "library(gprofiler2); cat(as.character(packageVersion('gprofiler2')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir gprofiler_output
    touch versions.yml
    """
}
