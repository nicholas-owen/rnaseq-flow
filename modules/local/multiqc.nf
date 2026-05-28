/*
 * MultiQC: Aggregate results from analysis
 */
process MULTIQC {
    label 'process_single'
    container 'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0'

    input:
    path multiqc_config
    path extra_multiqc_config
    path software_versions
    path workflow_summary
    path multiqc_files

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    multiqc -f . $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    touch versions.yml
    """
}
