/*
 * fastp: Ultra-fast all-in-one FASTQ preprocessor (QC + Trimming)
 */
process FASTP {
    tag "$meta.id"
    label 'process_medium'
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastp.fastq.gz") , emit: reads
    tuple val(meta), path("*.json")           , emit: json
    tuple val(meta), path("*.html")           , emit: html
    path "versions.yml"                       , emit: versions

    script:
    def args = task.ext.args ?: ''
    /*
     * fastp supports both single-end and paired-end.
     * We determine the input type by checking the number of read files.
     *
     * The reads keep their own filenames. MultiQC's fastp parser takes the
     * sample name from the --in1 path recorded in the JSON's command string
     * rather than from the .json filename, so this reports under whatever the
     * samplesheet pointed at -- typically a run accession. That is intentional:
     * it preserves the link back to the source data. MULTIQC is given a
     * sample-name mapping built from the samplesheet so the report still shows
     * one row per sample; see modules/local/multiqc.nf.
     */
    if (meta.single_end) {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --out1 ${meta.id}.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        """
        fastp \\
            --in1 ${reads[0]} \\
            --in2 ${reads[1]} \\
            --out1 ${meta.id}_1.fastp.fastq.gz \\
            --out2 ${meta.id}_2.fastp.fastq.gz \\
            --thread $task.cpus \\
            --json ${meta.id}.fastp.json \\
            --html ${meta.id}.fastp.html \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    }

    stub:
    if (meta.single_end) {
        """
        touch ${meta.id}.fastp.fastq.gz
        touch ${meta.id}.fastp.json
        touch ${meta.id}.fastp.html
        touch versions.yml
        """
    } else {
        """
        touch ${meta.id}_1.fastp.fastq.gz
        touch ${meta.id}_2.fastp.fastq.gz
        touch ${meta.id}.fastp.json
        touch ${meta.id}.fastp.html
        touch versions.yml
        """
    }
}
