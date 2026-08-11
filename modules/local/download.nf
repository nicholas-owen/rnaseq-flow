process DOWNLOAD_REFS {
    label 'process_low'
    // Use a container with python requests. 
    // quay.io/biocontainers/python-requests is not standard.
    // quay.io/biocontainers/multiqc has python/requests?
    // quay.io/biocontainers/genomepy includes python/requests.
    container 'quay.io/biocontainers/genomepy:0.16.1--pyh7cba7a3_0'

    input:
    val species
    val source
    val release   // Ensembl release to pin (e.g. 102), or 'current'
    path script

    output:
    // download_refs.py exits non-zero unless all three reference files were
    // fetched and verified, so genome/gtf/transcriptome are not optional.
    // The genome glob anchors on '.dna.' so the soft- and repeat-masked
    // variants (dna_sm / dna_rm) could never be confused for it.
    path "references/*"                , emit: files
    path "references/*.dna.*.fa.gz"    , emit: fasta
    path "references/*.cdna.all.fa.gz" , emit: transcript_fasta
    path "references/*.gtf.gz"         , emit: gtf
    path "references/download_log.txt" , emit: log
    // Machine-readable provenance for the set: source, resolved release,
    // assembly, annotation version, and each file's URL, size and SHA-256.
    // Written last by the script, so its presence means everything downloaded
    // and verified.
    path "references/reference_metadata.json", emit: metadata
    path "versions.yml"                , emit: versions

    script:
    """
    python3 ${script} "${species}" "${source}" references "${release}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    // Names mirror the real Ensembl layout so the output globs above are
    // actually exercised by -stub-run.
    """
    mkdir references
    touch references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    touch references/Homo_sapiens.GRCh38.115.gtf.gz
    touch references/Homo_sapiens.GRCh38.cdna.all.fa.gz
    touch references/download_log.txt
    touch references/reference_metadata.json
    touch versions.yml
    """
}
