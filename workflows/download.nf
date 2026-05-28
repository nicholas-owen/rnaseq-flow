include { DOWNLOAD_REFS } from '../modules/local/download'
include { DOWNLOAD_GMT  } from '../modules/local/download_gmt'

workflow DOWNLOAD {
    take:
    species
    source

    main:
    DOWNLOAD_REFS (
        species,
        source,
        params.download_release ?: 'current',
        file("${projectDir}/assets/download_refs.py")
    )

    if (params.download_gmt) {
        DOWNLOAD_GMT (
            params.organism, // Use the unified 'organism' param
            file("${projectDir}/assets/download_gmt.R")
        )
    }
}
