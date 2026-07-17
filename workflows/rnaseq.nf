/*
 * Import modules
 */
include { FASTQC        } from '../modules/local/fastqc'
include { FASTP         } from '../modules/local/fastp'
include { MULTIQC       } from '../modules/local/multiqc'
include { QUARTO_REPORT } from '../modules/local/quarto'
include { STAR_ALIGN    } from '../modules/local/star'
include { HISAT2_ALIGN  } from '../modules/local/hisat2'
include { SALMON_QUANT  } from '../modules/local/salmon'
include { KALLISTO_QUANT} from '../modules/local/kallisto'
include { FEATURECOUNTS } from '../modules/local/featurecounts'
include { DESEQ2        } from '../modules/local/deseq2'
include { EDGER         } from '../modules/local/edger'
include { GSEA          } from '../modules/local/gsea'
include { GPROFILER     } from '../modules/local/gprofiler'
include { GTF2BED       } from '../modules/local/gtf2bed'
include { RSEQC         } from '../modules/local/rseqc'
include { SAMTOOLS_INDEX} from '../modules/local/samtools_index'
include { BAM_TO_BIGWIG } from '../modules/local/bam_to_bigwig'
include { RMATS         } from '../modules/local/rmats'
include { STAR_FUSION   } from '../modules/local/star_fusion'
include { ISOFORM_SWITCH} from '../modules/local/isoform_switch'
include { TXIMPORT      } from '../modules/local/tximport'
include { DTU           } from '../modules/local/dtu'
include { FEATURECOUNTS_EXON } from '../modules/local/featurecounts_exon'
include { DIFFSPLICE    } from '../modules/local/diffsplice'
include { GTF2GENEINFO  } from '../modules/local/gtf2geneinfo'

/*
 * Resolve the DESeq2/edgeR model formula:
 *   - an explicit --design wins (a leading '~' is added if missing);
 *   - else '~ batch + condition' when the samplesheet has a 'batch' column;
 *   - else the plain '~ condition'.
 */
def resolveDesign(samplesheet) {
    if (params.design) {
        def d = params.design.toString().trim()
        return d.startsWith('~') ? d : "~ ${d}"
    }
    def cols = file(samplesheet).splitCsv(header: true).first().keySet()
    return cols.contains('batch') ? '~ batch + condition' : '~ condition'
}

/*
 * Main workflow
 */
workflow RNASEQ {
    take:
    ch_reads // channel: [ val(meta), [ reads ] ]

    main:
    ch_versions = channel.empty()

    // Result directories fed to the Quarto analysis report. Each defaults to an
    // empty list and is reassigned to the real channel if that stage runs, so
    // the report renders whichever sections have data.
    ch_quarto_deseq2    = channel.value([])
    ch_quarto_edger     = channel.value([])
    ch_quarto_gsea      = channel.value([])
    ch_quarto_gprofiler = channel.value([])

    // Count-matrix entry: when --counts is set the pipeline skips QC and
    // alignment entirely and enters at differential expression (Level 3),
    // consuming the pre-computed matrix instead of FASTQs.
    def from_counts = params.counts != null

    // Define execution levels
    def run_level = 100 // Default: run all
    if (params.stop_at == 'preQC')       run_level = 1
    else if (params.stop_at == 'postQC') run_level = 2
    else if (params.stop_at == 'DE')     run_level = 3
    else if (params.stop_at == 'GSEA')   run_level = 4
    
    //
    // MODULE: FastQC (Raw) (Level >= 1)
    //
    if (run_level >= 1 && !from_counts) {
        FASTQC ( ch_reads )
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        
        //
        // MODULE: fastp (Trimming)
        //
        FASTP ( ch_reads )
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
    }

    //
    // MODULE: Alignment / Quantification (Level >= 2)
    //
    ch_align_results = channel.empty()
    ch_bams          = channel.empty()

    // Initialize downstream channels to empty in case skipped.
    // ch_multiqc_files is declared here (not later) so it can be safely
    // referenced from inside conditional blocks below.
    ch_featurecounts_results = channel.empty()  // count files -> DESeq2 / edgeR
    ch_featurecounts_mqc     = channel.empty()  // .summary files -> MultiQC
    ch_featurecounts_exon    = channel.empty()  // exon counts -> diffSplice
    ch_rseqc_results         = channel.empty()
    ch_starfusion_results    = channel.empty()
    ch_multiqc_files         = channel.empty()
    ch_tx_quant              = channel.empty()  // Salmon/Kallisto quant dirs -> tximport

    if (run_level >= 2 && !from_counts) {
        if (params.aligner == 'star') {
            if (!params.star_index) error "STAR index not provided via --star_index"
            STAR_ALIGN( FASTP.out.reads, file(params.star_index), file(params.gtf) ) 
            ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())
            ch_align_results = ch_align_results.mix(STAR_ALIGN.out.log_final.collect { row -> row[1] })
            ch_bams = STAR_ALIGN.out.bam
            
            // STAR-Fusion (runs if Fusion is requested AND we are at Level 4+ OR just run it if CTAT provided? 
            // Usually fusion is advanced analysis. Let's put it in Level 4 as per plan).
            if (run_level >= 4 && params.ctat_lib) {
                STAR_FUSION (
                    STAR_ALIGN.out.chimeric_junction,
                    file(params.ctat_lib)
                )
                ch_versions = ch_versions.mix(STAR_FUSION.out.versions.first())
                ch_starfusion_results = ch_starfusion_results.mix(STAR_FUSION.out.fusions_abridged)
            }
        }
        else if (params.aligner == 'hisat2') {
            if (!params.hisat2_index) error "HISAT2 index not provided via --hisat2_index"
            HISAT2_ALIGN( FASTP.out.reads, file(params.hisat2_index) )
            ch_versions = ch_versions.mix(HISAT2_ALIGN.out.versions.first())
            ch_align_results = ch_align_results.mix(HISAT2_ALIGN.out.summary.collect { row -> row[1] })
            ch_bams = HISAT2_ALIGN.out.bam
        }
        else if (params.aligner == 'salmon') {
            if (!params.salmon_index) error "Salmon index not provided via --salmon_index"
            SALMON_QUANT( FASTP.out.reads, file(params.salmon_index) )
            ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
            ch_align_results = ch_align_results.mix(SALMON_QUANT.out.results.collect { row -> row[1] })
            ch_tx_quant = SALMON_QUANT.out.results.collect { row -> row[1] }
            
            // Isoform Switch (Level 4+)
            if (run_level >= 4 && params.isoform_switch && params.transcript_fasta) {
                ISOFORM_SWITCH (
                    file(params.input),
                    file(params.transcript_fasta),
                    file(params.gtf),
                    SALMON_QUANT.out.results.collect { row -> row[1] },
                    file("${projectDir}/assets/isoform_switch.R")
                )
                ch_versions = ch_versions.mix(ISOFORM_SWITCH.out.versions)
            }
        }
        else if (params.aligner == 'kallisto') {
            if (!params.kallisto_index) error "Kallisto index not provided via --kallisto_index"
            if (params.strandedness == 'auto') {
                log.warn "kallisto cannot infer strandedness (it produces no BAM for RSeQC). " +
                         "Running library-type-agnostic; pass --strandedness forward or reverse " +
                         "explicitly if your library is stranded."
            }
            KALLISTO_QUANT( FASTP.out.reads, file(params.kallisto_index) )
            ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
            ch_align_results = ch_align_results.mix(KALLISTO_QUANT.out.results.collect { row -> row[1] })
            ch_tx_quant = KALLISTO_QUANT.out.results.collect { row -> row[1] }
        }
        
        //
        // Post Alignment QC / Processing (Star/Hisat)
        //
        ch_bam_bai = channel.empty()
        
        if (params.aligner == 'star' || params.aligner == 'hisat2') {
             SAMTOOLS_INDEX ( ch_bams )
             ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
             ch_bam_bai = SAMTOOLS_INDEX.out.bam_bai
             
             // RSeQC (PostQC)
             if (params.gtf) {
                 GTF2BED ( file(params.gtf), file("${projectDir}/assets/gtf2bed12.py") )
                 ch_versions = ch_versions.mix(GTF2BED.out.versions)
                 
                 RSEQC ( ch_bam_bai, GTF2BED.out.bed, file("${projectDir}/assets/parse_strandedness.py") )
                 ch_versions = ch_versions.mix(RSEQC.out.versions.first())
                 ch_rseqc_results = ch_rseqc_results.mix(RSEQC.out.infer_experiment.collect { row -> row[1] })
                 ch_rseqc_results = ch_rseqc_results.mix(RSEQC.out.read_distribution.collect { row -> row[1] })
                 ch_rseqc_results = ch_rseqc_results.mix(RSEQC.out.genebody_coverage.collect { row -> row[1] })
             }
             
             // BigWig
             BAM_TO_BIGWIG( ch_bam_bai )
             ch_versions = ch_versions.mix(BAM_TO_BIGWIG.out.versions.first())
             
             // rMATS (Level 4+) - collapse all per-sample BAMs into one task.
             if (run_level >= 4 && params.input.endsWith('.csv')) {
                ch_rmats_input = ch_bam_bai
                    .map { meta, bam, bai -> [ bam, bai, meta.single_end ] }
                    // flat:false keeps each [bam, bai, single_end] as a sub-list
                    // instead of flattening every element into one long list.
                    .collect(flat: false)
                    .map { rows ->
                         def bams = rows.collect { row -> row[0] }
                         def bais = rows.collect { row -> row[1] }
                         [ [ id:'all_samples', single_end: rows[0][2] ], bams, bais ]
                    }
                RMATS ( file(params.input), ch_rmats_input, file(params.gtf), file("${projectDir}/assets/run_rmats.py") )
                ch_versions = ch_versions.mix(RMATS.out.versions)
             }
         
             // FeatureCounts gene counts (input to DESeq2 / edgeR at Level >= 3).
             if (params.gtf) {
                 def gtf_file = file(params.gtf)

                 // Per-sample strandedness for featureCounts: the RSeQC
                 // inference when --strandedness is 'auto', otherwise the
                 // user-supplied --strandedness value applied to every sample.
                 ch_strand = (params.strandedness == 'auto')
                     ? RSEQC.out.strandedness.map { meta, txt -> [ meta, txt.text.trim() ] }
                     : ch_bams.map { meta, _bam -> [ meta, params.strandedness ] }

                 ch_fc_input = ch_bams.join(ch_strand)
                     .map { meta, bam, strand -> [ meta, bam, gtf_file, strand ] }

                 FEATURECOUNTS ( ch_fc_input )
                 ch_versions = ch_versions.mix(FEATURECOUNTS.out.versions.first())
                 ch_featurecounts_results = FEATURECOUNTS.out.counts.collect { row -> row[1] }
                 ch_featurecounts_mqc     = FEATURECOUNTS.out.summary.collect { row -> row[1] }

                 // Exon-level counts for the optional edgeR diffSplice test.
                 if (run_level >= 3 && params.diffsplice && params.input.endsWith('.csv')) {
                     FEATURECOUNTS_EXON ( ch_fc_input )
                     ch_versions = ch_versions.mix(FEATURECOUNTS_EXON.out.versions.first())
                     ch_featurecounts_exon = FEATURECOUNTS_EXON.out.counts.collect { row -> row[1] }
                 }
             }
        }
    }
    
    //
    // MODULE: Differential Expression (Level >= 3)
    //
    // STAR/HISAT2 feed DESeq2 + edgeR via the featureCounts gene matrix.
    // Salmon/Kallisto feed them via TXIMPORT, which summarises transcript-level
    // quantification to gene level (this is what enables transcript-based
    // aligners to run gene-level differential expression).
    //
    if (run_level >= 3 && params.input.endsWith('.csv')) {
        ch_de_counts = channel.empty()
        ch_gene_info = channel.empty()
        def run_de   = false

        // Gene-symbol / biotype table parsed once from the GTF; used to
        // annotate the DESeq2 / edgeR result tables with readable symbols.
        if (params.gtf) {
            GTF2GENEINFO ( file(params.gtf), file("${projectDir}/assets/gtf2geneinfo.py") )
            ch_versions  = ch_versions.mix(GTF2GENEINFO.out.versions)
            ch_gene_info = GTF2GENEINFO.out.gene_info
        }

        if (from_counts) {
            // Pre-computed matrix: feed it straight to DESeq2/edgeR (both R
            // scripts read it in --matrix mode). No aligner ran.
            ch_de_counts = channel.fromPath(params.counts)
            run_de = true
        }
        else if (params.aligner == 'star' || params.aligner == 'hisat2') {
            ch_de_counts = ch_featurecounts_results
            run_de = true
            if (params.dtu) {
                log.warn "--dtu (differential transcript usage) needs a pseudo-aligner; " +
                         "it has no effect with --aligner ${params.aligner}."
            }

            // Optional: edgeR diffSplice exon-usage test (genome aligners).
            if (params.diffsplice) {
                if (params.gtf) {
                    DIFFSPLICE (
                        file(params.input),
                        ch_featurecounts_exon,
                        file(params.gtf),
                        'exon',
                        params.aligner,
                        file("${projectDir}/assets/diffsplice_edger.R")
                    )
                    ch_versions = ch_versions.mix(DIFFSPLICE.out.versions)
                } else {
                    log.warn "--diffsplice needs --gtf for exon-level counting; skipping."
                }
            }
        }
        else if (params.aligner == 'salmon' || params.aligner == 'kallisto') {
            if (params.gtf) {
                TXIMPORT (
                    file(params.input),
                    ch_tx_quant,
                    file(params.gtf),
                    params.aligner,
                    file("${projectDir}/assets/tximport.R")
                )
                ch_versions  = ch_versions.mix(TXIMPORT.out.versions)
                ch_de_counts = TXIMPORT.out.txi
                run_de = true

                // Optional: differential transcript usage (DEXSeq).
                if (params.dtu) {
                    DTU (
                        file(params.input),
                        ch_tx_quant,
                        file(params.gtf),
                        params.aligner,
                        file("${projectDir}/assets/dtu_dexseq.R")
                    )
                    ch_versions = ch_versions.mix(DTU.out.versions)
                }

                // Optional: edgeR diffSplice transcript-usage test.
                if (params.diffsplice) {
                    DIFFSPLICE (
                        file(params.input),
                        ch_tx_quant,
                        file(params.gtf),
                        'transcript',
                        params.aligner,
                        file("${projectDir}/assets/diffsplice_edger.R")
                    )
                    ch_versions = ch_versions.mix(DIFFSPLICE.out.versions)
                }
            } else {
                log.warn "Differential expression for --aligner ${params.aligner} needs --gtf " +
                         "(transcript-to-gene map); skipping DESeq2/edgeR."
            }
        }

        if (run_de) {
            def de_design = resolveDesign(params.input)
            log.info "Differential-expression design: ${de_design}"

            DESEQ2 ( file(params.input), de_design, ch_gene_info, ch_de_counts, file("${projectDir}/assets/deseq2.R") )
            ch_versions = ch_versions.mix(DESEQ2.out.versions)
            ch_quarto_deseq2 = DESEQ2.out.results

            EDGER ( file(params.input), de_design, ch_gene_info, ch_de_counts, file("${projectDir}/assets/edger.R") )
            ch_versions = ch_versions.mix(EDGER.out.versions)
            ch_quarto_edger = EDGER.out.results

            //
            // MODULE: Enrichment (Level >= 4)
            //
            if (run_level >= 4) {
                if (params.gmt) {
                     GSEA ( file(params.gmt), DESEQ2.out.results, file("${projectDir}/assets/gsea.R") )
                     ch_versions = ch_versions.mix(GSEA.out.versions)
                     ch_quarto_gsea = GSEA.out.results
                }
                GPROFILER ( params.organism, DESEQ2.out.results, file("${projectDir}/assets/gprofiler.R") )
                ch_versions = ch_versions.mix(GPROFILER.out.versions)
                ch_quarto_gprofiler = GPROFILER.out.results
            }
        }
    }
    else if (run_level >= 3 && !params.input.endsWith('.csv')) {
        log.warn "Differential Expression steps skipped because input is not a CSV samplesheet."
    }

    //
    // MODULE: MultiQC
    //
    // Skipped entirely for a --counts run: there are no QC / alignment logs to
    // aggregate, so `multiqc .` would find no modules and fail. The Quarto
    // report then receives an empty MultiQC directory and renders without its
    // QC section (the .qmd already guards every QC chunk on the data existing).
    //
    // ch_multiqc_files was initialised to channel.empty() at the top of the
    // workflow, so the conditional mixes below are always safe.
    ch_multiqc_data = channel.value([])
    if (!from_counts) {
        if (run_level >= 1) {
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { row -> row[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.html.collect { row -> row[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { row -> row[1] })
            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect { row -> row[1] })
        }

        ch_multiqc_files = ch_multiqc_files.mix(ch_align_results)
        ch_multiqc_files = ch_multiqc_files.mix(ch_featurecounts_mqc)
        ch_multiqc_files = ch_multiqc_files.mix(ch_rseqc_results)
        ch_multiqc_files = ch_multiqc_files.mix(ch_starfusion_results)

        // .collect() so MultiQC runs ONCE over all staged files, not once per file.
        // The custom config (auto-detected by `multiqc .`) styles the report; the
        // CSS + logo it references are staged alongside it so MultiQC finds them.
        MULTIQC (
            file("${projectDir}/assets/multiqc_config.yml"),
            [ file("${projectDir}/assets/multiqc_custom.css"),
              file("${projectDir}/assets/multiqc_logo.png") ],
            ch_versions.unique().collectFile(name: 'software_versions.yml'),
            [],
            ch_multiqc_files.collect()
        )
        ch_versions = ch_versions.mix(MULTIQC.out.versions.first())
        ch_multiqc_data = MULTIQC.out.data
    }

    //
    // MODULE: Quarto analysis report
    //
    // Always runs; the DE/enrichment result directories are passed when those
    // stages ran, else an empty list, and the report renders whichever sections
    // have data. multiqc_data is empty for a --counts run (see above).
    QUARTO_REPORT (
        ch_multiqc_data,
        file("${projectDir}/assets/analysis_report.qmd"),
        ch_quarto_deseq2,
        ch_quarto_edger,
        ch_quarto_gsea,
        ch_quarto_gprofiler
    )
    // Single unnamed emit (the strict parser wants the name omitted when there
    // is only one): the accumulated version channel, plus the reporting steps'.
    emit:
    ch_versions.mix(QUARTO_REPORT.out.versions)
}
