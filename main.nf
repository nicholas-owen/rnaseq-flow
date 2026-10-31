#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Import the main workflow
 */
include { RNASEQ; resolveContributors; resolveProject } from './workflows/rnaseq'
include { DOWNLOAD } from './workflows/download'
include { BUILD_INDICES } from './workflows/build_indices'

/*
 * Parameter help (--help), typo detection and required-input checks, driven by
 * nextflow_schema.json. Called at the start of the entry workflow's main:
 * section.
 *
 * Nextflow's strict script parser (default from 25.10 / 26.x) forbids top-level
 * statements in a script, so this logic lives in a function rather than running
 * at the top of the file. The schema is loaded defensively: if it is missing or
 * unreadable the pipeline still runs, with help / parameter checking skipped.
 */

/*
 * A parameter as a number, or null when it is not one.
 *
 * Needed because a value from the command line arrives as a String while the
 * same parameter declared in nextflow.config is a number, so neither a cast nor
 * a type check works for both. `as Double` on a non-numeric String throws, and
 * the message names a Groovy internal rather than the parameter, so the caller
 * gets null and writes its own error.
 *
 * A top-level function rather than a local closure, because the strict parser
 * does not resolve a closure variable invoked as a function. That is the same
 * constraint recorded for the run-summary helpers in workflows/rnaseq.nf, and
 * the first version of this check hit it: every value failed, valid ones
 * included.
 */
def asNumber(value) {
    if (value == null) return null
    try {
        return value.toString().toDouble()
    }
    catch (Exception e) {
        return null
    }
}

def checkParameters() {
    def schema_json = null
    try {
        schema_json = new groovy.json.JsonSlurper().parseText(
            file("${projectDir}/nextflow_schema.json").text)
    } catch (Exception e) {
        log.warn "Could not read nextflow_schema.json (${e.message}); --help and parameter checking disabled."
    }

    if (schema_json) {
        def schema_groups = [:]   // group title -> [ param names ]
        def schema_params = [:]   // param name  -> definition map
        (schema_json['$defs'] ?: [:]).each { gkey, group ->
            def names = []
            (group['properties'] ?: [:]).each { pname, pdef ->
                schema_params[pname] = pdef
                names << pname
            }
            schema_groups[group['title'] ?: gkey] = names
        }

        if (schema_params.isEmpty()) {
            log.warn "nextflow_schema.json defines no parameters; --help and parameter checking disabled."
        }
        else {
            // --help : print the grouped parameter list and exit.
            if (params.help) {
                def sb = new StringBuilder("\nrnaseq-flow - pipeline parameters\n")
                schema_groups.each { title, names ->
                    // Schema entries marked "hidden" are set by execution
                    // profiles rather than on the command line; they are still
                    // accepted by the typo check below, just not advertised. A
                    // group whose parameters are all hidden is skipped entirely.
                    def shown = names.findAll { pn -> !(schema_params[pn]?.get('hidden')) }
                    if (!shown) return
                    sb << "\n${title}\n"
                    shown.each { pn ->
                        def d   = schema_params[pn] ?: [:]
                        def dft = d.containsKey('default') ? "  (default: ${d['default']})" : ""
                        sb << String.format("  --%-18s %s%s%n", pn, (d['description'] ?: ''), dft)
                    }
                }
                sb << "\nExample:\n"
                sb << "  nextflow run main.nf --input samplesheet.csv --aligner star --gtf genes.gtf -profile docker\n"
                log.info sb.toString()
                System.exit(0)
            }

            // Typo detection : any --parameter not declared in the schema is
            // almost certainly a typo, and would otherwise be silently ignored.
            def known   = schema_params.keySet() as Set
            def unknown = params.keySet().findAll { p -> !known.contains(p.toString()) }
            if (unknown) {
                error "Unknown parameter(s): " + unknown.collect { p -> "--${p}" }.sort().join(', ') +
                      "\nThis is usually a typo. Run 'nextflow run main.nf --help' for the valid parameter list."
            }
        }
    }

    // Required-input check.
    if (params.input == null && !params.download_refs && !params.build_indices) {
        error "Please specify an input samplesheet or glob pattern with --input"
    }

    // An execution engine has to be chosen. Each of -profile docker/singularity/
    // conda enables exactly one and disables the others, but there is no default:
    // with no profile, every container directive is ignored and each script runs
    // against whatever happens to be on PATH. That is a silently uncontainerised
    // run for a pipeline whose reproducibility rests on containers, so fail
    // rather than produce results nobody can trace to a software version.
    // workflow.containerEngine covers docker and singularity and is null when
    // neither is on. Conda has no equivalent property, so it is read from the
    // active profile list, which is why enabling conda through a bare -c file
    // rather than -profile conda would still trip this. That is rare, and the
    // message says which profiles to use.
    //
    // Skipped for -stub-run, which executes no command and so needs no engine.
    // That also keeps the stub usable as a CI smoke test without pulling
    // anything.
    def engine_on = workflow.containerEngine ||
                    (workflow.profile?.toString()?.tokenize(',')?.contains('conda') ?: false)
    if (!workflow.stubRun && !engine_on) {
        error "No execution engine selected. Add one of -profile docker, " +
              "-profile singularity or -profile conda (combine with others using " +
              "a comma, e.g. -profile test_yeast,docker). Without one, container " +
              "directives are ignored and every tool is taken from PATH."
    }

    // Significance thresholds, checked here rather than left to R. The analysis
    // scripts validate them too, but they run after alignment and
    // quantification, so a typo would cost the whole run before being noticed.
    //
    // asNumber() is a top-level function (see above it) rather than a local
    // closure: the strict parser does not resolve a closure variable invoked
    // as a function, which is the same constraint recorded in
    // workflows/rnaseq.nf for its run-summary helpers.
    def padj_cut = asNumber(params.padj_cutoff)
    if (padj_cut == null) {
        error "--padj_cutoff must be a number; got '${params.padj_cutoff}'"
    }
    if (padj_cut <= 0 || padj_cut > 1) {
        error "--padj_cutoff must be greater than 0 and at most 1; got ${padj_cut}. " +
              "It is an adjusted p-value, so a value above 1 selects every gene."
    }
    def lfc_cut = asNumber(params.lfc_cutoff)
    if (lfc_cut == null) {
        error "--lfc_cutoff must be a number; got '${params.lfc_cutoff}'"
    }
    if (lfc_cut < 0) {
        error "--lfc_cutoff must be 0 or more; got ${lfc_cut}. It is compared against " +
              "an absolute log2 fold change, so a negative value selects every gene. " +
              "Use 0 to select on the adjusted p-value alone."
    }

    // The STAR route cannot run without an annotation: STAR_ALIGN passes
    // --sjdbGTFfile unconditionally, and the workflow hands it file(params.gtf).
    // Without --gtf that call failed with "Argument of `file()` function cannot
    // be null", which names neither the parameter nor the aligner. HISAT2 takes
    // no GTF at alignment time, so a GTF-less run there is legitimate and only
    // skips the annotation-dependent stages.
    //
    // Exempt when alignment never happens: --stop_at preQC stops at read QC,
    // and --counts enters below alignment and already demands --gtf of its own
    // accord just below.
    if (params.aligner == 'star' && !params.gtf && !params.counts &&
        params.stop_at != 'preQC' && !params.build_indices && !params.download_refs) {
        error "--aligner star requires --gtf: STAR is given the annotation as " +
              "--sjdbGTFfile. Pass --gtf, or use --aligner hisat2, which does not " +
              "need one at alignment time."
    }

    // Count-matrix entry (--counts): starts at differential expression, so it
    // needs the samplesheet (for the design) and the GTF (for gene annotation),
    // and it is incompatible with the pre-DE stop points.
    if (params.counts) {
        if (params.input == null) {
            error "--counts requires --input: the samplesheet supplies the condition/design."
        }
        if (!params.gtf) {
            error "--counts requires --gtf: gene-level annotation is mandatory (results must " +
                  "carry gene symbols and biotypes, not bare gene IDs)."
        }
        if (params.stop_at in ['preQC', 'postQC']) {
            error "--counts enters at differential expression; --stop_at ${params.stop_at} " +
                  "would leave nothing to run."
        }
        if (params.diffsplice || params.dtu || params.isoform_switch || params.ctat_lib) {
            log.warn "--counts: splicing / DTU / isoform-switch / fusion analyses need reads or " +
                     "BAMs and cannot run from a count matrix; those options are ignored."
        }
        log.info "Count-matrix input: skipping QC and alignment, entering at differential expression."
    }

    // Fail-fast on file/directory parameters that do not exist. Nextflow stages
    // a missing path as a dangling symlink rather than refusing to start, so the
    // error surfaces only when the process that needs it runs: a mistyped --gmt
    // failed at GSEA, after alignment and differential expression had already
    // taken a quarter of an hour. Checking here costs nothing and reports every
    // bad path at once.
    //
    // Skipped: values containing a glob (resolved later, and legitimately may
    // match nothing yet) and remote URIs (s3://, https://, ...), which cannot be
    // checked cheaply or offline.
    def path_params = ['gmt', 'gtf', 'transcript_fasta', 'counts',
                       'genome_fasta', 'star_index', 'hisat2_index',
                       'salmon_index', 'kallisto_index', 'ctat_lib']
    def bad_paths = path_params.findAll { p ->
        def v = params[p]
        if (!v) return false
        def s = v.toString()
        if (s.contains('://') || s.contains('*') || s.contains('?')) return false
        !file(s).exists()
    }
    if (bad_paths) {
        error "Path not found:\n" +
              bad_paths.collect { p -> "  --${p} ${params[p]}" }.join('\n') +
              "\nCheck each path exists and is readable from the launch directory."
    }
}

/*
 * Fail-fast samplesheet validation.
 *
 * Runs before any process is scheduled (no containers pulled, no work done) so
 * mistakes surface immediately. All problems are collected and reported at once.
 * Checks:
 *   - required columns (sample, R1, condition) are present
 *   - sample ids are non-empty and unique
 *   - every R1 (and R2 when given) FASTQ file exists
 *   - there are >= 2 conditions and >= 2 replicates per condition
 * The condition/replicate rules are hard errors when differential expression
 * will run (DESeq2/edgeR need them), and warnings otherwise (e.g. --stop_at preQC).
 */
/*
 * Is this string already a syntactically valid R name, i.e. does make.names()
 * leave it alone?
 *
 * It matters because DESeq2 builds resultsNames() from make.names()-sanitised
 * level names, and deseq2.R looks a coefficient up by constructing that name
 * itself: paste0("condition_", condA, "_vs_", condB). A condition containing a
 * hyphen, a space or a leading digit therefore produces a name DESeq2 does not
 * have, and results(name=) errors out after alignment and quantification have
 * already run.
 *
 * "TGF-beta", "6h" and "Ctrl 1" are all ordinary labels and all fail. The
 * pipeline's own example uses "TGFb", which passes, so this has been luck
 * rather than protection.
 *
 * R's rule: letters, digits, '.' and '_' only, starting with a letter or with
 * a dot not followed by a digit, and not a reserved word. make.names() would
 * rewrite anything else.
 */
def isValidRName(String s) {
    if (!(s ==~ /^(?:[a-zA-Z]|\.(?![0-9]))[a-zA-Z0-9._]*$/)) return false
    // Reserved words get a '.' appended by make.names(), so they fail too.
    return !(s in ['if', 'else', 'repeat', 'while', 'function', 'for', 'next',
                   'break', 'TRUE', 'FALSE', 'NULL', 'Inf', 'NaN', 'NA',
                   'NA_integer_', 'NA_real_', 'NA_character_'])
}

def validateSamplesheet(samplesheet_path) {
    def sheet = file(samplesheet_path)
    if (!sheet.exists()) {
        error "Samplesheet not found: ${samplesheet_path}"
    }

    def rows = sheet.splitCsv(header: true)
    if (!rows) {
        error "Samplesheet '${samplesheet_path}' has a header but no data rows."
    }

    def errors   = []
    def warnings = []

    // Required columns. R1 is only needed when the pipeline reads FASTQs; a
    // --counts run supplies expression directly, so it needs sample + condition
    // (plus any design covariates) but no read files.
    def cols = rows[0].keySet()
    def required_cols = (params.counts == null) ? ['sample', 'R1', 'condition']
                                                : ['sample', 'condition']
    required_cols.each { c ->
        if (!cols.contains(c)) errors << "missing required column '${c}'"
    }
    // --design covariate columns must exist; 'condition' must be in the design.
    if (params.design) {
        def dvars = (params.design =~ /[A-Za-z_]\w*/).collect { v -> v } as Set
        dvars.findAll { d -> d != 'condition' }.each { v ->
            if (!cols.contains(v)) errors << "--design variable '${v}' is not a samplesheet column"
        }
        if (!dvars.contains('condition')) {
            errors << "--design must include 'condition' (the variable DESeq2/edgeR contrasts)"
        }
    }
    if (errors) {
        error "Samplesheet validation failed:\n  - " + errors.join("\n  - ")
    }

    def has_batch        = cols.contains('batch')
    def seen_samples     = [] as Set
    def condition_counts = [:].withDefault { 0 }
    def batch_values     = [] as Set

    // Any other --design covariate column. 'condition' and 'batch' are checked
    // by name above; a third variable introduced through
    // --design '~ timepoint + condition' had its column checked for existence
    // and its values checked for nothing.
    def design_cols = [] as Set
    if (params.design) {
        design_cols = ((params.design =~ /[A-Za-z_]\w*/).collect { v -> v } as Set)
            .findAll { v -> v != 'condition' && v != 'batch' && cols.contains(v) }
    }
    def design_values = [:]
    design_cols.each { c -> design_values[c] = [] as Set }

    rows.eachWithIndex { row, i ->
        def lineno = i + 2   // +1 for 0-index, +1 for the header line
        def s = row.sample?.trim()
        if (!s) {
            errors << "row ${lineno}: empty 'sample' value"
        } else if (seen_samples.contains(s)) {
            errors << "row ${lineno}: duplicate sample id '${s}'"
        } else {
            seen_samples << s
        }

        def cond = row.condition?.trim()
        if (!cond) {
            errors << "row ${lineno} (sample '${s}'): empty 'condition' value"
        } else {
            condition_counts[cond] += 1
        }

        if (has_batch) {
            def b = row.batch?.trim()
            if (!b) errors << "row ${lineno} (sample '${s}'): empty 'batch' value"
            else    batch_values << b
        }

        design_cols.each { c ->
            def v = row[c]?.toString()?.trim()
            if (!v) errors << "row ${lineno} (sample '${s}'): empty '${c}' value"
            else    design_values[c] << v
        }

        // FASTQ existence (skipped during -stub-run, which never reads the data,
        // and for --counts, which has no read files to check)
        if (!workflow.stubRun && params.counts == null) {
            def r1 = row.R1?.trim()
            if (!r1) {
                errors << "row ${lineno} (sample '${s}'): empty 'R1' path"
            } else if (!file(r1).exists()) {
                errors << "row ${lineno} (sample '${s}'): R1 file not found: ${r1}"
            }
            def r2 = row.R2?.trim()
            if (r2 && !file(r2).exists()) {
                errors << "row ${lineno} (sample '${s}'): R2 file not found: ${r2}"
            }
        }
    }

    // The baseline condition must exist, when one was named explicitly.
    //
    // It is the denominator of every contrast, so it decides which way every
    // fold change and every PSI difference points. Asking for a level that is
    // not in the samplesheet is a mistake worth failing on: the analysis scripts
    // would otherwise fall back to alphabetical order and produce a complete,
    // plausible set of results oriented against the wrong baseline.
    //
    // Only checked when the user named one. Left at the default, an absent
    // 'REF' keeps the existing behaviour -- the scripts warn and fall back --
    // because a samplesheet with no REF condition has always been allowed.
    def ref_level = params.reference_level?.toString()
    if (ref_level && ref_level != 'REF' && !condition_counts.containsKey(ref_level)) {
        errors << "--reference_level '${ref_level}' is not a condition in the samplesheet. " +
                  "Found: ${condition_counts.keySet().sort().join(', ')}"
    }

    // Condition, batch and --design covariate values have to survive
    // make.names() unchanged. See isValidRName() above for why: DESeq2 names
    // its coefficients from the sanitised levels, and deseq2.R reconstructs
    // that name by hand, so a mismatch fails only once alignment and
    // quantification are done.
    //
    // A covariate fails more gently than a condition, because nothing rebuilds
    // its coefficient name by hand: a timepoint of '6h' reaches resultsNames()
    // as 'timepointX6h' instead of 'timepoint6h', a name that is in no input
    // file, and the model is otherwise fine. Checked all the same, since it is
    // one loop here against explaining that discrepancy to someone later.
    //
    // Checked here rather than in the R scripts so it lands with every other
    // samplesheet problem, before any container is pulled.
    def bad_names = []
    condition_counts.keySet().each { c ->
        if (!isValidRName(c)) bad_names << "condition '${c}'"
    }
    batch_values.each { b ->
        if (!isValidRName(b)) bad_names << "batch '${b}'"
    }
    design_values.each { col, vals ->
        vals.each { v ->
            if (!isValidRName(v)) bad_names << "${col} '${v}'"
        }
    }
    if (bad_names) {
        errors << "these values are not valid R names, which DESeq2 requires for " +
                  "its model levels: ${bad_names.join(', ')}. Use letters, digits, " +
                  "'.' and '_' only, starting with a letter. For example 'TGF-beta' " +
                  "becomes 'TGFbeta' or 'TGF_beta', and '6h' becomes 'h6' or 'T6h'."
    }

    // Differential-expression design checks
    def de_msgs = []
    if (condition_counts.size() < 2) {
        de_msgs << "found ${condition_counts.size()} condition(s) [${condition_counts.keySet().join(', ')}]; " +
                   "differential expression needs at least 2"
    }
    condition_counts.each { cond, n ->
        if (n < 2) {
            de_msgs << "condition '${cond}' has only ${n} replicate(s); DESeq2/edgeR need at least 2 per condition"
        }
    }

    // DE runs unless the pipeline is stopped before it.
    def de_will_run = !(params.stop_at in ['preQC', 'postQC'])
    if (de_will_run) {
        errors += de_msgs
    } else {
        warnings += de_msgs
    }

    if (has_batch && batch_values.size() < 2) {
        warnings << "the 'batch' column has only ${batch_values.size()} level(s); it adds nothing to the model"
    }

    warnings.each { w -> log.warn "Samplesheet: ${w}" }

    if (errors) {
        error "Samplesheet validation failed (${errors.size()} problem(s)):\n  - " +
              errors.join("\n  - ")
    }

    log.info "Samplesheet OK: ${seen_samples.size()} sample(s), " +
             "${condition_counts.size()} condition(s) " +
             "[${condition_counts.collect { k, v -> "${k}:${v}" }.join(', ')}]"
}

/*
 * Run-summary formatting helpers, used by the entry workflow's onComplete:
 * section. These are top-level functions rather than local closures because the
 * strict parser (Nextflow 25.10 / 26.x) does not resolve a closure variable
 * that is called from inside another closure (e.g. inside a .collect { } that
 * builds the HTML), whereas a top-level function resolves everywhere.
 */
def esc(s) {
    s == null ? '' : s.toString().replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
}
def num(v) {
    (v == null || !v.toString().trim().isNumber()) ? 0d : v.toString().trim().toDouble()
}
def fmtTime(ms) {
    double t = ms as double
    if (t <= 0) return '0s'
    long s = Math.round(t / 1000d)
    long h = s.intdiv(3600L)
    long m = (s % 3600L).intdiv(60L)
    long sec = s % 60L
    def p = []
    if (h)         p << "${h}h"
    if (m)         p << "${m}m"
    if (sec || !p) p << "${sec}s"
    p.join(' ')
}
def fmtMem(b) {
    double v = b as double
    if (v <= 0) return '—'
    def u = ['B', 'KB', 'MB', 'GB', 'TB']
    int i = 0
    // Range .each replaces a while loop (removed in the strict language); the
    // guard makes the extra iterations no-ops once v has been scaled down.
    (0..<4).each { _k ->
        if (v >= 1024d) {
            v /= 1024d
            i += 1
        }
    }
    i == 0 ? "${Math.round(v)} B" : String.format('%.1f %s', v, u[i])
}

/*
 * Main entry point
 */
workflow {

    main:
    // Help / typo-detection / required-input checks (may print help and exit).
    checkParameters()

    // Attribution and project metadata, checked here so a misspelled CRediT
    // role or a malformed RAiD fails before any process is scheduled rather
    // than at the end of a run that has already cost an hour. The results are
    // discarded: buildRunManifest() resolves them again when it writes the
    // manifest, and both calls are pure.
    resolveContributors()
    resolveProject()

    // --validate_only concerns the samplesheet, which the helper workflows do
    // not take. Warn rather than fail: the flag is harmless here, but silently
    // ignoring it would let someone believe a reference download had been
    // checked over when nothing was.
    if (params.validate_only && (params.download_refs || params.build_indices)) {
        log.warn "--validate_only has no effect with --download_refs or --build_indices; " +
                 "neither takes a samplesheet."
    }

    if (params.download_refs) {
        DOWNLOAD( params.download_species, params.download_source )
    }
    else if (params.build_indices) {
        BUILD_INDICES( params.aligner )
    }
    else {
        //
        // Create input channel
        //
        if (params.counts) {
            // Count-matrix entry: validate the (relaxed) samplesheet, but build
            // no read channel: RNASEQ enters directly at differential
            // expression and consumes params.counts, not FASTQs.
            validateSamplesheet( params.input )
            channel.empty().set { ch_reads }
        }
        else if (params.input.endsWith('.csv')) {
            // Fail fast on a bad samplesheet before scheduling any work.
            validateSamplesheet( params.input )

            channel
                .fromPath(params.input)
                .splitCsv(header: true)
                .map { row ->
                    def meta = [:]
                    meta.id = row.sample
                    meta.condition = row.condition

                    // Check if R2 exists and is not empty to determine single/paired
                    def fastq_1 = file(row.R1)
                    def fastq_2 = row.R2 ? file(row.R2) : null

                    meta.single_end = (fastq_2 == null)

                    if (meta.single_end) {
                        return [ meta, [ fastq_1 ] ]
                    } else {
                        return [ meta, [ fastq_1, fastq_2 ] ]
                    }
                }
                .set { ch_reads }
        } else {
            // Glob input: no samplesheet, so there is no condition column and
            // differential expression cannot run. Read layout is taken from what
            // fromFilePairs actually grouped, not guessed from the pattern.
            //
            // `size: -1` means "however many files share a key". It replaced
            // `size: params.input.count('*') > 1 ? 2 : 1`, which counted
            // asterisks in the pattern: the commonest paired-end glob there is,
            // /data/*_R{1,2}.fastq.gz, has exactly one, so it was read as
            // single-end and each mate was aligned as a separate sample with
            // half the evidence and a doubled sample count (m7). Counting
            // asterisks could not work, because the alternation group rather
            // than the wildcard is what makes a pattern paired.
            channel
                .fromFilePairs( params.input, size: -1 )
                .ifEmpty { error "Cannot find any reads matching: ${params.input}" }
                .map { name, reads ->
                    // Always a list under size: -1, but normalised anyway so the
                    // size test below cannot be fooled by a bare path.
                    def files = (reads instanceof List) ? reads : [ reads ]
                    if (files.size() > 2) {
                        error "Read pattern '${params.input}' grouped ${files.size()} files " +
                              "under sample '${name}': ${files.collect { it.name }.join(', ')}. " +
                              "Expected one file for single-end or two for paired-end. Narrow the " +
                              "pattern, or use a samplesheet with --input, which names each file."
                    }
                    // Derived once, from the files in hand. The previous version
                    // assigned this from `reads instanceof Path` and then
                    // overwrote it in the very next statement.
                    def meta = [ id: name, single_end: files.size() == 1 ]
                    return [ meta, files ]
                }
                .set { ch_reads }
        }

        //
        // Run workflow
        //
        if (params.validate_only) {
            // Validation has already run above, and exits non-zero on any
            // problem, so reaching this point means the inputs are usable.
            // Stop here rather than scheduling work: the point of the flag is
            // to answer "is my samplesheet right?" without waiting for a run.
            //
            // A glob input never reaches validateSamplesheet -- there is no
            // sample or condition column to check -- so say so plainly instead
            // of implying it passed something.
            if (!params.counts && !params.input.toString().endsWith('.csv')) {
                log.warn "--validate_only: the input is a file glob, not a samplesheet, " +
                         "so there is nothing to validate beyond the read-pair match."
            }
            log.info "--validate_only: inputs are valid. No processes were run."
        }
        else {
            RNASEQ ( ch_reads )
        }
    }

    /*
     * Run-completion summary.
     *
     * After every run (success or failure) this writes
     *   <outdir>/pipeline_info/run_summary.html
     * A self-contained page that links the MultiQC report and the key result
     * directories, and tabulates job time, peak memory and CPU usage per process.
     * The per-process figures are aggregated from the Nextflow execution trace
     * (trace.raw = true, set in nextflow.config, keeps that file machine-readable).
     * A concise version is also printed to the console.
     *
     * This is the entry workflow's onComplete: section (Nextflow strict syntax);
     * it replaces the former top-level `workflow.onComplete { }` handler, which
     * the strict parser does not allow.
     */
    onComplete:
    try {
        // --- aggregate per-process resources from the execution trace --------
        def traceHits = []
        try { traceHits = files("${params.tracedir}/execution_trace_*.txt") }
        catch (_ignored) { traceHits = [] }
        def traceFile = traceHits ? traceHits.sort { a -> a.name }.last() : null
        def procStats = [:]
        boolean traceOk = false
        if (traceFile && traceFile.exists()) {
            def lines = traceFile.readLines()
            if (lines.size() > 1) {
                traceOk = true
                def hdr = (lines[0].split('\t') as List).collect { c -> c.trim() }
                int iName = hdr.indexOf('name')
                int iStat = hdr.indexOf('status')
                int iReal = hdr.indexOf('realtime')
                int iCpu  = hdr.indexOf('%cpu')
                int iRss  = hdr.indexOf('peak_rss')
                lines.drop(1).each { ln ->
                    def f = ln.split('\t')
                    if (iName < 0 || iName >= f.size()) return
                    def proc = f[iName].replaceFirst(/ \(.*\)\s*$/, '').trim()
                    if (proc.contains(':')) proc = proc.tokenize(':').last()
                    if (!proc) return
                    def st = procStats[proc]
                    if (st == null) {
                        st = [tasks: 0, real: 0d, cpuSum: 0d, cpuN: 0, rss: 0d, failed: 0]
                        procStats[proc] = st
                    }
                    st.tasks += 1
                    if (iReal >= 0 && iReal < f.size()) st.real += num(f[iReal])
                    if (iCpu >= 0 && iCpu < f.size()) {
                        double c = num(f[iCpu])
                        if (c > 0) {
                            st.cpuSum += c
                            st.cpuN += 1
                        }
                    }
                    if (iRss >= 0 && iRss < f.size()) {
                        double r = num(f[iRss])
                        if (r > st.rss) st.rss = r
                    }
                    if (iStat >= 0 && iStat < f.size() && f[iStat] == 'FAILED') st.failed += 1
                }
            }
        }
        def procRows = procStats.collect { k, v -> [name: k] + v }.sort { r -> -r.real }
        int    totTasks = (procRows.sum { r -> r.tasks } ?: 0) as int
        double totReal  = (procRows.sum { r -> r.real } ?: 0d) as double

        // --- per-process resource table --------------------------------------
        def procTable
        if (procRows) {
            procTable = procRows.collect { p ->
                double meanCpu = p.cpuN ? (p.cpuSum / p.cpuN) : 0d
                double avg     = p.tasks ? (p.real / p.tasks) : 0d
                """        <tr${p.failed ? ' class="fail"' : ''}>
          <td class="proc">${esc(p.name)}</td>
          <td class="n">${p.tasks}${p.failed ? " (${p.failed} failed)" : ''}</td>
          <td class="n">${fmtTime(p.real)}</td>
          <td class="n">${fmtTime(avg)}</td>
          <td class="n">${fmtMem(p.rss)}</td>
          <td class="n">${String.format('%.0f%%', meanCpu)}</td>
        </tr>"""
            }.join('\n')
            procTable += """
        <tr class="tot">
          <td class="proc">All processes (${procRows.size()})</td>
          <td class="n">${totTasks}</td>
          <td class="n">${fmtTime(totReal)}</td>
          <td class="n"></td><td class="n"></td><td class="n"></td>
        </tr>"""
        } else {
            procTable = '        <tr><td colspan="6" class="proc">' +
                (traceOk ? 'No task records were found in the execution trace.'
                         : 'Execution trace not available - per-process figures could not be computed.') +
                '</td></tr>'
        }

        // --- key outputs ------------------------------------------------------
        def candidates = [
            ['MultiQC report',                 'multiqc/multiqc_report.html'],
            ['Quarto analysis report',         'quarto_report/analysis_report.html'],
            ['DESeq2 differential expression', 'deseq2_output'],
            ['edgeR differential expression',  'edger_output'],
            ['GSEA enrichment',                'gsea_output'],
            ['gProfiler enrichment',           'gprofiler_output'],
            ['Alternative splicing (rMATS)',   'rmats_output'],
            ['Differential transcript usage',  'dtu_output'],
            ['Differential splicing',          'diffsplice_output'],
            ['Gene fusions (STAR-Fusion)',     'star_fusion'],
            ['Isoform switching',              'isoform_switch'],
            ['Gene counts (featureCounts)',    'featurecounts'],
            ['Transcript-to-gene counts',      'tximport'],
            ['Coverage tracks (BigWig)',       'bam_to_bigwig'],
            ['STAR index',                     'star_index'],
            ['HISAT2 index',                   'hisat2_index'],
            ['Salmon index',                   'salmon_index'],
            ['Kallisto index',                 'kallisto_index'],
        ]
        def outItems = []
        candidates.each { c ->
            def p = file("${params.outdir}/${c[1]}")
            if (p.exists()) {
                def dir = p.isDirectory()
                outItems << ('          <li><a href="../' + c[1] + (dir ? '/' : '') + '">' +
                    esc(c[0]) + '</a> <span class="path">' + esc(c[1]) + (dir ? '/' : '') +
                    '</span></li>')
            }
        }
        def outRows = outItems ? outItems.join('\n')
            : ('          <li class="path">No published result directories were found under ' +
               esc(params.outdir) + '/.</li>')

        // Nextflow's own detailed execution report (same folder as this file).
        def repHits = []
        try { repHits = files("${params.tracedir}/execution_report_*.html") }
        catch (_ignored) { repHits = [] }
        def execLine = ''
        if (repHits) {
            def rn = repHits.sort { a -> a.name }.last().name
            execLine = '<p class="note">Detailed Nextflow execution report: ' +
                '<a href="' + esc(rn) + '">' + esc(rn) + '</a></p>'
        }

        // --- run details ------------------------------------------------------
        boolean ok = workflow.success
        def details = [
            ['Pipeline',       "${workflow.manifest.name ?: 'rnaseq-flow'} v${workflow.manifest.version ?: '?'}"],
            ['Run name',       workflow.runName],
            ['Profile',        workflow.profile],
            ['Command',        workflow.commandLine],
            ['Started',        workflow.start],
            ['Completed',      workflow.complete],
            ['Duration',       workflow.duration],
            ['Exit status',    workflow.exitStatus],
            ['Work directory', workflow.workDir],
            ['Nextflow',       workflow.nextflow?.version],
        ]
        def detailRows = details.collect { d ->
            '        <tr><td class="k">' + esc(d[0]) + '</td><td class="v">' +
            esc(d[1]) + '</td></tr>'
        }.join('\n')

        def stubNote   = workflow.stubRun ? ' — stub run (resource figures are placeholders)' : ''
        def errorBlock = (!ok && workflow.errorMessage)
            ? ('\n      <div class="err">' + esc(workflow.errorMessage) + '</div>') : ''

        def html = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>rnaseq-flow - run summary</title>
<style>
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:-apple-system,'Segoe UI',Roboto,Helvetica,Arial,sans-serif;background:#f1f3f5;color:#212529;line-height:1.5;padding:32px 16px}
.wrap{max-width:920px;margin:0 auto;background:#fff;border-radius:10px;overflow:hidden;box-shadow:0 1px 5px rgba(0,0,0,.12)}
header{background:#1F4E79;color:#fff;padding:22px 28px}
header h1{font-size:22px;font-weight:700;letter-spacing:.5px}
header .sub{font-size:13px;opacity:.82;margin-top:2px}
.status{padding:13px 28px;font-weight:600;font-size:15px;color:#fff}
.status.ok{background:#2f9e44}
.status.fail{background:#c92a2a}
section{padding:20px 28px;border-top:1px solid #e9ecef}
section h2{font-size:13px;color:#1F4E79;margin-bottom:12px;text-transform:uppercase;letter-spacing:.6px}
table{width:100%;border-collapse:collapse;font-size:13px}
table.kv td{padding:5px 8px;vertical-align:top}
table.kv td.k{color:#868e96;width:150px;white-space:nowrap}
table.kv td.v{font-family:Consolas,Menlo,monospace;word-break:break-all}
table.res th{background:#1F4E79;color:#fff;text-align:left;padding:7px 10px;font-size:12px}
table.res td{padding:6px 10px;border-bottom:1px solid #e9ecef}
table.res td.n{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
table.res td.proc{font-weight:600}
table.res tr.tot td{font-weight:700;border-top:2px solid #1F4E79;background:#f8f9fa}
table.res tr.fail td{background:#fff0f0}
ul.outputs{list-style:none}
ul.outputs li{padding:6px 0;border-bottom:1px solid #f1f3f5}
ul.outputs a{color:#2E75B6;font-weight:600;text-decoration:none}
ul.outputs a:hover{text-decoration:underline}
.path{color:#868e96;font-family:Consolas,monospace;font-size:12px;margin-left:6px}
.note{color:#868e96;font-size:12px;margin-top:10px}
.err{background:#fff0f0;border-left:3px solid #c92a2a;padding:10px 12px;font-family:Consolas,monospace;font-size:12px;white-space:pre-wrap;margin-top:6px}
footer{padding:14px 28px;color:#868e96;font-size:12px;background:#f8f9fa}
</style>
</head>
<body>
<div class="wrap">
  <header>
    <h1>rnaseq-flow</h1>
    <div class="sub">Run completion summary</div>
  </header>
  <div class="status ${ok ? 'ok' : 'fail'}">${ok ? 'Completed successfully' : 'Run failed'}${stubNote}</div>
  <section>
    <h2>Run details</h2>
    <table class="kv">
${detailRows}
    </table>${errorBlock}
  </section>
  <section>
    <h2>Key outputs</h2>
    <ul class="outputs">
${outRows}
    </ul>
    ${execLine}
  </section>
  <section>
    <h2>Resource usage by process</h2>
    <table class="res">
      <thead><tr><th>Process</th><th>Tasks</th><th>Total time</th><th>Avg / task</th><th>Peak memory</th><th>Mean CPU</th></tr></thead>
      <tbody>
${procTable}
      </tbody>
    </table>
    <p class="note">Total time is the summed task run-time for the process; peak memory is the largest peak_rss of any of its tasks; mean CPU averages %cpu across its tasks (above 100% means multiple cores were used). Figures are aggregated from the Nextflow execution trace.</p>
  </section>
  <footer>Generated ${esc(new Date().format('yyyy-MM-dd HH:mm:ss'))} &middot; rnaseq-flow ${esc(workflow.manifest.version ?: '')}</footer>
</div>
</body>
</html>
"""

        def summaryFile = file("${params.tracedir}/run_summary.html")
        summaryFile.parent.mkdirs()
        summaryFile.text = html

        // --- concise console summary -----------------------------------------
        def heaviest = procRows.take(3).collect { r -> "${r.name} (${fmtTime(r.real)})" }
        log.info(
            '\n' + ('-' * 66) + '\n' +
            " rnaseq-flow run ${ok ? 'COMPLETE' : 'FAILED'}  |  duration ${workflow.duration}\n" +
            " ${procRows.size()} process(es), ${totTasks} task(s), ${fmtTime(totReal)} total compute time\n" +
            (heaviest ? " heaviest: ${heaviest.join(', ')}\n" : '') +
            " summary: ${summaryFile.toUriString()}\n" +
            // The RO-Crate cannot be finished from inside the pipeline: nf-prov
            // writes ro-crate-metadata.json after this handler has already run,
            // measured at about 1.6 s later, so there is no hook that sees the
            // file. The command is printed here instead, every run, because a
            // crate without checksums or credit looks complete and is not.
            ('-' * 66))
    }
    catch (Exception e) {
        log.warn "Could not write the run-completion summary: ${e.message}"
    }

    // --- finish the RO-Crate ------------------------------------------------
    //
    // In a JVM shutdown hook rather than here, because nf-prov writes
    // ro-crate-metadata.json about 1.6 s after this handler returns: no
    // Nextflow hook ever sees the file. A shutdown hook runs later than every
    // hook and does see it, which was measured before relying on it.
    //
    // The work is in lib/CrateEnricher.groovy. It is Groovy rather than a
    // script because the launching machine is only guaranteed to have Java;
    // every tool the pipeline uses lives in a container.
    if (params.ro_crate.toString() != 'false') {
        def crateOutdir = params.outdir.toString()
        def crateManifest = "${params.tracedir}/run_manifest.json".toString()
        // Copied so the hook does not read Nextflow state during shutdown, and
        // so the enricher can resolve "Input file specified by params.<name>"
        // back to the file it was copied from.
        def crateParams = new LinkedHashMap(params)
        Runtime.runtime.addShutdownHook(new Thread({
            try {
                // Said before the work, not only after it. Every published file
                // is read and hashed with SHA-256 here, which on a cohort with
                // BAMs and BigWigs is minutes of real time starting after the
                // console has already printed "run COMPLETE". Unannounced that
                // reads as a hang, and a Ctrl-C part way through leaves
                // ro-crate-metadata.json half enriched with nothing recording
                // which half (m12). The checksums are what make a crate whose
                // inputs are references rather than copies verifiable at all,
                // so they stay; the wait is stated instead.
                if (new File(crateOutdir, CrateEnricher.CRATE).exists()) {
                    println " Completing the RO-Crate: reading and checksumming every published file. Let it finish."
                }
                def started = System.currentTimeMillis()
                def summary = CrateEnricher.enrich(crateOutdir, crateManifest, crateParams)
                // println, not log: the logger is being torn down by now.
                def took = String.format('%.1f', (System.currentTimeMillis() - started) / 1000d)
                if (summary) println " ${summary}, ${took}s"
            }
            catch (Exception e) {
                println " RO-Crate could not be completed: ${e}"
            }
        }))
    }
}
