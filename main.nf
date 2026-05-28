#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * Parameter help (--help) and typo detection, driven by nextflow_schema.json.
 * The schema is loaded defensively: if it is missing or unreadable the pipeline
 * still runs, with help / parameter checking simply skipped.
 */
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
                sb << "\n${title}\n"
                names.each { pn ->
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

        // Typo detection : any --parameter not declared in the schema is almost
        // certainly a typo, and would otherwise be silently ignored.
        def known   = schema_params.keySet() as Set
        def unknown = params.keySet().findAll { !known.contains(it.toString()) }
        if (unknown) {
            error "Unknown parameter(s): " + unknown.collect { "--${it}" }.sort().join(', ') +
                  "\nThis is usually a typo. Run 'nextflow run main.nf --help' for the valid parameter list."
        }
    }
}

/*
 * Workflow checks
 */
if (params.input == null && !params.download_refs && !params.build_indices) {
    error "Please specify an input samplesheet or glob pattern with --input"
}

/*
 * Import the main workflow
 */
include { RNASEQ } from './workflows/rnaseq'
include { DOWNLOAD } from './workflows/download'
include { BUILD_INDICES } from './workflows/build_indices'

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

    // Required columns
    def cols = rows[0].keySet()
    ['sample', 'R1', 'condition'].each { c ->
        if (!cols.contains(c)) errors << "missing required column '${c}'"
    }
    // --design covariate columns must exist; 'condition' must be in the design.
    if (params.design) {
        def dvars = (params.design =~ /[A-Za-z_]\w*/).collect { it } as Set
        dvars.findAll { it != 'condition' }.each { v ->
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

        // FASTQ existence (skipped during -stub-run, which never reads the data)
        if (!workflow.stubRun) {
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

    warnings.each { log.warn "Samplesheet: ${it}" }

    if (errors) {
        error "Samplesheet validation failed (${errors.size()} problem(s)):\n  - " +
              errors.join("\n  - ")
    }

    log.info "Samplesheet OK: ${seen_samples.size()} sample(s), " +
             "${condition_counts.size()} condition(s) " +
             "[${condition_counts.collect { k, v -> "${k}:${v}" }.join(', ')}]"
}

/*
 * Main entry point
 */
workflow {

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
        if (params.input.endsWith('.csv')) {
            // Fail fast on a bad samplesheet before scheduling any work.
            validateSamplesheet( params.input )

            Channel
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
            // Fallback to original Glob pattern logic
            Channel
                .fromFilePairs( params.input, size: params.input.count('*') > 1 ? 2 : 1 )
                .ifEmpty { error "Cannot find any reads matching: ${params.input}" }
                .map { name, reads ->
                    def meta = [:]
                    meta.id = name
                    meta.single_end = reads instanceof Path
                    if (reads instanceof List && reads.size() == 2) {
                        meta.single_end = false
                    } else {
                         meta.single_end = true
                    }
                    return [ meta, reads ]
                }
                .set { ch_reads }
        }

        //
        // Run workflow
        //
        RNASEQ ( ch_reads )
    }
}

/*
 * Run-completion summary.
 *
 * After every run (success or failure) this writes
 *   <outdir>/pipeline_info/run_summary.html
 * — a self-contained page that links the MultiQC report and the key result
 * directories, and tabulates job time, peak memory and CPU usage per process.
 * The per-process figures are aggregated from the Nextflow execution trace
 * (trace.raw = true, set in nextflow.config, keeps that file machine-readable).
 * A concise version is also printed to the console.
 */
workflow.onComplete {
    try {
        def esc = { s -> (s == null ? '' : s.toString())
            .replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;') }
        def num = { v -> (v == null || !v.toString().trim().isNumber())
            ? 0d : v.toString().trim().toDouble() }
        def fmtTime = { ms ->
            double t = ms as double
            if (t <= 0) return '0s'
            long s = Math.round(t / 1000d)
            long h = s.intdiv(3600L), m = (s % 3600L).intdiv(60L), sec = s % 60L
            def p = []
            if (h)         p << "${h}h"
            if (m)         p << "${m}m"
            if (sec || !p) p << "${sec}s"
            p.join(' ')
        }
        def fmtMem = { b ->
            double v = b as double
            if (v <= 0) return '—'
            def u = ['B', 'KB', 'MB', 'GB', 'TB']
            int i = 0
            while (v >= 1024d && i < 4) { v /= 1024d; i++ }
            i == 0 ? "${Math.round(v)} B" : String.format('%.1f %s', v, u[i])
        }

        // --- aggregate per-process resources from the execution trace --------
        def traceHits = []
        try { traceHits = files("${params.tracedir}/execution_trace_*.txt") }
        catch (ignored) { traceHits = [] }
        def traceFile = traceHits ? traceHits.sort { it.name }.last() : null
        def procStats = [:]
        boolean traceOk = false
        if (traceFile && traceFile.exists()) {
            def lines = traceFile.readLines()
            if (lines.size() > 1) {
                traceOk = true
                def hdr = (lines[0].split('\t') as List).collect { it.trim() }
                int iName = hdr.indexOf('name'),     iStat = hdr.indexOf('status')
                int iReal = hdr.indexOf('realtime'), iCpu  = hdr.indexOf('%cpu')
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
                        double c = num(f[iCpu]); if (c > 0) { st.cpuSum += c; st.cpuN += 1 }
                    }
                    if (iRss >= 0 && iRss < f.size()) {
                        double r = num(f[iRss]); if (r > st.rss) st.rss = r
                    }
                    if (iStat >= 0 && iStat < f.size() && f[iStat] == 'FAILED') st.failed += 1
                }
            }
        }
        def procRows = procStats.collect { k, v -> [name: k] + v }.sort { -it.real }
        int    totTasks = (procRows.sum { it.tasks } ?: 0) as int
        double totReal  = (procRows.sum { it.real } ?: 0d) as double

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
        catch (ignored) { repHits = [] }
        def execLine = ''
        if (repHits) {
            def rn = repHits.sort { it.name }.last().name
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
        def heaviest = procRows.take(3).collect { "${it.name} (${fmtTime(it.real)})" }
        log.info(
            '\n' + ('-' * 66) + '\n' +
            " rnaseq-flow run ${ok ? 'COMPLETE' : 'FAILED'}  |  duration ${workflow.duration}\n" +
            " ${procRows.size()} process(es), ${totTasks} task(s), ${fmtTime(totReal)} total compute time\n" +
            (heaviest ? " heaviest: ${heaviest.join(', ')}\n" : '') +
            " summary: ${summaryFile.toUriString()}\n" +
            ('-' * 66))
    }
    catch (Exception e) {
        log.warn "Could not write the run-completion summary: ${e.message}"
    }
}
