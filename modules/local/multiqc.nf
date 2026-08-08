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
    path samplesheet

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    path "*_plots"                  , optional:true, emit: plots
    // The aggregated tool versions, as a report section and as a plain file.
    // Previously software_versions.yml was a process *input* only, so a run's
    // versions were never published and could be recovered only from work/.
    path "software_versions_mqc.yml", optional:true, emit: software_versions
    path "versions.yml"             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def sheet = samplesheet ? "${samplesheet}" : ''
    // Sample-name reconciliation.
    //
    // Each tool names samples after whatever it was handed: FastQC after the
    // FASTQ filename ('DRR392077_1'), fastp after the --in1 path recorded in its
    // JSON ('DRR392077'), STAR and featureCounts after meta.id ('ypd_rep1').
    // MultiQC keys General Statistics on the sample name, so nothing merged --
    // six samples appeared as 24 rows, each mostly empty, and no row described
    // any sample completely.
    //
    // The fix is applied here rather than by renaming the reads upstream. The
    // published fastqc/<accession>_fastqc.html files and the paths recorded in
    // multiqc_sources.txt are the run's provenance: for public data the
    // accession is what ties results back to the archive, and rewriting the
    // filenames to normalise a *report* would throw that away. --replace-names
    // changes only what MultiQC displays; the source paths it records still
    // name the original files.
    //
    // Matching is on a substring by default, so one 'DRR392077 -> ypd_rep1' line
    // also covers 'DRR392077_1' and 'DRR392077_2'; the mate suffix survives, and
    // should, since R1 and R2 are separate measurements.
    // NOTE: this script block is written flush-left on purpose. The Python
    // heredoc below must start at column 0, which defeats Nextflow's
    // stripIndent(), so every other line has to be unindented too or the
    // versions.yml heredoc terminator stops matching.
    """
if [ -n "${sheet}" ] && [ -f "${sheet}" ]; then
    python3 - "${sheet}" <<'PY'
import csv, os, re, sys

# Two files are produced from the samplesheet:
#
#   replace_names.tsv          - the display mapping passed to --replace-names.
#   sample_provenance_mqc.tsv  - MultiQC custom content that adds a "Source
#                                FASTQ" column to General Statistics, so the
#                                original name is visible in the report and not
#                                merely recoverable from a filename. The Quarto
#                                report builds its QC table from MultiQC's
#                                general-stats TSV, so it inherits the column.
#
# Columns are looked up by name, not position.
EXT  = re.compile(r'\\.(fastq|fq)(\\.gz)?\$', re.I)
MATE = re.compile(r'[._-]?R?[12]\$')

rename, prov, seen = [], [], set()
with open(sys.argv[1], newline='') as fh:
    for row in csv.DictReader(fh):
        row = { (k or '').strip(): (v or '').strip() for k, v in row.items() }
        sample, r1, r2 = row.get('sample', ''), row.get('R1', ''), row.get('R2', '')
        if not sample or not r1:
            continue
        stem1 = EXT.sub('', os.path.basename(r1))
        base  = MATE.sub('', stem1)
        # Whatever mate suffix this dataset uses ('_1', '_R1', ...) survives the
        # substring replacement, so the per-mate rows are keyed the same way.
        suffix1 = stem1[len(base):]
        if base and base != sample and base not in seen:
            seen.add(base)
            rename.append(f"{base}\\t{sample}")
        prov.append(f"{sample}\\t{base}")
        if r2:
            stem2 = EXT.sub('', os.path.basename(r2))
            suffix2 = stem2[len(MATE.sub('', stem2)):]
            prov.append(f"{sample}{suffix1}\\t{stem1}")
            prov.append(f"{sample}{suffix2}\\t{stem2}")

with open('replace_names.tsv', 'w') as fh:
    fh.write("\\n".join(rename) + ("\\n" if rename else ""))

# MultiQC custom content: the '#'-prefixed YAML header declares this as extra
# General Statistics columns; the body is a plain Sample/value table.
with open('sample_provenance_mqc.tsv', 'w') as fh:
    fh.write("# id: 'sample_provenance'\\n")
    fh.write("# plot_type: 'generalstats'\\n")
    fh.write("# pconfig:\\n")
    fh.write("#     - source_fastq:\\n")
    fh.write("#         title: 'Source FASTQ'\\n")
    fh.write("#         description: 'Original FASTQ name from the samplesheet, before sample-name mapping'\\n")
    fh.write("#         scale: false\\n")
    fh.write("Sample\\tsource_fastq\\n")
    fh.write("\\n".join(prov) + ("\\n" if prov else ""))
PY
fi

# Software versions.
#
# The pipeline aggregates every process's versions.yml into software_versions.yml
# and stages it here, but MultiQC has no parser for a bare
# "PROCESS": {tool: version} mapping, so it was silently ignored: the Software
# Versions section showed only FastQC and fastp, the two tools whose own output
# embeds a version MultiQC parses for itself. STAR, samtools, subread, DESeq2,
# edgeR, fgsea, rMATS, RSeQC and the rest were collected and then dropped.
#
# Wrapping the same data in MultiQC's custom-content header turns it into a real
# report section. The file is also published, so a run's tool versions are
# recoverable from the results rather than only from work/.
if [ -f software_versions.yml ]; then
    python3 - <<'PY'
import yaml

with open('software_versions.yml') as fh:
    versions = yaml.safe_load(fh) or {}

rows = []
for process in sorted(versions):
    tools = versions[process] or {}
    # Processes are keyed '<WORKFLOW>:<PROCESS>'; the workflow prefix is the
    # same for every entry and only adds noise to the table.
    name = str(process).split(':')[-1]
    for tool in sorted(tools):
        rows.append((name, str(tool), str(tools[tool]).strip()))

html = ['<table class="table" style="width:auto;">',
        '<thead><tr><th>Process</th><th>Software</th><th>Version</th></tr></thead>',
        '<tbody>']
for name, tool, version in rows:
    html.append(f"<tr><td><samp>{name}</samp></td>"
                f"<td>{tool}</td><td>{version}</td></tr>")
html += ['</tbody></table>']

with open('software_versions_mqc.yml', 'w') as fh:
    # A distinct id: MultiQC has its own built-in "Software Versions" section
    # (populated only from tools whose output embeds a version, i.e. FastQC and
    # fastp here), and reusing 'software_versions' risks colliding with it.
    fh.write("id: 'rnaseq_flow_software_versions'\n")
    fh.write("section_name: 'rnaseq-flow Software Versions'\n")
    fh.write("description: 'Collected at run time from the tools each process actually invoked.'\n")
    fh.write("plot_type: 'html'\n")
    fh.write("data: |\n")
    for line in html:
        fh.write("    " + line + "\n")

print(f"Software versions: {len(rows)} tool(s) across "
      f"{len({r[0] for r in rows})} process(es)")
PY
fi

# An empty mapping means the sample ids already match the filenames, in which
# case MultiQC is left to its own devices.
RENAME=""
if [ -s replace_names.tsv ]; then
    RENAME="--replace-names replace_names.tsv"
    echo "MultiQC sample-name mapping (\$(wc -l < replace_names.tsv) entries):"
    cat replace_names.tsv
fi

multiqc -f . \$RENAME $args

# Flush-left, and <<END_VERSIONS rather than <<-END_VERSIONS. Nextflow strips
# the *common* leading whitespace from a script block; the Python heredoc above
# has to sit at column 0, so nothing is stripped from the rest, and a
# space-indented terminator would never match (<<- strips tabs, not spaces).
# That is exactly what happened: the heredoc ran to EOF and versions.yml came
# out indented with a literal END_VERSIONS line in it.
cat <<END_VERSIONS > versions.yml
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
