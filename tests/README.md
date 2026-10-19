# Developer checks

Static checks that run without containers or data. The runnable yeast example
for users is in [examples/yeast/](../examples/yeast/README.md).

## Static checks (`check_report.py`)

Run before a pipeline run, not after. No containers, no data, a second or two:

```bash
python3 tests/check_report.py
```

It catches five classes of mistake, each of which has reached a real run at
least once. All of them share a property that makes them expensive: they surface
at the *end* of a pipeline, or not at all.

| Check | Why |
|---|---|
| Non-ASCII in emitted R code | The rendering container runs in the C locale, so R escapes any non-ASCII byte to `<U+XXXX>` on output. Not confined to raw HTML: it reached a `cat()` call, a DT table caption and a plotly facet label |
| R chunk syntax | At R's top level a statement ends at the newline after an `if` body, so an `else` starting the next line never parses. Broke report rendering twice |
| Raw-HTML balance | An unbalanced brace or `<style>` tag mangles the page silently rather than erroring |
| Chunk `eval=` gates | A gate naming an undefined variable makes its chunk silently never run |
| Module heredocs | Nextflow strips a script block's *common* leading whitespace, so one line at column 0 leaves `<<-END_VERSIONS` unable to match its terminator. `versions.yml` comes out malformed, the aggregated `software_versions.yml` becomes invalid YAML, and MULTIQC fails a long way from the cause. Shipped twice |

It covers `assets/analysis_report.qmd`, every `assets/*.R`, and all the modules
in `modules/local/`. Exits non-zero on failure, so it works as a commit or CI
gate. `Rscript` is needed for the chunk-parsing check; its absence is reported
as a failure rather than skipped, so a green run means all five checks ran.

Pair it with `nextflow lint .`: between them they cover the static failures
that would otherwise cost a full run to discover.
