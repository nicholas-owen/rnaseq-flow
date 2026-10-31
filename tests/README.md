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

## Channel topology (`validate_topology.py`)

```bash
python3 tests/validate_topology.py
```

Not a Nextflow parser. It checks five mechanical things a reading eye misses:
every process call gets as many positional arguments as the process declares
inputs, every `<PROCESS>.out.<name>` names something the process emits, every
process has a `stub:` block, every label a module references exists in
`conf/base.config`, and every `includeConfig` target is on disk. Non-zero exit
on any error.

It takes the repository root as an optional argument and otherwise works it out
from its own location, so it runs from any directory.

This file was absent from the repository until 2026-10-30. `.gitignore` carried
an unanchored `validate*.py`, which matched it at any depth, so the one check
that would have caught a wrong call arity was in nobody's clone. The rule is now
`/validate*.py`, which ignores scratch validators at the root and nothing else.

## Unit test (`unit_rmats_libtype.nf`)

```bash
nextflow run ./tests/unit_rmats_libtype.nf
```

Seven assertions against `rmatsLibType()`, which collapses the per-sample RSeQC
verdicts to the one `--libType` rMATS takes for a whole cohort. Nextflow is
needed to resolve the include; nothing else is, since the function is pure and
takes a plain list.

Keep the `./`. A script path containing a slash is read as a `<owner>/<repo>`
project name, so `nextflow run tests/unit_rmats_libtype.nf` tries to pull
`nextflow-io/tests` from GitHub.

Six `no process matching config selector` warnings are expected. The selectors
are in `conf/modules.config` and this script includes no process at all.

The mapping is worth asserting because getting it the wrong way round is
silent: rMATS accepts any `--libType` and produces a full set of PSI values
either way. The three error paths are worth asserting because no dataset here
reaches them, the yeast cohort agreeing on `reverse` and GSE151251 being
unstranded.

## Stub run (`stub_fixtures.sh` and CI)

```bash
bash tests/stub_fixtures.sh /tmp/fixtures
nextflow run main.nf -profile test -stub-run \
    --input /tmp/fixtures/samplesheet.csv --aligner star \
    --star_index /tmp/fixtures/star_index \
    --gtf /tmp/fixtures/genome.gtf.gz --outdir /tmp/out
```

A stub run executes each process's `stub:` block, so it exercises the whole
channel graph in seconds with no tool, no container and no data. The fixture
script builds the empty inputs Nextflow still insists on staging.

`.github/workflows/ci.yml` runs this on both the genome and pseudo-aligner
routes and then asserts each result directory exists. The assertions carry the
value: a silently skipped analysis is an absent directory and nothing in the
log, which is how a release shipped with `--aligner star --gtf`-less runs
completing green and empty. Five negative cases sit beside them: STAR refusing
a GTF-less run by name, a run with no `--gtf` saying so and skipping, an
analysis dropped by `--stop_at` saying so, and a condition or a `--design`
covariate value `make.names()` would rewrite being refused before any work
begins.

The job sets `defaults.run.shell: bash` for `pipefail`. Without it a step's
`nextflow ... | tee log` reports `tee`'s exit status, which is always zero, so
a negative test fires its failure branch on a pipeline that refused exactly as
intended.

Deliberately no container engine in that job. A stub run still containerises
its tasks when an engine is enabled, so it would pull around thirty images to
run `mkdir` and `touch`; `checkParameters()` exempts stub runs from the engine
requirement for this reason.
