#!/usr/bin/env python3
"""Static checks for the Quarto analysis report and the R analysis scripts.

    python3 assets/test/check_report.py

Exits non-zero if anything fails, so it can gate a commit or a CI job.

Why this exists: the report is rendered by the very last process of a pipeline
run, inside a container. A mistake in it is not visible until every expensive
stage has already completed, and some mistakes do not fail the render at all --
they quietly produce a broken document. Each check below corresponds to a bug
that actually shipped.

  1. Non-ASCII in emitted code.  The rendering container runs in the C locale,
     where R escapes any non-ASCII byte to "<U+XXXX>" on output. This is not
     specific to raw HTML: it hit a cat() call, a DT table caption and a plotly
     facet label, so anything R emits has to stay ASCII. Use an HTML entity
     (&mdash;) where the output is HTML, and a plain ASCII character where it is
     text drawn inside a plot.

  2. R chunk syntax.  At R's top level a statement ends at the newline after an
     `if` body, so an `else` that starts the next line never parses. This broke
     rendering twice, and only at the end of a full run.

  3. Raw-HTML block balance.  The report carries one ```{=html} block holding
     its CSS and a small script. An unbalanced brace or <style> tag silently
     mangles the page rather than erroring.

  4. Chunk gates.  Every `eval=` option must name a variable the document
     defines, or the chunk silently never runs.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
QMD = os.path.join(ROOT, "assets", "analysis_report.qmd")
R_SCRIPTS = os.path.join(ROOT, "assets")

failures = []


def report(label, ok, detail=""):
    print(f"  {'PASS' if ok else 'FAIL'}  {label}" + (f"  {detail}" if detail and not ok else ""))
    if not ok:
        failures.append(label)


def emitted_code(line):
    """The part of a line that R will execute, i.e. not a whole-line comment.

    Deliberately crude: a '#' inside a string literal would truncate the line
    early, which can only ever *miss* a problem, never invent one. The check is
    a guard, not a parser.
    """
    stripped = line.lstrip()
    if stripped.startswith("#"):
        return ""
    return line.split("#", 1)[0]


def check_non_ascii(path, chunks=None):
    """Non-ASCII inside code that R will emit. Comments are fine -- they only
    ever appear in the folded source view."""
    bad = []
    if chunks is None:
        blocks = [("", open(path, encoding="utf-8").read())]
    else:
        blocks = chunks
    for name, body in blocks:
        for i, line in enumerate(body.split("\n"), 1):
            chars = sorted({c for c in emitted_code(line) if ord(c) >= 128})
            if chars:
                where = f"{name}:{i}" if name else f"line {i}"
                codes = " ".join(f"U+{ord(c):04X} ({c})" for c in chars)
                bad.append((where, codes, line.strip()[:88]))
    return bad


def main():
    src = open(QMD, encoding="utf-8").read()
    chunks = [((o.strip().split(",")[0] or f"chunk{n}").strip(), b)
              for n, (o, b) in enumerate(
                  re.findall(r"^```\{r([^}]*)\}\n(.*?)^```", src, re.S | re.M), 1)]
    raw = re.findall(r"^```\{=html\}\n(.*?)^```", src, re.S | re.M)

    print(f"analysis_report.qmd: {len(chunks)} R chunk(s), {len(raw)} raw-html block(s)\n")

    # ---- 1. non-ASCII ------------------------------------------------------
    print("non-ASCII in emitted code")
    bad = check_non_ascii(QMD, chunks)
    for f in sorted(os.listdir(R_SCRIPTS)):
        if f.endswith(".R"):
            p = os.path.join(R_SCRIPTS, f)
            bad += [(f"{f}:{w.split()[-1]}", c, l)
                    for w, c, l in check_non_ascii(p)]
    report("no non-ASCII reaches R output", not bad)
    for where, codes, line in bad:
        print(f"        {where}  {codes}")
        print(f"          {line}")

    # ---- 2. chunk syntax ---------------------------------------------------
    print("\nR chunk syntax")
    rscript = shutil.which("Rscript")
    if not rscript:
        report("Rscript available to parse chunks", False,
               "-- install R, or run this where R is present")
    else:
        broken = []
        for name, body in chunks:
            with tempfile.NamedTemporaryFile("w", suffix=".R", delete=False,
                                             encoding="utf-8") as fh:
                fh.write(body)
                tmp = fh.name
            r = subprocess.run([rscript, "-e", f'invisible(parse("{tmp}"))'],
                               capture_output=True, text=True)
            os.unlink(tmp)
            if r.returncode != 0:
                broken.append((name, r.stderr.strip().splitlines()[:3]))
        report(f"all {len(chunks)} chunks parse", not broken)
        for name, err in broken:
            print(f"        {name}")
            for line in err:
                print(f"          {line}")

    # ---- 3. raw-html balance ----------------------------------------------
    print("\nraw HTML block")
    html = "".join(raw)
    report("braces balanced", html.count("{") == html.count("}"),
           f"-- {html.count('{')} open / {html.count('}')} close")
    report("<style> tags balanced", html.count("<style>") == html.count("</style>"))
    report("<script> tags balanced", html.count("<script>") == html.count("</script>"))

    # ---- 4. chunk gates ----------------------------------------------------
    print("\nchunk eval= gates")
    gates = set(re.findall(r"eval\s*=\s*!?\(?([A-Za-z_][A-Za-z0-9_.]*)", src))
    defined = set(re.findall(r"^\s*([A-Za-z_][A-Za-z0-9_.]*)\s*<-", src, re.M))
    missing = sorted(g for g in gates if g not in defined and g not in {"TRUE", "FALSE"})
    report(f"all {len(gates)} gates resolve to a defined variable", not missing,
           f"-- undefined: {missing}")

    print()
    if failures:
        print(f"{len(failures)} CHECK(S) FAILED: {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
