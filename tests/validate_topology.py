#!/usr/bin/env python3
"""
Static topology validator for the rnaseq-flow Nextflow DSL2 pipeline.

It is NOT a Nextflow parser. It mechanically checks the things most likely to
break a run that a human review can miss:

  1. Every process invoked in a workflow gets the same number of positional
     arguments as the process declares input channels.
  2. Every `<PROCESS>.out.<name>` reference uses a name the process `emit:`s.
  3. Every process has a `stub:` block (needed for -stub-run validation).
  4. Process labels referenced by modules exist in conf/base.config.
  5. conf files referenced by includeConfig exist.

Exit code is non-zero if any ERROR is found.
"""
import os, re, sys

# Defaults to the repository root, the parent of tests/, so the checker works
# from any working directory rather than only from the root. An explicit path
# still wins, which is how it was called when it lived at the root itself.
ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MOD_DIR = os.path.join(ROOT, "modules", "local")
WF_DIR  = os.path.join(ROOT, "workflows")

errors, warnings = [], []
def err(m):  errors.append(m)
def warn(m): warnings.append(m)

# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def strip_comments(text):
    """Remove // line comments and /* */ block comments (good enough here)."""
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.S)
    out = []
    for line in text.split('\n'):
        # drop // comments but keep '//' that are inside quotes (rare here)
        in_s = None
        res = ''
        i = 0
        while i < len(line):
            c = line[i]
            if in_s:
                res += c
                if c == in_s:
                    in_s = None
            elif c in ('"', "'"):
                in_s = c
                res += c
            elif c == '/' and i + 1 < len(line) and line[i+1] == '/':
                break
            else:
                res += c
            i += 1
        out.append(res)
    return '\n'.join(out)

def split_top_args(s):
    """Split a process-call argument string on top-level commas."""
    args, depth, cur = [], 0, ''
    for c in s:
        if c in '([{':
            depth += 1; cur += c
        elif c in ')]}':
            depth -= 1; cur += c
        elif c == ',' and depth == 0:
            args.append(cur.strip()); cur = ''
        else:
            cur += c
    if cur.strip():
        args.append(cur.strip())
    return args

def match_paren(text, open_idx):
    """Given index of '(', return index of the matching ')'."""
    depth = 0
    for i in range(open_idx, len(text)):
        if text[i] == '(':
            depth += 1
        elif text[i] == ')':
            depth -= 1
            if depth == 0:
                return i
    return -1

# ----------------------------------------------------------------------
# Parse modules
# ----------------------------------------------------------------------
INPUT_KW = ('tuple', 'path', 'val', 'each', 'env', 'stdin')
processes = {}   # name -> {'inputs': int, 'emits': set, 'stub': bool, 'file': str}

for fn in sorted(os.listdir(MOD_DIR)):
    if not fn.endswith('.nf'):
        continue
    raw = open(os.path.join(MOD_DIR, fn)).read()
    text = strip_comments(raw)
    m = re.search(r'process\s+([A-Z][A-Z0-9_]*)\s*\{', text)
    if not m:
        warn(f"{fn}: no process block found")
        continue
    name = m.group(1)

    # input block: from 'input:' to the next section keyword
    inp = re.search(r'\binput:\s*(.*?)\n\s*(output:|script:|shell:|exec:|stub:)',
                     text, flags=re.S)
    n_inputs = 0
    if inp:
        for line in inp.group(1).split('\n'):
            t = line.strip()
            if t and t.split('(')[0].split()[0] if t else False:
                kw = t.split('(')[0].split()[0]
                if kw in INPUT_KW:
                    n_inputs += 1

    # output block
    out = re.search(r'\boutput:\s*(.*?)\n\s*(script:|shell:|exec:|stub:)',
                    text, flags=re.S)
    emits = set()
    if out:
        for em in re.finditer(r'emit:\s*([A-Za-z_][A-Za-z0-9_]*)', out.group(1)):
            emits.add(em.group(1))

    has_stub = bool(re.search(r'\n\s*stub:\s*\n', text))
    label    = re.findall(r"label\s+'([^']+)'", text)

    processes[name] = {'inputs': n_inputs, 'emits': emits, 'stub': has_stub,
                       'file': fn, 'labels': label}

print(f"Parsed {len(processes)} process modules.\n")

# ----------------------------------------------------------------------
# Check 3: stub blocks present
# ----------------------------------------------------------------------
for name, info in sorted(processes.items()):
    if not info['stub']:
        err(f"{info['file']}: process {name} has no stub: block")

# ----------------------------------------------------------------------
# Check 4: labels exist in conf/base.config
# ----------------------------------------------------------------------
base_cfg = os.path.join(ROOT, 'conf', 'base.config')
known_labels = set()
if os.path.exists(base_cfg):
    btext = open(base_cfg).read()
    known_labels = set(re.findall(r'withLabel:\s*([A-Za-z0-9_]+)', btext))
else:
    err("conf/base.config not found")
for name, info in processes.items():
    for lab in info['labels']:
        if lab not in known_labels:
            err(f"{info['file']}: process {name} uses label '{lab}' "
                f"not defined in conf/base.config")

# ----------------------------------------------------------------------
# Check 5: includeConfig targets exist
# ----------------------------------------------------------------------
nf_cfg = os.path.join(ROOT, 'nextflow.config')
if os.path.exists(nf_cfg):
    for inc in re.findall(r"includeConfig\s+'([^']+)'", open(nf_cfg).read()):
        if not os.path.exists(os.path.join(ROOT, inc)):
            err(f"nextflow.config: includeConfig '{inc}' -> file missing")

# ----------------------------------------------------------------------
# Checks 1 & 2: process calls and .out references in workflow files
# ----------------------------------------------------------------------
wf_files = [os.path.join(WF_DIR, f) for f in os.listdir(WF_DIR) if f.endswith('.nf')]
wf_files.append(os.path.join(ROOT, 'main.nf'))

for wf in wf_files:
    if not os.path.exists(wf):
        continue
    text = strip_comments(open(wf).read())
    rel = os.path.relpath(wf, ROOT)

    # process invocations: NAME ( ... )
    for m in re.finditer(r'\b([A-Z][A-Z0-9_]+)\s*\(', text):
        name = m.group(1)
        if name not in processes:
            continue   # workflow call (RNASEQ/DOWNLOAD/BUILD_INDICES) or Groovy
        open_idx = text.index('(', m.start())
        close_idx = match_paren(text, open_idx)
        if close_idx < 0:
            err(f"{rel}: unbalanced parentheses in call to {name}")
            continue
        arg_str = text[open_idx+1:close_idx]
        n_args = len(split_top_args(arg_str))
        expected = processes[name]['inputs']
        if n_args != expected:
            err(f"{rel}: {name}() called with {n_args} args but module "
                f"declares {expected} input channel(s)")
        else:
            print(f"  OK  {rel}: {name}({n_args} args) == {expected} inputs")

    # .out.<emit> references
    for m in re.finditer(r'\b([A-Z][A-Z0-9_]+)\.out\.([A-Za-z0-9_]+)', text):
        name, emit = m.group(1), m.group(2)
        if name not in processes:
            continue
        if emit not in processes[name]['emits']:
            err(f"{rel}: {name}.out.{emit} -> '{emit}' is not an emit of "
                f"{name} (emits: {sorted(processes[name]['emits'])})")

# ----------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------
print()
for w in warnings:
    print(f"WARNING: {w}")
for e in errors:
    print(f"ERROR:   {e}")
print()
if errors:
    print(f"FAILED: {len(errors)} error(s), {len(warnings)} warning(s)")
    sys.exit(1)
print(f"PASSED: 0 errors, {len(warnings)} warning(s)")
