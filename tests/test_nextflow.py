#!/usr/bin/env python3
"""
Static checks on the Nextflow layer.

Nextflow cannot be installed where this pipeline was developed, so the
workflow is the one part of it that has never been executed. That is a
real gap and this does not close it - but most of what goes wrong in an
unrun workflow is not subtle. It is a module passing `--foo` to a tool
whose flag is `--bar`, an `out.thing` that no process emits, a `params.x`
nobody declared, or an include pointing at a file that was renamed. Every
one of those is decidable without running anything, and every one of them
would otherwise surface as a failed job on the farm an hour into a run.

So each module's command line is checked against the tool's real
argparse, taken from `--help` rather than from a list kept in step by
hand. Each `.out.name` is checked against the emits of the process it
names. Each parameter is checked against nextflow.config.

One cardinality question IS decidable statically and is checked, because
it is the only failure here that corrupts results instead of stopping
the run: pairing a per-sample channel with a single-item queue channel
silently reduces the fan-out to one task. See
test_fanout_is_not_capped.

What this still cannot tell you: whether a `collect()` blocks, whether
the DAG deadlocks, or whether a join drops keys. Those need Nextflow.
The checklist in the README says what to verify on the first real run.

Run: python3 tests/test_nextflow.py
"""

from __future__ import annotations

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BIN = os.path.join(ROOT, "bin")

PASS, FAIL = [], []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def nf_files() -> list[str]:
    out = []
    for base, _dirs, files in os.walk(ROOT):
        parts = base.split(os.sep)
        # Nextflow's scratch is a directory literally named "work"; the
        # substring test that used to be here also excluded
        # "subworkflows", so every subworkflow went unchecked.
        if "work" in parts or ".nextflow" in parts:
            continue
        out += [os.path.join(base, f) for f in files if f.endswith(".nf")]
    return sorted(out)


def tool_flags(tool: str) -> set[str]:
    """The flags a tool really accepts, asked of the tool itself."""
    try:
        r = subprocess.run([sys.executable, os.path.join(BIN, tool), "--help"],
                           capture_output=True, text=True, timeout=60)
    except Exception:
        return set()
    return set(re.findall(r"(--[a-z0-9][a-z0-9-]*)", r.stdout))


def strip_comments(text: str) -> str:
    """Comments mention flags too, and a comment cannot break a run.

    Scanning them produced a false alarm on the first pass: a note about
    `--input_dir` matched as `--input`, because the underscore ends the
    flag pattern. A checker that cries wolf gets switched off.
    """
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"//[^\n]*", "", text)


def processes() -> dict[str, dict]:
    """Every process, with its emits and the script it runs."""
    out: dict[str, dict] = {}
    for path in nf_files():
        text = open(path, encoding="utf-8").read()
        for m in re.finditer(r"process\s+(\w+)\s*\{", text):
            name = m.group(1)
            body = text[m.end():]
            depth, end = 1, 0
            for i, c in enumerate(body):
                depth += (c == "{") - (c == "}")
                if depth == 0:
                    end = i
                    break
            body = body[:end]
            out[name] = {
                "file": path,
                "emits": set(re.findall(r"emit:\s*(\w+)", body)),
                "outputs": re.findall(r'path\s+"([^"]+)"', body),
                "body": strip_comments(body),
            }
    return out


def input_arity(body: str) -> int:
    """How many things a process declares in its input block.

    One declaration per line, and a tuple counts once however many
    elements it holds - which is how Nextflow counts them at the call
    site too.
    """
    m = re.search(r"\n\s*input:\s*\n(.*?)\n\s*(output:|script:|stub:|exec:)",
                  body, re.S)
    if not m:
        return 0
    n = 0
    for line in m.group(1).splitlines():
        line = line.strip()
        if line and not line.startswith("//"):
            n += 1
    return n


def call_arity(text: str, name: str) -> list[int]:
    """How many arguments each call site passes to `name`.

    Counts commas at bracket depth zero, so a nested call or a closure
    containing commas is one argument, not several.
    """
    out = []
    for m in re.finditer(rf"(?<![\w.]){name}\s*\(", text):
        i = m.end()
        depth, args = 1, 0
        # Groovy permits a trailing comma before the closing bracket, so
        # commas cannot simply be counted: an argument exists only if
        # something non-blank followed the comma before it.
        current = False
        while i < len(text) and depth:
            c = text[i]
            if c in "([{":
                depth += 1
                current = True
            elif c in ")]}":
                depth -= 1
                if depth == 0:
                    break
                current = True
            elif c == "," and depth == 1:
                if current:
                    args += 1
                current = False
            elif not c.isspace():
                current = True
            i += 1
        if current:
            args += 1
        out.append(args)
    return out


def test_call_arity() -> None:
    """A process must be called with as many arguments as it declares.

    This is the failure that adding an input causes, and the one static
    check that would have caught it: the pipeline parses fine, the flags
    all exist, the channels all resolve, and it dies at runtime with a
    message about the wrong number of arguments. Since the Nextflow layer
    cannot be executed in development here, arity is checked rather than
    discovered.
    """
    print("\nprocess call arity")
    procs = processes()
    text = ""
    for path in nf_files():
        with open(path, encoding="utf-8") as fh:
            body = fh.read()
        # A process's own definition is not a call of it.
        body = re.sub(r"process\s+\w+\s*\{", "", body)
        text += strip_comments(body) + "\n"

    bad = []
    for name, info in sorted(procs.items()):
        want = input_arity(info["body"])
        for got in call_arity(text, name):
            if got != want:
                bad.append(f"{name}: declares {want}, called with {got}")
    check("every process is called with the arity it declares",
          not bad, "; ".join(bad))

    # And the check itself must be able to fail, or it proves nothing.
    fake = "SYNTHESISE(a, b)"
    check("the arity check would notice a missing argument",
          call_arity(fake, "SYNTHESISE") == [2]
          and input_arity(procs["SYNTHESISE"]["body"]) != 2,
          f"{call_arity(fake, 'SYNTHESISE')} vs "
          f"{input_arity(procs['SYNTHESISE']['body'])}")


def test_placeholders_are_distinct() -> None:
    """Two absent optional inputs must not stage under one filename.

    Nextflow stages every input of a process into a single directory, so
    if two optional inputs both fall back to `assets/NO_FILE` the run
    dies before the first tool starts:

        Process `PREFLIGHT` input file name collision -- There are
        multiple input files for each of the following file names: NO_FILE

    It only fires when two optionals are absent AT ONCE, which is why it
    survived every static check and the first stub run attempt. Each
    optional input therefore gets its own placeholder, and this asserts
    they stay distinct.
    """
    print("\noptional-input placeholders")
    text = ""
    for path in nf_files():
        with open(path, encoding="utf-8") as fh:
            text += strip_comments(fh.read()) + "\n"

    used = re.findall(r"assets/(NO_[A-Z0-9_]*)", text)
    dupes = sorted({p for p in used if used.count(p) > 1})
    # NO_RUNS legitimately appears twice: the manifest is read once for
    # pre-flight and once for the nanopore branch, and they are the same
    # file, so no collision is possible.
    dupes = [d for d in dupes if d != "NO_RUNS"]
    check("no placeholder is reused across different optional inputs",
          not dupes,
          f"{dupes} appear more than once - two absent inputs would collide")

    missing = [p for p in sorted(set(used))
               if not os.path.exists(os.path.join(ROOT, "assets", p))]
    check("every placeholder referenced actually exists",
          not missing, str(missing))

    # A module that tests for one exact placeholder name would silently
    # start passing the wrong thing the moment a new one is introduced.
    exact = re.findall(r"\.name\s*!=\s*'NO_FILE'", text)
    check("modules recognise placeholders by prefix, not one exact name",
          not exact, f"{len(exact)} module(s) still compare against NO_FILE")


def test_parsed_stubs_have_content() -> None:
    """A stub output that gets parsed must be structurally valid.

    `-stub-run` exists to prove the DAG connects, and it cannot do that
    if a placeholder breaks the first thing that reads it. An empty
    file_resolution.tsv passes `touch`, passes every other check here,
    and then fails at runtime with `Missing 'header' in CSV file` - and
    because the same file feeds a `.first()`, an empty one would instead
    leave a later process waiting on a channel that never emits.

    So: any emit consumed by splitCsv must be written with real content
    in the stub, not merely touched.
    """
    print("\nstub outputs that are parsed downstream")
    procs = processes()
    text = ""
    for path in nf_files():
        with open(path, encoding="utf-8") as fh:
            text += strip_comments(fh.read()) + "\n"

    # Which emits are fed into splitCsv anywhere?
    parsed = set(re.findall(r"\.out\.(\w+)\s*\n?\s*\.splitCsv", text))
    parsed |= set(re.findall(r"\.out\.(\w+)\.splitCsv", text))
    check("at least one parsed emit was found, so this test is live",
          bool(parsed), "no .out.X.splitCsv found - has the wiring changed?")

    bad = []
    for name, info in sorted(procs.items()):
        body = info["body"]
        out_block = re.search(r"\n\s*output:(.*?)\n\s*(script:|stub:|exec:)",
                              body, re.S)
        stub = re.search(r"\n\s*stub:\s*\n(.*)$", body, re.S)
        if not out_block or not stub:
            continue
        for emit in parsed:
            m = re.search(rf'path\s+"([^"]+)"\s*,?\s*emit:\s*{emit}\b',
                          out_block.group(1))
            if not m:
                continue
            fn = m.group(1)
            # A redirect into the file, wherever it sits: the command
            # may be split across lines with a continuation, which an
            # anchored single-line pattern misses - and a check with a
            # false positive gets switched off as fast as one that
            # cries wolf.
            wrote = re.search(rf'>>?\s*\S*{re.escape(fn)}', stub.group(1))
            # A stub that runs the REAL tool needs no redirect: the
            # content is real by construction. PREFLIGHT does exactly
            # that, because a stub cannot read its inputs and this
            # workflow routes on a column of them - a hardcoded
            # file_resolution.tsv decides the platform instead of
            # reflecting it, so -profile test_ont -stub-run ran the
            # Illumina branch. Pre-flight is also cheap enough that
            # faking it saves nothing.
            script = re.search(r"\n\s*script:(.*?)\n\s*stub:", body, re.S)
            tools = set(re.findall(r"(\w+\.py)", script.group(1))) if script else set()
            real = any(re.search(rf"(?<![\w/]){re.escape(t)}\b", stub.group(1))
                       for t in tools)
            if not wrote and not real:
                bad.append(f"{name}.{emit} ({fn}) is only touched")
    check("every stub output that is parsed downstream has real content",
          not bad, "; ".join(bad))

    # The "runs the real tool" exemption is narrow on purpose: a stub
    # that merely mentions a tool in a comment must not qualify, and a
    # stub that touches a parsed file must still fail.
    fake_body = ('\n  output:\n  path "x/f.tsv", emit: file_resolution\n'
                 '  script:\n  """\n  thing.py --out x\n  """\n'
                 '  stub:\n  """\n  touch x/f.tsv\n  """\n')
    m = re.search(r"\n\s*stub:\s*\n(.*)$", fake_body, re.S)
    sc = re.search(r"\n\s*script:(.*?)\n\s*stub:", fake_body, re.S)
    tools = set(re.findall(r"(\w+\.py)", sc.group(1)))
    check("a stub that only touches a parsed file is still rejected",
          not any(re.search(rf"(?<![\w/]){re.escape(t)}\b", m.group(1))
                  for t in tools)
          and not re.search(r">>?\s*\S*f\.tsv", m.group(1)),
          "the exemption is too broad")


def test_publishdir_does_not_double_directories() -> None:
    """publishDir must not repeat a directory the outputs already carry.

    Every tool here writes into a named subdirectory and declares its
    outputs with that path - `synthesis/calls.tsv`, `demux/summary.json`.
    publishDir preserves that relative path, so a publishDir of
    `${params.outdir}/synthesis` produces

        results/synthesis/synthesis/calls.tsv

    The run succeeds, every file exists, and nothing says the layout is
    wrong. It was found by looking for a file and not finding it.
    """
    print("\npublishDir against declared output paths")
    with open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8") as fh:
        cfg = strip_comments(fh.read())

    # Directory prefixes the tools already write into.
    dirs = set()
    for info in processes().values():
        for out in info["outputs"]:
            if "/" in out:
                dirs.add(out.split("/")[0])
    check("output paths with a directory prefix were found",
          bool(dirs), "none - has the output style changed?")

    # Any publishDir path ending in one of those names doubles it.
    bad = []
    for m in re.finditer(r'publishDir\s*=\s*\[[^\]]*?path\s*:\s*([^,\]]+)', cfg):
        expr = m.group(1)
        for d in dirs:
            if re.search(rf'/{re.escape(d)}"\s*\}}?\s*$', expr.strip()):
                bad.append(f"{expr.strip()} duplicates '{d}/'")
    check("no publishDir path repeats a directory the outputs supply",
          not bad, "; ".join(bad))


def test_publishdir_paths_are_lazy() -> None:
    """A publishDir path that reads params must be a closure.

    Nextflow evaluates a plain string when it parses the config, which
    is BEFORE a profile has set params.outdir. So `-profile test`, which
    sets outdir to ./results-test, published most of the run to
    ./results instead - while the base publishDir, already a closure,
    used the right directory. One run, two destinations, no warning.

    A closure is evaluated per task, after profiles are applied.
    """
    print("\npublishDir resolution")
    with open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8") as fh:
        cfg = strip_comments(fh.read())

    # path: "...${params...}..."  - a string, evaluated too early.
    eager = re.findall(r'path\s*:\s*"[^"]*\$\{params\.[^"]*"', cfg)
    check("every publishDir path using params is a closure",
          not eager,
          f"{len(eager)} eager path(s), e.g. {eager[:1]}")

    # And the path must be a closure at all. Asserting on GStrings
    # specifically was wrong: after the publishDir was simplified to
    # `{ params.outdir }` there are no GStrings left, and a test that
    # counts them starts failing for a reason that is not a fault.
    paths = re.findall(r'publishDir\s*=\s*\[[^\]]*?path\s*:\s*([^,\]]+)', cfg)
    check("a publishDir path was found, so this test is live", bool(paths))
    not_closures = [p for p in paths if not p.strip().startswith("{")]
    check("and every publishDir path is a closure",
          not not_closures, str(not_closures))


def test_tools_are_executable() -> None:
    """Nextflow puts bin/ on PATH and runs the tools directly.

    So a tool without the executable bit fails with exit 126 and
    `Permission denied` - after pre-flight has passed and the first
    sample has been staged. Nothing else notices: the shell drivers call
    `python3 tool.py`, the tests import the module, and only Nextflow
    invokes the file as a command.

    A shebang matters for the same reason: without one the kernel has no
    interpreter to hand the file to.
    """
    print("\ntools are executable")
    bad_x, bad_s = [], []
    for f in sorted(os.listdir(BIN)):
        if not (f.endswith(".py") or f.endswith(".sh")):
            continue
        p = os.path.join(BIN, f)
        if not os.access(p, os.X_OK):
            bad_x.append(f)
        with open(p, encoding="utf-8") as fh:
            if not fh.readline().startswith("#!"):
                bad_s.append(f)
    check("every tool in bin/ is executable", not bad_x, str(bad_x))
    check("and every one starts with a shebang", not bad_s, str(bad_s))


def test_flags_exist() -> None:
    """A flag a tool does not have is a job that fails an hour in."""
    print("\nmodule command lines against the tools they call")
    procs = processes()
    for name, p in sorted(procs.items()):
        called = re.findall(r"^\s*(\w+\.py)\s*\\\\\s*$", p["body"], re.M)
        if not called:
            continue
        tool = called[0]
        if not os.path.exists(os.path.join(BIN, tool)):
            check(f"{name}: calls a tool that exists", False, tool)
            continue
        real = tool_flags(tool)
        # only the literal flags; interpolated ones are checked by eye
        used = set(re.findall(r"(--[a-z0-9][a-z0-9-]*)", p["body"]))
        # flags built inside groovy string vars are in the body too, which
        # is what we want - they still have to be real
        bogus = sorted(f for f in used - real
                       if not f.startswith("--version"))
        check(f"{name}: every flag exists in {tool}", not bogus,
              f"not accepted: {', '.join(bogus)}")


def test_emits_resolve() -> None:
    """`.out.thing` where nothing emits `thing` fails at parse time on the
    farm and is silent here."""
    print("\nchannel references resolve to something emitted")
    procs = processes()
    # subworkflow emits too
    for path in nf_files():
        text = open(path, encoding="utf-8").read()
        for m in re.finditer(r"^workflow\s+(\w+)\s*\{(.*?)^\}",
                             text, re.S | re.M):
            wf, body = m.group(1), m.group(2)
            em = re.search(r"emit:(.*?)$", body, re.S)
            if em:
                names = set(re.findall(r"^\s*(\w+)\s*=", em.group(1), re.M))
                names |= set(re.findall(r"^\s*(\w+)\s*$", em.group(1), re.M))
                procs[wf] = {"emits": names, "file": path,
                             "outputs": [], "body": body}

    bad = []
    for path in nf_files():
        text = open(path, encoding="utf-8").read()
        for m in re.finditer(r"\b([A-Z][A-Z0-9_]*)\.out\.(\w+)", text):
            proc, out_name = m.group(1), m.group(2)
            if proc not in procs:
                bad.append(f"{os.path.basename(path)}: {proc} is not defined")
            elif out_name not in procs[proc]["emits"]:
                bad.append(
                    f"{os.path.basename(path)}: {proc}.out.{out_name} "
                    f"(emits: {', '.join(sorted(procs[proc]['emits'])) or 'none'})"
                )
    check("every .out reference names a real emit", not bad,
          "; ".join(bad[:4]))


def test_includes_resolve() -> None:
    print("\nincludes point at files that exist")
    bad = []
    for path in nf_files():
        base = os.path.dirname(path)
        for m in re.finditer(r"include\s*\{[^}]*\}\s*from\s*'([^']+)'",
                             open(path, encoding="utf-8").read()):
            target = m.group(1)
            resolved = os.path.normpath(os.path.join(base, target + ".nf"))
            if not os.path.exists(resolved):
                bad.append(f"{os.path.basename(path)} -> {target}")
    check("every include resolves", not bad, "; ".join(bad))


def test_params_declared() -> None:
    """A params typo reads as null and the flag silently vanishes from the
    command line - the worst kind of failure, because the job succeeds."""
    print("\nparameters used are parameters declared")
    cfg = open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8").read()
    block = re.search(r"params\s*\{(.*?)\n\}", cfg, re.S)
    declared = set(re.findall(r"^\s*(\w+)\s*=", block.group(1), re.M)) \
        if block else set()
    declared |= {"outdir", "publish_dir_mode", "help"}

    used = set()
    for path in nf_files():
        used |= set(re.findall(
            r"params\.(\w+)",
            strip_comments(open(path, encoding="utf-8").read())))
    missing = sorted(used - declared)
    check("every params.x used is declared in nextflow.config",
          not missing, f"undeclared: {', '.join(missing)}")


def test_outputs_are_written() -> None:
    """A declared output the tool never writes fails the task after the
    work is already done."""
    print("\ndeclared outputs against what the tools write")
    procs = processes()
    sources = {f: open(os.path.join(BIN, f), encoding="utf-8").read()
               for f in os.listdir(BIN) if f.endswith(".py")}
    bad = []
    for name, p in sorted(procs.items()):
        called = re.findall(r"^\s*(\w+\.py)\s*\\\\\s*$", p["body"], re.M)
        if not called or called[0] not in sources:
            continue
        src = sources[called[0]]
        for out in p["outputs"]:
            if out in ("versions.yml",) or "*" in out or "${" in out:
                continue
            leaf = os.path.basename(out)
            # A path handed to the tool as an argument is one the tool
            # writes; the name lives in the module, not in the source.
            if out in p["body"] and leaf not in src and "--" in p["body"]:
                if re.search(re.escape(out) + r"\s*(\\|$)", p["body"], re.M):
                    continue
            if leaf not in src:
                bad.append(f"{name}: {out} (not written by {called[0]})")
    check("every fixed output path is one the tool writes", not bad,
          "; ".join(bad[:4]))


def call_args(text: str, name: str) -> list[list[str]]:
    """The argument expressions at each call site of `name`.

    Same bracket walk as call_arity, but keeping the text.
    """
    out = []
    for m in re.finditer(rf"(?<![\w.]){name}\s*\(", text):
        i = m.end()
        depth, start, args = 1, i, []
        while i < len(text) and depth:
            c = text[i]
            if c in "([{":
                depth += 1
            elif c in ")]}":
                depth -= 1
                if depth == 0:
                    break
            elif c == "," and depth == 1:
                args.append(text[start:i].strip())
                start = i + 1
            i += 1
        tail = text[start:i].strip()
        if tail:
            args.append(tail)
        out.append([a for a in args if a])
    return out


# Operators that turn a queue channel into a value channel. A value
# channel is inexhaustible: it is re-read for every item of whatever it
# is paired with.
VALUE_OP = re.compile(r"\.(first|collect|toList|count|last|reduce|sum|min|max)"
                      r"\s*[({]")
LITERAL = re.compile(r"^(\[|\"|'|file\s*\(|params\.|Channel\.value\s*\()")


def assignments(text: str) -> dict[str, str]:
    """Every `name = expr` in a workflow body, expr as written."""
    out: dict[str, str] = {}
    for m in re.finditer(r"^\s{4}(\w+)\s*=\s*", text, re.M):
        i, depth = m.end(), 0
        while i < len(text):
            c = text[i]
            if c in "([{":
                depth += 1
            elif c in ")]}":
                depth -= 1
            elif c == "\n" and depth == 0:
                # A statement continues only while a bracket is open or
                # the next line starts with an operator or a dot.
                nxt = text[i + 1:i + 200].lstrip()
                if not nxt.startswith((".", "?", ":", "+")):
                    break
            i += 1
        out[m.group(1)] = text[m.end():i].strip()
    return out


def is_value(expr: str, defs: dict[str, str], seen: frozenset = frozenset(),
             value_out: frozenset = frozenset()) -> bool:
    """Would this expression give a value channel (or a plain literal)?

    value_out names processes whose outputs are value channels. A
    process invoked with nothing but value channels emits value
    channels - Nextflow's rule, and the reason `.first()` on such an
    output draws a "the operator first is useless" warning. The rule is
    encoded here so the workflow does not have to carry a defensive
    .first() that Nextflow objects to.
    """
    expr = expr.strip()
    if LITERAL.match(expr):
        return True
    m = re.match(r"^(\w+)\.out\b", expr)
    if m and m.group(1) in value_out and not re.search(
            r"\.(flatten|flatMap|transpose|splitCsv|splitText|mix)\s*[({]", expr):
        return True
    if VALUE_OP.search(expr):
        # ...and nothing after it turns it back into a queue. flatMap,
        # splitCsv, transpose and friends all do.
        tail = expr[VALUE_OP.search(expr).end():]
        return not re.search(r"\.(flatten|flatMap|transpose|splitCsv|splitText|mix)\s*[({]",
                             tail)
    if expr in defs and expr not in seen:
        return is_value(defs[expr], defs, seen | {expr}, value_out)
    return False


def workflow_bodies() -> dict[str, str]:
    """Each workflow (main and subworkflows) with its body and take list."""
    out = {}
    for path in nf_files():
        text = strip_comments(open(path, encoding="utf-8").read())
        for m in re.finditer(r"workflow\s*(\w*)\s*\{", text):
            body = text[m.end():]
            depth, end = 1, len(body)
            for i, c in enumerate(body):
                depth += (c == "{") - (c == "}")
                if depth == 0:
                    end = i
                    break
            out[m.group(1) or "MAIN"] = body[:end]
    return out


def test_fanout_is_not_capped() -> None:
    """A per-sample process must not be paired with a single-item queue.

    Nextflow pairs a process's inputs positionally and stops at the
    SHORTEST channel. Channel.fromPath on one file, and any process
    output, are queue channels holding one item - so pairing 8 samples
    with one reference file yields one task, not eight, and Nextflow
    reports no error because that is its documented behaviour. The
    missing samples simply never appear, and synthesis reports
    confidently on the ones that did.

    This is the only bug class found so far that corrupts results
    rather than stopping the run, which is why it is checked
    statically: `[0%] 0 of 1` in the log is the sole symptom, and it
    looks like a small test run.
    """
    print("\nno per-sample process is capped by a single-item queue channel")
    procs = processes()

    # A process fans out over samples when its first input is a tuple
    # carrying meta. Everything after that must be a value channel.
    fanout = {}
    for name, p in procs.items():
        m = re.search(r"\n\s*input:\s*\n(.*?)\n\s*(output:|script:|stub:|exec:)",
                      p["body"], re.S)
        if not m:
            continue
        lines = [l.strip() for l in m.group(1).splitlines() if l.strip()]
        if lines and re.match(r"tuple\s+val\s*\(\s*meta", lines[0]):
            fanout[name] = len(lines)

    check("processes that fan out over samples were found", bool(fanout),
          str(sorted(fanout)))

    bodies = workflow_bodies()
    # A subworkflow parameter inherits its requirement from how the
    # subworkflow uses it: passed to a fan-out process at any position
    # after the first, it must be a value channel in every caller.
    sub_scalar: dict[str, set[int]] = {}
    sub_takes: dict[str, list[str]] = {}
    for wf, body in bodies.items():
        m = re.search(r"take:\s*\n(.*?)\n\s*main:", body, re.S)
        if not m:
            continue
        takes = [l.strip().split()[0] for l in m.group(1).splitlines() if l.strip()]
        sub_takes[wf] = takes
        need = set()
        for callee in fanout:
            for args in call_args(body, callee):
                for pos, a in enumerate(args):
                    if pos and a in takes:
                        need.add(takes.index(a))
        if need:
            sub_scalar[wf] = need

    # Which processes emit value channels. Computed to a fixed point,
    # because one all-value process feeding another makes the second
    # all-value too. Conservative: an argument that cannot be PROVEN to
    # be a value channel disqualifies the process, so the worst case is
    # demanding a .first() that was not strictly needed - never
    # excusing one that was.
    all_bodies = {w: (b, assignments(b)) for w, b in bodies.items()}
    value_out: set[str] = set()
    for _ in range(4):
        grew = False
        for name in procs:
            if name in value_out:
                continue
            calls = [a for _w, (b, _d) in all_bodies.items()
                     for a in call_args(b, name)]
            if not calls:
                continue
            ok = all(is_value(arg, d, frozenset(), frozenset(value_out))
                     for _w, (b, d) in all_bodies.items()
                     for args in call_args(b, name) for arg in args)
            if ok:
                value_out.add(name)
                grew = True
        if not grew:
            break
    check("processes emitting value channels were identified",
          "PREFLIGHT" in value_out,
          f"value_out={sorted(value_out)}")

    bad = []
    verified: set[str] = set()
    for wf, body in bodies.items():
        defs = assignments(body)
        # A subworkflow's own take parameters are not judged here: the
        # obligation travels outward to every caller, which is checked
        # below via sub_scalar. Judging them in place would demand a
        # .first() inside the subworkflow that the caller has already
        # applied.
        exempt = set(sub_takes.get(wf, []))
        targets = {n: set(range(1, k)) for n, k in fanout.items()}
        targets.update(sub_scalar)
        for callee, positions in targets.items():
            if callee == wf:
                continue
            for args in call_args(body, callee):
                if callee in sub_scalar:
                    verified.add(callee)
                for pos in positions:
                    if pos >= len(args):
                        continue
                    if args[pos] in exempt:
                        continue
                    if not is_value(args[pos], defs,
                                    value_out=frozenset(value_out)):
                        bad.append(f"{wf}: {callee} arg {pos + 1} "
                                   f"`{args[pos][:40]}` is a queue channel")

    check("every scalar input to a fan-out process is a value channel",
          not bad, "; ".join(bad[:6]))

    # The exemption above is only sound if every subworkflow carrying a
    # deferred obligation is actually called somewhere. An uncalled one
    # would have its take parameters excused and never checked.
    unchecked = sorted(set(sub_scalar) - verified)
    check("every subworkflow with a deferred obligation has a caller",
          not unchecked, ", ".join(unchecked))

    # Liveness: the check must actually reject the shape it exists to
    # catch. A vacuous assertion here would have let the original bug
    # through unchanged.
    probe = {"ch_q": "Channel.fromPath(params.reference)",
             "ch_v": "Channel.fromPath(params.reference).first()"}
    check("and it rejects a bare fromPath", not is_value("ch_q", probe))
    check("and accepts one with .first()", is_value("ch_v", probe))
    check("and rejects a bare process output",
          not is_value("PREFLIGHT.out.insert_fasta", {}))
    check("and sees through splitCsv after a value op",
          not is_value("X.collect().splitCsv()", {}))


def test_labels_are_sized_everywhere() -> None:
    """Every resource label must be sized in every profile that resizes.

    The `local` and `stub` profiles exist to shrink jobs to a laptop or
    to nothing. They do that with withLabel blocks, and a withLabel
    selector takes precedence over the generic process block - so a
    label they forget keeps the FARM sizing. That has already happened
    twice here, once as "CPU 6 > 4 available" and once as a stub task
    asking for 24 GB. Both surfaced only at run time.
    """
    print("\nresource labels are sized in every profile")
    cfg = open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8").read()

    used = set()
    for path in nf_files():
        text = strip_comments(open(path, encoding="utf-8").read())
        used |= set(re.findall(r"label\s+'([a-z_]+)'", text))
    check("modules declare resource labels", bool(used), str(sorted(used)))

    def block(name: str) -> str:
        """The body of a named profile, by brace matching."""
        m = re.search(rf"\n    {name}\s*\{{", cfg)
        if not m:
            return ""
        i, depth = m.end(), 1
        while i < len(cfg) and depth:
            depth += (cfg[i] == "{") - (cfg[i] == "}")
            i += 1
        return cfg[m.end():i]

    # The top-level process block is everything before `profiles {`.
    top = cfg[:cfg.index("\nprofiles {")]
    scopes = {"process (default)": top}
    for prof in ("local", "stub"):
        body = block(prof)
        # Only profiles that resize anything are held to this.
        if "withLabel" in body:
            scopes[f"profile {prof}"] = body

    check("the profiles that resize were found", len(scopes) == 3,
          str(sorted(scopes)))

    bad = []
    for scope, body in scopes.items():
        have = set(re.findall(r"withLabel:\s*([a-z_]+)", body))
        for lab in sorted(used - have):
            bad.append(f"{scope} does not size {lab}")
    check("every label used is sized in every resizing scope", not bad,
          "; ".join(bad))

    # And nothing is sized that no module uses - a label renamed in the
    # modules but left behind here is dead config that looks live.
    stale = []
    for scope, body in scopes.items():
        for lab in sorted(set(re.findall(r"withLabel:\s*([a-z_]+)", body)) - used):
            stale.append(f"{scope} sizes unused {lab}")
    check("and no label is sized that no module declares", not stale,
          "; ".join(stale))


def test_usage_manual_is_current() -> None:
    """Every parameter must be documented, and vice versa.

    A usage manual is the one artefact that rots without anything
    failing: add a param, forget the doc, and the next person reads a
    file that is confidently wrong. Nothing else in this repo would
    notice. So the two are compared - both directions, because a
    documented parameter that no longer exists misleads just as
    effectively as an undocumented one.
    """
    print("\nthe usage manual matches the config")
    doc_path = os.path.join(ROOT, "docs", "USAGE.md")
    check("the usage manual exists", os.path.exists(doc_path), doc_path)
    if not os.path.exists(doc_path):
        return
    doc = open(doc_path, encoding="utf-8").read()
    cfg = open(os.path.join(ROOT, "nextflow.config"), encoding="utf-8").read()

    # The params block is everything up to the first closing brace at
    # column 0.
    block = cfg[cfg.index("params {"):]
    block = block[:block.index("\n}")]
    declared = set(re.findall(r"^\s{4}([a-z_][a-z0-9_]*)\s*=", block, re.M))
    # Internal plumbing, not a knob anyone sets on the command line.
    declared -= {"help"}
    check("parameters were found in the config", len(declared) > 20,
          str(len(declared)))

    documented = set(re.findall(r"`--([a-z_][a-z0-9_]*)`", doc))

    missing = sorted(declared - documented)
    check("every parameter in nextflow.config is in the manual",
          not missing, ", ".join(missing))

    # The reverse: something documented that no longer exists. Allow
    # Nextflow's own flags, which are not params.
    nextflow_own = {"resume", "profile", "stub_run", "bg", "with_report",
                    "allow_missing"}
    stale = sorted(documented - declared - nextflow_own)
    check("the manual documents nothing that has been removed",
          not stale, ", ".join(stale))

    # And the defaults quoted in the manual must be the real ones.
    wrong = []
    for name in sorted(declared):
        m = re.search(rf"^\s{{4}}{name}\s*=\s*(.+?)\s*(?://.*)?$",
                      block, re.M)
        if not m:
            continue
        real = m.group(1).strip().strip("'\"")
        # Only check rows of the parameter tables: | `--x` | `default` |
        row = re.search(rf"\|\s*`--{name}`\s*\|\s*([^|]+?)\s*\|", doc)
        if not row:
            continue
        shown = row.group(1).strip().strip("`")
        if shown in ("—", "-", ""):
            # documented as "no default"; config must agree
            if real != "null":
                wrong.append(f"{name}: manual says none, config says {real}")
            continue
        if shown.replace("'", "") != real.replace("'", ""):
            wrong.append(f"{name}: manual {shown!r} vs config {real!r}")
    check("every default quoted in the manual is the real one",
          not wrong, "; ".join(wrong[:5]))


def test_stubs_present() -> None:
    print("\nevery process has a stub")
    procs = processes()
    missing = [n for n, p in procs.items() if "stub:" not in p["body"]]
    check("a stub exists for -stub-run", not missing,
          ", ".join(missing))
    no_ver = [n for n, p in procs.items() if "versions" not in p["emits"]]
    check("every process emits its versions", not no_ver, ", ".join(no_ver))


def main() -> int:
    test_includes_resolve()
    test_params_declared()
    test_publishdir_does_not_double_directories()
    test_publishdir_paths_are_lazy()
    test_tools_are_executable()
    test_parsed_stubs_have_content()
    test_placeholders_are_distinct()
    test_call_arity()
    test_fanout_is_not_capped()
    test_labels_are_sized_everywhere()
    test_usage_manual_is_current()
    test_flags_exist()
    test_emits_resolve()
    test_outputs_are_written()
    test_stubs_present()

    print(f"\n{'='*62}")
    print(f"{len(PASS)} passed, {len(FAIL)} failed")
    if FAIL:
        print("\nfailures:")
        for name, detail in FAIL:
            print(f"  {name}: {detail}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
