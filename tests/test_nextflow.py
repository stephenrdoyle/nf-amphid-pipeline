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

What this cannot tell you: whether the channel cardinalities are right,
whether a `collect()` blocks, or whether the DAG deadlocks. Those need
Nextflow. The checklist in the README says what to verify on the first
real run.

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
