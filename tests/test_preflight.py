#!/usr/bin/env python3
"""
Tests for pre-flight validation.

Every hard error must actually fire, and the happy path must stay clean.
A validator that silently passes bad input is worse than no validator, so
each case below breaks exactly one thing and asserts the specific complaint.

Run: python3 tests/test_preflight.py
"""

from __future__ import annotations

import copy
import importlib.util
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

spec = importlib.util.spec_from_file_location(
    "preflight", os.path.join(ROOT, "bin", "preflight.py")
)
pf = importlib.util.module_from_spec(spec)
# register before exec: dataclasses resolve their module from sys.modules
sys.modules["preflight"] = pf
spec.loader.exec_module(pf)

sys.path.insert(0, HERE)
import make_fixtures as mf  # noqa: E402

INDEX_MAP = os.path.join(ROOT, "assets", "index_map.tsv")

PASS, FAIL = [], []


def check(name: str, condition: bool, detail: str = "") -> None:
    if condition:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


# --------------------------------------------------------------------------
# unit tests: the pieces most likely to be subtly wrong
# --------------------------------------------------------------------------


def test_tag_derivation() -> None:
    print("\ntag derivation")
    tag = "TGGTCAATTTCAAGAGAT"

    # clean insertion, no shared bases at the boundary
    parent = "AAAACCCCGGGG" + "ACACACACAC" + "TTTTGGGGCCCC"
    control = parent[:12] + tag + parent[12:]
    r = pf.derive_tag(control, parent, tag)
    check("clean insertion recovers tag", r.tag == tag, r.tag)
    check("clean insertion recovers position", r.start == 12, str(r.start))

    # the hard case: parent ends in TG exactly where the tag begins TG, so
    # the naive common-prefix runs into the tag and returns a rotation
    parent2 = "AAAACCCCGGTG" + "ACACACACAC" + "TTTTGGGGCCCC"
    control2 = parent2[:12] + tag + parent2[12:]
    r2 = pf.derive_tag(control2, parent2, tag)
    check("ambiguous boundary still recovers the exact tag", r2.tag == tag, r2.tag)
    check(
        "ambiguous boundary tag occurs at the reported position",
        control2[r2.start : r2.start + r2.length] == r2.tag,
    )

    # left-alignment must be deterministic even with no expected tag given
    r3 = pf.derive_tag(control2, parent2, "")
    check(
        "left-aligned without an expected tag is self-consistent",
        control2[r3.start : r3.start + r3.length] == r3.tag,
    )

    # substitution-only decoy
    parent4 = "AAAACCCCGGGGTTTT"
    control4 = "AAAACCCCGGGGTTTA"
    r4 = pf.derive_tag(control4, parent4)
    check("substitution-only decoy detected", r4.kind == "substitutions", r4.kind)
    check("substitution position found", r4.substitutions == [15], str(r4.substitutions))

    # identical
    r5 = pf.derive_tag(parent4, parent4)
    check("identical control has no substitutions", r5.substitutions == [])


def test_distance_helpers() -> None:
    print("\ndistance helpers")
    check("hamming", pf.hamming("ACGT", "ACCT") == 1)
    check("levenshtein insertion", pf.levenshtein("ACGT", "ACCGT") == 1)
    check("approx substring exact", pf.approx_substring_distance("CGT", "AACGTAA") == 0)
    check("approx substring 1 edit", pf.approx_substring_distance("CGT", "AACCTAA") == 1)
    check(
        "approx substring absent",
        pf.approx_substring_distance("GGGGGGGG", "AAAAAAAAAAAA") == 8,
    )
    check("rotations count", len(pf.rotations("ACGT")) == 3)

    # the numbers quoted in the design spec
    p = pf.detection_probability(18, 0.02, 0)
    check("18bp exact match at 2% error ~69.5%", abs(p - 0.695) < 0.005, f"{p:.4f}")
    p3 = pf.detection_probability(18, 0.05, 3)
    check("18bp <=3 edits at 5% error ~98.9%", abs(p3 - 0.989) < 0.005, f"{p3:.4f}")


def test_plate_geometry() -> None:
    print("\nplate geometry")
    check(
        "odd rows and odd columns is Q1",
        pf.quadrant_of({"A", "C", "E"}, {1, 3, 5}).startswith("Q1"),
    )
    check(
        "even rows and even columns is Q4",
        pf.quadrant_of({"B", "D"}, {2, 4}).startswith("Q4"),
    )
    check(
        "mixed parity spans quadrants",
        "spans" in pf.quadrant_of({"A", "B"}, {1, 2}),
    )


def test_sanitisation() -> None:
    print("\nidentifier sanitisation")
    check("space becomes underscore", pf.sanitise_id("iCCH 001") == "iCCH_001")
    check("slash becomes underscore", pf.sanitise_id("a/b") == "a_b")
    check("hyphens and dots survive", pf.sanitise_id("iCCH-001.r1") == "iCCH-001.r1")
    check(
        "distinct ids can collide after sanitising",
        pf.sanitise_id("iCCH 001") == pf.sanitise_id("iCCH/001"),
    )


# --------------------------------------------------------------------------
# end-to-end: the happy path, then one break at a time
# --------------------------------------------------------------------------


def build_case(tmp: str, mutate_ref=None, mutate_ss=None, mutate_runs=None) -> tuple[str, str, str]:
    """Write a valid fixture set into `tmp`, optionally mutated."""
    import random

    rng = random.Random(20260918)
    ref = mf.build_reference(rng)
    i7, i5 = mf.load_index_map(INDEX_MAP)
    samples = mf.build_samplesheet(i7, i5, 88, 8, "RUN01", "illumina")
    runs_path = os.path.join(tmp, "runs.tsv")
    run_dir = os.path.join(tmp, "run01")

    if mutate_ref:
        ref = mutate_ref(copy.deepcopy(ref))
    if mutate_ss:
        samples = mutate_ss(copy.deepcopy(samples))

    ref_path = os.path.join(tmp, "reference.tsv")
    ss_path = os.path.join(tmp, "samplesheet.csv")
    mf.write_reference(ref, ref_path)
    mf.write_samplesheet(samples, ss_path)
    mf.write_runs(runs_path, "RUN01", "illumina", run_dir)
    mf.write_fake_fastqs(samples, run_dir)

    if mutate_runs:
        mutate_runs(runs_path)

    return ref_path, ss_path, runs_path


def run_preflight(tmp: str, ref: str, ss: str, runs: str) -> tuple[int, pf.Diagnostics]:
    """Run via the library rather than the CLI so the diagnostics can be
    inspected directly."""
    out = os.path.join(tmp, "out")
    code = pf.main(
        [
            "--reference", ref,
            "--samplesheet", ss,
            "--runs", runs,
            "--index-map", INDEX_MAP,
            "--outdir", out,
        ]
    )
    # re-read the structured summary for assertions
    import json

    with open(os.path.join(out, "preflight_summary.json"), encoding="utf-8") as fh:
        summary = json.load(fh)
    diag = pf.Diagnostics(errors=summary["errors"], warnings=summary["warnings"])
    return code, diag


def expect_error(name: str, substring: str, **mutators) -> None:
    tmp = tempfile.mkdtemp()
    try:
        ref, ss, runs = build_case(tmp, **mutators)
        code, diag = run_preflight(tmp, ref, ss, runs)
        hit = any(substring.lower() in e.lower() for e in diag.errors)
        check(
            name,
            code == 1 and hit,
            f"exit={code}, errors={diag.errors[:3]}" if not hit else f"exit={code}",
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_happy_path() -> None:
    print("\nhappy path")
    tmp = tempfile.mkdtemp()
    try:
        ref, ss, runs = build_case(tmp)
        code, diag = run_preflight(tmp, ref, ss, runs)
        check("valid inputs pass", code == 0, f"errors={diag.errors}")
        out = os.path.join(tmp, "out")
        for f in (
            "reference.fasta",
            "tag_report.tsv",
            "tag_report.txt",
            "file_resolution.tsv",
            "preflight_summary.json",
            "validation.log",
        ):
            check(f"writes {f}", os.path.exists(os.path.join(out, f)))

        # every sample resolved two files. Read by header, not by index:
        # a column inserted later should not silently break the assertion.
        import csv as _csv

        with open(os.path.join(out, "file_resolution.tsv"), encoding="utf-8") as fh:
            rows = list(_csv.DictReader(fh, delimiter="\t"))
        check("all samples resolved R1+R2",
              all(r["n_files"] == "2" for r in rows),
              str({r["sample_id"]: r["n_files"] for r in rows if r["n_files"] != "2"}))
        check("wells derived for all samples", all(r["well"] for r in rows))
        check("plate row and column emitted separately",
              all(r["plate_row"] and r["plate_col"] for r in rows))
        check("well label is unpadded, matching the lab plate map",
              any(r["well"].endswith(("1", "3", "5")) and len(r["well"]) == 2
                  for r in rows),
              str(sorted({r["well"] for r in rows})[:5]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_reference_errors() -> None:
    print("\nreference validation")

    def dup_id(ref):
        ref.append(dict(ref[0]))
        return ref

    expect_error("duplicate seq_id", "appears more than once", mutate_ref=dup_id)

    def identical_seq(ref):
        ref[1]["sequence"] = ref[0]["sequence"]
        return ref

    expect_error("identical sequences", "identical", mutate_ref=identical_seq)

    def bad_type(ref):
        ref[0]["type"] = "decoy"
        return ref

    expect_error("unknown type", "not in", mutate_ref=bad_type)

    def orphan_parent(ref):
        for r in ref:
            if r["type"] == "control":
                r["parent_id"] = "does_not_exist"
                break
        return ref

    expect_error("control parent unresolved", "does not match any seq_id",
                 mutate_ref=orphan_parent)

    def no_parent(ref):
        for r in ref:
            if r["type"] == "control":
                r["parent_id"] = ""
                break
        return ref

    expect_error("control without parent", "parent_id is required",
                 mutate_ref=no_parent)

    def parent_wrong_locus(ref):
        for r in ref:
            if r["type"] == "control":
                r["locus"] = "cox1" if r["locus"] == "ITS2" else "ITS2"
                break
        return ref

    expect_error("control locus differs from parent", "differs from parent",
                 mutate_ref=parent_wrong_locus)

    def quant_no_units(ref):
        for r in ref:
            if r["type"] == "quant":
                r["conc_units"] = ""
                break
        return ref

    expect_error("quant without units", "conc_units is required",
                 mutate_ref=quant_no_units)

    def quant_no_conc(ref):
        for r in ref:
            if r["type"] == "quant":
                r["default_conc"] = ""
                break
        return ref

    expect_error("quant without concentration", "default_conc is required",
                 mutate_ref=quant_no_conc)

    def bad_bases(ref):
        ref[0]["sequence"] = ref[0]["sequence"][:-3] + "XQZ"
        return ref

    expect_error("non-IUPAC sequence", "non-IUPAC", mutate_ref=bad_bases)

    def gapped(ref):
        ref[0]["sequence"] = ref[0]["sequence"][:-3] + "A-A"
        return ref

    expect_error("gapped sequence", "gaps", mutate_ref=gapped)

    def no_targets(ref):
        return [r for r in ref if r["type"] != "target"]

    expect_error("no targets", "nothing to quantify", mutate_ref=no_targets)

    def tag_in_target(ref):
        # plant the tag inside a target: read-level detection would then call
        # target reads control
        for r in ref:
            if r["type"] == "target":
                s = r["sequence"]
                r["sequence"] = s[:100] + mf.TAG + s[100 + len(mf.TAG):]
                break
        return ref

    expect_error("tag present in a target", "occurs in", mutate_ref=tag_in_target)

    def tag_in_vector(ref):
        # A chance tag-like window in the backbone is a warning, not an
        # error: the vector carries no primer sites, so it barely amplifies,
        # and the read-level search is anchored where its reads never sit.
        for r in ref:
            if r["type"] == "vector":
                s0 = r["sequence"]
                r["sequence"] = s0[:200] + mf.TAG[:-2] + "CC" + s0[200 + len(mf.TAG):]
                break
        return ref

    tmp = tempfile.mkdtemp()
    try:
        ref_p, ss_p, runs_p = build_case(tmp, mutate_ref=tag_in_vector)
        code, diag = run_preflight(tmp, ref_p, ss_p, runs_p)
        check(
            "near-tag in the vector warns rather than blocking",
            code == 0 and any("vector" in w and "edit match" in w for w in diag.warnings),
            f"exit={code}, errors={diag.errors[:2]}",
        )
        check(
            "one collision is reported once, not once per control",
            sum(1 for w in diag.warnings if "edit match inside" in w) <= 1,
            str([w[:50] for w in diag.warnings if "edit match inside" in w]),
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    def exact_tag_in_vector(ref):
        for r in ref:
            if r["type"] == "vector":
                s0 = r["sequence"]
                r["sequence"] = s0[:200] + mf.TAG + s0[200 + len(mf.TAG):]
                break
        return ref

    expect_error("exact tag in the vector is still an error", "occurs in",
                 mutate_ref=exact_tag_in_vector)

    def wrong_tag(ref):
        for r in ref:
            if r["type"] == "control":
                s = r["sequence"]
                i = s.find(mf.TAG)
                r["sequence"] = s[:i] + "ACGTACGTACGTACGTAC" + s[i + len(mf.TAG):]
                break
        return ref

    expect_error("inserted sequence is not the expected tag",
                 "does not match the expected tag", mutate_ref=wrong_tag)


def test_samplesheet_errors() -> None:
    print("\nsamplesheet validation")

    def dup_sample(ss):
        ss.append(dict(ss[5]))
        ss[-1]["index_i7"] = ss[10]["index_i7"]
        return ss

    expect_error("duplicate sample_id in a run", "appears 2 times",
                 mutate_ss=dup_sample)

    def mixed_platform(ss):
        ss[4]["platform"] = "nanopore"
        return ss

    expect_error("mixed platforms in one run", "mixed platforms",
                 mutate_ss=mixed_platform)

    def bad_sample_type(ss):
        ss[6]["sample_type"] = "negative"
        return ss

    expect_error("unknown sample_type", "sample_type", mutate_ss=bad_sample_type)

    def unknown_index(ss):
        ss[7]["index_i5"] = "ACGTACGT"
        return ss

    expect_error("index not in the map", "not in the index map",
                 mutate_ss=unknown_index)

    def collide_well(ss):
        ss[8]["index_i5"] = ss[9]["index_i5"]
        ss[8]["index_i7"] = ss[9]["index_i7"]
        return ss

    expect_error("two samples in one well", "resolves to well",
                 mutate_ss=collide_well)

    def sanitise_collision(ss):
        ss[10]["sample_id"] = "iCCH 007"
        ss[11]["sample_id"] = "iCCH/007"
        return ss

    expect_error("sample_ids collide after sanitising", "sanitise to",
                 mutate_ss=sanitise_collision)

    # The fixture writes fastqs from whatever sheet it is given, so a
    # resolution failure is simulated by deleting a file afterwards rather
    # than by mutating the sheet.
    tmp = tempfile.mkdtemp()
    try:
        ref, ss_path, runs = build_case(tmp)
        # remove one sample's R2 to simulate a resolution failure
        run_dir = os.path.join(tmp, "run01")
        victim = sorted(f for f in os.listdir(run_dir) if "_R2_" in f)[0]
        os.remove(os.path.join(run_dir, victim))
        code, diag = run_preflight(tmp, ref, ss_path, runs)
        hit = any("no R2 file matched" in e for e in diag.errors)
        check("missing fastq hard-fails", code == 1 and hit,
              f"exit={code}, errors={diag.errors[:2]}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_warnings_not_errors() -> None:
    print("\nwarnings must not block")

    def drop_blank(ss):
        return [s for s in ss if s["sample_type"] != "neg_extraction"]

    tmp = tempfile.mkdtemp()
    try:
        ref, ss_path, runs = build_case(tmp, mutate_ss=drop_blank)
        code, diag = run_preflight(tmp, ref, ss_path, runs)
        warned = any("no neg_extraction" in w for w in diag.warnings)
        check("missing blank warns but passes", code == 0 and warned,
              f"exit={code}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    def drop_vector(ref):
        return [r for r in ref if r["type"] != "vector"]

    tmp = tempfile.mkdtemp()
    try:
        ref_p, ss_p, runs_p = build_case(tmp, mutate_ref=drop_vector)
        code, diag = run_preflight(tmp, ref_p, ss_p, runs_p)
        warned = any("vector" in w for w in diag.warnings)
        check("missing vector warns but passes", code == 0 and warned, f"exit={code}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    def odd_corr(ref):
        ref[0]["corr_factor"] = "42"
        return ref

    tmp = tempfile.mkdtemp()
    try:
        ref_p, ss_p, runs_p = build_case(tmp, mutate_ref=odd_corr)
        code, diag = run_preflight(tmp, ref_p, ss_p, runs_p)
        warned = any("corr_factor" in w for w in diag.warnings)
        check("implausible corr_factor warns but passes", code == 0 and warned,
              f"exit={code}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_sentinel_detection() -> None:
    print("\nsentinel and observability")
    import json
    import random

    for n, unused, expect_obs in ((88, 8, 8), (96, 0, 0)):
        tmp = tempfile.mkdtemp()
        try:
            rng = random.Random(20260918)
            ref = mf.build_reference(rng)
            i7, i5 = mf.load_index_map(INDEX_MAP)
            samples = mf.build_samplesheet(i7, i5, n, unused, "RUN01", "illumina")
            ref_p = os.path.join(tmp, "reference.tsv")
            ss_p = os.path.join(tmp, "samplesheet.csv")
            runs_p = os.path.join(tmp, "runs.tsv")
            run_dir = os.path.join(tmp, "run01")
            mf.write_reference(ref, ref_p)
            mf.write_samplesheet(samples, ss_p)
            mf.write_runs(runs_p, "RUN01", "illumina", run_dir)
            mf.write_fake_fastqs(samples, run_dir)
            code, diag = run_preflight(tmp, ref_p, ss_p, runs_p)
            with open(os.path.join(tmp, "out", "preflight_summary.json"), encoding="utf-8") as fh:
                plates = json.load(fh)["plates"]["RUN01"]
            check(
                f"{n} samples + {unused} unused gives {expect_obs} sentinels",
                plates["sentinels"] == expect_obs,
                str(plates),
            )
            if expect_obs == 0:
                check(
                    "a filled plate warns that misassignment is unmeasurable",
                    any("UNMEASURABLE" in w for w in diag.warnings),
                )
            check(f"{n}-sample quadrant identified", plates["quadrant"].startswith("Q1"),
                  plates["quadrant"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


def test_reference_only_mode() -> None:
    """The panel must be checkable before a samplesheet exists - and that
    mode must not quietly skip the panel checks themselves."""
    print("\nreference-only mode")
    import json
    import random

    def run_ref_only(ref_path: str, outdir: str):
        code = pf.main(["--reference", ref_path, "--outdir", outdir])
        with open(os.path.join(outdir, "preflight_summary.json"), encoding="utf-8") as fh:
            return code, json.load(fh)

    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(20260918)
        ref = mf.build_reference(rng)
        good = os.path.join(tmp, "reference.tsv")
        mf.write_reference(ref, good)

        code, summary = run_ref_only(good, os.path.join(tmp, "a"))
        check("valid panel passes without a samplesheet", code == 0,
              str(summary["errors"]))
        check("mode recorded", summary["mode"] == "reference_only")
        check("tags still derived", len(summary["tags"]) == 5, str(len(summary["tags"])))
        check("fasta still written",
              os.path.exists(os.path.join(tmp, "a", "reference.fasta")))
        check("no samplesheet-dependent output",
              not os.path.exists(os.path.join(tmp, "a", "file_resolution.tsv")))
        check("no plate geometry claimed", summary["plates"] == {})

        # the panel checks must still bite
        import copy as _copy

        bad = _copy.deepcopy(ref)
        bad[1]["sequence"] = bad[0]["sequence"]
        dup = os.path.join(tmp, "dup.tsv")
        mf.write_reference(bad, dup)
        code, summary = run_ref_only(dup, os.path.join(tmp, "b"))
        check("identical sequences still rejected",
              code == 1 and any("identical" in e for e in summary["errors"]))

        bad2 = _copy.deepcopy(ref)
        for r in bad2:
            if r["type"] == "quant":
                r["conc_units"] = ""
                break
        noq = os.path.join(tmp, "noq.tsv")
        mf.write_reference(bad2, noq)
        code, summary = run_ref_only(noq, os.path.join(tmp, "c"))
        check("quant units still required",
              code == 1 and any("conc_units" in e for e in summary["errors"]))

        bad3 = _copy.deepcopy(ref)
        for r in bad3:
            if r["type"] == "target":
                s = r["sequence"]
                r["sequence"] = s[:80] + mf.TAG + s[80 + len(mf.TAG):]
                break
        tagged = os.path.join(tmp, "tagged.tsv")
        mf.write_reference(bad3, tagged)
        code, summary = run_ref_only(tagged, os.path.join(tmp, "d"))
        check("tag inside a target still rejected",
              code == 1 and any("occurs in" in e for e in summary["errors"]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # the shipped template must itself be well formed apart from placeholders
    tmpl = os.path.join(ROOT, "assets", "reference_template.tsv")
    check("template exists", os.path.exists(tmpl))
    with open(tmpl, encoding="utf-8") as fh:
        rows = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    header = rows[0].rstrip("\n").split("\t")
    expected = ["seq_id", "locus", "type", "purpose", "taxon", "parent_id",
                "default_conc", "conc_units", "corr_factor", "sequence"]
    check("template header matches the schema", header == expected, str(header))
    check(f"every template row has {len(expected)} fields",
          all(len(r.rstrip("\n").split("\t")) == len(expected) for r in rows),
          str([len(r.rstrip('\n').split('\t')) for r in rows]))
    # the template must itself validate, or it teaches the wrong shape
    tmp2 = tempfile.mkdtemp()
    try:
        import json
        real = os.path.join(tmp2, "t.tsv")
        with open(tmpl, encoding="utf-8") as fh, open(real, "w", encoding="utf-8") as out:
            import random as _r
            rng = _r.Random(1)
            for ln in fh:
                if ln.startswith("#") or not ln.strip():
                    continue
                if "REPLACE_WITH_SEQUENCE" in ln:
                    ln = ln.replace("REPLACE_WITH_SEQUENCE",
                                    "".join(rng.choice("ACGT") for _ in range(300)))
                out.write(ln)
        code = pf.main(["--reference", real, "--outdir", os.path.join(tmp2, "o")])
        with open(os.path.join(tmp2, "o", "preflight_summary.json"), encoding="utf-8") as fh:
            summ = json.load(fh)
        check("template parses with the real validator",
              not any("column" in e for e in summ["errors"]),
              str([e for e in summ["errors"] if "column" in e]))
    finally:
        shutil.rmtree(tmp2, ignore_errors=True)


def test_purpose() -> None:
    """`purpose` exists to keep a second locus out of the composition
    denominator. Without it, adding beta-tubulin counts every species twice
    with a denominator mixing two loci at different copy numbers - numbers
    that stay plausible and are quietly wrong."""
    print("\nreference purpose")
    import json
    import random

    REF_COLS_P = ["seq_id", "locus", "type", "purpose", "taxon", "parent_id",
                  "default_conc", "conc_units", "corr_factor", "sequence"]

    def write(rows, path):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\t".join(REF_COLS_P) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in REF_COLS_P) + "\n")

    def run(rows):
        tmp = tempfile.mkdtemp()
        try:
            rp = os.path.join(tmp, "reference.tsv")
            write(rows, rp)
            out = os.path.join(tmp, "out")
            code = pf.main(["--reference", rp, "--outdir", out])
            with open(os.path.join(out, "preflight_summary.json"), encoding="utf-8") as fh:
                summary = json.load(fh)
            return code, summary, pf.Diagnostics(errors=summary["errors"],
                                                 warnings=summary["warnings"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    base = mf.build_reference(random.Random(20260918))

    # backward compatible: no purpose set anywhere behaves as before
    code, summary, diag = run([dict(r, purpose="") for r in base])
    n_targets = sum(1 for r in base if r["type"] == "target")
    check("blank purpose passes", code == 0, str(diag.errors))
    check("blank purpose means speciation",
          len(summary["reference"]["composition_denominator"]) == n_targets,
          str(len(summary["reference"]["composition_denominator"])))
    check("nothing marked for variants by default",
          summary["reference"]["variant_references"] == [])

    # a second locus marked variants must leave the denominator alone
    rows = []
    for r in base:
        p = "variants" if r["locus"] == "cox1" and r["type"] == "target" else ""
        rows.append(dict(r, purpose=p))
    code, summary, diag = run(rows)
    n_cox1 = sum(1 for r in base if r["locus"] == "cox1" and r["type"] == "target")
    check("variants excluded from the composition denominator",
          len(summary["reference"]["composition_denominator"]) == n_targets - n_cox1,
          str(len(summary["reference"]["composition_denominator"])))
    check("variant references listed",
          any("cox1" in x for x in summary["reference"]["variant_references"]))
    check("controls inherit purpose from their parent",
          any(x.endswith("_ctrl") for x in summary["reference"]["variant_references"]),
          str(summary["reference"]["variant_references"]))

    # both purposes at once is legitimate
    rows = [dict(r, purpose="speciation,variants"
                 if r["type"] == "target" else "") for r in base]
    code, summary, diag = run(rows)
    check("speciation,variants counts for both",
          code == 0
          and len(summary["reference"]["composition_denominator"]) == n_targets
          and len(summary["reference"]["variant_references"]) >= n_targets,
          str(diag.errors))

    # a bad value is an error, not a silent ignore
    rows = [dict(r, purpose="genotyping" if r["seq_id"] == "Nam_ITS2_h1" else "")
            for r in base]
    code, summary, diag = run(rows)
    check("unknown purpose value is rejected",
          code == 1 and any("purpose" in e and "not in" in e for e in diag.errors),
          str(diag.errors[:2]))

    # a control disagreeing with its parent is a mistake worth flagging
    rows = [dict(r, purpose="variants" if r["seq_id"] == "Nam_ITS2_h1_ctrl" else "")
            for r in base]
    code, summary, diag = run(rows)
    check("control disagreeing with its parent warns",
          code == 0 and any("differs from its parent" in w for w in diag.warnings),
          str([w[:60] for w in diag.warnings]))

    # one entry flagged differently from its locus-mates
    rows = [dict(r, purpose="variants" if r["seq_id"] == "Nam_ITS2_h2" else "")
            for r in base]
    code, summary, diag = run(rows)
    check("mixed purpose within a locus warns",
          code == 0 and any("mixes purposes" in w for w in diag.warnings),
          str([w[:60] for w in diag.warnings]))

    # purpose on a quant or vector is meaningless
    rows = [dict(r, purpose="speciation" if r["type"] == "quant" else "")
            for r in base]
    code, summary, diag = run(rows)
    check("purpose on a quant is flagged as ignored",
          any("is ignored" in w for w in diag.warnings),
          str([w[:60] for w in diag.warnings]))


def test_primer_extent() -> None:
    """Extent is the one thing primers are needed for. A panel where some
    entries carry flanking sequence outside the primers and others do not
    biases assignment toward the longer, and hides inside genuine length
    variation where no proxy finds it."""
    print("\nprimer extent check")
    import random

    def write_primers(tmp, fwd, rev):
        p = os.path.join(tmp, "primers.tsv")
        with open(p, "w", encoding="utf-8") as fh:
            fh.write("name\tdirection\tsequence\n")
            for i, s0 in enumerate(fwd, 1):
                fh.write(f"F{i}\tfwd\t{s0}\n")
            for i, s0 in enumerate(rev, 1):
                fh.write(f"R{i}\trev\t{s0}\n")
        return p

    def run(tmp, ref_path, primers_path):
        import json

        out = os.path.join(tmp, "out")
        code = pf.main(["--reference", ref_path, "--primers", primers_path,
                        "--outdir", out])
        with open(os.path.join(out, "preflight_summary.json"), encoding="utf-8") as fh:
            summary = json.load(fh)
        return code, pf.Diagnostics(errors=summary["errors"],
                                    warnings=summary["warnings"]), summary

    # the fixture panel is uniform: every entry runs primer to primer
    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(20260918)
        ref = mf.build_reference(rng)
        ref_p = os.path.join(tmp, "reference.tsv")
        mf.write_reference(ref, ref_p)

        # recover the primers the fixture used, from a clean target pair
        its2 = [r for r in ref if r["type"] == "target" and r["locus"] == "ITS2"]
        cox1 = [r for r in ref if r["type"] == "target" and r["locus"] == "cox1"]
        fwd = [its2[0]["sequence"][:20], cox1[0]["sequence"][:20]]
        rev = [pf.revcomp(its2[0]["sequence"][-20:]),
               pf.revcomp(cox1[0]["sequence"][-20:])]
        pp = write_primers(tmp, fwd, rev)

        code, diag, summary = run(tmp, ref_p, pp)
        check("uniform panel passes", code == 0, str(diag.errors))
        check("no extent warning on a uniform panel",
              not any("extent is not uniform" in w for w in diag.warnings),
              str([w[:60] for w in diag.warnings]))
        ext = summary["reference"]["extent"]
        # 5 targets + 5 controls + 2 quants; the vector does not amplify
        check("extent recorded for every amplifiable entry", len(ext) == 12, str(len(ext)))
        check("all entries measured as primer-to-primer",
              all(e["extra_5p"] == 0 and e["extra_3p"] == 0 for e in ext))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # now give two entries flanking sequence, as the real panel had
    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(20260918)
        ref = mf.build_reference(rng)
        its2 = [r for r in ref if r["type"] == "target" and r["locus"] == "ITS2"]
        cox1 = [r for r in ref if r["type"] == "target" and r["locus"] == "cox1"]
        fwd = [its2[0]["sequence"][:20], cox1[0]["sequence"][:20]]
        rev = [pf.revcomp(its2[0]["sequence"][-20:]),
               pf.revcomp(cox1[0]["sequence"][-20:])]
        for r in ref:
            if r["seq_id"] in ("Tri_ITS2_h1", "Tri_ITS2_h1_ctrl"):
                r["sequence"] = "ATGTCGACGCT" + r["sequence"] + "CCGCTGCTACCACTT"
        ref_p = os.path.join(tmp, "reference.tsv")
        mf.write_reference(ref, ref_p)
        pp = write_primers(tmp, fwd, rev)

        code, diag, summary = run(tmp, ref_p, pp)
        check("mixed extent warns but does not block", code == 0, str(diag.errors))
        check("mixed extent is reported",
              any("extent is not uniform" in w for w in diag.warnings))
        check("the offending entries are named",
              any("Tri_ITS2_h1" in w for w in diag.warnings
                  if "extent is not uniform" in w))
        ext = {e["seq_id"]: e for e in summary["reference"]["extent"]}
        check("5' excess measured exactly", ext["Tri_ITS2_h1"]["extra_5p"] == 11,
              str(ext["Tri_ITS2_h1"]["extra_5p"]))
        check("3' excess measured exactly", ext["Tri_ITS2_h1"]["extra_3p"] == 15,
              str(ext["Tri_ITS2_h1"]["extra_3p"]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # a trimmed-to-insert entry has no primer at all
    tmp = tempfile.mkdtemp()
    try:
        rng = random.Random(20260918)
        ref = mf.build_reference(rng)
        its2 = [r for r in ref if r["type"] == "target" and r["locus"] == "ITS2"]
        cox1 = [r for r in ref if r["type"] == "target" and r["locus"] == "cox1"]
        fwd = [its2[0]["sequence"][:20], cox1[0]["sequence"][:20]]
        rev = [pf.revcomp(its2[0]["sequence"][-20:]),
               pf.revcomp(cox1[0]["sequence"][-20:])]
        for r in ref:
            if r["seq_id"] == "Asc_ITS2_h1":
                r["sequence"] = r["sequence"][25:-25]
        ref_p = os.path.join(tmp, "reference.tsv")
        mf.write_reference(ref, ref_p)
        pp = write_primers(tmp, fwd, rev)
        code, diag, _ = run(tmp, ref_p, pp)
        check("insert-only entry is flagged",
              any("no forward primer found" in w or "no reverse primer found" in w
                  for w in diag.warnings),
              str([w[:70] for w in diag.warnings]))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_assay_config() -> None:
    """The assay definition must be swappable, override-able, forward
    compatible, and recorded in the summary - otherwise a result cannot be
    traced back to the assay that produced it."""
    print("\nassay configuration")
    import json

    real = os.path.join(ROOT, "assets", "assay", "ampseq_helminth.json")

    diag = pf.Diagnostics()
    assay = pf.load_assay(real, diag)
    check("shipped assay config loads", diag.ok, str(diag.errors))
    check("name read", assay.name == "ampseq_helminth", assay.name)
    check("version read", assay.version == "1.0.0", assay.version)
    check(
        "tag read from config",
        assay.section("control_tag")["sequence"] == "TGGTCAATTTCAAGAGAT",
    )
    check(
        "error profiles sorted by error rate",
        [n for n, _ in assay.error_profiles][0] == "illumina",
        str(assay.error_profiles),
    )
    check(
        "comment keys are not treated as settings",
        not any(k.startswith("_") for k in assay.section("control_tag")),
    )

    # forward compatibility: a config missing a key added later still works
    tmp = tempfile.mkdtemp()
    try:
        partial = os.path.join(tmp, "partial.json")
        with open(partial, "w", encoding="utf-8") as fh:
            json.dump({"name": "minimal", "version": "0.1"}, fh)
        d2 = pf.Diagnostics()
        a2 = pf.load_assay(partial, d2)
        check("sparse config loads", d2.ok, str(d2.errors))
        check(
            "missing sections fall back to defaults",
            a2.section("tolerances")["corr_factor_max"] == 10.0,
        )
        check("its own values still win", a2.name == "minimal")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # bad input is reported, not raised
    tmp = tempfile.mkdtemp()
    try:
        broken = os.path.join(tmp, "broken.json")
        open(broken, "w").write("{not json")
        d3 = pf.Diagnostics()
        pf.load_assay(broken, d3)
        check("malformed JSON is an error", any("not valid JSON" in e for e in d3.errors))

        d4 = pf.Diagnostics()
        pf.load_assay(os.path.join(tmp, "absent.json"), d4)
        check("missing file is an error", any("not found" in e for e in d4.errors))

        anon = os.path.join(tmp, "anon.json")
        with open(anon, "w", encoding="utf-8") as fh:
            json.dump({"control_tag": {"sequence": "ACGT"}}, fh)
        d5 = pf.Diagnostics()
        pf.load_assay(anon, d5)
        check(
            "unnamed assay warns about traceability",
            any("traced back" in w for w in d5.warnings),
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # a genuinely different assay: different tag, tighter tolerances
    tmp = tempfile.mkdtemp()
    try:
        other_tag = "CCAGTTGACCTAAGTCGA"
        alt = os.path.join(tmp, "other.json")
        with open(alt, "w", encoding="utf-8") as fh:
            json.dump(
                {
                    "name": "other_assay",
                    "version": "2.0",
                    "control_tag": {"sequence": other_tag, "max_edits": 1},
                    "tolerances": {"corr_factor_min": 0.9, "corr_factor_max": 1.1},
                },
                fh,
            )

        # panel built with the DEFAULT tag must now fail against this assay
        ref, ss, runs = build_case(tmp)
        out = os.path.join(tmp, "out")
        code = pf.main(
            ["--reference", ref, "--samplesheet", ss, "--runs", runs,
             "--index-map", INDEX_MAP, "--assay", alt, "--outdir", out]
        )
        with open(os.path.join(out, "preflight_summary.json"), encoding="utf-8") as fh:
            summary = json.load(fh)
        check(
            "wrong tag for the declared assay is caught",
            code == 1
            and any("does not match the expected tag" in e for e in summary["errors"]),
            str(summary["errors"][:1]),
        )
        check("assay name recorded in the summary",
              summary["assay"]["name"] == "other_assay")
        check("assay source recorded", summary["assay"]["source"] == alt)
        check("tolerance came from the config",
              summary["assay"]["tag_max_edits"] == 1)

        # tightened corr_factor range must now flag the panel's 1.42
        check(
            "tightened tolerance changes what warns",
            any("corr_factor" in w and "outside" in w for w in summary["warnings"]),
        )

        # an explicit flag overrides the config
        out2 = os.path.join(tmp, "out2")
        pf.main(
            ["--reference", ref, "--samplesheet", ss, "--runs", runs,
             "--index-map", INDEX_MAP, "--assay", alt,
             "--expected-tag", "TGGTCAATTTCAAGAGAT", "--outdir", out2]
        )
        with open(os.path.join(out2, "preflight_summary.json"), encoding="utf-8") as fh:
            s2 = json.load(fh)
        check(
            "explicit --expected-tag overrides the config",
            s2["assay"]["control_tag"] == "TGGTCAATTTCAAGAGAT",
            s2["assay"]["control_tag"],
        )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_barcode_repeats() -> None:
    """Dual index with and without the tandem duplication. Both are valid;
    they differ in whether each index carries its own internal replicate,
    and therefore in whether a frame slip can be detected."""
    print("\nbarcode repeats")
    import json
    import random

    def run_with(reps: int, platform: str):
        tmp = tempfile.mkdtemp()
        try:
            rng = random.Random(20260918)
            ref = mf.build_reference(rng)
            i7, i5 = mf.load_index_map(INDEX_MAP)
            samples = mf.build_samplesheet(i7, i5, 24, 4, "RUN01", platform)
            ref_p = os.path.join(tmp, "reference.tsv")
            ss_p = os.path.join(tmp, "samplesheet.csv")
            runs_p = os.path.join(tmp, "runs.tsv")
            run_dir = os.path.join(tmp, "run01")
            mf.write_reference(ref, ref_p)
            mf.write_samplesheet(samples, ss_p)
            mf.write_runs(runs_p, "RUN01", platform, run_dir)
            if platform == "illumina":
                mf.write_fake_fastqs(samples, run_dir)
            else:
                os.makedirs(run_dir, exist_ok=True)

            lines = open(runs_p).read().splitlines()
            head = lines[0].split("\t")
            row = lines[1].split("\t")
            row[head.index("barcode_repeats")] = str(reps)
            open(runs_p, "w").write("\t".join(head) + "\n" + "\t".join(row) + "\n")

            code, diag = run_preflight(tmp, ref_p, ss_p, runs_p)
            with open(os.path.join(tmp, "out", "preflight_summary.json"), encoding="utf-8") as fh:
                strat = json.load(fh)["barcode_strategy"]["RUN01"]
            return code, diag, strat
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    code, diag, strat = run_with(2, "nanopore")
    check("repeats=2 passes", code == 0, str(diag.errors))
    check("repeats=2 requires copy agreement", strat["require_copy_agreement"] is True)
    check(
        "repeats=2 protects against frame slip",
        strat["frame_slip_protection"] == "copy agreement",
    )
    check(
        "repeats=2 raises no single-copy warning",
        not any("barcode_repeats=1" in w for w in diag.warnings),
    )

    code, diag, strat = run_with(1, "nanopore")
    check("repeats=1 passes", code == 0, str(diag.errors))
    check("repeats=1 drops copy agreement", strat["require_copy_agreement"] is False)
    check("repeats=1 has no frame-slip protection",
          strat["frame_slip_protection"] == "none")
    check(
        "repeats=1 warns about the lost cross-check",
        any("no second copy to cross-check" in w for w in diag.warnings),
    )
    check(
        "repeats=1 escalates the rotational collisions",
        any("live misassignment routes" in w for w in diag.warnings),
    )
    check(
        "indels stay disabled in both modes",
        strat["allow_indels"] is False,
    )

    # None of this applies to Illumina, which the sequencer demultiplexes.
    code, diag, strat = run_with(1, "illumina")
    check("repeats=1 on Illumina is not flagged", code == 0
          and not any("no second copy" in w for w in diag.warnings))
    check("Illumina marked as not applicable", strat["applies"] is False)

    # Default must be 2 when the manifest omits the column entirely.
    tmp = tempfile.mkdtemp()
    try:
        ref, ss, runs = build_case(tmp)
        lines = open(runs).read().splitlines()
        head = lines[0].split("\t")
        row = lines[1].split("\t")
        i = head.index("barcode_repeats")
        del head[i]
        del row[i]
        open(runs, "w").write("\t".join(head) + "\n" + "\t".join(row) + "\n")
        code, diag = run_preflight(tmp, ref, ss, runs)
        with open(os.path.join(tmp, "out", "preflight_summary.json"), encoding="utf-8") as fh:
            reps = json.load(fh)["runs"]["RUN01"]["barcode_repeats"]
        check("omitted column defaults to 2", reps == 2, str(reps))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_index_set_numbers() -> None:
    print("\nindex set characterisation")
    diag = pf.Diagnostics()
    im = pf.parse_index_map(INDEX_MAP, diag)
    check("index map parses cleanly", diag.ok, str(diag.errors))
    check("24 i7 indices", len(im.i7) == 24, str(len(im.i7)))
    check("16 i5 indices", len(im.i5) == 16, str(len(im.i5)))
    # Unpadded, matching the lab's own 384 plate map (A1 not A01), so well
    # labels join against hand-written plate records.
    check("well A1 derived", im.well("TAAGGCGA", "CTCTCTAT") == "A1",
          str(im.well("TAAGGCGA", "CTCTCTAT")))
    check("well P24 derived", im.well("TCGACGTC", "TTATGCGA") == "P24")
    check("well C10 derived", im.well("GTAGAGGA", "GTAAGGAG") == "C10")
    check("row and column available separately",
          im.row_col("GTAGAGGA", "GTAAGGAG") == ("C", 10))

    rep = pf.index_set_report(im, pf.Diagnostics())
    check("i7 min Hamming is 4", rep["i7"]["min_hamming"] == 4, str(rep["i7"]))
    check("i7 min edit distance is 2", rep["i7"]["min_levenshtein"] == 2)
    check("i5 min Hamming is 4", rep["i5"]["min_hamming"] == 4)
    check("cross-set distance is 1", rep["cross_min_distance"] == 1)
    check(
        "i7 tolerates 1 substitution but 0 indels",
        rep["i7"]["safe_substitutions"] == 1 and rep["i7"]["safe_edits"] == 0,
    )
    check("rotational collisions found", len(rep["rotational_collisions"]) == 4,
          str(len(rep["rotational_collisions"])))


def main() -> int:
    test_tag_derivation()
    test_distance_helpers()
    test_plate_geometry()
    test_sanitisation()
    test_index_set_numbers()
    test_reference_only_mode()
    test_purpose()
    test_primer_extent()
    test_assay_config()
    test_barcode_repeats()
    test_happy_path()
    test_reference_errors()
    test_samplesheet_errors()
    test_warnings_not_errors()
    test_sentinel_detection()

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
