#!/usr/bin/env python3
"""
Tests for the index misassignment model (ampseqlib.fit_misassignment)
and the plate-position tests shared by both syntheses.

Simulated plates with a known per-pair rate along rows and columns
(simulate.simulate_plate_counts), at 42, 96 and 384 wells, with and
without unused index combinations: the fitted rates must recover the
truth, and the unused combinations, where there are any, must agree.
Then the cases where a rate must NOT be fitted - no misassignment at
all, too small a plate, samples laid out by group down the columns,
contamination that does not follow the indexes - each of which says so
rather than fitting noise.

The read-level simulator's per-end misreads are tested in
test_simulate.py; the calls built on the rates in test_synthesis_db.py.

Run: python3 tests/test_misassignment.py        (prints a results table)
"""

from __future__ import annotations

import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "bin"))

import ampseqlib as L  # noqa: E402
import simulate as sim  # noqa: E402

PASS, FAIL = [], []
TABLE: list[str] = []


def check(name: str, cond: bool, detail: str = "") -> None:
    if cond:
        PASS.append(name)
        print(f"  ok   {name}")
    else:
        FAIL.append((name, detail))
        print(f"  FAIL {name}" + (f" - {detail}" if detail else ""))


def plate(n_rows, n_cols, rate_row, rate_col, seed=1, **kw):
    p = sim.simulate_plate_counts(random.Random(seed), n_rows, n_cols,
                                  rate_row=rate_row, rate_col=rate_col, **kw)
    fit = L.fit_misassignment(p["geom"], p["counts"], random.Random(seed + 7))
    unused = L.sentinel_rates(p["matrix"], fit.get("for_calls"))
    return p, fit, unused


def rate_cell(fit, unused, ax, truth) -> str:
    a = fit[ax]
    if a["basis"] != "fitted":
        return f"{a['basis']} | |"
    lo, hi = a["ci95"]
    u = (unused or {}).get(ax)
    return (f"{a['rate']:.2e} ({100 * (a['rate'] / truth - 1):+.0f}%) | "
            f"{lo:.1e}-{hi:.1e} | "
            + (f"{u['rate']:.2e}" if u else "-"))


def test_recovery() -> None:
    """The fitted rates recover the truth within 20%, on combinatorial
    plates of 42 (8 x 6), 96 and 384 wells, full and with unused
    combinations; and the unused combinations, where they exist and
    measure an axis, agree within 25%."""
    print("\nrecovery of known rates")
    TABLE.append("| wells | unused | true row / col | row fit (error) | "
                 "row 95% CI | row unused | col fit (error) | col 95% CI "
                 "| col unused |")
    TABLE.append("|---|---|---|---|---|---|---|---|---|")
    for n_rows, n_cols, unused in ((8, 6, 6), (8, 6, 0), (8, 12, 0),
                                   (8, 12, 8), (16, 24, 0), (16, 24, 16)):
        for rr, rc in ((4e-4, 7e-4), (1e-4, 2e-3)):
            n = n_rows * n_cols - unused
            p, fit, su = plate(n_rows, n_cols, rr, rc, unused=unused)
            for ax, truth in (("row", rr), ("col", rc)):
                a = fit[ax]
                ok = a["basis"] == "fitted" and \
                    abs(a["rate"] / truth - 1) <= 0.2
                # A plate of 48 at 1e-4 predicts 10-20 leaked reads a
                # pair: falling back, and saying so, is right there too.
                thin = n < 60 and truth <= 1e-4
                check(f"{n} wells, {unused} unused: {ax} {truth:.0e} "
                      "recovered" + (" (or an explicit fallback)"
                                     if thin else ""),
                      ok or (thin and a["basis"] == "default"
                             and bool(a["reason"])),
                      f"{a['basis']} {a['rate']} {a['reason']}")
                if su and su.get(ax):
                    check(f"{n} wells: unused combinations agree on {ax}",
                          abs(su[ax]["rate"] / truth - 1) <= 0.25,
                          f"{su[ax]['rate']:.2e} vs {truth:.0e}")
            TABLE.append(
                f"| {n} | {unused} | {rr:.0e} / {rc:.0e} | "
                + rate_cell(fit, su, "row", rr) + " | "
                + rate_cell(fit, su, "col", rc) + " |")


def test_no_misassignment() -> None:
    """No misassignment at all: the wells next to heavy sources show
    nothing, and the rate is a bound - with nothing checked against it
    and nothing taken off the blanks."""
    print("\nno misassignment")
    for n_rows, n_cols in ((8, 12), (16, 24)):
        p, fit, _su = plate(n_rows, n_cols, 0.0, 0.0)
        for ax in ("row", "col"):
            check(f"{n_rows * n_cols} wells, rate 0: {ax} is a bound, not a "
                  "fit", fit[ax]["basis"] == "upper_bound",
                  f"{fit[ax]['basis']} {fit[ax]['reason']}")
        check("nothing is checked or taken off at a bound",
              fit["for_calls"] == {"row": 0.0, "col": 0.0}
              and fit["for_blanks"] == {"row": 0.0, "col": 0.0},
              str(fit["for_calls"]))
    TABLE.append("")
    TABLE.append("| case | row | col |")
    TABLE.append("|---|---|---|")
    TABLE.append(f"| 96 wells, no misassignment | {fit['row']['basis']} | "
                 f"{fit['col']['basis']} |")


def test_fallback() -> None:
    """Where the run cannot measure a rate it says why, checks calls at
    the default, and takes nothing off the blanks: a 3 x 4 plate with
    no heavy source, and a plate too shallow for any leak to show."""
    print("\nthe fallback engages when the contrast is not there")
    cases = (("12 wells, no heavy source", (3, 4, 4e-4, 7e-4),
              {"heavy_column": False, "seed": 2}),
             ("96 wells at 5,000 reads", (8, 12, 4e-4, 7e-4),
              {"depth": 5000}))
    for label, args, kw in cases:
        p, fit, _su = plate(*args, **kw)
        bases = {ax: fit[ax]["basis"] for ax in ("row", "col")}
        fell = [ax for ax in ("row", "col") if bases[ax] == "default"]
        check(f"{label}: at least one axis falls back", bool(fell),
              str(bases))
        for ax in fell:
            check(f"{label}: {ax} says why", bool(fit[ax]["reason"]))
            check(f"{label}: {ax} calls checked at the default",
                  fit["for_calls"][ax] == L.MISASSIGN_DEFAULT_RATE)
            check(f"{label}: {ax} takes nothing off the blanks",
                  fit["for_blanks"][ax] == 0.0)
        TABLE.append(f"| {label} | {bases['row']} | {bases['col']} |")


def test_biology_down_the_columns() -> None:
    """Samples laid out by group - each column its own host - share
    taxa down the columns for real. That must not be read as
    misassignment; the row rate is still measured."""
    print("\nsamples grouped by column are not misassignment")
    p, fit, _su = plate(8, 12, 4e-4, 0.0, grouped_by_column=True,
                        heavy_column=False)
    check("column: not fitted from shared biology",
          fit["col"]["basis"] != "fitted"
          or fit["col"]["rate"] < 1e-4,
          f"{fit['col']['basis']} {fit['col']['rate']} "
          f"{fit['col']['reason']}")
    check("row: still measured",
          fit["row"]["basis"] == "fitted"
          and abs(fit["row"]["rate"] / 4e-4 - 1) <= 0.2,
          f"{fit['row']['basis']} {fit['row']['rate']}")
    TABLE.append(f"| 96 wells, samples grouped by column, no column "
                 f"misassignment | {fit['row']['basis']} "
                 f"{fit['row']['rate'] or 0:.1e} | {fit['col']['basis']} |")


def test_not_following_the_indexes() -> None:
    """Contamination that does not follow the indexes - every well
    receiving the same small share of the whole plate's reads, as from a
    shared tube - is not misassignment, and the placebo (wells shuffled
    over the plate) says so. A test at p 0.02 is wrong now and then by
    design, so over eight plates, with and without a heavy column: at
    most one axis in eight may be fitted."""
    print("\ncontamination that ignores the indexes")
    fitted = tests = 0
    for heavy in (True, False):
        for seed in range(3, 7):
            p = sim.simulate_plate_counts(random.Random(seed), 8, 12,
                                          rate_row=0.0, rate_col=0.0,
                                          heavy_column=heavy)
            tot: dict = {}
            for cs in p["true"].values():
                for k, v in cs.items():
                    tot[k] = tot.get(k, 0) + v
            rng = random.Random(seed + 100)
            counts = {}
            for sid, cs in p["true"].items():
                got = {}
                for k in tot:
                    n = cs.get(k, 0) + sim.poisson(
                        rng, 4e-4 * (tot[k] - cs.get(k, 0)))
                    if n >= 5:
                        got[k] = n
                counts[sid] = got
            fit = L.fit_misassignment(p["geom"], counts,
                                      random.Random(seed))
            for ax in ("row", "col"):
                tests += 1
                fitted += fit[ax]["basis"] == "fitted"
    check(f"fitted on at most 1 axis in 8 ({fitted} of {tests})",
          fitted <= tests / 8, f"{fitted} of {tests}")
    TABLE.append(f"| 96 wells x 8 plates, contamination from the whole "
                 f"plate, not along indexes | fitted on {fitted} of {tests} "
                 "axes | |")


def test_layout_needed() -> None:
    print("\nno layout, no model")
    fit = L.fit_misassignment({}, {}, random.Random(1))
    check("no positions: not applicable, and said", not fit["applicable"]
          and "no plate positions" in fit["note"])
    fit = L.fit_misassignment({"A1": ("A", 1), "A2": ("A", 2)},
                              {"A1": {"x": 9}}, random.Random(1))
    check("one row: not applicable", not fit["applicable"])


def test_expected_reads() -> None:
    """The expectation itself: rate x reads in the other wells of the
    row, plus rate x the column's; nothing from wells sharing neither."""
    print("\nexpected misassigned reads")
    geom = {"A1": ("A", 1), "A2": ("A", 2), "B1": ("B", 1),
            "B2": ("B", 2)}
    counts = {"A1": {"t": 1000}, "A2": {"t": 100}, "B1": {"t": 10},
              "B2": {}}
    sh = L.IndexSharing(geom, counts)
    check("sources: row, column, neither",
          sh.sources("B2", "t") == (10, 100, 1000), str(sh.sources("B2", "t")))
    check("expected = row rate x row + column rate x column",
          abs(sh.expected("B2", "t", 0.01, 0.1) - (0.01 * 10 + 0.1 * 100))
          < 1e-9)


def test_unused_combinations() -> None:
    print("\nunused combinations measure the rate directly")
    m = [{"index_i5": "a", "index_i7": "x", "reads": "100000",
          "assigned": "yes"},
         {"index_i5": "a", "index_i7": "y", "reads": "50", "assigned": "no"},
         {"index_i5": "b", "index_i7": "x", "reads": "300000",
          "assigned": "yes"},
         {"index_i5": "b", "index_i7": "y", "reads": "150",
          "assigned": "no"}]
    r = L.sentinel_rates(m)
    check("row rate from cells with an empty column",
          r["row"] and abs(r["row"]["rate"] - 5e-4) < 1e-12, str(r))
    check("no column measurement where no cell's column has reads",
          r["col"] is None)
    check("a full plate has none", L.sentinel_rates(
        [x for x in m if x["assigned"] == "yes"]) is None)


def test_plate_effects() -> None:
    """The shared position tests: a real row effect is found, noise is
    not, and too few wells are not tested at all."""
    print("\nplate position tests")
    rng = random.Random(1)
    rows = [r for r in "ABCDEFGH" for _c in range(12)]
    cols = [c for _r in "ABCDEFGH" for c in range(1, 13)]
    noise = [rng.random() for _ in rows]
    rowfx = [rng.random() + (2 if r in "AB" else 0) for r in rows]
    f = L.plate_effects(rows, cols, [("noise", noise), ("rowfx", rowfx),
                                     ("flat", [1.0] * len(rows))],
                        random.Random(2), 500)
    by = {x["quantity"]: x for x in f}
    check("a row effect is found", by["rowfx"]["row_p"] < 0.01,
          str(by["rowfx"]))
    check("noise is not", by["noise"]["row_p"] > 0.01
          and by["noise"]["col_p"] > 0.01, str(by["noise"]))
    check("a constant is not tested", "flat" not in by)
    check("under 12 wells nothing is tested",
          L.plate_effects(rows[:11], cols[:11], [("noise", noise[:11])],
                          random.Random(2), 100) == [])


def main() -> int:
    test_expected_reads()
    test_unused_combinations()
    test_layout_needed()
    test_plate_effects()
    test_recovery()
    test_no_misassignment()
    test_fallback()
    test_biology_down_the_columns()
    test_not_following_the_indexes()
    print("\n" + "\n".join(TABLE))
    print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
    for n, d in FAIL:
        print(f"  FAIL {n}: {d}")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
