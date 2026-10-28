#!/usr/bin/env python3
"""
Database-mode assignment for nanopore: one sample's reads against a
taxonomically annotated reference database release.

Panel mode asks "which of these few targets is this read?", and a read
can only ever be called something on the panel. Database mode asks "what
is in this sample?" against tens of thousands of 18S V4 records, and the
answer has to be honest about two limits that a panel never meets: the
database does not contain everything that will be amplified, and V4 does
not separate every pair of sister species. So the unit of classification
here is not the read, and the answer is not always a species.

Why not classify reads. Benchmarked with minimap2 against simulated reads
from genera and species held OUT of the index, an LCA over equally best
hits named 79-85% of reads from missing species or genera at species
level - a confident wrong name. No identity floor can fix that at read
level: at 3% error a read from a species that IS present sits at ~97%
identity, and a read from its missing sister at ~96%, because sister taxa
differ by about 1% over V4. The two distributions overlap. A consensus of
many reads has ~0.2% error, and there a species floor of ~0.995 kept 90%
of correct species calls while wrong calls on missing taxa fell to 5-8%.

So, per sample:

  1. Trim each read at the primers, which orients it and gives it the
     extent of a database insert.
  2. Map every read to the database (minimap2) and group reads by their
     primary hit.
  3. Build a consensus for each group, guided by its reference: pile the
     reads up on the reference through minimap2's cs tags, take the
     majority base at each position, re-align to that draft and repeat
     until it stops changing.
  4. Split hidden haplotypes by linkage. A group whose reads come from two
     sequences gives a consensus of the majority, and the minority is lost
     inside it. Their differences show as alleles that are too frequent to
     be read errors AND that occur together in the same reads: an error is
     independent of every other error, a haplotype's differences travel
     together. The largest linked set is the second haplotype's
     signature, its reads form a new group, and both are resolved again.
  5. Merge consensus sequences within one edit of each other: groups that
     were one taxon split across near-identical references converge on
     the same sequence.
  6. Count reads per consensus by mapping every read to the final set, so
     reads from groups too small for a consensus of their own are counted
     too - or left unclustered, if they fit their own database hit better
     than any consensus (a taxon too rare here to build one).
  7. Classify each consensus: LCA over the database hits within 0.5% of
     the best alignment score, then cut back to the deepest rank whose
     identity floor the best hit clears. The identity is reported beside
     the rank, so a near-identical relative is visible, not hidden, and a
     consensus that clears no floor down to family is reported as
     "unclassified, NN% to X" rather than given X's name.

The idea in step 4 - that diversity within a group is evidence of a
second sequence when it is LINKED - is Stephen Doyle's. On simulated
samples it raised hidden-haplotype recovery from 0% to 37% and exact
consensus for present taxa from 87% to 98%, for about three false splits
per 30 samples; merging cut spurious consensus sequences from 33 to 5-10.

Output is relative 18S abundance, never organisms or copies: 18S copy
number varies widely between taxa.

Standard-library Python plus minimap2, the pipeline's one external
dependency, used in database mode only.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter, defaultdict
from itertools import combinations

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ampseqlib import (  # noqa: E402
    RANKS, bounded_edit, load_primers, load_release_manifest, load_taxonomy,
    open_reads, read_fasta, read_table, release_version, writer,
)

# minimap2's index parameters for the database. k=19, w=10 gave the same
# accuracy as the map-ont default k=15 on the leave-out benchmark and ran
# 1.6x faster; the index builds in 1-2 s, so it is built once per run.
INDEX_ARGS = ["-x", "map-ont", "-k", "19", "-w", "10"]

# Floors on the best hit's identity, per rank. A consensus is named at a
# rank only if its best hit is at least this close; otherwise it stops at
# the next rank up. From the leave-out benchmark at consensus error rates:
# species 0.995 kept 90% of correct species calls and cut wrong calls on
# missing species to 5-8%. Genus was 0.98 in the prototype, which let
# 27% of consensus from missing genera take a wrong genus on a 150-sample
# leave-out run of v0.4.0; 0.99 cut that to 20% at no cost to present
# taxa, the rest moving honestly to family (docs/REFERENCE_DATABASE.md
# section 11.3, tests/bench_db_floors.py).
DEFAULT_FLOORS = {"species": 0.995, "genus": 0.99, "family": 0.95}

CS = re.compile(r"(:\d+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+)")


# --------------------------------------------------------------------------
# minimap2
# --------------------------------------------------------------------------


class Minimap2:
    """Runs minimap2 and parses its PAF. One instance per sample, so the
    executable and thread count are decided once."""

    def __init__(self, exe: str = "minimap2", threads: int = 1):
        path = shutil.which(exe) or (exe if os.path.isfile(exe) else None)
        if not path:
            raise SystemExit(
                f"minimap2 not found ('{exe}'). Database mode needs it; see "
                "docs/USAGE.md, field setup")
        self.exe, self.threads = path, max(1, threads)
        self.calls = 0

    def version(self) -> str:
        out = subprocess.run([self.exe, "--version"], capture_output=True,
                             text=True)
        return out.stdout.strip()

    def map(self, target: str, query: str, extra=(), preset=("-x", "map-ont"),
            cs: bool = False) -> list[dict]:
        args = [self.exe, *preset, "-c", "-t", str(self.threads),
                *(["--cs"] if cs else []), *extra, target, query]
        proc = subprocess.run(args, capture_output=True, text=True)
        self.calls += 1
        if proc.returncode != 0:
            raise SystemExit(f"minimap2 failed ({' '.join(args)}):\n"
                             f"{proc.stderr[-2000:]}")
        rows = []
        for line in proc.stdout.splitlines():
            f = line.split("\t")
            tags = {t[:2]: t[5:] for t in f[12:]}
            rows.append({
                "q": f[0], "qlen": int(f[1]), "qs": int(f[2]),
                "qe": int(f[3]), "strand": f[4], "t": f[5],
                "tlen": int(f[6]), "ts": int(f[7]), "te": int(f[8]),
                "match": int(f[9]), "alen": int(f[10]),
                "AS": int(tags.get("AS", 0)), "cs": tags.get("cs", ""),
                "tp": tags.get("tp", "P"),
            })
        return rows


def write_fasta(path: str, recs) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in recs:
            fh.write(f">{name}\n{seq}\n")


def identity(hit: dict) -> float:
    """Matching bases over the alignment block PLUS any of the query left
    unaligned at either end.

    minimap2's own ratio (column 10 over column 11) ignores the ends it
    clipped, so a consensus aligned over only part of its length would
    look as close as one aligned end to end. A consensus and a record are
    both primer-to-primer inserts, so an unaligned end is a real
    difference and is charged as one."""
    clip = hit["qs"] + (hit["qlen"] - hit["qe"])
    return hit["match"] / max(1, hit["alen"] + clip)


# --------------------------------------------------------------------------
# reference-guided consensus
# --------------------------------------------------------------------------


class Pileup:
    """The reads of one group laid over a draft sequence.

    base[p]   Counter of what the reads have at draft position p: a base,
              or '-' for a deletion;
    ins[p]    Counter of insertions placed before position p;
    cover[p]  how many reads span p;
    alleles   per read, its differences from the draft - substitutions
              and deletions only, keyed by draft position - which is what
              the linkage test compares between reads;
    span      per read, the draft interval it covers.
    """

    def __init__(self, ref: str, alns: list[dict]):
        n = len(ref)
        self.ref = ref
        self.base = [Counter() for _ in range(n)]
        self.ins = [Counter() for _ in range(n + 1)]
        self.cover = [0] * (n + 1)
        self.alleles: dict[str, dict[int, str]] = {}
        self.span: dict[str, tuple[int, int]] = {}
        for a in alns:
            p, al = a["ts"], {}
            start = p
            for op in CS.findall(a["cs"]):
                kind = op[0]
                if kind == ":":
                    for _ in range(int(op[1:])):
                        self.base[p][ref[p]] += 1
                        self.cover[p] += 1
                        p += 1
                elif kind == "*":
                    b = op[2].upper()
                    self.base[p][b] += 1
                    self.cover[p] += 1
                    al[p] = b
                    p += 1
                elif kind == "-":
                    for _ in range(len(op) - 1):
                        self.base[p]["-"] += 1
                        self.cover[p] += 1
                        al[p] = "-"
                        p += 1
                else:
                    self.ins[p][op[1:].upper()] += 1
            self.alleles[a["q"]] = al
            self.span[a["q"]] = (start, p)

    def consensus(self, n_reads: int) -> str:
        """Majority base at every position at least half the reads span,
        with an insertion wherever most reads spanning it carry one."""
        out = []
        need = max(1, n_reads // 2)
        for p in range(len(self.ref)):
            if self.cover[p] < need:
                continue
            if self.ins[p] and sum(self.ins[p].values()) > self.cover[p] / 2:
                out.append(self.ins[p].most_common(1)[0][0])
            b = self.base[p].most_common(1)[0][0]
            if b != "-":
                out.append(b)
        end = len(self.ref)
        if self.ins[end] and sum(self.ins[end].values()) > \
                self.cover[end - 1] / 2:
            out.append(self.ins[end].most_common(1)[0][0])
        return "".join(out)

    def error_rates(self) -> tuple[float, float]:
        """(substitution, deletion) rates of the reads against the draft.

        These set what counts as too frequent to be an error in the
        linkage test, and they are measured, not assumed: the run's
        chemistry, basecaller and quality are all in them. A second
        haplotype inflates them slightly, which only makes the test more
        conservative."""
        sub = dele = tot = 0
        for p, c in enumerate(self.base):
            tot += self.cover[p]
            for b, k in c.items():
                if b == "-":
                    dele += k
                elif b != self.ref[p]:
                    sub += k
        if not tot:
            return 0.02, 0.01
        return (min(0.2, max(0.002, sub / tot)),
                min(0.2, max(0.001, dele / tot)))


def binom_sf(k: int, n: int, p: float) -> float:
    """P(X >= k) for X ~ Binomial(n, p), summed in log space."""
    if k <= 0:
        return 1.0
    if p <= 0:
        return 0.0
    lp, lq = math.log(p), math.log1p(-p)
    tot = 0.0
    for i in range(k, n + 1):
        tot += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1)
                        - math.lgamma(n - i + 1) + i * lp + (n - i) * lq)
    return tot


def in_homopolymer(ref: str, p: int) -> bool:
    return (p > 0 and ref[p - 1] == ref[p]) or \
        (p + 1 < len(ref) and ref[p + 1] == ref[p])


class Settings:
    """The consensus and splitting parameters, in one place, so the
    summary can record exactly what a sample was resolved with."""

    def __init__(self, **kw):
        self.min_reads = kw.get("min_reads", 5)
        self.max_reads = kw.get("max_reads", 300)
        self.hap_min_frac = kw.get("hap_min_frac", 0.10)
        self.hap_p = kw.get("hap_p", 1e-6)
        self.min_cooccur = kw.get("min_cooccur", 3)
        self.cooccur_ratio = kw.get("cooccur_ratio", 3.0)
        self.max_split_depth = kw.get("max_split_depth", 2)
        self.max_iter = kw.get("max_iter", 4)
        self.merge_edits = kw.get("merge_edits", 1)
        self.split = kw.get("split", True)

    def as_dict(self) -> dict:
        return dict(vars(self))


def signature(pu: Pileup, st: Settings) -> list[tuple[int, str]] | None:
    """The linked set of non-consensus alleles marking a second
    haplotype, or None.

    A candidate allele must be at least hap_min_frac of the reads at its
    position, and improbable as read error under the measured rate
    (binomial tail below hap_p). Deletions in homopolymers are passed
    over: that is where nanopore indels concentrate and where their
    placement smears between neighbouring positions.

    Two candidates are linked when at least min_cooccur reads spanning
    both carry both, and that is at least cooccur_ratio times what
    independence would give. Read errors are independent of each other;
    a haplotype's differences are not. The largest connected set of
    linked alleles is the signature - two or more alleles, because one
    frequent allele on its own cannot be told from a systematic error.
    """
    sub_rate, del_rate = pu.error_rates()
    ref = pu.ref
    floor = max(st.min_reads * 2, 10)
    cand = []
    for p in range(len(ref)):
        c = pu.cover[p]
        if c < floor:
            continue
        for b, k in pu.base[p].items():
            if b == ref[p] or k < max(st.min_cooccur, st.hap_min_frac * c):
                continue
            if b == "-" and in_homopolymer(ref, p):
                continue
            if binom_sf(k, c, sub_rate / 3 if b != "-" else del_rate) \
                    > st.hap_p:
                continue
            cand.append((p, b))
    if len(cand) < 2:
        return None
    carriers = {x: {q for q, al in pu.alleles.items() if al.get(x[0]) == x[1]}
                for x in cand}
    adj: dict = defaultdict(set)
    for x, y in combinations(cand, 2):
        lo, hi = min(x[0], y[0]), max(x[0], y[0])
        cov = {q for q, (s, e) in pu.span.items() if s <= lo and e > hi}
        if len(cov) < 10:
            continue
        a, b = carriers[x] & cov, carriers[y] & cov
        both = len(a & b)
        expect = len(a) * len(b) / len(cov)
        if both >= st.min_cooccur and both >= st.cooccur_ratio * max(expect, 0.5):
            adj[x].add(y)
            adj[y].add(x)
    best: list = []
    seen: set = set()
    for v in sorted(adj):
        if v in seen:
            continue
        comp, stack = set(), [v]
        while stack:
            u = stack.pop()
            if u in comp:
                continue
            comp.add(u)
            stack.extend(adj[u] - comp)
        seen |= comp
        if len(comp) > len(best):
            best = sorted(comp)
    return best if len(best) >= 2 else None


class ConsensusBuilder:
    """Reference-guided consensus for the groups of one sample."""

    def __init__(self, mm: Minimap2, reads: dict[str, str], tmp: str,
                 st: Settings):
        self.mm, self.reads, self.tmp, self.st = mm, reads, tmp, st
        self.splits = 0

    def build(self, qs: list[str], start: str) -> tuple[str, Pileup | None]:
        """Iterate pileup consensus from `start` until it is stable;
        return it and the reads' pileup against it."""
        draft = start
        qfa = os.path.join(self.tmp, "group.fa")
        tfa = os.path.join(self.tmp, "draft.fa")
        write_fasta(qfa, ((q, self.reads[q]) for q in qs))
        pu = None
        for _ in range(self.st.max_iter + 1):
            write_fasta(tfa, [("draft", draft)])
            alns = [a for a in self.mm.map(tfa, qfa, cs=True,
                                           extra=("--secondary=no",))
                    if a["tp"] == "P"]
            if not alns:
                return draft, None
            pu = Pileup(draft, alns)
            new = pu.consensus(len(alns))
            if new == draft or not new:
                return draft, pu
            draft = new
        return draft, pu

    def resolve(self, qs: list[str], start: str, depth: int = 0,
                origin: str = "") -> list[dict]:
        cons, pu = self.build(qs, start)
        if pu is None:
            return []
        if self.st.split and depth < self.st.max_split_depth:
            sig = signature(pu, self.st)
            if sig:
                need = (len(sig) + 1) // 2
                hap = [q for q in qs
                       if sum(1 for p, b in sig
                              if pu.alleles.get(q, {}).get(p) == b) >= need]
                hs = set(hap)
                rest = [q for q in qs if q not in hs]
                if (len(hap) >= max(self.st.min_reads,
                                    self.st.hap_min_frac * len(qs))
                        and len(rest) >= self.st.min_reads):
                    self.splits += 1
                    tag = f"{origin}split{depth + 1}"
                    return (self.resolve(hap, cons, depth + 1, tag + "a.")
                            + self.resolve(rest, cons, depth + 1, tag + "b."))
        return [{"consensus": cons, "reads": list(qs),
                 "origin": origin.rstrip(".") or "group"}]


def merge_close(clusters: list[dict], max_edits: int) -> list[dict]:
    """Merge consensus sequences within max_edits of each other, largest
    first, so the survivor is the better-supported sequence."""
    merged: list[dict] = []
    for c in sorted(clusters, key=lambda x: (-len(x["reads"]), x["consensus"])):
        for m in merged:
            if abs(len(m["consensus"]) - len(c["consensus"])) <= max_edits and \
                    bounded_edit(c["consensus"], m["consensus"], max_edits) \
                    <= max_edits:
                m["reads"].extend(c["reads"])
                m["merged"] += 1
                break
        else:
            merged.append({**c, "reads": list(c["reads"]), "merged": 0})
    return merged


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------


def lca_depth(names, tax: dict) -> int:
    """Index into RANKS of the deepest rank at which every record agrees;
    -1 when they disagree even at domain. A blank rank (a PR2
    placeholder) breaks agreement: no record can vouch for a name it
    does not have."""
    depth = -1
    for d, r in enumerate(RANKS):
        vals = {tax.get(t, {}).get(r, "") for t in names}
        if len(vals) == 1 and "" not in vals:
            depth = d
        else:
            break
    return depth


def apply_floors(depth: int, ident: float, floors: dict) -> int:
    """Cut a rank back until the identity clears its floor. Ranks
    without a floor are kept as the LCA gives them."""
    for r in ("species", "genus", "family", "order", "class"):
        if r in floors and depth >= RANKS.index(r) and ident < floors[r]:
            depth = RANKS.index(r) - 1
    return depth


def classify_hits(hits: list[dict], tax: dict, floors: dict,
                  margin: float) -> dict:
    """One consensus's call from its database hits."""
    if not hits:
        return {"rank": "", "depth": -1, "identity": 0.0, "coverage": 0.0,
                "best_hit": "", "near_hits": 0, "lca_depth": -1}
    top = max(h["AS"] for h in hits)
    near = [h for h in hits if h["AS"] >= top * (1 - margin)]
    best = max(near, key=lambda h: (h["AS"], identity(h), h["t"]))
    ident = identity(best)
    names = sorted({h["t"] for h in near})
    # A record with no taxonomy at all - v0.4.0 has nine, conflicts its
    # build left unresolved - can vouch for nothing, and as a member of
    # the LCA it vetoed everything: human 18S tied with one stopped at
    # no rank. It is left out of the LCA, unless it is all there is.
    known = [t for t in names if tax.get(t, {}).get("domain")]
    lca = lca_depth(known or names, tax)
    depth = apply_floors(lca, ident, floors)
    return {
        "depth": depth, "rank": RANKS[depth] if depth >= 0 else "",
        "identity": ident,
        "coverage": (best["qe"] - best["qs"]) / max(1, best["qlen"]),
        "best_hit": best["t"], "near_hits": len(names),
        "lca_depth": lca,
    }


def deepest_named(lineage: dict) -> str:
    for r in reversed(RANKS):
        if lineage.get(r):
            return lineage[r]
    return ""


# PR2 fills ranks that have no name with its parent's name and an _X per
# level (Chromadorea_X, Craniata_XX), and leaves a species unnamed as
# Genus_sp.
_PLACEHOLDER = re.compile(r"_X+$")


def is_placeholder(name: str) -> bool:
    return bool(_PLACEHOLDER.search(name or ""))


def label_for(call: dict, tax: dict, label_depth: int) -> tuple[str, str]:
    """(taxon at the rank reached, label for the report).

    Down to label_depth (family) a call is the taxon's name. Above it
    the consensus cleared no floor that would let it be named: it is
    "unclassified, NN% to X", where X is the best hit's genus (or its
    deepest named rank) - so a gap in the database reads as a gap, and
    how far off the nearest record is stays on the page.

    Where the name at the rank reached is one of PR2's placeholders, the
    label says what it means instead of repeating it. "Chromadorea_X" in
    a report reads as a taxon; it is "no family name recorded, somewhere
    in Chromadorea" - on a real VESPA run it was mostly Haemonchus at
    100%, stopped at family by a tie with a second genus, and the label
    said nothing a reader could use. Now: "Chromadorea, unnamed family
    (nearest Haemonchus)". The taxon column keeps PR2's name, so tables
    still join on it.

    And a call stopped above family by a TIE between records - the floors
    let it through, the LCA did not - is "X, unresolved below <rank>
    (NN% to Y)", not "unclassified": the database holds the sequence and
    cannot choose between records, which is not the same as a gap.""" 
    lin = tax.get(call["best_hit"], {})
    d = call["depth"]
    taxon = lin.get(RANKS[d], "") if d >= 0 else ""
    near = lin.get("genus") or deepest_named(lin) or call["best_hit"] or "nothing"
    if d >= label_depth:
        if not is_placeholder(taxon):
            return taxon, taxon
        named = next((lin.get(RANKS[k]) for k in range(d - 1, -1, -1)
                      if lin.get(RANKS[k]) and not is_placeholder(lin.get(RANKS[k]))),
                     "")
        nearest = "" if not near or is_placeholder(near) or near == named \
            else f" (nearest {near})"
        return taxon, f"{named or taxon}, unnamed {RANKS[d]}{nearest}"
    if not call["best_hit"]:
        return taxon, "unclassified, no database hit"
    if call.get("lca_depth", -2) == d and d >= 0:
        # Not stopped by an identity floor but by records disagreeing
        # (or the two halves of a pair pointing at different ones):
        # the database has the sequence, it cannot say which it is.
        # "Unclassified" would claim a gap that is not there.
        named = next((lin.get(RANKS[k]) for k in range(d, -1, -1)
                      if lin.get(RANKS[k]) and not is_placeholder(lin.get(RANKS[k]))),
                     taxon)
        return taxon, (f"{named}, unresolved below {RANKS[d]} "
                       f"({100 * call['identity']:.0f}% to {near})")
    return taxon, f"unclassified, {100 * call['identity']:.0f}% to {near}"


# --------------------------------------------------------------------------
# the sample
# --------------------------------------------------------------------------


def parse_floors(text: str) -> dict:
    out = {}
    for part in (text or "").split(","):
        if not part.strip():
            continue
        k, _, v = part.partition("=")
        k = k.strip()
        if k not in RANKS:
            raise SystemExit(f"--floors: '{k}' is not a rank ({', '.join(RANKS)})")
        out[k] = float(v)
    return out


def trim_reads(path: str, primers) -> tuple[dict[str, str], Counter]:
    """Oriented, primer-trimmed inserts by read id, and what was lost.

    The release's records are primer-trimmed inserts, so reads must be
    too: identity is measured over the insert, and a primer left on one
    side and not the other would be charged as a difference."""
    out: dict[str, str] = {}
    fate: Counter = Counter()
    from ampseqlib import revcomp
    for rid, seq, _q in open_reads(path):
        fate["reads_in"] += 1
        seq = seq.upper()
        hit = None
        for s in (seq, revcomp(seq)):
            f = primers.match(s, "fwd")
            r = primers.match(revcomp(s), "rev")
            if f is not None and r is not None:
                hit = s[f[1]:len(s) - r[1]]
                break
        if hit is None:
            fate["no_primer"] += 1
            continue
        if len(hit) < 50:
            fate["too_short"] += 1
            continue
        # Duplicate ids happen when a split read keeps its parent's name;
        # minimap2 would merge them, so each gets its own.
        key = rid if rid not in out else f"{rid}_{fate['reads_in']}"
        out[key] = hit
    return out, fate


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--reads", required=True, help="the sample's FASTQ")
    ap.add_argument("--database", required=True,
                    help="reference database release directory")
    ap.add_argument("--index", help="minimap2 index of the release's "
                    "refdb.fasta, built with " + " ".join(INDEX_ARGS) +
                    " (default: build from refdb.fasta on the fly)")
    ap.add_argument("--primers", required=True)
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--outdir", default="assign_db")
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--min-cluster-reads", type=int, default=5,
                    help="reads a group needs for a consensus of its own")
    ap.add_argument("--max-consensus-reads", type=int, default=300,
                    help="reads drawn from a larger group to build its "
                    "consensus; every read is still counted")
    ap.add_argument("--haplotype-min-frac", type=float, default=0.10)
    ap.add_argument("--no-split", action="store_true",
                    help="do not split groups into haplotypes")
    ap.add_argument("--merge-edits", type=int, default=1)
    ap.add_argument("--lca-margin", type=float, default=0.005,
                    help="hits within this fraction of the best alignment "
                    "score take part in the LCA")
    ap.add_argument("--floors", default=",".join(
        f"{k}={v}" for k, v in DEFAULT_FLOORS.items()),
        help="identity floor per rank, e.g. species=0.995,genus=0.99")
    ap.add_argument("--recruit-margin", type=int, default=10,
                    help="a read is counted to a consensus unless its own "
                    "database hit scores more than this much better "
                    "(minimap2 alignment score)")
    ap.add_argument("--primer-max-edits", type=int, default=3,
                    help="edits allowed when finding each primer (panel mode "
                    "uses 2). 3, because the degenerate VESPA reverse primer "
                    "meets most bacterial 16S at three: at 2 a fifth of a "
                    "real run's reads - nearly all off-target 16S - were "
                    "dropped as no_primer instead of counted as bacteria")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args(argv)
    t0 = time.time()

    sid = args.sample_id
    os.makedirs(args.outdir, exist_ok=True)
    man = load_release_manifest(args.database)
    version = release_version(man)
    tax = load_taxonomy(os.path.join(args.database, "taxonomy.tsv"))
    resolution = {r["id"]: r.get("deepest_rank", "")
                  for r in read_table(os.path.join(args.database,
                                                   "resolution_nanopore.tsv"))}
    db_fasta = os.path.join(args.database, "refdb.fasta")
    index = args.index or db_fasta
    floors = parse_floors(args.floors)
    st = Settings(min_reads=args.min_cluster_reads,
                  max_reads=args.max_consensus_reads,
                  hap_min_frac=args.haplotype_min_frac,
                  merge_edits=args.merge_edits, split=not args.no_split)
    mm = Minimap2(args.minimap2, args.threads)

    primers = load_primers(args.primers)
    primers.max_edits = args.primer_max_edits
    reads, fate = trim_reads(args.reads, primers)
    tmp = tempfile.mkdtemp(prefix=f"assign_db_{sid}_", dir=args.outdir)
    clusters: list[dict] = []
    db_primary: dict[str, dict] = {}
    try:
        rfa = os.path.join(tmp, "reads.fa")
        write_fasta(rfa, reads.items())
        groups: dict[str, list[str]] = defaultdict(list)
        if reads:
            for a in mm.map(index, rfa, preset=INDEX_ARGS,
                            extra=("--secondary=no",)):
                if a["tp"] == "P" and a["q"] not in db_primary:
                    db_primary[a["q"]] = a
                    groups[a["t"]].append(a["q"])
        fate["unmapped"] = len(reads) - len(db_primary)

        # Consensus per group large enough to have one. The reference
        # sequence is the starting draft; the database is read once, for
        # just the records that head a group.
        heads = {t for t, qs in groups.items() if len(qs) >= st.min_reads}
        refseq = {i: s for i, s in read_fasta(db_fasta) if i in heads}
        cb = ConsensusBuilder(mm, reads, tmp, st)
        rng = random.Random(f"{args.seed}:{sid}")
        built = []
        for t in sorted(heads, key=lambda t: (-len(groups[t]), t)):
            qs = sorted(groups[t])
            if len(qs) > st.max_reads:
                qs = sorted(rng.sample(qs, st.max_reads))
            for c in cb.resolve(qs, refseq[t]):
                c["group"] = t
                built.append(c)
        clusters = merge_close(built, st.merge_edits)

        # Count every read against the final set.
        counts: Counter = Counter()
        member: dict[str, int] = {}
        if clusters:
            cfa = os.path.join(tmp, "consensus.fa")
            write_fasta(cfa, ((f"c{k}", c["consensus"])
                              for k, c in enumerate(clusters)))
            best: dict[str, dict] = {}
            for a in mm.map(cfa, rfa, extra=("--secondary=no",)):
                if a["tp"] == "P" and a["q"] not in best:
                    best[a["q"]] = a
            for q, a in best.items():
                own = db_primary.get(q)
                if own is not None and own["AS"] - a["AS"] > args.recruit_margin:
                    continue
                k = int(a["t"][1:])
                counts[k] += 1
                member[q] = k
        fate["clustered"] = len(member)
        fate["unclustered"] = len(reads) - len(member)

        # Classify the consensus sequences that kept any reads.
        keep = [k for k in range(len(clusters)) if counts[k] > 0]
        calls: dict[int, dict] = {}
        if keep:
            kfa = os.path.join(tmp, "classify.fa")
            write_fasta(kfa, ((f"c{k}", clusters[k]["consensus"]) for k in keep))
            hits: dict[str, list] = defaultdict(list)
            for a in mm.map(index, kfa, preset=INDEX_ARGS,
                            extra=("--secondary=yes", "-N", "1000",
                                   "-p", "0.9")):
                hits[a["q"]].append(a)
            for k in keep:
                calls[k] = classify_hits(hits.get(f"c{k}", []), tax, floors,
                                         args.lca_margin)
        mm_version = mm.version()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # ---- outputs ---------------------------------------------------------
    label_depth = RANKS.index("family")
    total = sum(counts.values())
    order = sorted(keep, key=lambda k: (-counts[k], clusters[k]["consensus"]))
    rows = []
    for n, k in enumerate(order, 1):
        c, call = clusters[k], calls[k]
        taxon, label = label_for(call, tax, label_depth)
        lin = tax.get(call["best_hit"], {})
        d = call["depth"]
        rows.append({
            "sample_id": sid, "cluster": f"{sid}_c{n}",
            "reads": counts[k],
            "rel_abundance": f"{counts[k] / total:.6f}" if total else "",
            "rank": call["rank"], "taxon": taxon, "label": label,
            "identity": f"{call['identity']:.4f}",
            "coverage": f"{call['coverage']:.3f}",
            "lineage": ";".join(lin.get(r, "") for r in RANKS[:d + 1]),
            "best_hit": call["best_hit"],
            "best_hit_lineage": ";".join(lin.get(r, "") for r in RANKS),
            "near_hits": call["near_hits"],
            "lca_rank": RANKS[call["lca_depth"]] if call["lca_depth"] >= 0 else "",
            "db_resolution": resolution.get(call["best_hit"], ""),
            "origin": c["origin"], "merged": c["merged"],
            "length": len(c["consensus"]),
            "database_version": version,
            "consensus": c["consensus"],
        })
    cols = ["sample_id", "cluster", "reads", "rel_abundance", "rank", "taxon",
            "label", "identity", "coverage", "lineage", "best_hit",
            "best_hit_lineage", "near_hits", "lca_rank", "db_resolution",
            "origin", "merged", "length", "database_version", "consensus"]
    with open(os.path.join(args.outdir, f"{sid}.clusters.tsv"), "w",
              encoding="utf-8") as fh:
        w = writer(fh, cols)
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(args.outdir, f"{sid}.consensus.fasta"), "w",
              encoding="utf-8") as fh:
        for r in rows:
            fh.write(f">{r['cluster']} reads={r['reads']} rank={r['rank'] or 'none'} "
                     f"taxon={r['taxon'] or 'none'} identity={r['identity']}\n"
                     f"{r['consensus']}\n")

    # Per taxon at the rank each call reached. Two clusters of one
    # species (haplotypes) are one row; a species call and a genus call
    # in the same genus are NOT merged - they are different statements.
    agg: dict = {}
    for r in rows:
        key = (r["rank"], r["taxon"], r["label"])
        a = agg.setdefault(key, {"reads": 0, "clusters": 0, "identity": 0.0,
                                 "lineage": r["lineage"]})
        a["reads"] += r["reads"]
        a["clusters"] += 1
        a["identity"] = max(a["identity"], float(r["identity"]))
    with open(os.path.join(args.outdir, f"{sid}.taxa.tsv"), "w",
              encoding="utf-8") as fh:
        w = writer(fh, ["sample_id", "rank", "taxon", "label", "reads",
                        "rel_abundance", "clusters", "best_identity",
                        "lineage", "database_version"])
        w.writeheader()
        for (rank, taxon, label), a in sorted(agg.items(),
                                              key=lambda kv: -kv[1]["reads"]):
            w.writerow({"sample_id": sid, "rank": rank, "taxon": taxon,
                        "label": label, "reads": a["reads"],
                        "rel_abundance": f"{a['reads'] / total:.6f}" if total else "",
                        "clusters": a["clusters"],
                        "best_identity": f"{a['identity']:.4f}",
                        "lineage": a["lineage"], "database_version": version})

    summary = {
        "sample_id": sid,
        "mode": "database",
        "database": {"version": version, "name": man.get("name", ""),
                     "path": os.path.abspath(args.database)},
        "minimap2": mm_version,
        "elapsed_s": round(time.time() - t0, 2),
        "reads": {k: fate.get(k, 0) for k in
                  ("reads_in", "no_primer", "too_short", "unmapped",
                   "clustered", "unclustered")},
        "reads_note": (
            "reads_in = no_primer + too_short + clustered + unclustered. "
            "unmapped reads (no database hit) are counted within "
            "clustered or unclustered, depending on whether a consensus "
            "took them. rel_abundance in the cluster and taxon tables is "
            "over clustered reads: relative 18S abundance, never organisms "
            "or copies."),
        "clusters": len(rows),
        "groups": len(heads),
        "haplotype_splits": cb.splits,
        "minimap2_calls": mm.calls,
        "settings": {**st.as_dict(), "floors": floors,
                     "lca_margin": args.lca_margin,
                     "recruit_margin": args.recruit_margin,
                     "index_args": INDEX_ARGS},
        "rank_reached": dict(Counter(r["rank"] or "none" for r in rows)),
    }
    with open(os.path.join(args.outdir, f"{sid}.assign_db.json"), "w",
              encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"{sid}: {fate['reads_in']} reads, {len(rows)} clusters, "
          f"{fate['clustered']} clustered, {fate['unclustered']} unclustered, "
          f"{cb.splits} haplotype split(s), {summary['elapsed_s']} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
