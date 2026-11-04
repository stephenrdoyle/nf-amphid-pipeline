# Index redesign — 18bp composite indexes

Design note. Of the implementation steps in §7, step 1 (separate anchor
and junction tolerances) is in the pipeline and step 2 (measure the anchor
split) has been run: see §5.1. Nothing else here is implemented.

Every number below is either **measured** on the 36.9M-read reference
MinION run, **simulated** against the real index sequences, or
**unverified** and flagged as such. Nothing is asserted from the
literature or from memory.

---

## 1. The decision

Extend each existing 8bp Illumina index with a 10bp UDP index, giving an
18bp composite read as **one word** (`barcode_repeats = 1`).

- **Oligo cost: unchanged.** Still 16 i5 + 24 i7 = 40 primers, just
  longer. The alternative that would fix more — one unique index per
  sample — needs ~385 primers and is not affordable.
- **Illumina compatibility: unchanged.** The index read is 8 cycles and
  stops inside X, so the sequencer never sees the extension.
- **Indel tolerance: 0 → 3.**

## 2. Why, in one table

Sequence-Levenshtein is the metric, not Hamming and not plain
Levenshtein. A barcode is read *in context* — a deletion inside it pulls
a base in from the adapter — and seq-Lev is the only one of the three
that models that. It is what `design_indices.py` should be scoring and
currently is not.

| | 8bp X alone | 16bp X+X today | **18bp X+Y** |
|---|---|---|---|
| i5 min seq-Lev | 2 → 0 indels | 4 → 1 indel | **8 → 3 indels** |
| i7 min seq-Lev | 2 → 0 indels | 3 → 1 indel | **8 → 3 indels** |
| derived substitution tolerance | 1 | 1 | **3–4** |
| rotational collisions | — | 4 (i7) | **0** |
| max homopolymer | 2 | 2 | **2** |

The tolerance and the rotational-collision count both come out of the
set geometry automatically — `BarcodeSet.build` derives them, so no
parameter needs changing.

## 3. Layout — UNVERIFIED, confirm before ordering

The Illumina index read starts at the adapter and reads inward, so X
must be **adapter-proximal at both ends**. That makes the composite
order differ between i5 and i7:

```
i5:  P5 ---- [X 8bp][Y 10bp] ---- [i5 anchor] ---- amplicon
i7:  amplicon ---- [i7 anchor] ---- [Y 10bp][X 8bp] ---- rc(P7)
```

So the word the demultiplexer reads outward from the anchor is `X+Y` at
the i5 end and `Y+X` at the i7 end. **This is inferred from the
construct in `assets/assay/ampseq_helminth.json`, not read off an oligo
spec.** If the i7 index read actually starts from the other side, the
i7 table below needs flipping — same sequences, different order. The
distances were optimised for the orientation shown, so a flip means
re-running the optimisation.

## 4. The indexes

X is unchanged and stays in its existing plate position. Only Y is new.

### i5 — plate row, composite = X + Y

| name | row | X (Illumina, 8bp) | Y (UDP, 10bp) | from | composite 18bp |
|---|---|---|---|---|---|
| S502 | A | CTCTCTAT | AACACGTGGA | UDP0073 | CTCTCTATAACACGTGGA |
| S503 | B | TATCCTCT | GTGCGTCCTT | UDP0052 | TATCCTCTGTGCGTCCTT |
| S505 | C | GTAAGGAG | GCTTACGGAC | UDP0089 | GTAAGGAGGCTTACGGAC |
| S506 | D | ACTGCATA | TTGGCTCCGC | UDP0063 | ACTGCATATTGGCTCCGC |
| S507 | E | AAGGAGTA | GTTGATAGTG | UDP0009 | AAGGAGTAGTTGATAGTG |
| S508 | F | CTAAGCCT | GAACCATGAA | UDP0090 | CTAAGCCTGAACCATGAA |
| S510 | G | CGTCTAAT | TCAATCCATT | UDP0027 | CGTCTAATTCAATCCATT |
| S511 | H | TCTCTCCG | AGACACATTA | UDP0083 | TCTCTCCGAGACACATTA |
| S513 | I | TCGACTAG | CGGTGGCGAA | UDP0047 | TCGACTAGCGGTGGCGAA |
| S515 | J | TTCTAGCT | CGTTGCTTAC | UDP0077 | TTCTAGCTCGTTGCTTAC |
| S516 | K | CCTAGAGT | GATTAAGGTG | UDP0094 | CCTAGAGTGATTAAGGTG |
| S517 | L | GCGTAAGA | TCGCCGCTAG | UDP0084 | GCGTAAGATCGCCGCTAG |
| S518 | M | CTATTAAG | GATACCTCCT | UDP0087 | CTATTAAGGATACCTCCT |
| S520 | N | AAGGCTAT | ACTAGCCGTG | UDP0034 | AAGGCTATACTAGCCGTG |
| S521 | O | GAGCCTTA | GAGTCTCTCC | UDP0054 | GAGCCTTAGAGTCTCTCC |
| S522 | P | TTATGCGA | GATCACCGCG | UDP0037 | TTATGCGAGATCACCGCG |

### i7 — plate column, composite = Y + X

| name | col | X (Illumina, 8bp) | Y (UDP, 10bp) | from | composite 18bp |
|---|---|---|---|---|---|
| N701 | 1 | TAAGGCGA | TATCATGAGA | UDP0073 | TATCATGAGATAAGGCGA |
| N702 | 2 | CGTACTAG | AGCACATCCT | UDP0067 | AGCACATCCTCGTACTAG |
| N703 | 3 | AGGCAGAA | GCCGTAACCG | UDP0016 | GCCGTAACCGAGGCAGAA |
| N704 | 4 | TCCTGAGC | GGATACCAGA | UDP0077 | GGATACCAGATCCTGAGC |
| N705 | 5 | GGACTCCT | TACGGCCGGT | UDP0036 | TACGGCCGGTGGACTCCT |
| N706 | 6 | TAGGCATG | CTGGATATGT | UDP0093 | CTGGATATGTTAGGCATG |
| N707 | 7 | CTCTCTAC | CAATCTATGA | UDP0094 | CAATCTATGACTCTCTAC |
| N710 | 8 | CGAGGCTG | GTAAGGCATA | UDP0083 | GTAAGGCATACGAGGCTG |
| N711 | 9 | AAGAGGCA | CGCTTGAAGT | UDP0090 | CGCTTGAAGTAAGAGGCA |
| N712 | 10 | GTAGAGGA | ATTGCCGAGT | UDP0041 | ATTGCCGAGTGTAGAGGA |
| N714 | 11 | GCTCATGA | TGCGTGTCAC | UDP0031 | TGCGTGTCACGCTCATGA |
| N715 | 12 | ATCTCAGG | GTCGATTACA | UDP0037 | GTCGATTACAATCTCAGG |
| N716 | 13 | ACTCGCTA | AGAGTGCGGC | UDP0029 | AGAGTGCGGCACTCGCTA |
| N718 | 14 | GGAGCTAC | CCATCCACGC | UDP0076 | CCATCCACGCGGAGCTAC |
| N719 | 15 | GCGTAGTA | TGCCTTGATC | UDP0008 | TGCCTTGATCGCGTAGTA |
| N720 | 16 | CGGAGCCT | TGGACTCGTA | UDP0072 | TGGACTCGTACGGAGCCT |
| N721 | 17 | TACGCTGC | ACCGCGCAAT | UDP0063 | ACCGCGCAATTACGCTGC |
| N722 | 18 | ATGCGCAG | ACAATAGAGT | UDP0056 | ACAATAGAGTATGCGCAG |
| N723 | 19 | TAGCGCTC | AACAAGGCGT | UDP0023 | AACAAGGCGTTAGCGCTC |
| N724 | 20 | ACTGAGCG | TCCTATTGTG | UDP0005 | TCCTATTGTGACTGAGCG |
| N726 | 21 | CCTAAGAC | GCGCGATGTT | UDP0006 | GCGCGATGTTCCTAAGAC |
| N727 | 22 | CGATCAGT | GTATGTAGAA | UDP0024 | GTATGTAGAACGATCAGT |
| N728 | 23 | TGCAGCTA | TTGTAACGGT | UDP0082 | TTGTAACGGTTGCAGCTA |
| N729 | 24 | TCGACGTC | GGTGGAATAC | UDP0095 | GGTGGAATACTCGACGTC |

Y values are drawn from the IDT/Illumina UD 96-plex set, adapter
orientation. A Y may repeat between the i5 and i7 tables — they are
separate sets and never compared to each other.

## 5. Evidence

### Measured on the reference run (36,877,899 reads)

Per-end failure budget, as % of all ends:

| | i5 | i7 |
|---|---|---|
| resolved (ok + ok_single_copy) | 59.88% | 51.15% |
| **matching failure** (no_copy_matched + copies_disagreed) | 19.71% | 10.28% |
| **no barcode attempted** (anchor_not_found + barcode_off_read) | 20.40% | **38.57%** |

Tier A 26.32%, tier B 51.67%, C 0.25%, D 18.41%, E 3.35%.

The index set can only address the *matching failure* row. On i7 the
anchor row is **3.8× larger** than the row this redesign fixes.

### Simulated, on the real sequences

Exactly k indels, no substitutions, window extracted from the anchor the
way `call_end` does it:

| k | kind | A 8bp×2 ok | A misassigned | B 18bp ok | B misassigned |
|---|---|---|---|---|---|
| 1 | del | 57.3% | 0.10% | **100.0%** | **0.00%** |
| 2 | del | 33.5% | 1.77% | 37.7% | 0.00% |
| 3 | del | 18.8% | 1.10% | 6.4% | 0.00% |

**One indel is the common case and B recovers all of it.** At 3+ indels
B is worse, for two reasons worth knowing: the matcher's `cap=3`
declines anything further out regardless of geometry (raising it to 4 is
licensed by seq-Lev 8 and untested), and the tandem design has a real
advantage there — damage concentrated in one 8bp copy leaves the other
clean, which a single word cannot do.

Misassignment is **0.00% in every condition tested**. Scheme A runs at
0.1–3%. For a diagnostic that matters more than the k=3 column.

Under a mixed error model, scheme A reproduced the observed i7 matching
failure (16.7%) to within half a point at 10% error/base without being
tuned to. That is the main reason to trust the projections.

### Projection

Tier A 26.32% → **~34%** from the composite alone. With anchorless
scanning recovering part of the anchor bucket:

| anchor bucket recovered | tier A | vs today |
|---|---|---|
| 0% (composite alone) | 34.3% | 1.3× |
| 50% | 46.4% | 1.8× |
| 75% | 53.0% | 2.0× |

**Not included:** the window grows 16bp → 18bp, so `barcode_off_read`
(3.94% i5 / 9.51% i7) rises by an unquantified amount. Needs the
read-length distribution near the truncation boundary.

**Rows 2–4 now look out of reach.** §5.1 measured the anchor bucket: at
most 1–4% of its ends carry a readable index. The composite-alone row is
the realistic target.

### 5.1 Measured: what the anchor bucket holds (2026-11-04)

Run on 854,515 split TITR01 reads (one of the four laptop-slice batch
files, `0000_FBH87365_pass_dd4877b8_8e846a4f_0`), with today's tandem 8bp
indexes. Here `anchor_not_found` is 21.2% of i5 ends and 24.7% of i7 ends,
lower than the reference run's 38.57% i7.

**Anchor tolerance sweep** (`--anchor-max-edits`, after the step 1 fix):

| edits | i5 no anchor | i7 no anchor | i5 index read | i7 index read | i7 no_copy_matched | tier A |
|---|---|---|---|---|---|---|
| 2 (default) | 21.2% | 24.7% | 62.6% | 51.6% | 10.7% | 32.0% |
| 3 | 18.1% | 21.4% | 63.6% | 53.1% | 11.9% | 33.4% |
| 4 | 7.0% | 16.6% | 65.3% | 53.7% | 15.8% | 34.7% |
| 6 | 0.3% | 1.1% | 65.4% | 55.2% | 28.2% | 35.7% |

`anchor_not_found` collapses, but readable indexes barely move: at 6
edits the i7 bucket falls 23.6 points and only 3.6 become an index; 17.5
move to `no_copy_matched`. A 14–15bp anchor at 6 edits matches chance
sequence in almost any 120bp window, so most of those "anchors" are
spurious and the window beside them holds no index.

The same sweep on a simulated TITR01 plate (200 reads/well, scored
against truth; clean reads only):

| edits | HAC: right well | HAC: wrong well | fast: right well | fast: wrong well |
|---|---|---|---|---|
| 2 | 89.4% | 18 | 80.3% | 123 |
| 3 | 89.9% | 19 | 84.9% | 132 |
| 4 | 89.9% | 20 | 85.5% | 136 |
| 6 | 89.7% | 56 | 85.3% | 196 |

Wrong-well calls roughly track extra assignments up to 4 edits and
triple at 6. The simulator damages anchors only by base errors, never by
losing an end, so it prices the looser search and does not measure its
benefit. The real run cannot price it: TITR01 uses 94 of 96 wells, so a
wrong pairing almost always names a real well and tier C stays at 0.2%.
A real wrong-well check needs a plate with unused combinations (VONT01,
50 of 96).

**Outer landmark.** The index copies sit between the inner anchor and the
outer P5 tail (29bp) or rc(P7) tail (24bp). For each end whose inner
anchor was missing, the outer tail was searched (edits ≤ length/5) and
the copies read from that side:

| | i5 | i7 |
|---|---|---|
| inner anchor missing | 181,423 (21.2%) | 211,194 (24.7%) |
| … no outer tail either | 94.1% | 93.3% |
| … outer tail, both copies agree on a declared index | 3.8% (0.80% of ends) | 4.7% (1.17% of ends) |
| … outer tail, one copy matches | 1.4% | 1.2% |
| … outer tail, no index | 0.6% | 0.8% |

Validated on ends the inner anchor resolves: where the outer method read
an index it was the same one 99.7% of the time ("both copies agree,
different index": 30 ends at each end).

Caveat: even resolved ends lack the outer tail 28% (i5) and 37% (i7) of
the time, so these reads often do not keep it (basecaller trimming or
concatemer splitting, not checked). Its absence alone does not prove an
end is missing. But missing-anchor ends lack it 94% of the time, against
28–37% on good ends.

**Conclusion.** The anchor bucket is mostly ends that are not in the read:
truncated molecules, or PCR2 never adding the index. Neither a looser
anchor nor anchorless scanning recovers more than about 1–4% of it. The
outer-tail route is accurate and could add about one point of tier A.
Telling "never indexed" from "read stops short" is a library question: a
TapeStation or gel trace of the final library.

## 6. Why this unlocks anchorless scanning

The demultiplexer anchors first because the current sets cannot be
scanned for — measured:

| | min i5-vs-i7 distance | random window matching some index |
|---|---|---|
| current 8bp | **1** — indistinguishable | 4.115% at edit ≤1 |
| new 18bp | **6** — unambiguous | **0.005%** at edit ≤3 |

At 18bp the anchor stops being a prerequisite and becomes an
optimisation: use it to place the window cheaply when present, **fall
back to scanning when it is not**. That is the only route to the 38.57%
of i7 ends currently never attempted.

Caveat: 0.005% is against *random* sequence. Real ITS2 is GC-structured
and the true rate could be higher — testable by scanning real tier-D
reads for composites and counting hits where no index should be.

## 7. Implementation, in order

**1. Fix `--anchor-max-edits` (blocking). Done 2026-11-04:**
`--anchor-max-edits` now reaches the end anchors and `--junction-max-edits`
governs the tier E scan; both default to 2, and output at the defaults is
byte-identical (VONT01, 683k reads). `assign_read` takes it and
calls `call_end` twice without passing it, so anchor finding is pinned
at 2 and the flag only reaches `has_internal_adapter`. Raising it today
destroys tier A via the tier E check. Split into two parameters —
anchor tolerance and junction tolerance — since their cost asymmetries
are opposite: a missed anchor loses a read, a false junction hit
destroys a good one.

**2. Measure the anchor split. Done 2026-11-04, see §5.1: the anchor is
mostly absent, not too strictly matched.** Sweep anchor tolerance 2/3/4/6 on
existing reads. If `anchor_not_found` collapses, the anchor is present
and too strictly matched; if it barely moves, PCR2 never added the index
and no index design helps. **This decides whether rows 2–4 of the
projection are reachable, and it should happen before any oligo order.**

**3. Index map schema.** `assets/index_map.tsv` needs both
representations — the 18bp composite for nanopore and the 8bp X for
Illumina. The samplesheet already accepts index *names*, so carrying
`S502` and resolving per platform avoids a second samplesheet.

**4. `barcode_repeats = 1`** in the run manifest. Existing code handles
an 18bp single index without modification — `span = length × repeats`
and the window sits adjacent to the anchor, verified on synthetic reads.

**5. Pre-flight's `repeats=1` warning needs to key on geometry, not on
the repeat count.** It currently warns that there is no second copy to
cross-check and that rotational collisions become live misassignment
routes. For these sets that is misleading: the 18bp word has **zero**
rotational collisions and is stronger than today's two 8bp copies on
every axis.

**6. Anchorless scanning as a fallback** — the large one, and only worth
building if step 2 says the anchor bucket contains readable indexes.

## 8. Open questions

- **Which copy does Illumina read?** Determines the composite order
  (§3). Inferred, not confirmed. Load-bearing.
- **How much of `anchor_not_found` carries a readable index?** Answered
  in §5.1: 1–4% with today's indexes. Whether the 18bp composite reads
  more of the damaged windows is still untested.
- **Which Illumina index read is primed from which side?** As understood
  (unverified against Illumina documentation): the i7 index read starts
  from the insert side, so X would need to sit next to the i7 anchor, and
  the i5 direction differs between forward-strand (MiSeq) and
  reverse-complement (NextSeq, NovaSeq v1.5/X) instruments.
- **Does `barcode_off_read` rise enough to matter at 18bp?** Needs the
  read-length distribution.
- **Does `cap=4` recover the k=3 indel row?** One parameter, untested.
- **Is seq-Lev 8 optimal?** Local search converged there from 14
  restarts with many pairs at the minimum, which suggests a structural
  ceiling given the 96-member UDP pool. Not proven; an exact solver
  would settle it.

## 9. Rejected, and why

**384 unique indexes** (one per sample, 14bp at seq-Lev 5). Would make
tier B assignable and take usable reads from 26% to ~78% — by far the
largest available gain. Needs ~385 primers against the current 40. Not
affordable.

**Putting the whole plate address in one end** as 24 X × 16 Y = 384
composites. Collapses: two wells sharing an X differ only in Y, so a
product structure inherits the worse half. Measured min seq-Lev **2** —
no better than today.

```
TACGCTGCAACACGTGGA
TAGCGCTCAACACGTGGA
```

Uniqueness and distance are not obtainable from the same concatenation.

**Dropping the tandem duplication without extending** — ρ = 0.006/0.008
measured, so the copies fail independently and the second one is doing
real work. It is only worth replacing with something stronger.
