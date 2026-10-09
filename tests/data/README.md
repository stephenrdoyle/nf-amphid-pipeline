# Fixture data for `-profile test`

Generated, not curated — rebuild at any time with:

```bash
python3 tests/make_fixtures.py --outdir tests/data --platform illumina \
    --n-samples 8 --leave-unused 2
python3 bin/simulate.py --reference tests/data/reference.tsv \
    --samplesheet tests/data/samplesheet.tsv --runs tests/data/runs.tsv \
    --assay assets/assay/ampseq_helminth.json --outdir tests/data/sim \
    --reads-per-sample 300
```

Then point `runs.tsv`'s `path` column at `tests/data/sim/RUN01`.

Eight wells: one extraction blank, one NTC, one mock, one positive control
and four field samples, with two index combinations deliberately left
unused as hopping sentinels. `sim/RUN01.truth.tsv` records what every read
really was and `sim/RUN01.composition.tsv` the input copies per sample, so
the fixture can be used to check accuracy and not only that the code runs.
