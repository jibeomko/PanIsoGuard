# A worked locus: one toy gene from input files to verdicts

`data/` holds a 1,700 bp toy chromosome `chrT` with one gene, GENE_T (+ strand): three reference
transcripts (TX1 uses all five exons, TX2 skips exon 2, TX3 skips exon 4) and nine isoforms from a
caller (`iso_known` = TX1, `iso_ism` = a truncated TX1, `iso_A`–`iso_G` novel), plus the inputs a run
can take: a SQANTI3 classification, a second caller, STAR short-read junctions, a long-read BAM and
two haplotype FASTAs. `data/build_toy.sh` shows how every file was made (`make_toy.py`, SQANTI3
6.0.1, STAR); the outputs are committed, so nothing here needs those tools.

Adding the inputs one at a time changes the verdict of `iso_A`, the isoform with one novel junction:

| input added | `iso_A` |
|---|---|
| SQANTI3 classification + caller GTF | `UNRESOLVED` |
| + reference GTF | `UNRESOLVED` (one novel junction, nothing to test it with) |
| + short-read `SJ.out.tab` | `ARTIFACT` (unsupported, non-canonical) |
| + long-read BAM | `ARTIFACT` (no mapping sign either) |
| + haplotype FASTAs, provenance `unknown` | `UNRESOLVED` (the rescue is held by the circularity firewall) |
| + provenance `wgs` | `REFERENCE_BIAS` |

Two scripts keep this folder honest, and both run in CI:

- `check_worked_locus.py` recomputes every verdict of the nine isoforms under nine input
  combinations with the Python standard library only (readers, junction classification, short-read
  support, BAM features, splice motifs, fingerprints, the rules) and compares them with the binary
  (CTest `integration_worked_locus`).
- `rules_replay.py` is the rulebook of `RuleEngine::evaluate()` in 42 lines of Python, replayed on the
  six runs above against the binary's `attribution.jsonl` (CTest `integration_rules_replay`).

```
panisoguard example worked_locus      # copies this folder next to you
python3 worked_locus/check_worked_locus.py $(command -v panisoguard)
python3 worked_locus/rules_replay.py $(command -v panisoguard)
```
