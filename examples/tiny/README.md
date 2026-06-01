# Tiny PanIsoGuard example

This is a minimal text-only dataset for a 5-minute smoke run. It demonstrates
three expected outcomes without requiring aligners, samtools, or large inputs:

| isoform | evidence pattern | expected class |
|---------|------------------|----------------|
| `iso_known` | matches the reference transcript | `HIGH_CONF_KNOWN` |
| `iso_novel_supported` | one novel junction supported by `SJ.tab` | `HIGH_CONF_NOVEL` |
| `iso_novel_artifact` | unsupported novel junction with SQANTI3 non-canonical/RTS/degradation priors | `ARTIFACT` |

## Run

From the repository root after building PanIsoGuard:

```bash
cmake --build build -j
cd examples/tiny
./run.sh
```

Or point to an installed binary:

```bash
PANISOGUARD=panisoguard ./run.sh
```

The script writes `output/sample.{adjudicated.tsv,attribution.jsonl,provenance.log}`
and compares them with the checked-in files under `expected/`.

## Files

| Path | Purpose |
|------|---------|
| `data/classification.tsv` | SQANTI3-like classification table consumed as QC priors |
| `data/caller.gtf` | caller isoform structures; transcript IDs match `classification.tsv` |
| `data/reference.gtf` | minimal reference catalog with one known intron |
| `data/short_read.SJ.tab` | STAR-like junction support for the supported novel isoform |
| `expected/` | deterministic expected output from `panisoguard adjudicate` |
