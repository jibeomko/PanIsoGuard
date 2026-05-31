# PanIsoGuard orchestration (Snakemake)

"Raw-to-report" convenience wrapper: long-read RNA FASTQ → one shared alignment →
isoform callers **in parallel** → SQANTI3 → PanIsoGuard `combine` + `adjudicate`.

This is a thin, optional wrapper around the external callers — **not** the core
contribution. The core C++ tool (`panisoguard`) is the adjudication layer; this
just runs the upstream tools reproducibly and fast.

## Why it is faster than a naive pipeline

1. **One alignment, reused.** `minimap2` runs once per sample (rule `align`); FLAIR,
   IsoQuant, Bambu, and PanIsoGuard's mapping axis all consume that single sorted,
   indexed BAM instead of each re-aligning. (FLAIR here reuses the BAM via
   `bam2Bed12` rather than `flair align` re-running minimap2.)
2. **Callers run concurrently.** flair / isoquant / bambu are independent rules, so
   `snakemake --cores N` schedules them in parallel (wall-clock ≈ slowest caller,
   not the sum) — also across samples.
3. **Threaded tools.** Each rule passes `threads` to minimap2/caller/SQANTI3.
4. **Pick fewer callers.** Set `callers:` to the 2–3 you actually want; PanIsoGuard
   integrates whatever you run.

The remaining wall-clock is dominated by the callers themselves (minutes–hours per
sample) — external tools PanIsoGuard cannot make intrinsically faster. The wrapper
removes the *redundant* cost (re-alignment, serial execution).

## Run

```bash
# edit config.yaml: samples, reference, callers, threads, panisoguard path, sqanti3_dir
cd workflow
snakemake --use-conda --cores 16            # full run
snakemake -n                                 # dry-run: inspect the DAG first
snakemake --use-conda --cores 16 results/adjudicate/sample1.adjudicated.tsv   # one target
```

Each tool runs in its own conda env (`envs/*.yaml`) via `--use-conda`, so nothing
needs to be on the base PATH (mirrors the manuscript pipeline's pattern).

## Outputs (per sample)

- `results/align/{s}.bam(.bai)` — shared alignment
- `results/callers/{caller}/{s}.{caller}.gtf` — each enabled caller
- `results/sqanti3/{s}_classification.txt` — SQANTI3 QC on the primary caller
- `results/combine/{s}.caller_support_matrix.tsv` — multi-caller integration
- `results/adjudicate/{s}.adjudicated.tsv` (+ `.attribution.jsonl`, `.provenance.log`)

## Status / caveats

- **Validated:** the DAG shape and the `align` / `combine` / `adjudicate` rules
  (PanIsoGuard's own commands are exercised by the test suite and on real data).
- **Adapt to your versions:** the `flair` / `isoquant` / `bambu` / `sqanti3` rule
  command lines and output paths are version-sensitive templates — check them
  against your installed tools before a production run. The FLAIR rule and the
  reference paths mirror the project's existing `scripts/pipeline/*.sh`.
- The variant axis (`haplotype_fastas`) needs personalized FASTAs (e.g. from
  matched WGS via `bcftools consensus`); leave empty to skip it.
