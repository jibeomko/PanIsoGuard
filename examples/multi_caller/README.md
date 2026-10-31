# Multi-caller quickstart (self-contained, offline, <1 s)

The shortest end-to-end run of PanIsoGuard's headline flow — **integrate several callers,
then let cross-caller agreement decide which novel isoforms to trust** — on tiny committed
fixtures (no downloads, no real callers, deterministic).

```bash
./run.sh            # combine -> adjudicate --caller-support
```

## What it shows

Three toy callers report novel isoforms that overlap to different degrees (chains keyed by
splice structure, not by name):

| isoform | recovered by | verdict |
|---------|--------------|---------|
| `isoA` | flair + isoquant + bambu (**3 callers**) | **NOVEL_PROBABLE** |
| `isoB` | flair + isoquant (**2 callers**) | **NOVEL_PROBABLE** |
| `isoC` | flair only (**1 caller**) | **UNRESOLVED** |

With no short-read evidence, a single-caller novel is held `UNRESOLVED`, while a novel chain
recovered by two or more callers is `NOVEL_PROBABLE`. The consensus axis never reaches
`NOVEL_CONFIRMED`: caller agreement is methodological, not experimental, corroboration.
`run.sh` diffs its output against [expected/](expected/).

## The pipeline (same three commands on real data)

```
combine  --gtf flair:… --gtf isoquant:… --gtf bambu:… --ref-gtf … --out matrix.tsv
adjudicate --classification <sqanti> --isoforms-gtf <caller> --ref-gtf … --caller-support matrix.tsv --out-prefix sample
```

For the full real-data, truth-scored version (5 callers, gffcompare/TAMA head-to-head) see
[../../benchmark/multicaller](../../benchmark/multicaller) and
[../../benchmark/merge_comparison](../../benchmark/merge_comparison). For a single-caller
example see [../tiny](../tiny).
