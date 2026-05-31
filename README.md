# PanIsoGuard

**Caller-agnostic adjudication of long-read RNA-seq novel isoforms.**

PanIsoGuard is a post-processing / decision layer that ingests the *novel* isoform
calls produced by long-read isoform callers (FLAIR, IsoQuant, Bambu, ESPRESSO,
TALON, …) and/or [SQANTI3](https://github.com/ConesaLab/SQANTI3) output, together
with a BAM and **optional** evidence inputs (short-read `SJ.tab`, personalized
haplotype FASTA, pangenome `GFA` / `rpvg` output), and re-classifies each novel
call into a confidence class with a **machine-readable mechanistic attribution**
and a **provenance / circularity flag**.

It does not recompute SQANTI3's QC features and does not claim to beat the SQANTI3
random-forest filter or FLAIR2. Its contribution is the **adjudication logic** —
how orthogonal evidence axes are integrated into a transparent, auditable verdict
(every threshold and the full label-projection table live in a runtime
[`config/rules.default.toml`](config/rules.default.toml), and each verdict carries
a `rule_trace`).

> **Status: pre-alpha (M0 scaffold).** Only the build skeleton and a CLI stub
> exist. APIs, formats, and rule thresholds are unstable and uncalibrated.

## Confidence classes

`HIGH_CONF_KNOWN` · `HIGH_CONF_NOVEL` · `MEDIUM_CONF_NOVEL` · `LOW_CONF_PARTIAL` ·
`PAN_REF_RESCUED_FALSE_NOVEL` · `AMBIGUOUS` · `ARTIFACT` — emitted as a
deterministic projection of a 2-axis evidence grid (novelty-support ×
artifact-mechanism). (A graph/haplotype-path rescue class is reserved for the
future pangenome tier and is not emitted yet.)

## Evidence tiers

| Tier | Input | Required? | Mechanism |
|------|-------|-----------|-----------|
| 0 | SQANTI3 classification (priors) + STAR `SJ.tab` | recommended | short-read junction corroboration |
| 1 | BAM (HiFi/ONT) | recommended | read-level degradation / mapping / chimera features |
| 2 | personalized haplotype FASTA (`--reference-haplotype`) | optional | variant-created/destroyed splice-site motif |
| 3 | pangenome `GFA` / `rpvg` text | optional | haplotype-path rescue (subprocess/file only in v1.0) |

Tiers 0–2 build and run from a plain `cmake .. && make` against htslib alone.

## Build

Requires a C++17 compiler, CMake ≥ 3.20, and **htslib ≥ 1.18**
(`sam_itr_regarray`, `hts_set_thread_pool`).

```bash
git clone <repo> && cd PanIsoGuard
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build            # version smoke test
./build/panisoguard --version
```

htslib is discovered from `$CONDA_PREFIX` automatically; override with
`-DCMAKE_PREFIX_PATH=/path/to/prefix`.

The experimental in-process pangenome (GBZ) backend is **off by default**; enable
with `-DWITH_PANGENOME_LIB=ON` (and `--experimental-gbz` at runtime). It is
unsupported in v1.0.

## License

MIT — see [LICENSE](LICENSE).
