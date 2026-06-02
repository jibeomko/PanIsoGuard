#!/usr/bin/env python3
"""SQANTI-SIM threshold sweep harness (the release-blocking AUPRC sweep).

The expensive half of the SQANTI-SIM protocol is producing the truth set
(SQANTI-SIM -> PBSIM3 -> minimap2 -> FLAIR -> SQANTI3); see run.sh. Once those
artifacts exist, the *sweep* is cheap: re-run `panisoguard adjudicate` under a grid
of config/rules.default.toml threshold settings and score each against the same
truth with score.py. This harness automates that and emits:

  * a sweep.tsv row per config (auprc, precision, specificity, nnc_recall),
  * the operating-point row flagged (the shipped default-0.0.1 config),
  * optionally, a benchmark/results/ tracked metrics.json from a chosen row.

Note on sweepable knobs: this protocol's truth SJ.tab is truth-derived with a fixed
n_uniq=50 and canonical motif (prep_truth_sj.py), so the short-read stringency knobs
(sj_min_uniq_reads below ~50, sj_require_canonical_motif) are near-inert here and the
grid includes a value above 50 to expose the support cliff; the knob that actually
moves verdicts is perc_A_degradation_threshold (degradation mechanism -> HIGH vs MED).

Usage:
  sweep.py PIG --work DIR [--out sweep.tsv] [--emit-metrics results/sqanti_sim/metrics.json]
           [--scope gencode_vNN_chr22] [--full]

DIR is the run.sh work directory; defaults follow run.sh's output names (override
with --classification/--isoforms-gtf/--ref-gtf/--sj-tab/--truth as needed).
"""
import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCORE = HERE / "score.py"

# Default config grid. (sj_min_uniq_reads, sj_require_canonical_motif, perc_A_threshold)
GRID_DEFAULT = [
    (3, True, 60.0),    # <- default-0.0.1 (operating point)
    (3, True, 50.0),
    (3, True, 70.0),
    (3, False, 60.0),
    (60, True, 60.0),   # above truth n_uniq=50 -> exposes the short-read support cliff
]
GRID_FULL = [(s, c, p)
             for s in (1, 3, 5, 10, 60)
             for c in (True, False)
             for p in (50.0, 60.0, 70.0)]

DEFAULT_CONFIG = (3, True, 60.0)


def write_toml(path: Path, sj_uniq: int, sj_canon: bool, perc_a: float) -> None:
    path.write_text(
        "[axis_novelty_support]\n"
        f"sj_min_uniq_reads = {sj_uniq}\n"
        f"sj_require_canonical_motif = {'true' if sj_canon else 'false'}\n"
        "[axis_artifact.degradation]\n"
        f"max_perc_A_downstream_TTS = {perc_a}\n"
    )


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"FAILED: {' '.join(map(str, cmd))}\n{r.stdout}\n{r.stderr}\n")
        sys.exit(1)
    return r.stdout


def parse_score(out: str) -> dict:
    m_prs = re.search(r"precision=([\d.]+)\s+recall=([\d.]+)\s+specificity\(false rejected\)=([\d.]+)", out)
    m_auprc = re.search(r"AUPRC = ([\d.]+)\s+\(baseline = ([\d.]+)\)", out)
    m_nnc = re.search(r"novel_not_in_catalog\s+\d+/\d+\s+recall=([\d.]+)", out)
    if not (m_prs and m_auprc):
        sys.stderr.write("could not parse score.py output:\n" + out + "\n")
        sys.exit(1)
    return {
        "precision_genuine": float(m_prs.group(1)),
        "recall_genuine": float(m_prs.group(2)),
        "specificity_false": float(m_prs.group(3)),
        "auprc": float(m_auprc.group(1)),
        "auprc_baseline": float(m_auprc.group(2)),
        "nnc_recall": float(m_nnc.group(1)) if m_nnc else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("pig")
    ap.add_argument("--work", required=True)
    ap.add_argument("--classification")
    ap.add_argument("--isoforms-gtf")
    ap.add_argument("--ref-gtf")
    ap.add_argument("--sj-tab")
    ap.add_argument("--truth")
    ap.add_argument("--out", default=str(HERE.parent / "results" / "sqanti_sim" / "sweep.tsv"))
    ap.add_argument("--emit-metrics", help="write a tracked metrics.json from the operating point")
    ap.add_argument("--scope", default="gencode_chr22")
    ap.add_argument("--full", action="store_true", help="use the full 5x2x3 grid")
    a = ap.parse_args()

    w = Path(a.work)
    cls = a.classification or str(w / "sqanti_out" / "flair_classification.txt")
    iso = a.isoforms_gtf or str(w / "flair.isoforms.gtf")
    ref = a.ref_gtf or str(w / "chr22_modified.gtf")
    sj = a.sj_tab or str(w / "truth.SJ.tab")
    truth = a.truth or str(w / "truth.truth.tsv")
    for p in (cls, iso, ref, sj, truth):
        if not Path(p).exists():
            sys.stderr.write(f"missing truth artifact: {p}\n(run benchmark/sqanti_sim/run.sh first)\n")
            return 2

    grid = GRID_FULL if a.full else GRID_DEFAULT
    rows = []
    with tempfile.TemporaryDirectory(prefix="pig_sweep_") as tmp:
        t = Path(tmp)
        for (sj_uniq, sj_canon, perc_a) in grid:
            cfg = t / f"cfg_{sj_uniq}_{int(sj_canon)}_{int(perc_a)}.toml"
            write_toml(cfg, sj_uniq, sj_canon, perc_a)
            pfx = t / cfg.stem
            run([a.pig, "adjudicate", "--config", str(cfg),
                 "--classification", cls, "--isoforms-gtf", iso,
                 "--ref-gtf", ref, "--sj-tab", sj, "--out-prefix", str(pfx)])
            sc = parse_score(run([sys.executable, str(SCORE), iso, truth, str(pfx) + ".adjudicated.tsv"]))
            is_op = (sj_uniq, sj_canon, perc_a) == DEFAULT_CONFIG
            rows.append({
                "config_id": "default-0.0.1" if is_op else f"sj{sj_uniq}_canon{int(sj_canon)}_pA{int(perc_a)}",
                "sj_min_uniq_reads": sj_uniq, "sj_require_canonical_motif": str(sj_canon).lower(),
                "perc_A_degradation_threshold": perc_a,
                "bam_max_low_mapq_frac": 0.5, "bam_max_supplementary_frac": 0.5,
                "operating_point": "yes" if is_op else "no", **sc,
            })
            print(f"  {rows[-1]['config_id']:24s} auprc={sc['auprc']:.4f} "
                  f"prec={sc['precision_genuine']:.3f} spec={sc['specificity_false']:.3f} "
                  f"nnc_recall={sc['nnc_recall']:.3f}")

    cols = ["config_id", "sj_min_uniq_reads", "sj_require_canonical_motif",
            "perc_A_degradation_threshold", "bam_max_low_mapq_frac", "bam_max_supplementary_frac",
            "auprc", "precision_genuine", "specificity_false", "nnc_recall", "operating_point", "source"]
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("# SQANTI-SIM threshold sweep produced by benchmark/sqanti_sim/sweep.py "
                 f"(scope={a.scope}). operating_point=yes is config/rules.default.toml.\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            r["source"] = f"sweep.py ({a.scope})"
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"wrote {len(rows)} rows -> {out}")

    op = next(r for r in rows if r["operating_point"] == "yes")
    best = max(rows, key=lambda r: r["auprc"])
    print(f"operating point (default-0.0.1): auprc={op['auprc']:.4f} prec={op['precision_genuine']:.3f}")
    print(f"best auprc in grid: {best['config_id']} auprc={best['auprc']:.4f}")

    if a.emit_metrics:
        import json
        ver = subprocess.run([a.pig, "version"], capture_output=True, text=True).stdout.splitlines()
        env = {
            "command": "bash benchmark/sqanti_sim/run.sh && python3 benchmark/sqanti_sim/sweep.py "
                       "$PIG --work benchmark/sqanti_sim/work --emit-metrics benchmark/results/sqanti_sim/metrics.json",
            "data_provenance": f"SQANTI-SIM ({a.scope}); truth via run.sh (SQANTI-SIM/PBSIM3/FLAIR/SQANTI3)",
            "generated_utc": None,
            "metrics": {
                "auprc_genuine": op["auprc"], "auprc_baseline": op["auprc_baseline"],
                "precision_genuine": op["precision_genuine"], "specificity_false": op["specificity_false"],
                "recall_genuine": op["recall_genuine"], "nnc_recall": op["nnc_recall"],
            },
            "notes": "Operating point = config/rules.default.toml (default-0.0.1). Full grid in sweep.tsv.",
            "protocol": "sqanti_sim", "ruleset_version": "builtin-0.0.1", "scope": a.scope,
            "source": None, "status": "tracked",
            "tool_version": ver[0].split()[-1] if ver else None,
        }
        Path(a.emit_metrics).parent.mkdir(parents=True, exist_ok=True)
        Path(a.emit_metrics).write_text(json.dumps(env, indent=2, sort_keys=True) + "\n")
        print(f"wrote tracked metrics -> {a.emit_metrics}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
