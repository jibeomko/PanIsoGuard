#!/usr/bin/env python3
"""Envelope for the held-out chromosomes: the chr22 head-to-head scoring (score_h2h / emit_metrics of
../sqanti3_filter_h2h, imported unchanged) applied to each chromosome's run, settings none and sr_1000000,
with the same bootstrap CIs. Thresholds were never touched after chr22.

Usage: emit_heldout.py CHROM=ROOT [CHROM=ROOT ...] --emit-metrics out.json     (ROOT has work/ and h2h/)
"""
import datetime, json, os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels, load, metrics
from emit_metrics import risk_coverage, methods_for

SETTINGS = ["none", "sr_1000000"]


def one(root):
    work, h2h = f"{root}/work", f"{root}/h2h"
    base_cls = f"{work}/sqanti_out/flair_classification.txt"
    label, cat = build_labels(f"{work}/flair.isoforms.gtf", f"{work}/truth.truth.tsv", base_cls)
    ids = sorted(t for t in label if cat[t] != "full-splice_match")
    rng = random.Random(1); resamples = [[ids[rng.randrange(len(ids))] for _ in ids] for _ in range(2000)]
    ci = lambda a, b: (lambda d: [round(d[int(0.025 * len(d))], 3), round(d[int(0.975 * len(d)) - 1], 3)])(sorted(metrics(a, r, label)["auprc"] - metrics(b, r, label)["auprc"] for r in resamples))
    out = {"n_non_fsm": len(ids), "n_genuine": sum(label[t] == "genuine_novel" for t in ids), "n_false": sum(label[t] == "false_novel" for t in ids),
           "bootstrap_auprc_diff_vs_allSJ3_95ci": {}, "bootstrap_auprc_diff_vs_SQ3rules_95ci_no_short_reads": {}}
    for s in SETTINGS:
        cls = base_cls if s == "none" else f"{h2h}/{s}/sq/flair_classification.txt"
        specs = methods_for(h2h, s, cls); preds = {n: load(k, p) for n, (k, p) in specs.items()}
        out[s] = {}
        for n, pred in preds.items():
            m = metrics(pred, ids, label)
            out[s][n] = {"TP": m["tp"], "FP": m["fp"], **{k: (None if m[k] != m[k] else round(m[k], 3)) for k in ("prec", "rec", "spec", "f1", "auprc")},
                         **risk_coverage(specs[n][0], specs[n][1], pred, ids, label)}
        ref = "SQ3rules" if s == "none" else "allSJ3"; cis = {n: ci(preds[n], preds[ref]) for n in ("PIG", "PIG+bam")}
        if s == "none": out["bootstrap_auprc_diff_vs_SQ3rules_95ci_no_short_reads"] = cis
        else: out["bootstrap_auprc_diff_vs_allSJ3_95ci"][s] = cis
        print(os.path.basename(root), s, {n: out[s][n]["auprc"] for n in preds}, cis, file=sys.stderr)
    prov = {}
    for line in open(f"{h2h}/none/pig.provenance.log"):
        if "\t" in line and not line.startswith("#"):
            k, v = line.rstrip("\n").split("\t", 1); prov.setdefault(k, v)
    return out, prov


def main():
    specs = [a.split("=", 1) for a in sys.argv[1:] if "=" in a and not a.startswith("--")]
    out = sys.argv[sys.argv.index("--emit-metrics") + 1]
    res, prov = {}, {}
    for chrom, root in specs:
        res[chrom], prov = one(root)
    env = {"protocol": "sqanti_sim_heldout", "scope": "gencode_v49_heldout_chromosomes_sqanti_sim", "status": "tracked",
           "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
           "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
           "data_provenance": "SQANTI-SIM on GENCODE v49, one chromosome at a time (" + ", ".join(c for c, _ in specs) + "), the chr22 protocol of "
                              "benchmark/sqanti_sim + benchmark/sqanti3_filter_h2h repeated unchanged (design equal -nt 1500, NIC 300, NNC 300, ISM 200, "
                              "100,000 PBSIM3 reads, seed 1; FLAIR; SQANTI3 6.0.1; short reads none or 1 M wgsim pairs -> STAR 2-pass). Public/simulated.",
           "command": "bash benchmark/sqanti_sim_heldout/run.sh <chrom>; benchmark/sqanti_sim_heldout/emit_heldout.py CHROM=ROOT ... --emit-metrics ...",
           "metrics": res,
           "notes": "Held-out test of the thresholds set on chr22 (chain_min_full_reads 3, mapping gates, support axis): nothing was re-tuned. "
                    "Per chromosome and short-read setting, non-FSM genuine-vs-false metrics per method as in sqanti3_filter_h2h, with the same bootstrap CIs."}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(env, open(out, "w"), indent=2, sort_keys=True); open(out, "a").write("\n"); print("wrote", out, file=sys.stderr)


if __name__ == "__main__":
    main()
