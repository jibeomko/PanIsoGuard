#!/usr/bin/env python3
"""Write the tracked envelope for the head-to-head (../results/sqanti3_filter_h2h/metrics.json):
per short-read setting, non-FSM metrics per method, plus bootstrap 95 % CIs (2000 resamples
of the non-FSM isoforms, seed 1) of the AUPRC difference PIG / PIG+bam minus allSJ3 (and
minus SQ3rules when there are no short reads).

Usage: emit_metrics.py <sqanti_sim work> <h2h work (run.sh OUT)> --emit-metrics out.json
"""
import datetime, json, os, random, sys

from score_h2h import build_labels, load, metrics

SETTINGS = ["none", "oracle", "sr_250000", "sr_1000000", "sr_4000000"]


def methods_for(h2h, s, cls):
    m = {"SQ3rules": ("rules", f"{h2h}/{s}/rules/rules_RulesFilter_classification.txt"),
         "SQ3ml": ("ml", f"{h2h}/{s}/ml/ml_ML_classification.txt")}
    if s != "none":
        m["SQ3strict"] = ("rules", f"{h2h}/{s}/strict/strict_RulesFilter_classification.txt")
        m["allSJ3"] = ("mincov", cls)
    m["PIG"] = ("pig", f"{h2h}/{s}/pig.adjudicated.tsv")
    m["PIG+bam"] = ("pig", f"{h2h}/{s}/pigbam.adjudicated.tsv")
    return m


def main():
    work, h2h = sys.argv[1], sys.argv[2]
    out = sys.argv[sys.argv.index("--emit-metrics") + 1]
    base_cls = f"{work}/sqanti_out/flair_classification.txt"
    label, cat = build_labels(f"{work}/flair.isoforms.gtf", f"{work}/truth.truth.tsv", base_cls)
    ids = sorted(t for t in label if cat[t] != "full-splice_match")
    rng = random.Random(1)
    resamples = [[ids[rng.randrange(len(ids))] for _ in ids] for _ in range(2000)]

    def ci(pred_a, pred_b):
        d = sorted(metrics(pred_a, r, label)["auprc"] - metrics(pred_b, r, label)["auprc"] for r in resamples)
        return [round(d[int(0.025 * len(d))], 3), round(d[int(0.975 * len(d)) - 1], 3)]

    env_metrics = {"bootstrap_auprc_diff_vs_allSJ3_95ci": {}, "bootstrap_auprc_diff_vs_SQ3rules_95ci_no_short_reads": {}}
    for s in SETTINGS:
        cls = base_cls if s == "none" else f"{h2h}/{s}/sq/flair_classification.txt"
        preds = {n: load(k, p) for n, (k, p) in methods_for(h2h, s, cls).items()}
        env_metrics[s] = {}
        for n, pred in preds.items():
            m = metrics(pred, ids, label)
            env_metrics[s][n] = {"TP": m["tp"], "FP": m["fp"],
                                 **{k: (None if m[k] != m[k] else round(m[k], 3)) for k in ("prec", "rec", "spec", "f1", "auprc")}}
        ref = "SQ3rules" if s == "none" else "allSJ3"
        cis = {n: ci(preds[n], preds[ref]) for n in ("PIG", "PIG+bam")}
        if s == "none":
            env_metrics["bootstrap_auprc_diff_vs_SQ3rules_95ci_no_short_reads"] = cis
        else:
            env_metrics["bootstrap_auprc_diff_vs_allSJ3_95ci"][s] = cis
        print(s, {n: env_metrics[s][n]["auprc"] for n in preds}, cis, file=sys.stderr)

    prov = {}
    for line in open(f"{h2h}/none/pig.provenance.log"):
        if "\t" in line and not line.startswith("#"):
            k, v = line.rstrip("\n").split("\t", 1); prov.setdefault(k, v)
    env = {
        "protocol": "sqanti3_filter_h2h", "scope": "gencode_v49_chr22_sqanti_sim", "status": "tracked",
        "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
        "generated_utc": datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
        "data_provenance": "SQANTI-SIM GENCODE v49 chr22 truth set of benchmark/sqanti_sim (FLAIR + SQANTI3 6.0.1); "
                           "non-FSM multi-exon isoforms, 730 genuine / 140 false. Short reads: none, the oracle truth "
                           "SJ.tab, or wgsim 2x100 pairs (250k/1M/4M) from the simulated transcripts weighted by "
                           "long-read counts -> STAR 2.7.4a 2-pass. Public/simulated.",
        "command": "bash benchmark/sqanti3_filter_h2h/run.sh (after benchmark/sqanti_sim/run.sh); "
                   "benchmark/sqanti3_filter_h2h/emit_metrics.py <work> <h2h> --emit-metrics ...; bootstrap CIs: "
                   "2000 resamples of the non-FSM isoforms, seed 1",
        "metrics": env_metrics,
        "notes": "Per setting: non-FSM genuine-vs-false metrics per method. allSJ3 = keep an isoform iff every "
                 "junction has >= 3 unique short reads (SQANTI3 min_cov). SQ3strict = SQANTI3 rules filter with "
                 "min_cov >= 3 required. SQ3ml uses SQANTI3's default TP/TN (FSM vs NNC) which is mismatched to "
                 "this simulation (most NNC are genuine), so it is reported but not used for conclusions. "
                 "Conclusion: PanIsoGuard is not more accurate than a one-line short-read rule; PIG+bam's small "
                 "AUPRC gain is within the bootstrap CI at every depth since 0.0.6 (the BAM gate needs >= 3 "
                 "spanning reads); with no short reads it is worse than the SQANTI3 rules filter.",
    }
    json.dump(env, open(out, "w"), indent=2, sort_keys=True); open(out, "a").write("\n")
    print("wrote", out, file=sys.stderr)


if __name__ == "__main__":
    main()
