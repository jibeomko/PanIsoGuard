#!/usr/bin/env python3
"""Envelope for the ONT re-simulation of the chr22 truth set: the chr22 head-to-head scoring
(score_h2h / emit_metrics of ../sqanti3_filter_h2h, imported unchanged) on the NanoSim run, settings
none and sr_1000000, with a third PanIsoGuard method (BAM, --platform ont) next to the default preset,
the number of isoforms whose class differs between the two presets, and the share of gated novel
junctions whose raw indel-near fraction, and whose excess over the reads' own indel rate, pass the gate.

Usage: emit_ont.py ROOT --emit-metrics out.json     (ROOT has work/ and h2h/)
"""
import datetime, json, os, random, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels, load, metrics, table
from emit_metrics import risk_coverage, methods_for

SETTINGS = ["none", "sr_1000000"]
PIGS = ("PIG", "PIG+bam", "PIG+bam(ont)")


def indel_gate_rates(junctions_tsv, gate=0.5, min_spanning=3):
    """Over novel junctions with >= min_spanning reads: share whose raw indel-near fraction exceeds the gate,
    and share whose excess over the background does."""
    rows = [r for r in table(junctions_tsv) if r.get("known", "") in ("false", "0", "no", "novel") or r.get("junction_known", "") in ("false", "0")]
    if not rows:   # column naming differs between versions: fall back to every junction with a background value
        rows = table(junctions_tsv)
    gated = [r for r in rows if r.get("bam_n_spanning", "NA") not in ("NA", "") and int(r["bam_n_spanning"]) >= min_spanning]
    raw = sum(1 for r in gated if float(r["bam_frac_indel_near"]) > gate)
    excess = sum(1 for r in gated if float(r["bam_frac_indel_near"]) - float(r["bam_indel_near_background"]) > gate)
    bg = [float(r["bam_indel_near_background"]) for r in gated]
    return {"gated_junctions": len(gated), "raw_over_gate": raw, "excess_over_gate": excess,
            "mean_background": round(sum(bg) / len(bg), 3) if bg else None}


def one(root):
    work, h2h = f"{root}/work", f"{root}/h2h"
    base_cls = f"{work}/sqanti_out/flair_classification.txt"
    label, cat = build_labels(f"{work}/flair.isoforms.gtf", f"{work}/truth.truth.tsv", base_cls)
    ids = sorted(t for t in label if cat[t] != "full-splice_match")
    rng = random.Random(1); resamples = [[ids[rng.randrange(len(ids))] for _ in ids] for _ in range(2000)]
    ci = lambda a, b: (lambda d: [round(d[int(0.025 * len(d))], 3), round(d[int(0.975 * len(d)) - 1], 3)])(sorted(metrics(a, r, label)["auprc"] - metrics(b, r, label)["auprc"] for r in resamples))
    out = {"n_non_fsm": len(ids), "n_genuine": sum(label[t] == "genuine_novel" for t in ids), "n_false": sum(label[t] == "false_novel" for t in ids),
           "n_flair_isoforms": len(cat),
           "bootstrap_auprc_diff_vs_allSJ3_95ci": {}, "bootstrap_auprc_diff_vs_SQ3rules_95ci_no_short_reads": {},
           "hifi_vs_ont_preset": {}, "indel_near_gate": {}}
    for s in SETTINGS:
        cls = base_cls if s == "none" else f"{h2h}/{s}/sq/flair_classification.txt"
        specs = methods_for(h2h, s, cls); specs["PIG+bam(ont)"] = ("pig", f"{h2h}/{s}/pigbam_ont.adjudicated.tsv")
        preds = {n: load(k, p) for n, (k, p) in specs.items()}
        out[s] = {}
        for n, pred in preds.items():
            m = metrics(pred, ids, label)
            out[s][n] = {"TP": m["tp"], "FP": m["fp"], **{k: (None if m[k] != m[k] else round(m[k], 3)) for k in ("prec", "rec", "spec", "f1", "auprc")},
                         **risk_coverage(specs[n][0], specs[n][1], pred, ids, label)}
        ref = "SQ3rules" if s == "none" else "allSJ3"; cis = {n: ci(preds[n], preds[ref]) for n in PIGS}
        if s == "none": out["bootstrap_auprc_diff_vs_SQ3rules_95ci_no_short_reads"] = cis
        else: out["bootstrap_auprc_diff_vs_allSJ3_95ci"][s] = cis
        hifi = {r["isoform_id"]: r["confidence_class"] for r in table(f"{h2h}/{s}/pigbam.adjudicated.tsv")}
        ont = {r["isoform_id"]: r["confidence_class"] for r in table(f"{h2h}/{s}/pigbam_ont.adjudicated.tsv")}
        out["hifi_vs_ont_preset"][s] = {"isoforms_with_a_different_class": sum(1 for t in hifi if hifi[t] != ont.get(t))}
        out["indel_near_gate"][s] = indel_gate_rates(f"{h2h}/{s}/pigbam.junctions.tsv")
        print(s, {n: out[s][n]["auprc"] for n in preds}, cis, out["hifi_vs_ont_preset"][s], out["indel_near_gate"][s], file=sys.stderr)
    prov = {}
    for line in open(f"{h2h}/none/pig.provenance.log"):
        if "\t" in line and not line.startswith("#"):
            k, v = line.rstrip("\n").split("\t", 1); prov.setdefault(k, v)
    return out, prov


def main():
    root = [a for a in sys.argv[1:] if not a.startswith("--") and a != sys.argv[sys.argv.index("--emit-metrics") + 1]][0]
    out = sys.argv[sys.argv.index("--emit-metrics") + 1]
    res, prov = one(root)
    env = {"protocol": "sqanti_sim_ont", "scope": "gencode_v49_chr22_ont_nanosim_sqanti_sim", "status": "tracked",
           "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
           "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
           "data_provenance": "SQANTI-SIM on GENCODE v49 chr22, the design of benchmark/sqanti_sim reused unchanged (same reduced annotation, "
                              "hidden transcripts and requested counts); reads re-simulated as ONT cDNA with NanoSim (SQANTI-SIM's pre-trained "
                              "human WTC11 cDNA guppy model, 100,000 reads, seed 1); minimap2 -ax splice; FLAIR; SQANTI3 6.0.1; short reads none or "
                              "1 M wgsim pairs -> STAR 2-pass; SQANTI3 rules / ML / strict filters as in benchmark/sqanti3_filter_h2h. Public/simulated.",
           "command": "bash benchmark/sqanti_sim_ont/run.sh; benchmark/sqanti_sim_ont/emit_ont.py ROOT --emit-metrics ...",
           "metrics": res,
           "notes": "ONT counterpart of the chr22 truth set with the thresholds set on PacBio HiFi left untouched. Per short-read setting, "
                    "non-FSM genuine-vs-false metrics per method as in sqanti3_filter_h2h (same bootstrap CIs), PanIsoGuard with the BAM under the "
                    "default (hifi) preset and under --platform ont, the number of isoforms whose class differs between the presets, and, over "
                    "the gated novel junctions, how many the raw indel-near fraction would flag against how many the excess over the reads' own "
                    "indel rate flags (ruleset 0.0.7)."}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(env, open(out, "w"), indent=2, sort_keys=True); open(out, "a").write("\n"); print("wrote", out, file=sys.stderr)


if __name__ == "__main__":
    main()
