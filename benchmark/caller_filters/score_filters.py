#!/usr/bin/env python3
"""Post-caller filtering options on the same caller output: nothing, SQANTI3's rules filter, or
PanIsoGuard (+ BAM, with and without 1 M short-read pairs), for each of the five chr22 SQANTI-SIM
callers; and, where a caller has its own filter knob, that knob tightened (Bambu NDR 0.5 -> 0.3 /
0.2 / 0.1; TALON before and after talon_filter_transcripts) against PanIsoGuard on the permissive
output. Labels as in sqanti3_filter_h2h (exact intron-chain match to the simulated truth, non-FSM
multi-exon isoforms).

Usage: score_filters.py --work W --pig-dir CHAIN_SUPPORT_OUT
         --caller NAME=GTF,CLS,RULES_CLS ...            (the tracked output of each caller)
         --family NAME=LABEL=GTF,CLS[,RULES_CLS] ...    (alternative settings of one caller; the
                                                         first LABEL listed is the permissive base)
         [--emit-metrics ../results/caller_filters/metrics.json]
"""
import argparse, datetime, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels, load, table  # noqa: E402

POSITIVE = {"NOVEL_CONFIRMED", "NOVEL_PROBABLE", "REFERENCE_BIAS"}
DECISIVE = POSITIVE | {"ARTIFACT"}


def chain_keys(gtf):
    ex, strand, chrom = {}, {}, {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        t = f[8].split('transcript_id "')[1].split('"')[0]
        ex.setdefault(t, []).append((int(f[3]), int(f[4]))); strand[t] = f[6]; chrom[t] = f[0]
    out = {}
    for t, e in ex.items():
        e = sorted(e)
        if len(e) < 2: continue
        out[t] = f"{chrom[t]}|{strand[t]}|" + ",".join(f"{e[i][1] + 1}-{e[i + 1][0] - 1}" for i in range(len(e) - 1))
    return out


def pr(tp, fp, denom):
    return {"TP": tp, "FP": fp, "prec": round(tp / (tp + fp), 3) if tp + fp else None, "rec": round(tp / denom, 3) if denom else None}


def pig_point(path, ids, label, denom):
    cls = {r["isoform_id"]: r["confidence_class"] for r in table(path)}
    pos = [t for t in ids if cls.get(t) in POSITIVE]
    dec = [t for t in ids if cls.get(t) in DECISIVE]
    tp = sum(label[t] == "genuine_novel" for t in pos); fp = len(pos) - tp
    wrong = sum((cls[t] in POSITIVE) != (label[t] == "genuine_novel") for t in dec)
    out = pr(tp, fp, denom)
    out.update(coverage=round(len(dec) / len(ids), 3) if ids else None, selective_risk=round(wrong / len(dec), 3) if dec else None)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--pig-dir", required=True)
    ap.add_argument("--caller", action="append", default=[]); ap.add_argument("--family", action="append", default=[])
    ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    truth = f"{a.work}/truth.truth.tsv"
    res = {"post_caller_filters": {}, "caller_own_knob": {}}

    print("== three post-caller filters on the same output (recall = share of the caller's own genuine novel isoforms kept)")
    print(f"{'caller':9s} {'n':>5s} {'genuine':>7s} | {'none P/R':>14s} | {'SQ3 rules P/R':>14s} | {'PIG+BAM P/R  cov/risk':>26s} | {'PIG+BAM+1M P/R  cov/risk':>26s}")
    for spec in a.caller:
        name, rest = spec.split("=", 1); gtf, cls, rules = rest.split(",")
        label, cat = build_labels(gtf, truth, cls)
        ids = sorted(t for t in label if cat[t] != "full-splice_match")
        g = sum(label[t] == "genuine_novel" for t in ids)
        row = {"n_non_fsm": len(ids), "n_genuine": g, "none": pr(g, len(ids) - g, g)}
        rp = load("rules", rules)
        kept = [t for t in ids if rp.get(t, (False,))[0]]
        tpr = sum(label[t] == "genuine_novel" for t in kept)
        row["sqanti3_rules"] = pr(tpr, len(kept) - tpr, g)
        row["panisoguard_bam"] = pig_point(f"{a.pig_dir}/{name}.nosr.on.adjudicated.tsv", ids, label, g)
        row["panisoguard_bam_sr1M"] = pig_point(f"{a.pig_dir}/{name}.sr1M.on.adjudicated.tsv", ids, label, g)
        res["post_caller_filters"][name] = row
        f = lambda d: f"{d['prec']} / {d['rec']}"
        f2 = lambda d: f"{d['prec']} / {d['rec']}  {d['coverage']} / {d['selective_risk']}"
        print(f"{name:9s} {len(ids):5d} {g:7d} | {f(row['none']):>14s} | {f(row['sqanti3_rules']):>14s} | {f2(row['panisoguard_bam']):>26s} | {f2(row['panisoguard_bam_sr1M']):>26s}")

    fams = {}
    for spec in a.family:
        fam, rest = spec.split("=", 1); lab, files = rest.split("=", 1)
        fams.setdefault(fam, []).append((lab, files.split(",")))
    for fam, members in fams.items():
        base_lab, base_files = members[0]
        blabel, bcat = build_labels(base_files[0], truth, base_files[1])
        bkeys = chain_keys(base_files[0])
        universe = {bkeys[t] for t in blabel if bcat[t] != "full-splice_match" and blabel[t] == "genuine_novel" and t in bkeys}
        print(f"\n== {fam}: its own knob vs PanIsoGuard on the permissive output (recall = share of the {len(universe)} genuine novel chains of '{base_lab}')")
        out = {"base": base_lab, "n_genuine_chains_base": len(universe), "settings": {}}
        for lab, files in members:
            label, cat = build_labels(files[0], truth, files[1]); keys = chain_keys(files[0])
            ids = [t for t in label if cat[t] != "full-splice_match"]
            tp_keys = {keys[t] for t in ids if label[t] == "genuine_novel" and t in keys} & universe
            fp = sum(label[t] == "false_novel" for t in ids)
            entry = {"n_non_fsm": len(ids), "as_is": pr(len(tp_keys), fp, len(universe))}
            if len(files) > 2:
                rp = load("rules", files[2]); kept = [t for t in ids if rp.get(t, (False,))[0]]
                entry["sqanti3_rules"] = pr(len({keys[t] for t in kept if label[t] == "genuine_novel" and t in keys} & universe),
                                            sum(label[t] == "false_novel" for t in kept), len(universe))
            for setting in ("nosr", "sr1M"):
                p = f"{a.pig_dir}/{lab}.{setting}.on.adjudicated.tsv"
                if os.path.exists(p):
                    cls = {r["isoform_id"]: r["confidence_class"] for r in table(p)}
                    pos = [t for t in ids if cls.get(t) in POSITIVE]
                    entry[f"panisoguard_bam{'_sr1M' if setting == 'sr1M' else ''}"] = pr(
                        len({keys[t] for t in pos if label[t] == "genuine_novel" and t in keys} & universe),
                        sum(label[t] == "false_novel" for t in pos), len(universe))
            out["settings"][lab] = entry
            cells = "  ".join(f"{k}: {v['prec']} / {v['rec']}" for k, v in entry.items() if isinstance(v, dict))
            print(f"  {lab:<14s} n={len(ids):5d}  {cells}")
        res["caller_own_knob"][fam] = out

    if a.emit_metrics:
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{a.pig_dir}/flair.nosr.on.provenance.log") if "\t" in l and not l.startswith("#"))
        env = {"protocol": "caller_filters", "scope": "gencode_v49_chr22_5caller_sqanti_sim", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "chr22 SQANTI-SIM truth set (benchmark/sqanti_sim) called by FLAIR, IsoQuant, Bambu (NDR 0.5), "
                                  "ESPRESSO and TALON (filtered) on one alignment (benchmark/multicaller); SQANTI3 6.0.1 per caller; "
                                  "PanIsoGuard runs of benchmark/chain_support (ruleset 0.0.5); Bambu NDR 0.3/0.2/0.1 and TALON "
                                  "before talon_filter_transcripts re-run for the caller-knob comparison. Public/simulated.",
               "command": "benchmark/caller_filters/run.sh; score_filters.py --work ... --pig-dir ... --caller ... --family ... --emit-metrics",
               "metrics": res,
               "notes": "Three post-caller filters on the same output: none, SQANTI3 rules (no short-read coverage column, so its "
                        "perc_A + RTS rules), PanIsoGuard + BAM without and with 1 M short-read pairs. Recall is relative to the "
                        "caller's own genuine novel isoforms (the caller's recall is outside the comparison). The caller-knob "
                        "section tightens Bambu's NDR and applies TALON's own filter, with recall relative to the permissive "
                        "setting's genuine chains, against PanIsoGuard applied to the permissive output."}
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
