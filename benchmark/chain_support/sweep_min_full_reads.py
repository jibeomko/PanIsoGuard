#!/usr/bin/env python3
"""Why >= 3 whole-chain reads: adjudicate FLAIR and IsoQuant on the chr22 SQANTI-SIM set with
chain_min_full_reads = 2, 3, 5, 10 (with 1 M short-read pairs and without) and score each against
the truth; and, from the decoy tables already computed, the share of wrong recombinant chains the
same thresholds would accept. Writes `min_full_reads_sweep` into the chain_support envelope.

Usage: sweep_min_full_reads.py --pig PATH --work W --mc MC --sj SJ --genome FA --out DIR
         --decoys sqsim_reduced=PREFIX[,truth.tsv] --decoys hg03516=PREFIX ... [--emit-metrics ENV]
"""
import argparse, csv, datetime, json, os, subprocess, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels, load, metrics, table  # noqa: E402

POSITIVE = {"NOVEL_CONFIRMED", "NOVEL_PROBABLE", "REFERENCE_BIAS"}; DECISIVE = POSITIVE | {"ARTIFACT"}
THRESHOLDS = (2, 3, 5, 10)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pig", required=True); ap.add_argument("--work", required=True); ap.add_argument("--mc", required=True)
    ap.add_argument("--sj", required=True); ap.add_argument("--genome", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--cfg", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "config", "rules.default.toml"))
    ap.add_argument("--decoys", action="append", default=[]); ap.add_argument("--emit-metrics")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    callers = {"flair": (f"{a.work}/flair.isoforms.gtf", f"{a.work}/sqanti_out/flair_classification.txt"),
               "isoquant": (f"{a.mc}/isoquant_out/iq/iq.transcript_models.gtf", f"{a.mc}/sq_isoquant/isoquant_classification.txt")}
    base = open(a.cfg).read()
    res = {"thresholds": list(THRESHOLDS), "callers": {}, "decoys": {}}
    print(f"{'caller':9s}{'setting':7s}{'thr':>4s}{'prec':>7s}{'rec':>7s}{'auprc':>7s}{'TP':>5s}{'FP':>4s}{'cov':>6s}{'risk':>7s}")
    for c, (gtf, cls) in callers.items():
        label, cat = build_labels(gtf, f"{a.work}/truth.truth.tsv", cls)
        ids = sorted(t for t in label if cat[t] != "full-splice_match")
        res["callers"][c] = {}
        for setting in ("sr1M", "nosr"):
            res["callers"][c][setting] = {}
            for thr in THRESHOLDS:
                cfg = f"{a.out}/cfg_min{thr}.toml"
                open(cfg, "w").write(base.replace("min_full_reads = 3", f"min_full_reads = {thr}"))
                pre = f"{a.out}/{c}.{setting}.min{thr}"
                if not os.path.exists(pre + ".adjudicated.tsv"):
                    cmd = [a.pig, "adjudicate", "--config", cfg, "--classification", cls, "--isoforms-gtf", gtf, "--ref-gtf", f"{a.work}/chr22_modified.gtf",
                           "--bam", f"{a.work}/aln.bam", "--reference", a.genome, "--out-prefix", pre]
                    if setting == "sr1M": cmd += ["--sj-tab", a.sj]
                    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                pred = load("pig", pre + ".adjudicated.tsv"); m = metrics(pred, ids, label)
                clsmap = {r["isoform_id"]: r["confidence_class"] for r in table(pre + ".adjudicated.tsv")}
                dec = [t for t in ids if clsmap[t] in DECISIVE]
                wrong = sum((clsmap[t] in POSITIVE) != (label[t] == "genuine_novel") for t in dec)
                e = {"TP": m["tp"], "FP": m["fp"], "prec": round(m["prec"], 3), "rec": round(m["rec"], 3), "auprc": round(m["auprc"], 3),
                     "coverage": round(len(dec) / len(ids), 3), "selective_risk": round(wrong / len(dec), 3) if dec else None}
                res["callers"][c][setting][str(thr)] = e
                print(f"{c:9s}{setting:7s}{thr:4d}{e['prec']:7.3f}{e['rec']:7.3f}{e['auprc']:7.3f}{e['TP']:5d}{e['FP']:4d}{e['coverage']:6.2f}{e['selective_risk']:7.3f}")
    for spec in a.decoys:
        lab, rest = spec.split("=", 1); parts = rest.split(","); prefix = parts[0]; truth = dict(l.rstrip("\n").split("\t") for l in open(parts[1])) if len(parts) > 1 else None
        rows = list(csv.DictReader(open(prefix + ".decoys.tsv"), delimiter="\t"))
        both = [r for r in rows if int(r["parent_a_full_reads"]) >= 3 and int(r["parent_b_full_reads"]) >= 3]
        out = {}
        for thr in THRESHOLDS:
            e = {"all": {"n": len(rows), "lifted": sum(int(r["n_full_reads"]) >= thr for r in rows)},
                 "both_parents_ge3": {"n": len(both), "lifted": sum(int(r["n_full_reads"]) >= thr for r in both)}}
            if truth:
                g = [r for r in rows if truth.get(f"{r['chrom']}|{r['strand']}|{r['introns_1based']}") == "genuine_novel"]; f = [r for r in rows if r not in g]
                e["genuine"] = {"n": len(g), "lifted": sum(int(r["n_full_reads"]) >= thr for r in g)}
                e["false"] = {"n": len(f), "lifted": sum(int(r["n_full_reads"]) >= thr for r in f)}
            out[str(thr)] = e
        res["decoys"][lab] = out
        print(f"decoys {lab}: " + "  ".join(f"thr{t}: all {v['all']['lifted']}/{v['all']['n']}, both-parents {v['both_parents_ge3']['lifted']}/{v['both_parents_ge3']['n']}" + (f", genuine {v['genuine']['lifted']}/{v['genuine']['n']}, false {v['false']['lifted']}/{v['false']['n']}" if 'false' in v else "") for t, v in out.items()))
    if a.emit_metrics:
        env = json.load(open(a.emit_metrics)); env["metrics"]["min_full_reads_sweep"] = res
        env["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n"); print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
