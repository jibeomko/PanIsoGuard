#!/usr/bin/env python3
"""SQANTI3's ML filter trained on the truth. Its default TP / TN labels (FSM / NNC) are wrong for the
SQANTI-SIM set, where most NNC are genuine, so the SQ3ml row of the head-to-head is not a fair
ceiling. Here the random forest is trained with the simulated truth itself (genuine = TP, false =
TN) in k-fold cross-validation over the non-FSM isoforms and scored on the held-out fold. The
filter refuses training sets under 250 isoforms per class, which rules out FLAIR (140 false
isoforms in all); TALON's permissive output (697 genuine / 1,063 false) is the one caller that
qualifies.

Usage: ml_truth_cv.py --work W --gtf CALLER.gtf --cls CALLER_classification.txt --label talon --out DIR
         [--folds 2] [--emit-metrics ../results/caller_filters/metrics.json]
"""
import argparse, datetime, json, os, random, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_h2h import build_labels, metrics, table  # noqa: E402

MIN_PER_CLASS = 250  # sqanti3_filter.py ml: "user-defined TP and TN sets must have >=250 isoforms"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--gtf", required=True); ap.add_argument("--cls", required=True)
    ap.add_argument("--label", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--folds", type=int, default=2); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--emit-metrics")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    label, cat = build_labels(a.gtf, f"{a.work}/truth.truth.tsv", a.cls)
    ids = sorted(t for t in label if cat[t] != "full-splice_match")
    rng = random.Random(a.seed); order = ids[:]; rng.shuffle(order)
    folds = [order[k::a.folds] for k in range(a.folds)]
    pred = {}
    for k, test in enumerate(folds):
        tset = set(test); train = [t for t in ids if t not in tset]
        tp = [t for t in train if label[t] == "genuine_novel"]; tn = [t for t in train if label[t] == "false_novel"]
        if min(len(tp), len(tn)) < MIN_PER_CLASS:
            sys.exit(f"{a.label}: fold {k} has TP {len(tp)} / TN {len(tn)}; sqanti3_filter ml needs >= {MIN_PER_CLASS} of each")
        d = f"{a.out}/{a.label}/fold{k}"; os.makedirs(d, exist_ok=True)
        open(f"{d}/tp.txt", "w").write("\n".join(tp) + "\n"); open(f"{d}/tn.txt", "w").write("\n".join(tn) + "\n")
        out_cls = f"{d}/ml_ML_classification.txt"
        if not os.path.exists(out_cls):
            subprocess.run(["sqanti3_filter.py", "ml", "--sqanti_class", a.cls, "--TP", f"{d}/tp.txt", "--TN", f"{d}/tn.txt", "-o", "ml", "-d", d, "--skip_report"],
                           check=True, stdout=open(f"{d}/log.txt", "w"), stderr=subprocess.STDOUT)
        if "Skipping ML filter" in open(f"{d}/log.txt").read():
            sys.exit(f"{a.label}: sqanti3_filter ml skipped the model (see {d}/log.txt)")
        for r in table(out_cls):
            if r["isoform"] in tset:
                prob = float(r["POS_MLprob"]) if r["POS_MLprob"] not in ("NA", "") else 0.0
                pred[r["isoform"]] = (r["filter_result"] == "Isoform", prob)
    m = metrics(pred, ids, label)
    res = {"TP": m["tp"], "FP": m["fp"], **{k: round(m[k], 3) for k in ("prec", "rec", "spec", "f1", "auprc")},
           "folds": a.folds, "seed": a.seed, "n_scored": len(pred), "n_genuine": sum(label[t] == "genuine_novel" for t in ids)}
    print(a.label, res)
    if a.emit_metrics:
        env = json.load(open(a.emit_metrics))
        env["metrics"].setdefault("sqanti3_ml_truth_cv", {})[a.label] = res
        env["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        if "[sqanti3_ml_truth_cv]" not in env.get("notes", ""):
            env["notes"] = env.get("notes", "") + (" [sqanti3_ml_truth_cv] SQANTI3's ML filter trained on the simulated truth itself (genuine = TP, "
                                                    "false = TN) in two-fold cross-validation and scored on the held-out fold: the supervised ceiling of "
                                                    "SQANTI3's own features. Only TALON's output has the >= 250 false isoforms the filter requires per class.")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n"); print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
