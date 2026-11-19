#!/usr/bin/env python3
"""Is the whole-chain evidence a read-count filter? Without short reads, on the SQANTI-SIM chr22
truth set: the caller's own read count as a filter and as a ranking, SQANTI3's default rules with
that count required as well (the `FL` descriptor), the whole-chain read count on its own, and
PanIsoGuard with the BAM.

Two modes.
  --write-fl CLASSIFICATION COUNTS OUT
      copy a SQANTI3 classification with its FL column filled from COUNTS (isoform<TAB>count), as
      `sqanti3_qc.py --fl_count` would fill it, so that `sqanti3_filter.py rules` can read it.
  --work W --caller NAME=GTF,CLASSIFICATION,COUNTS,PIG_PREFIX,RULES[,FL3_RULES,FL10_RULES] ...
      score. PIG_PREFIX is the `adjudicate --bam` run without short reads; RULES is SQANTI3's
      default rules filter output, FL<k>_RULES the same filter with FL >= k added to every rule.
      With --emit-metrics the result goes into the envelope as `read_count_baseline`.

A caller's count assigns each read to one isoform by the caller's own compatibility rule; a
whole-chain read carries every intron of the chain at exact coordinates in one alignment.
"""
import argparse, csv, datetime, json, os, random, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from score_h2h import build_labels, load, metrics  # noqa: E402

THRESHOLDS = (3, 5, 10, 20)


def read_counts(path):
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        try: out[f[0]] = float(f[1])
        except (IndexError, ValueError): continue   # header and summary lines
    return out


def write_fl(cls, counts, out):
    cnt = read_counts(counts)
    with open(cls) as f, open(out, "w") as o:
        hdr = f.readline().rstrip("\n").split("\t"); i = hdr.index("FL"); o.write("\t".join(hdr) + "\n")
        for line in f:
            r = line.rstrip("\n").split("\t"); r[i] = str(int(round(cnt.get(r[0], 0)))); o.write("\t".join(r) + "\n")


def summary(m): return {"TP": m["tp"], "FP": m["fp"], "prec": round(m["prec"], 3), "rec": round(m["rec"], 3), "auprc": round(m["auprc"], 3)}


def precision_at_recall(score, ids, label, recall):
    """Precision of the ranking at the loosest threshold that still gives at most `recall`."""
    n_pos = sum(label[t] == "genuine_novel" for t in ids); best = None
    for thr in sorted({score[t] for t in ids}, reverse=True):
        kept = [t for t in ids if score[t] >= thr]; tp = sum(label[t] == "genuine_novel" for t in kept)
        if tp / n_pos > recall: break
        best = {"threshold": thr, "TP": tp, "FP": len(kept) - tp, "prec": round(tp / len(kept), 3), "rec": round(tp / n_pos, 3)}
    return best


def one(work, spec, resamples_of):
    f = spec.split(","); gtf, cls, counts, pig, rules = f[:5]
    label, cat = build_labels(gtf, f"{work}/truth.truth.tsv", cls)
    ids = sorted(t for t in label if cat[t] != "full-splice_match")
    cnt = read_counts(counts); c = {t: cnt.get(t, 0.0) for t in ids}
    full = {}
    for line in open(pig + ".attribution.jsonl"):
        d = json.loads(line); full[d["isoform"]] = d["evidence"].get("chain_n_full_reads") or 0
    w = {t: full.get(t, 0) for t in ids}
    gen = [t for t in ids if label[t] == "genuine_novel"]; fal = [t for t in ids if label[t] == "false_novel"]
    preds = {"caller_count": {t: (c[t] >= 3, c[t]) for t in ids}, "whole_chain": {t: (w[t] >= 3, float(w[t])) for t in ids},
             "SQ3rules": load("rules", rules), "PIG+bam": load("pig", pig + ".adjudicated.tsv")}
    res = {"n": len(ids), "genuine": len(gen), "false": len(fal), "isoforms_without_a_count": sum(t not in cnt for t in ids),
           "caller_count": {"auprc": summary(metrics(preds["caller_count"], ids, label))["auprc"],
                            **{f"ge{k}": {x: v for x, v in summary(metrics({t: (c[t] >= k, c[t]) for t in ids}, ids, label)).items() if x != "auprc"} for k in THRESHOLDS}},
           "whole_chain_ge3": summary(metrics(preds["whole_chain"], ids, label)),
           "SQ3rules": summary(metrics(preds["SQ3rules"], ids, label)), "PIG+bam": summary(metrics(preds["PIG+bam"], ids, label))}
    for k, path in zip((3, 10), f[5:7]):
        p = load("rules", path); preds[f"SQ3rules_FL{k}"] = p; res[f"SQ3rules_FL{k}"] = summary(metrics(p, ids, label))
    pig_rec = res["PIG+bam"]["rec"]
    res["at_the_recall_of_PIG+bam"] = {"caller_count": precision_at_recall(c, ids, label, pig_rec), "whole_chain": precision_at_recall(w, ids, label, pig_rec)}
    res["same_threshold_3"] = {"false_passed_by_caller_count": sum(c[t] >= 3 for t in fal), "false_passed_by_whole_chain": sum(w[t] >= 3 for t in fal),
                               "genuine_passed_by_caller_count": sum(c[t] >= 3 for t in gen), "genuine_passed_by_whole_chain": sum(w[t] >= 3 for t in gen)}
    resamples = resamples_of(ids)
    def ci(a, b):
        d = sorted(metrics(preds[a], r, label)["auprc"] - metrics(preds[b], r, label)["auprc"] for r in resamples)
        return [round(d[int(0.025 * len(d))], 3), round(d[int(0.975 * len(d)) - 1], 3)]
    res["bootstrap_auprc_diff_95ci"] = {"PIG+bam_minus_caller_count": ci("PIG+bam", "caller_count"), "whole_chain_minus_caller_count": ci("whole_chain", "caller_count"),
                                        "caller_count_minus_SQ3rules": ci("caller_count", "SQ3rules")}
    return res


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--write-fl":
        write_fl(*sys.argv[2:5]); return
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--caller", action="append", required=True); ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    def resamples_of(ids):
        rng = random.Random(1); return [[ids[rng.randrange(len(ids))] for _ in ids] for _ in range(2000)]
    out = {}
    for spec in a.caller:
        name, rest = spec.split("=", 1); r = out[name] = one(a.work, rest, resamples_of); cc = r["caller_count"]
        print(f"== {name}: {r['n']} non-FSM isoforms ({r['genuine']} genuine, {r['false']} false)")
        for k in THRESHOLDS: print(f"  caller count >= {k:<2d}            TP {cc[f'ge{k}']['TP']:4d} FP {cc[f'ge{k}']['FP']:4d}  prec {cc[f'ge{k}']['prec']:.3f}  rec {cc[f'ge{k}']['rec']:.3f}")
        print(f"  caller count as a ranking     AUPRC {cc['auprc']:.3f}")
        for m in ("SQ3rules", "SQ3rules_FL3", "SQ3rules_FL10", "whole_chain_ge3", "PIG+bam"):
            if m in r: print(f"  {m:28s}  TP {r[m]['TP']:4d} FP {r[m]['FP']:4d}  prec {r[m]['prec']:.3f}  rec {r[m]['rec']:.3f}  AUPRC {r[m]['auprc']:.3f}")
        print(f"  at the recall of PIG+bam: {r['at_the_recall_of_PIG+bam']}")
        print(f"  three reads: {r['same_threshold_3']}")
        print(f"  AUPRC difference, 95 % CI: {r['bootstrap_auprc_diff_95ci']}")
    if a.emit_metrics:
        env = json.load(open(a.emit_metrics)); env["metrics"]["read_count_baseline"] = out
        env["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n"); print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
