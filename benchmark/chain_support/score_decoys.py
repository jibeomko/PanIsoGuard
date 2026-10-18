#!/usr/bin/env python3
"""Summarize `panisoguard decoys` runs: per dataset, the fraction of recombined known chains
(novel-in-catalog decoys) that >= min_full_reads long reads carry end to end, next to the same
fraction for the annotated parents. Where a SQANTI-SIM truth table is given, decoys that equal a
deleted (genuine) transcript are separated from the rest, so the decoy rate can be compared with
the true false-lift rate.

Usage: score_decoys.py --run LABEL=PREFIX [--run ...] [--truth LABEL=truth.truth.tsv ...]
                       [--emit-metrics ../results/chain_support/metrics.json]
PREFIX is the --out of `panisoguard decoys` (PREFIX.decoys.tsv / .json exist).
"""
import argparse, csv, datetime, json, sys


def load_run(prefix):
    summary = json.load(open(prefix + ".decoys.json"))
    rows = list(csv.DictReader(open(prefix + ".decoys.tsv"), delimiter="\t"))
    return summary, rows


def with_truth(rows, truth_path, min_full, both_parents=False):
    truth = dict(l.rstrip("\n").split("\t") for l in open(truth_path))
    out = {"genuine": {"n": 0, "lifted": 0}, "false": {"n": 0, "lifted": 0}}
    for r in rows:
        if both_parents and not (int(r["parent_a_full_reads"]) >= min_full and int(r["parent_b_full_reads"]) >= min_full):
            continue
        key = f"{r['chrom']}|{r['strand']}|{r['introns_1based']}"
        lab = "genuine" if truth.get(key) == "genuine_novel" else "false"
        out[lab]["n"] += 1
        out[lab]["lifted"] += int(r["n_full_reads"]) >= min_full
    for v in out.values():
        v["rate"] = round(v["lifted"] / v["n"], 4) if v["n"] else None
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True, help="LABEL=PREFIX")
    ap.add_argument("--truth", action="append", default=[], help="LABEL=truth.truth.tsv (SQANTI-SIM)")
    ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    truths = dict(t.split("=", 1) for t in a.truth)
    res = {}
    print(f"{'dataset':<14}{'decoys':>8}{'lifted':>8}{'rate':>8} | {'both parents>=N':>16}{'lifted':>8}{'rate':>8} | {'parents':>8}{'lifted':>8}{'rate':>8}")
    for spec in a.run:
        label, prefix = spec.split("=", 1)
        s, rows = load_run(prefix)
        d, b, p = s["decoys"], s["decoys_both_parents_lifted"], s["parents"]
        entry = {"min_full_reads": s["min_full_reads"], "max_per_gene": s["max_per_gene"],
                 "decoys": {"n": d["n"], "lifted": d["lifted"], "rate": round(d["rate"], 4)},
                 "decoys_both_parents_lifted": {"n": b["n"], "lifted": b["lifted"], "rate": round(b["rate"], 4)},
                 "parents": {"n": p["n"], "lifted": p["lifted"], "rate": round(p["rate"], 4)}}
        print(f"{label:<14}{d['n']:>8}{d['lifted']:>8}{d['rate']:>8.4f} | {b['n']:>16}{b['lifted']:>8}{b['rate']:>8.4f} | {p['n']:>8}{p['lifted']:>8}{p['rate']:>8.4f}")
        if label in truths:
            entry["by_truth"] = with_truth(rows, truths[label], s["min_full_reads"])
            entry["by_truth_both_parents_lifted"] = with_truth(rows, truths[label], s["min_full_reads"], both_parents=True)
            g, f = entry["by_truth"]["genuine"], entry["by_truth"]["false"]
            g2, f2 = entry["by_truth_both_parents_lifted"]["genuine"], entry["by_truth_both_parents_lifted"]["false"]
            print(f"  truth: genuine decoys {g['lifted']}/{g['n']} lifted, false decoys {f['lifted']}/{f['n']} lifted (false-lift rate {f['rate']});"
                  f" both parents >= N: genuine {g2['lifted']}/{g2['n']}, false {f2['lifted']}/{f2['n']} ({f2['rate']})")
        res[label] = entry
    if a.emit_metrics:
        env = json.load(open(a.emit_metrics))
        env["metrics"]["decoys"] = res
        env["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        env["notes"] = env.get("notes", "") + (
            " [decoys] `panisoguard decoys` recombines each gene's annotated transcripts at a shared intron into "
            "novel-in-catalog chains that are neither annotated nor a piece of an annotated chain, counts the reads "
            "carrying each end to end, and reports the fraction >= min_full_reads: a per-dataset upper bound on the "
            "chain axis's false-lift rate (some decoys are real, unannotated isoforms; on SQANTI-SIM the truth "
            "separates them).")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
