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
    ap.add_argument("--cross", action="append", default=[], help="LABEL_A,LABEL_B: are the decoys lifted in A "
                    "also carried in B (same annotation + seed -> identical decoy sets), next to A's lifted parents")
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
        entry["_rows"] = rows
    for spec in a.cross:
        la, lb = spec.split(",")
        A = {f"{r['chrom']}|{r['strand']}|{r['introns_1based']}": r for r in res[la]["_rows"]}
        B = {f"{r['chrom']}|{r['strand']}|{r['introns_1based']}": r for r in res[lb]["_rows"]}
        if set(A) != set(B):
            print(f"cross {la},{lb}: decoy sets differ (different annotation or seed), skipped", file=sys.stderr); continue
        n = res[la]["min_full_reads"]
        full = lambda r: int(r["n_full_reads"])
        pmin = lambda r: min(int(r["parent_a_full_reads"]), int(r["parent_b_full_reads"]))
        lifted = [k for k, r in A.items() if full(r) >= n]
        pa = {}
        for r in A.values(): pa[r["parent_a"]] = int(r["parent_a_full_reads"]); pa[r["parent_b"]] = int(r["parent_b_full_reads"])
        pb = {}
        for r in B.values(): pb[r["parent_a"]] = int(r["parent_a_full_reads"]); pb[r["parent_b"]] = int(r["parent_b_full_reads"])
        lp = [t for t, c in pa.items() if c >= n]
        cond = [k for k in lifted if pmin(A[k]) >= n and pmin(B[k]) >= n]
        c = {"lifted_decoys_in_a": len(lifted), "also_lifted_in_b": sum(full(B[k]) >= n for k in lifted),
             "also_ge1_in_b": sum(full(B[k]) >= 1 for k in lifted),
             "lifted_parents_in_a": len(lp), "parents_also_lifted_in_b": sum(pb.get(t, 0) >= n for t in lp),
             "lifted_decoys_both_parents_expressed_in_both": len(cond),
             "of_those_also_lifted_in_b": sum(full(B[k]) >= n for k in cond)}
        res.setdefault("cross_library", {})[f"{la}->{lb}"] = c
        print(f"cross {la}->{lb}: lifted decoys {c['lifted_decoys_in_a']} -> >= {n} in {lb}: {c['also_lifted_in_b']} "
              f"({c['also_lifted_in_b'] / max(1, c['lifted_decoys_in_a']):.2f}); lifted parents {c['lifted_parents_in_a']} -> "
              f"{c['parents_also_lifted_in_b']} ({c['parents_also_lifted_in_b'] / max(1, c['lifted_parents_in_a']):.2f}); "
              f"both parents expressed in both: {c['of_those_also_lifted_in_b']} / {c['lifted_decoys_both_parents_expressed_in_both']}")
    for v in res.values():
        if isinstance(v, dict): v.pop("_rows", None)
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
