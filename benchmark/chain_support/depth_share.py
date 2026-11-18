#!/usr/bin/env python3
"""Would a floor on the share of spanning reads improve the whole-chain lift?

`evidence.chain_n_spanning_reads` is the depth a chain is seen at: the reads at its first donor
that reach past its last intron. `[axis_chain] min_fraction` asks that the whole-chain reads be at
least that share of them. This script reads runs made with the default (no floor) and counts, for
each floor, the chain-lifted calls that would remain (a lift remains when
chain_n_full_reads >= floor x chain_n_spanning_reads, the rule the engine applies), split by the
truth each data set has:
  SQANTI-SIM chr22 (FLAIR, IsoQuant; no short reads)   genuine / false novel chains
  decoy chains (panisoguard decoys)                    by the truth where there is one
  LRGASP amplicon targets run as would-be isoforms     supported / unsupported / likely
  TUSCO single-isoform genes                           novel calls there are presumed false
  whole real runs                                      every chain-lifted isoform (what a floor costs)
Writes `depth_share` into the chain_support envelope.

Usage: depth_share.py --work W --chain OUT --mc MC --lrgasp L --tusco tusco_human.tsv
         --decoys LABEL=PREFIX[,truth.tsv] ... --run LABEL=CLASSIFICATION,PREFIX ... [--emit-metrics ENV]
"""
import argparse, collections, csv, datetime, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels  # noqa: E402

FLOORS = (0.005, 0.01, 0.02, 0.05, 0.10)
NOVEL = ("novel_in_catalog", "novel_not_in_catalog")


def lifts(jsonl):
    """isoform -> (whole-chain reads, spanning reads) of every isoform the chain axis lifted."""
    out = {}
    for line in open(jsonl):
        d = json.loads(line)
        if d["confidence_class"] == "NOVEL_PROBABLE" and d["verdict_reason"].startswith("chain_support"):
            e = d["evidence"]; out[d["isoform"]] = (e["chain_n_full_reads"], e["chain_n_spanning_reads"])
    return out


def kept(pairs):
    pairs = list(pairs)
    share = sorted(f / s for f, s in pairs if s)
    return {"lifted": len(pairs), "lowest_share": round(share[0], 4) if share else None,
            "median_share": round(share[len(share) // 2], 4) if share else None,
            "kept_at_floor": {str(t): sum(1 for f, s in pairs if f >= t * s) for t in FLOORS}}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--chain", required=True); ap.add_argument("--mc", required=True)
    ap.add_argument("--lrgasp", required=True); ap.add_argument("--tusco", required=True)
    ap.add_argument("--decoys", action="append", default=[]); ap.add_argument("--run", action="append", default=[]); ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    res = {"floors": list(FLOORS), "sqanti_sim": {}, "decoys": {}, "lrgasp_targets": {}, "tusco": {}, "real_runs": {}}

    callers = {"flair": (f"{a.work}/flair.isoforms.gtf", f"{a.work}/sqanti_out/flair_classification.txt"),
               "isoquant": (f"{a.mc}/isoquant_out/iq/iq.transcript_models.gtf", f"{a.mc}/sq_isoquant/isoquant_classification.txt")}
    for c, (gtf, cls) in callers.items():
        label, _ = build_labels(gtf, f"{a.work}/truth.truth.tsv", cls)
        up = lifts(f"{a.chain}/{c}.nosr.on.attribution.jsonl")
        res["sqanti_sim"][c] = {l: kept(v for t, v in up.items() if label.get(t) == l) for l in ("genuine_novel", "false_novel")}
        print(f"sqanti_sim {c}: " + "; ".join(f"{l} {d['lifted']} (lowest share {d['lowest_share']}) kept {d['kept_at_floor']}" for l, d in res["sqanti_sim"][c].items()))

    for spec in a.decoys:
        lab, rest = spec.split("=", 1); parts = rest.split(",")
        truth = dict(l.rstrip("\n").split("\t") for l in open(parts[1])) if len(parts) > 1 else None
        rows = [r for r in csv.DictReader(open(parts[0] + ".decoys.tsv"), delimiter="\t") if r["lifted"] == "true"]
        pair = lambda r: (int(r["n_full_reads"]), int(r["n_spanning_reads"]))
        e = {"all": kept(map(pair, rows)),
             "both_parents_ge3": kept(pair(r) for r in rows if int(r["parent_a_full_reads"]) >= 3 and int(r["parent_b_full_reads"]) >= 3)}
        if truth:
            gen = lambda r: truth.get(f"{r['chrom']}|{r['strand']}|{r['introns_1based']}") == "genuine_novel"
            e["genuine"] = kept(pair(r) for r in rows if gen(r)); e["false"] = kept(pair(r) for r in rows if not gen(r))
        res["decoys"][lab] = e
        print(f"decoys {lab}: " + "; ".join(f"{k} {d['lifted']} kept {d['kept_at_floor']}" for k, d in e.items()))

    L = a.lrgasp; W = f"{L}/targets_run"; tid = {}; cat = {}
    for setname, g, tsv in (("wtc11", "targets.gtf", "targets.tsv"), ("gencode", "targets_gencode.gtf", "targets_gencode.tsv")):
        for l in open(f"{L}/targets/{g}"):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 9: tid[f[8].split('transcript_id "')[1].split('"')[0]] = (setname, f[8].split('gene_id "')[1].split('"')[0])
        for r in csv.DictReader(open(f"{L}/targets/{tsv}"), delimiter="\t"): cat[(setname, r["target"])] = r["category"]
    sq = {r["isoform"]: r["structural_category"] for r in csv.DictReader(open(f"{W}/sq/tg_classification.txt"), delimiter="\t")}
    for p in ("pb", "ont"):
        up = lifts(f"{W}/adj_{p}_nosr.attribution.jsonl")
        res["lrgasp_targets"][p] = {v: kept(x for t, x in up.items() if t in tid and cat.get(tid[t]) == v and sq.get(t) in NOVEL)
                                    for v in ("supported", "unsupported", "likely")}
        print(f"lrgasp targets {p}: " + "; ".join(f"{v} {d['lifted']} (lowest share {d['lowest_share']}) kept {d['kept_at_floor']}" for v, d in res["lrgasp_targets"][p].items()))

    rows = [l.rstrip("\n").split("\t") for l in open(a.tusco) if l.strip() and not l.startswith("#")]
    ids = {r[0].split(".")[0] for r in rows}; names = {r[2] for r in rows if len(r) > 2}
    in_tusco = lambda gene: any(x.split(".")[0] in ids or x in names for x in gene.split("_"))
    for spec in a.run:
        lab, files = spec.split("=", 1); cls_path, prefix = files.split(",")
        cls = {r["isoform"]: r for r in csv.DictReader(open(cls_path), delimiter="\t")}
        up = lifts(prefix + ".attribution.jsonl")
        res["real_runs"][lab] = kept(up.values())
        res["tusco"][lab] = kept(x for t, x in up.items() if t in cls and cls[t]["structural_category"] in NOVEL and in_tusco(cls[t]["associated_gene"]))
        print(f"{lab}: every chain lift {res['real_runs'][lab]['lifted']} kept {res['real_runs'][lab]['kept_at_floor']}; in TUSCO genes {res['tusco'][lab]['lifted']} kept {res['tusco'][lab]['kept_at_floor']}")

    if a.emit_metrics:
        env = json.load(open(a.emit_metrics)); env["metrics"]["depth_share"] = res
        env["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n"); print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
