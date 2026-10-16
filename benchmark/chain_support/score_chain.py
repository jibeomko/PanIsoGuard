#!/usr/bin/env python3
"""Score the whole-chain long-read support axis (ON vs OFF) per caller against SQANTI-SIM truth.

For each caller and short-read setting: confident-novel precision / recall / AUPRC (as in
../sqanti3_filter_h2h), the OFF->ON transitions split by truth label (what the axis lifted),
and the risk-coverage point (decisive = NOVEL_CONFIRMED / NOVEL_PROBABLE / ARTIFACT /
REFERENCE_BIAS; selective risk = wrong decisive calls / decisive calls).

Usage: score_chain.py --work W --out OUTDIR --gtf caller:gtf ... --cls caller:classification ... [--emit-metrics m.json]
"""
import argparse, collections, datetime, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels, load, metrics, table  # noqa: E402

DECISIVE = {"NOVEL_CONFIRMED", "NOVEL_PROBABLE", "ARTIFACT", "REFERENCE_BIAS"}
POSITIVE = {"NOVEL_CONFIRMED", "NOVEL_PROBABLE", "REFERENCE_BIAS"}


def risk_coverage(cls, ids, label):
    dec = [t for t in ids if cls[t] in DECISIVE]
    wrong = sum(1 for t in dec if (cls[t] in POSITIVE) != (label[t] == "genuine_novel"))
    return {"coverage": round(len(dec) / len(ids), 3) if ids else None,
            "selective_risk": round(wrong / len(dec), 3) if dec else None, "n_decisive": len(dec), "n_wrong": wrong}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--gtf", action="append", default=[]); ap.add_argument("--cls", action="append", default=[])
    ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    gtf = dict(x.split(":", 1) for x in a.gtf); cls_path = dict(x.split(":", 1) for x in a.cls)
    res = {}
    print(f"{'caller':9s} {'setting':6s} {'n':>5s} {'genuine':>7s}  {'prec OFF->ON':>14s}  {'rec OFF->ON':>13s}  {'AUPRC OFF->ON':>14s}  lifted genuine/false  risk@cov OFF -> ON")
    for c in gtf:
        label, cat = build_labels(gtf[c], f"{a.work}/truth.truth.tsv", cls_path[c])
        ids = sorted(t for t in label if cat[t] != "full-splice_match")
        res[c] = {"n_non_fsm": len(ids), "n_genuine": sum(label[t] == "genuine_novel" for t in ids)}
        for setting in ("sr1M", "nosr"):
            out = {}
            pred = {}; cls = {}
            for state in ("off", "on"):
                p = f"{a.out}/{c}.{setting}.{state}.adjudicated.tsv"
                pred[state] = load("pig", p)
                cls[state] = {r["isoform_id"]: r["confidence_class"] for r in table(p)}
                m = metrics(pred[state], ids, label)
                out[state] = {k: (None if m[k] != m[k] else round(m[k], 3)) for k in ("prec", "rec", "f1", "auprc")}
                out[state].update(TP=m["tp"], FP=m["fp"], **risk_coverage(cls[state], ids, label))
            trans = collections.Counter((label[t], cls["off"][t], cls["on"][t]) for t in ids if cls["off"][t] != cls["on"][t])
            out["transitions_off_to_on"] = {f"{l}:{x}->{y}": n for (l, x, y), n in sorted(trans.items())}
            lifted_g = sum(n for (l, x, y), n in trans.items() if l == "genuine_novel" and y in POSITIVE)
            lifted_f = sum(n for (l, x, y), n in trans.items() if l == "false_novel" and y in POSITIVE)
            out["lifted_to_confident"] = {"genuine": lifted_g, "false": lifted_f}
            res[c][setting] = out
            o, n = out["off"], out["on"]
            print(f"{c:9s} {setting:6s} {len(ids):5d} {res[c]['n_genuine']:7d}  {str(o['prec']):>6s}->{str(n['prec']):<6s}  {str(o['rec']):>6s}->{str(n['rec']):<6s}  {str(o['auprc']):>6s}->{str(n['auprc']):<6s}  {lifted_g:7d}/{lifted_f:<5d}       {o['selective_risk']}@{o['coverage']} -> {n['selective_risk']}@{n['coverage']}")
    if a.emit_metrics:
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{a.out}/flair.sr1M.on.provenance.log") if "\t" in l and not l.startswith("#"))
        env = {"protocol": "chain_support", "scope": "gencode_v49_chr22_5caller_sqanti_sim", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "SQANTI-SIM GENCODE v49 chr22 (PBSIM3 HiFi, one shared minimap2 alignment) called by FLAIR, "
                                  "IsoQuant, Bambu, ESPRESSO and TALON (the benchmark/multicaller set), each classified with "
                                  "SQANTI3 6.0.1; short reads = 1 M wgsim pairs -> STAR (sqanti3_filter_h2h) or none. "
                                  "Public/simulated, no private data.",
               "command": "bash benchmark/chain_support/run.sh (adjudicate per caller with the chain axis ON = default and OFF = "
                          "min_full_reads 1000000, with and without --sj-tab); benchmark/chain_support/score_chain.py",
               "metrics": res,
               "notes": "Whole-chain read support (>= 3 primary reads carrying the entire intron chain, same BAM the callers "
                        "used) lifts UNKNOWN-support isoforms to NOVEL_PROBABLE. It separates caller-level false chains "
                        "(FLAIR, IsoQuant: their false novels mostly lack full-chain reads) from alignment-level ones "
                        "(Bambu, ESPRESSO, TALON: their false novels are carried by reads too), so the lift is a recall gain "
                        "at little precision cost for the former and a precision cost for the latter. Non-independent "
                        "evidence by construction: it never reaches NOVEL_CONFIRMED and never overrides a short-read verdict."}
        os.makedirs(os.path.dirname(a.emit_metrics), exist_ok=True)
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
