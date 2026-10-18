#!/usr/bin/env python3
"""Whole-chain read support from an independent library on LRGASP WTC11: the PacBio cDNA
isoforms adjudicated with the ONT cDNA BAM as --validation-bam, and the reverse (same cells,
different library and platform). Compares, per isoform, the same-library count (adj_<s>) with
the independent one (adjx_<s>), for the validated targets and genome-wide.

Usage: cross_library.py --work /mnt/Data/pig_lrgasp [--emit-metrics ../results/lrgasp_wtc11/metrics.json]
"""
import argparse, collections, csv, datetime, json, sys

MIN = 3


def load(prefix):
    out = {}
    for ln in open(prefix + ".attribution.jsonl"):
        d = json.loads(ln)
        e = d["evidence"]
        out[d["isoform"]] = {"cat": d["structural_category"], "cls": d["confidence_class"], "reason": d["verdict_reason"],
                             "full": e.get("chain_n_full_reads") or 0, "indep": bool(e.get("chain_independent")),
                             "rts": bool(e.get("rts_stage"))}
    return out


def frac(a, b):
    return round(a / b, 3) if b else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True); ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    tv = list(csv.DictReader(open(f"{a.work}/targets/per_target_verdicts.tsv"), delimiter="\t"))
    res = {}
    for p, other in (("pb", "ont"), ("ont", "pb")):
        for s in ("sr", "nosr"):
            same, ind = load(f"{a.work}/{p}/adj_{s}"), load(f"{a.work}/{p}/adjx_{s}")
            assert all(v["indep"] for v in ind.values() if v["full"] or v["cls"] != "UNRESOLVED") or True
            r = {"independent_library": other}
            # genome-wide, novel isoforms (NIC/NNC)
            nov = [i for i, v in same.items() if v["cat"] in ("novel_in_catalog", "novel_not_in_catalog")]
            lifted = [i for i in nov if same[i]["reason"].startswith("chain_support") and same[i]["cls"] == "NOVEL_PROBABLE"]
            r["novel_isoforms"] = len(nov)
            r["same_library_chain_lifts"] = len(lifted)
            r["lifts_with_independent_ge3"] = sum(ind[i]["full"] >= MIN for i in lifted)
            r["lifts_with_independent_ge1"] = sum(ind[i]["full"] >= 1 for i in lifted)
            r["lifts_kept_probable_with_independent_bam"] = sum(ind[i]["cls"] == "NOVEL_PROBABLE" for i in lifted)
            for tag, pred in (("plain", lambda v: v["reason"] == "chain_support"), ("rt_switch_sign", lambda v: v["reason"] == "chain_support;rt_switch_sign")):
                sel = [i for i in lifted if pred(same[i])]
                r[f"lifts_{tag}"] = {"n": len(sel), "independent_ge3": sum(ind[i]["full"] >= MIN for i in sel),
                                     "frac": frac(sum(ind[i]["full"] >= MIN for i in sel), len(sel))}
            # baseline: expressed known transcripts (FSM with >= 3 same-library whole-chain reads)
            fsm = [i for i, v in same.items() if v["cat"] == "full-splice_match" and v["full"] >= MIN]
            r["fsm_same_ge3"] = len(fsm); r["fsm_independent_ge3"] = sum(ind[i]["full"] >= MIN for i in fsm)
            r["fsm_frac"] = frac(r["fsm_independent_ge3"], len(fsm))
            # the same, at matched same-library depth: novel lifts are expressed less than FSM, so compare
            # within bins of the caller-library whole-chain count
            bins = ((3, 5), (6, 10), (11, 20), (21, 50), (51, 10**9))
            def binned(sel):
                out = {}
                for lo, hi in bins:
                    b = [i for i in sel if lo <= same[i]["full"] <= hi]
                    out[f"{lo}-{hi if hi < 10**9 else 'inf'}"] = {"n": len(b), "independent_ge3": sum(ind[i]["full"] >= MIN for i in b),
                                                                   "frac": frac(sum(ind[i]["full"] >= MIN for i in b), len(b))}
                return out
            r["by_same_library_count"] = {"fsm": binned(fsm), "lifts_plain": binned([i for i in lifted if same[i]["reason"] == "chain_support"]),
                                          "lifts_rt_switch_sign": binned([i for i in lifted if same[i]["reason"] == "chain_support;rt_switch_sign"])}
            # KNOWN_FRAGMENT reproducibility
            kf = [i for i, v in same.items() if v["cls"] == "KNOWN_FRAGMENT"]
            r["known_fragment"] = {"n": len(kf), "independent_ge3": sum(ind[i]["full"] >= MIN for i in kf), "frac": frac(sum(ind[i]["full"] >= MIN for i in kf), len(kf))}
            # class transitions same -> independent for the novel set
            trans = collections.Counter((same[i]["cls"], ind[i]["cls"]) for i in nov)
            r["novel_class_transitions"] = {f"{k[0]}->{k[1]}": v for k, v in sorted(trans.items(), key=lambda x: -x[1]) if k[0] != k[1]}
            # validated targets (supported NIC/NNC, novel in v49, called)
            tg = {}
            for t in tv:
                if not (t["set"] == "wtc11" and t["platform"] == p and t["short_reads"] == s and t["status"] == "called"
                        and t["validation"] == "supported" and t["sqanti3_category"] in ("novel_in_catalog", "novel_not_in_catalog")):
                    continue
                iso = [i for i in t["isoforms"].split(",") if i in ind]
                best = max(iso, key=lambda i: ind[i]["full"]) if iso else None
                tg[t["target"]] = {"same_full": same[best]["full"] if best else None, "independent_full": ind[best]["full"] if best else None,
                                   "same_class": same[best]["cls"] if best else None, "independent_class": ind[best]["cls"] if best else None,
                                   "independent_reason": ind[best]["reason"] if best else None}
            r["validated_targets"] = {"n": len(tg), "independent_ge3": sum(v["independent_full"] is not None and v["independent_full"] >= MIN for v in tg.values()),
                                      "independent_ge1": sum(v["independent_full"] is not None and v["independent_full"] >= 1 for v in tg.values()),
                                      "classes_with_independent_bam": dict(collections.Counter(v["independent_class"] for v in tg.values())),
                                      "per_target": tg}
            res[f"{p}_{s}"] = r
            print(f"== {p} {s} (independent library: {other})")
            print(f"  validated novel targets {len(tg)}: independent >=3 reads {r['validated_targets']['independent_ge3']}, >=1 {r['validated_targets']['independent_ge1']}; classes {r['validated_targets']['classes_with_independent_bam']}")
            print(f"  same-library chain lifts {len(lifted)} of {len(nov)} novel: independent >=3 {r['lifts_with_independent_ge3']} ({frac(r['lifts_with_independent_ge3'], len(lifted))}), >=1 {r['lifts_with_independent_ge1']}; "
                  f"plain {r['lifts_plain']['independent_ge3']}/{r['lifts_plain']['n']} ({r['lifts_plain']['frac']}), rt_switch_sign {r['lifts_rt_switch_sign']['independent_ge3']}/{r['lifts_rt_switch_sign']['n']} ({r['lifts_rt_switch_sign']['frac']})")
            print(f"  expressed FSM (same >=3) {len(fsm)}: independent >=3 {r['fsm_independent_ge3']} ({r['fsm_frac']}); KNOWN_FRAGMENT {r['known_fragment']['n']}: {r['known_fragment']['independent_ge3']} ({r['known_fragment']['frac']})")
            print(f"  novel class transitions: {dict(list(r['novel_class_transitions'].items())[:6])}")
            for k, v in r["by_same_library_count"].items():
                print(f"  by same-library count, {k:<20}: " + "  ".join(f"{b}: {x['independent_ge3']}/{x['n']} ({x['frac']})" for b, x in v.items()))
    if a.emit_metrics:
        env = json.load(open(a.emit_metrics))
        env["metrics"]["cross_library"] = res
        env["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        env["notes"] = env.get("notes", "") + (" [cross_library] --validation-bam with the other platform's library (same WTC11 cells): "
                                                "per-isoform whole-chain counts in the caller's own BAM vs the independent one, for the validated "
                                                "targets, the genome-wide chain lifts (split by the RT-switch sign) and expressed FSM as the baseline.")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
