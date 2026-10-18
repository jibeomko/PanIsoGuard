#!/usr/bin/env python3
"""Sensitivity of the tracked results to the ruleset-0.0.5 RT-switch policy.

Since 0.0.5 SQANTI3's `RTS_stage` no longer demotes a caller-consensus / whole-chain lift to
UNCONFIRMED; it is noted in the reason (`chain_support;rt_switch_sign`). The alternative policy
("hold": those isoforms are UNCONFIRMED, as in 0.0.4) can be read off the same outputs, so this
script tabulates, per tracked run, how many calls the policy decides and, where a truth or an
external filter exists, what they are. The ISM counterpart (`ism_chain_support;rts_sign`, a
KNOWN_FRAGMENT with the sign) is listed too.

Usage: rts_sensitivity.py --sqsim-work W --h2h H2H --e2e E2E_WORK --lrgasp L --threeway hg03516=DIR hg02717=DIR
                          --gm12878 rb_ont.adjudicated.tsv [--emit-metrics ../results/rts_sensitivity/metrics.json]
"""
import argparse, collections, csv, datetime, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import build_labels, load, metrics, table  # noqa: E402

RT_LIFT = lambda r: r["verdict_reason"] in ("chain_support;rt_switch_sign", "caller_consensus;rt_switch_sign",
                                            "chain_support;rt_switch_sign;independent_library")
KF_RTS = lambda r: r["verdict_reason"].startswith("ism_chain_support;rts_sign")


def hold_policy(rows):
    """The same rows under 'hold': RT-switch-annotated lifts become UNCONFIRMED."""
    out = []
    for r in rows:
        r = dict(r)
        if RT_LIFT(r):
            r["confidence_class"] = "UNCONFIRMED"; r["verdict"] = "unconfirmed"
            r["verdict_reason"] = r["verdict_reason"].replace("chain_support;rt_switch_sign", "chain_support_with_rt_switch") \
                                                     .replace("caller_consensus;rt_switch_sign", "consensus_with_rt_switch")
        out.append(r)
    return out


def pig_pred(rows):
    from score_h2h import RANK
    return {r["isoform_id"]: (r["confidence_class"] in ("NOVEL_CONFIRMED", "NOVEL_PROBABLE"), RANK.get(r["confidence_class"], 0.0)) for r in rows}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sqsim-work", required=True); ap.add_argument("--h2h", required=True)
    ap.add_argument("--e2e"); ap.add_argument("--lrgasp"); ap.add_argument("--threeway", nargs="*", default=[])
    ap.add_argument("--gm12878"); ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    res = {}
    W = a.sqsim_work
    label, cat = build_labels(f"{W}/flair.isoforms.gtf", f"{W}/truth.truth.tsv", f"{W}/sqanti_out/flair_classification.txt")
    ids = sorted(t for t in label if cat[t] != "full-splice_match")
    print(f"{'run':<22}{'RT lifts':>9}{'genuine':>8}{'false':>6} | annotate P / R / AUPRC -> hold P / R / AUPRC")
    for s in ("sr_1000000", "none"):
        rows = table(f"{a.h2h}/{s}/pigbam.adjudicated.tsv")
        lifts = [r for r in rows if RT_LIFT(r)]
        g = sum(label.get(r["isoform_id"]) == "genuine_novel" for r in lifts)
        f = sum(label.get(r["isoform_id"]) == "false_novel" for r in lifts)
        m1 = metrics(pig_pred(rows), ids, label); m2 = metrics(pig_pred(hold_policy(rows)), ids, label)
        kf = sum(KF_RTS(r) for r in rows)
        res[f"sqanti_sim_{s}_bam"] = {"rt_lifts": len(lifts), "genuine": g, "false": f, "known_fragment_rts_sign": kf,
                                      "annotate": {k: round(m1[k], 3) for k in ("prec", "rec", "auprc")},
                                      "hold": {k: round(m2[k], 3) for k in ("prec", "rec", "auprc")}}
        print(f"{'sqsim '+s+' +BAM':<22}{len(lifts):>9}{g:>8}{f:>6} | {m1['prec']:.3f} / {m1['rec']:.3f} / {m1['auprc']:.3f} -> {m2['prec']:.3f} / {m2['rec']:.3f} / {m2['auprc']:.3f}")
    if a.e2e:
        # truth as in end2end/score.py: a FLAIR chain equal to a hidden transcript is genuine
        ex = collections.defaultdict(list)
        for ln in open(f"{a.e2e}/flair.isoforms.gtf"):
            fld = ln.rstrip("\n").split("\t")
            if len(fld) >= 9 and fld[2] == "exon":
                ex[fld[8].split('transcript_id "')[1].split('"')[0]].append((int(fld[3]), int(fld[4])))
        chain = lambda e: tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))
        hidden = set()
        for ln in open(f"{a.e2e}/truth_hidden_chains.tsv"):
            t, ch = ln.rstrip("\n").split("\t")
            if ch: hidden.add(tuple(tuple(int(x) for x in p.split("-")) for p in ch.split(";")))
        rows = table(f"{a.e2e}/e2e_adj.adjudicated.tsv")
        lifts = [r for r in rows if RT_LIFT(r)]
        g = sum(chain(sorted(ex[r["isoform_id"]])) in hidden for r in lifts)
        res["end2end"] = {"rt_lifts": len(lifts), "genuine": g, "false": len(lifts) - g, "known_fragment_rts_sign": sum(KF_RTS(r) for r in rows)}
        print(f"{'end2end':<22}{len(lifts):>9}{g:>8}{len(lifts)-g:>6} | (genuine novels PROBABLE 19/20 -> {19-g}/20 under hold)")
    if a.lrgasp:
        tv = list(csv.DictReader(open(f"{a.lrgasp}/targets/per_target_verdicts.tsv"), delimiter="\t"))
        for p in ("pb", "ont"):
            for s in ("sr", "nosr"):
                rows = table(f"{a.lrgasp}/{p}/adj_{s}.adjudicated.tsv")
                lifts = [r for r in rows if RT_LIFT(r)]
                lifted_ids = {r["isoform_id"] for r in lifts}
                val = sorted({t["target"] for t in tv if t["set"] == "wtc11" and t["platform"] == p and t["short_reads"] == s
                              and t["status"] == "called" and t["validation"] == "supported"
                              and any(i in lifted_ids for i in t["isoforms"].split(","))})
                res[f"lrgasp_{p}_{s}"] = {"rt_lifts": len(lifts), "validated_targets_among_them": val,
                                         "known_fragment_rts_sign": sum(KF_RTS(r) for r in rows)}
                print(f"{'lrgasp '+p+' '+s:<22}{len(lifts):>9}{'':>8}{'':>6} | validated targets among them: {', '.join(val) or '-'}")
    for spec in a.threeway:
        lab, d = spec.split("=", 1)
        filt = {r["isoform"]: r["filter_result"] for r in csv.DictReader(open(f"{d}/armA/sq/rules/rules_RulesFilter_classification.txt"), delimiter="\t")}
        rows = table(f"{d}/armC/wgs.adjudicated.tsv")
        lifts = [r for r in rows if RT_LIFT(r)]
        sq = collections.Counter(filt.get(r["isoform_id"], "?") for r in lifts)
        res[f"refbias_3way_{lab}"] = {"rt_lifts": len(lifts), "sqanti3_filter": dict(sq), "known_fragment_rts_sign": sum(KF_RTS(r) for r in rows),
                                     "novel_probable_total": sum(r["confidence_class"] == "NOVEL_PROBABLE" for r in rows)}
        print(f"{'3-way '+lab:<22}{len(lifts):>9}{'':>8}{'':>6} | SQANTI3 filter: {dict(sq)}; NOVEL_PROBABLE total {res[f'refbias_3way_{lab}']['novel_probable_total']}")
    if a.gm12878:
        rows = table(a.gm12878)
        res["gm12878_ont"] = {"rt_lifts": sum(RT_LIFT(r) for r in rows)}
        print(f"{'gm12878 ONT':<22}{res['gm12878_ont']['rt_lifts']:>9}")
    if a.emit_metrics:
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{a.h2h}/sr_1000000/pigbam.provenance.log")
                    if "\t" in l and not l.startswith("#"))
        env = {"protocol": "rts_sensitivity", "scope": "all_tracked_runs", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "the adjudicated.tsv files of the tracked runs (SQANTI-SIM h2h 1M/none + BAM, end2end, LRGASP WTC11 pb/ont sr/nosr, refbias_3way arm C wgs, gm12878 ONT preset)",
               "command": "benchmark/rts_sensitivity/rts_sensitivity.py --sqsim-work ... --h2h ... --e2e ... --lrgasp ... --threeway ... --gm12878 ... --emit-metrics",
               "metrics": res,
               "notes": "Ruleset 0.0.5 lets a whole-chain / consensus lift stand when SQANTI3's RTS_stage is TRUE (reason ;rt_switch_sign). "
                        "'hold' is the 0.0.4 alternative (UNCONFIRMED). On SQANTI-SIM every such lift is genuine, but the simulator "
                        "produces no RT switching, so this measures only the flag's false-sign side; the real-data columns show how "
                        "many calls the choice moves (all of them SQANTI3-filter Artifacts) and which validated LRGASP isoforms it keeps."}
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
