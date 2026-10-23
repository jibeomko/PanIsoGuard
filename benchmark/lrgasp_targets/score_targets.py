#!/usr/bin/env python3
"""Score the LRGASP validation targets run as would-be isoforms (run.sh): per target set, platform and
short-read setting, what PanIsoGuard and SQANTI3's default rules filter do with the amplicon-supported
chains (positives) and the unsupported ones (negatives), restricted to chains novel against GENCODE v49.

Usage: score_targets.py --lrgasp /mnt/Data/pig_lrgasp [--emit-metrics m.json]
"""
import argparse, collections, csv, datetime, json, os, sys

NOVEL = ("novel_in_catalog", "novel_not_in_catalog")
CONF = ("NOVEL_CONFIRMED", "NOVEL_PROBABLE")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lrgasp", required=True); ap.add_argument("--emit-metrics"); a = ap.parse_args()
    L = a.lrgasp; W = f"{L}/targets_run"
    # transcript_id -> (set, target) from the two GTFs; target -> validation category from the TSVs
    tid = {}
    for setname, g in (("wtc11", "targets.gtf"), ("gencode", "targets_gencode.gtf")):
        for l in open(f"{L}/targets/{g}"):
            f = l.rstrip("\n").split("\t")
            if len(f) < 9: continue
            t = f[8].split('transcript_id "')[1].split('"')[0]; tg = f[8].split('gene_id "')[1].split('"')[0]; tid[t] = (setname, tg)
    cat = {}
    for setname, tsv in (("wtc11", "targets.tsv"), ("gencode", "targets_gencode.tsv")):
        for r in csv.DictReader(open(f"{L}/targets/{tsv}"), delimiter="\t"): cat[(setname, r["target"])] = r["category"]
    sq = {r["isoform"]: r for r in csv.DictReader(open(f"{W}/sq/tg_classification.txt"), delimiter="\t")}
    rules = {r["isoform"]: r["filter_result"] for r in csv.DictReader(open(f"{W}/rules/rules_RulesFilter_classification.txt"), delimiter="\t")}
    res = {}; rows = []
    for p in ("pb", "ont"):
        for s in ("nosr", "sr"):
            adj = {r["isoform_id"]: r for r in csv.DictReader(open(f"{W}/adj_{p}_{s}.adjudicated.tsv"), delimiter="\t")}
            full = {}
            for l in open(f"{W}/adj_{p}_{s}.attribution.jsonl"):
                d = json.loads(l); full[d["isoform"]] = d["evidence"].get("chain_n_full_reads")
            summ = {}
            for val in ("supported", "likely", "unsupported"):
                sel = [t for t in adj if tid.get(t) and cat.get(tid[t]) == val and sq[t]["structural_category"] in NOVEL and adj[t]["confidence_class"] != "KNOWN"]
                summ[val] = {"novel_targets": len(sel), "class_counts": dict(collections.Counter(adj[t]["confidence_class"] for t in sel)),
                             "panisoguard_confident": sum(adj[t]["confidence_class"] in CONF for t in sel), "panisoguard_artifact": sum(adj[t]["confidence_class"] == "ARTIFACT" for t in sel),
                             "sqanti3_rules_keeps": sum(rules.get(t) == "Isoform" for t in sel), "whole_chain_reads_ge3": sum((full.get(t) or 0) >= 3 for t in sel)}
                for t in sel:
                    rows.append({"platform": p, "short_reads": s, "set": tid[t][0], "target": tid[t][1], "validation": val, "sqanti3_category": sq[t]["structural_category"],
                                 "class": adj[t]["confidence_class"], "reason": adj[t]["verdict_reason"], "whole_chain_reads": full.get(t), "sqanti3_rules": rules.get(t)})
            res[f"{p}_{s}"] = summ
            print(f"{p} {s}: " + "; ".join(f"{v}: n {d['novel_targets']}, PIG confident {d['panisoguard_confident']}, artifact {d['panisoguard_artifact']}, >=3 reads {d['whole_chain_reads_ge3']}, SQ3 keeps {d['sqanti3_rules_keeps']}" for v, d in summ.items()))
    with open(f"{W}/per_target.tsv", "w") as o:
        w = csv.DictWriter(o, delimiter="\t", fieldnames=list(rows[0].keys())); w.writeheader(); [w.writerow(r) for r in rows]
    if a.emit_metrics:
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{W}/adj_pb_nosr.provenance.log") if "\t" in l and not l.startswith("#"))
        env = {"protocol": "lrgasp_targets", "scope": "lrgasp_wtc11_validation_targets_as_would_be_isoforms", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "REAL public LRGASP WTC11 (ENCODE PacBio cDNA ENCSR507JOF, ONT cDNA ENCSR539ZXJ, Illumina ENCSR673UKZ) and LRGASP "
                                  "Supplementary Data 15 amplicon validation; every multi-exon target chain of both sheets written as a GTF "
                                  "(benchmark/lrgasp_wtc11/targets.py), SQANTI3 6.0.1 vs GENCODE v49, PanIsoGuard with the platform's BAM, "
                                  "with and without the Illumina SJ.tab. No private data.",
               "command": "bash benchmark/lrgasp_targets/run.sh; benchmark/lrgasp_targets/score_targets.py --lrgasp /mnt/Data/pig_lrgasp --emit-metrics ...",
               "metrics": res,
               "notes": "Targets novel against GENCODE v49 (SQANTI3 NIC/NNC, not KNOWN) by amplicon validation: supported (positives), unsupported "
                        "(negatives), likely. Independent of any caller's reconstruction: the chains are adjudicated as would-be isoforms, so the "
                        "negatives test whether whole-chain reads or short reads would make PanIsoGuard confirm a chain the amplicons did not support."}
        os.makedirs(os.path.dirname(a.emit_metrics), exist_ok=True)
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n"); print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
