#!/usr/bin/env python3
"""Real-data negatives from TUSCO: the 46 human universal single-isoform control genes (Pardo-Palacios
et al., Nat Commun 2026; list from ConesaLab/SQANTI3 src/utilities/report_qc/tusco_human.tsv). A
novel isoform (NIC / NNC) a caller reports in one of these genes is a presumed false positive, so
what a filter does with it is a real-data precision check, small as the set is.

Usage: score_tusco.py --tusco data/tusco_human.tsv --run LABEL=CLS,ADJ[,RULES] ... [--emit-metrics ...]
  CLS   SQANTI3 classification of the caller's isoforms
  ADJ   PanIsoGuard <prefix>.adjudicated.tsv on the same isoforms
  RULES SQANTI3 rules-filter classification on the same isoforms (optional)
"""
import argparse, collections, csv, datetime, json, sys

CONFIDENT = {"NOVEL_CONFIRMED", "NOVEL_PROBABLE", "REFERENCE_BIAS"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tusco", required=True); ap.add_argument("--run", action="append", required=True); ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    rows = [l.rstrip("\n").split("\t") for l in open(a.tusco) if l.strip() and not l.startswith("#")]
    ids = {r[0].split(".")[0] for r in rows}; names = {r[2] for r in rows if len(r) > 2}
    in_tusco = lambda gene: any(x.split(".")[0] in ids or x in names for x in gene.split("_"))
    res = {"n_tusco_genes": len(ids), "runs": {}}
    tot = collections.Counter()
    print(f"{'run':<16}{'TUSCO isoforms':>15}{'novel':>7}{'PIG confident':>15}{'PIG artifact':>13}{'SQ3 rules keeps':>16}")
    for spec in a.run:
        label, files = spec.split("=", 1); f = files.split(",")
        cls = {r["isoform"]: r for r in csv.DictReader(open(f[0]), delimiter="\t")}
        adj = {r["isoform_id"]: r for r in csv.DictReader(open(f[1]), delimiter="\t")}
        rules = {r["isoform"]: r["filter_result"] for r in csv.DictReader(open(f[2]), delimiter="\t")} if len(f) > 2 else None
        allt = [i for i, r in cls.items() if in_tusco(r["associated_gene"])]
        novel = [i for i in allt if cls[i]["structural_category"] in ("novel_in_catalog", "novel_not_in_catalog") and i in adj]
        per = [{"isoform": i, "gene": cls[i]["associated_gene"], "category": cls[i]["structural_category"],
                "class": adj[i]["confidence_class"], "reason": adj[i]["verdict_reason"],
                "sqanti3_rules": rules.get(i) if rules else None} for i in novel]
        r = {"isoforms_in_tusco_genes": len(allt), "by_category": dict(collections.Counter(cls[i]["structural_category"] for i in allt)),
             "novel": len(novel), "panisoguard_classes": dict(collections.Counter(p["class"] for p in per)),
             "panisoguard_confident": sum(p["class"] in CONFIDENT for p in per),
             "panisoguard_artifact": sum(p["class"] == "ARTIFACT" for p in per),
             "sqanti3_rules_keeps": sum(p["sqanti3_rules"] == "Isoform" for p in per) if rules else None, "per_isoform": per}
        res["runs"][label] = r
        if not label.endswith("_sr"):   # one run per data set: the no-short-read run (every data set has one)
            tot["novel"] += len(novel); tot["pig"] += r["panisoguard_confident"]; tot["art"] += r["panisoguard_artifact"]
            if rules: tot["sq3"] += r["sqanti3_rules_keeps"]; tot["sq3_n"] += len(novel)
        print(f"{label:<16}{len(allt):>15}{len(novel):>7}{r['panisoguard_confident']:>15}{r['panisoguard_artifact']:>13}{str(r['sqanti3_rules_keeps']):>16}")
    res["total_without_short_reads"] = {"novel": tot["novel"], "panisoguard_confident": tot["pig"], "panisoguard_artifact": tot["art"],
                                        "sqanti3_rules_keeps": tot["sq3"], "sqanti3_rules_n": tot["sq3_n"]}
    print(f"total without short reads (one run per data set): {tot['novel']} presumed-false novel isoforms; PanIsoGuard confident {tot['pig']}, artifact {tot['art']}; SQANTI3 rules keeps {tot['sq3']} / {tot['sq3_n']}")
    if a.emit_metrics:
        env = {"protocol": "tusco_negatives", "scope": "tusco_human_universal_46_genes_real_data", "status": "tracked",
               "tool_version": "0.0.7", "ruleset_version": "builtin-0.0.5",
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "TUSCO human universal single-isoform genes (46; Pardo-Palacios et al. 2026, list from ConesaLab/SQANTI3) "
                                  "applied to IsoQuant isoforms on LRGASP WTC11 PacBio / ONT cDNA (ENCODE) and HPRC R2 HG03516 / HG02717 "
                                  "Kinnex Iso-Seq; PanIsoGuard runs of benchmark/lrgasp_wtc11 (adj_sr / adj_nosr) and benchmark/refbias_3way "
                                  "(arm C, wgs); SQANTI3 6.0.1 rules filter on the same classifications. Public data.",
               "command": "benchmark/tusco_negatives/score_tusco.py --tusco data/tusco_human.tsv --run LABEL=CLS,ADJ,RULES ... --emit-metrics",
               "metrics": res,
               "notes": "A novel isoform in a TUSCO gene is a presumed false positive (TUSCO's own FP definition), so this is the only "
                        "real-data precision check available; the set is small (46 genes), and a TUSCO gene can still express an "
                        "unannotated isoform, so a confident call here is suspicious, not proven wrong."}
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
