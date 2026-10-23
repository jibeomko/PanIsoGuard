#!/usr/bin/env python3
"""What a real run returns: the verdict-class distribution of PanIsoGuard's adjudicated table.

Counts `confidence_class` over all isoforms and over the novel ones (SQANTI3 NIC + NNC) in the
adjudicated tables other protocols already wrote (lrgasp_wtc11, refbias_3way), with the share of
probable calls that carry the RT-switch annotation and the reasons behind the undecided and held
calls. No new runs: the inputs are those protocols' outputs, so the tool and ruleset versions are
read from their provenance logs.

Usage: score_distribution.py --run LABEL=ADJUDICATED_TSV ... [--emit-metrics m.json]
"""
import argparse, collections, csv, datetime, json, os, sys

CLASSES = ["KNOWN", "KNOWN_FRAGMENT", "NOVEL_CONFIRMED", "NOVEL_PROBABLE", "REFERENCE_BIAS", "UNCONFIRMED", "ARTIFACT", "UNRESOLVED"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", action="append", required=True, help="LABEL=path/to/x.adjudicated.tsv"); ap.add_argument("--emit-metrics")
    a = ap.parse_args(); res = {}; prov = {}
    for spec in a.run:
        label, path = spec.split("=", 1)
        rows = list(csv.DictReader(open(path), delimiter="\t"))
        novel = [r for r in rows if r["structural_category"] in ("novel_in_catalog", "novel_not_in_catalog")]
        c_all = collections.Counter(r["confidence_class"] for r in rows); c_nov = collections.Counter(r["confidence_class"] for r in novel)
        res[label] = {"n_isoforms": len(rows), "all": {k: c_all.get(k, 0) for k in CLASSES}, "novel_nic_nnc": len(novel),
                      "novel": {k: c_nov.get(k, 0) for k in CLASSES},
                      "novel_probable_with_rts_sign": sum("rt_switch_sign" in r["verdict_reason"] for r in novel if r["confidence_class"] == "NOVEL_PROBABLE"),
                      "novel_unconfirmed_reasons": dict(collections.Counter(r["verdict_reason"] for r in novel if r["confidence_class"] == "UNCONFIRMED")),
                      "novel_unresolved_reasons": dict(collections.Counter(r["verdict_reason"] for r in novel if r["confidence_class"] == "UNRESOLVED"))}
        plog = path.replace(".adjudicated.tsv", ".provenance.log")
        if os.path.exists(plog):
            d = dict(l.rstrip("\n").split("\t", 1) for l in open(plog) if "\t" in l and not l.startswith("#"))
            prov.setdefault("tool_version", d.get("tool_version")); prov.setdefault("ruleset_version", d.get("ruleset_version"))
            assert d.get("ruleset_version") == prov["ruleset_version"], f"{label}: ruleset {d.get('ruleset_version')} != {prov['ruleset_version']}"
        print(f"{label}: {len(rows)} isoforms, {len(novel)} novel; " + ", ".join(f"{k} {c_nov.get(k, 0)}" for k in CLASSES if c_nov.get(k)))
    if a.emit_metrics:
        env = {"protocol": "verdict_distribution", "scope": "real_data_verdict_distribution_lrgasp_wtc11_hprc_kinnex", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "REAL public data, adjudicated tables written by benchmark/lrgasp_wtc11 (LRGASP WTC11 PacBio and ONT cDNA, "
                                  "IsoQuant 3.13 vs GENCODE v49, with and without the Illumina SJ.tab) and benchmark/refbias_3way (HPRC R2 "
                                  "HG03516 and HG02717 Kinnex Iso-Seq, IsoQuant vs GENCODE v49, personal haplotypes from the HPRC assemblies, "
                                  "no short reads). No private data.",
               "command": "benchmark/verdict_distribution/score_distribution.py --run LABEL=x.adjudicated.tsv ... --emit-metrics benchmark/results/verdict_distribution/metrics.json",
               "metrics": res,
               "notes": "Per run: confidence_class counts over all isoforms and over SQANTI3 NIC + NNC isoforms; how many NOVEL_PROBABLE calls carry "
                        "the rt_switch_sign annotation; the verdict_reason behind UNCONFIRMED (held on a sign) and UNRESOLVED (undecided) novel calls. "
                        "Descriptive: what a user sees after one run, not an accuracy claim."}
        os.makedirs(os.path.dirname(a.emit_metrics), exist_ok=True)
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
