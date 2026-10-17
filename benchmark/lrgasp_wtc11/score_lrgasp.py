#!/usr/bin/env python3
"""Score PanIsoGuard's verdicts on the LRGASP WTC11 experimentally validated isoforms.

For every validation target (an intron chain), per platform (pb = PacBio cDNA, ont = ONT cDNA)
and short-read setting (sr = with the Illumina SJ.tab, nosr = without): did our caller
(IsoQuant) report the same intron chain, and if so what class / verdict / reason did
PanIsoGuard give? Targets are 'supported', 'likely' or 'unsupported' by amplicon sequencing
(Supplementary Data 15). Also: how many of our isoforms contain the target chain
contiguously (a longer model that uses it).

Usage: score_lrgasp.py --work W [--emit-metrics m.json]
"""
import argparse, collections, csv, datetime, json, os, sys


def gtf_chains(path):
    ex = collections.defaultdict(list); chrom = {}; strand = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or 'transcript_id "' not in f[8]: continue
        t = f[8].split('transcript_id "')[1].split('"')[0]
        ex[t].append((int(f[3]), int(f[4]))); chrom[t] = f[0]; strand[t] = f[6]
    keys = {}
    for t, e in ex.items():
        e = sorted(e)
        if len(e) < 2: continue
        keys[t] = (chrom[t], strand[t], tuple((e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)))
    return keys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--emit-metrics")
    a = ap.parse_args(); W = a.work
    targets = [t for t in csv.DictReader(open(f"{W}/targets/targets.tsv"), delimiter="\t") if t["chain_key"] not in ("NA", "monoexonic")]
    tkey = {}
    for t in targets:
        chrom, strand, rest = t["chain_key"].split("|")
        tkey[t["target"]] = (chrom, strand, tuple(tuple(int(x) for x in p.split("-")) for p in rest.split(",")))
    res = {}
    rows = []
    for p in ("pb", "ont"):
        ours = gtf_chains(f"{W}/{p}/isoquant/iq/iq.transcript_models.gtf")
        exact = collections.defaultdict(list); contain = collections.defaultdict(list)
        by_chrom = collections.defaultdict(list)
        for iso, k in ours.items(): by_chrom[(k[0], k[1])].append((iso, k[2]))
        for tg, (chrom, strand, chain) in tkey.items():
            n = len(chain)
            for iso, c in by_chrom[(chrom, strand)]:
                if c == chain: exact[tg].append(iso)
                elif len(c) > n and any(c[i:i + n] == chain for i in range(len(c) - n + 1)): contain[tg].append(iso)
        for setting in ("sr", "nosr"):
            adj = {r["isoform_id"]: r for r in csv.DictReader(open(f"{W}/{p}/adj_{setting}.adjudicated.tsv"), delimiter="\t")}
            tab = collections.Counter(); detail = []
            for t in targets:
                tg = t["target"]; isos = exact.get(tg, [])
                if not isos:
                    status = "not_called"; cls = reason = cat = ""
                else:
                    r = adj[isos[0]]; status = "called"; cls = r["confidence_class"]; reason = r["verdict_reason"]; cat = r["structural_category"]
                tab[(t["category"], t["type"], status, cls)] += 1
                detail.append({"platform": p, "short_reads": setting, "target": tg, "transcript": t["transcript"], "lrgasp_type": t["type"],
                               "validation": t["category"], "detail": t["detail"], "status": status, "isoforms": isos[:3],
                               "n_containing": len(contain.get(tg, [])), "sqanti3_category": cat, "class": cls, "reason": reason,
                               "pbPlatformCount": t["pbPlatformCount"], "ontPlatformCount": t["ontPlatformCount"], "mediaCpm": t["mediaCpm"]})
            rows += detail
            summ = {}
            for val in ("supported", "likely", "unsupported"):
                sel = [d for d in detail if d["validation"] == val]
                called = [d for d in sel if d["status"] == "called"]
                summ[val] = {"n": len(sel), "called_exact": len(called),
                             "class_counts": dict(collections.Counter(d["class"] for d in called)),
                             "confident": sum(d["class"] in ("NOVEL_CONFIRMED", "NOVEL_PROBABLE", "KNOWN") for d in called),
                             "artifact": sum(d["class"] == "ARTIFACT" for d in called),
                             "by_type": {ty: dict(collections.Counter(f'{d["status"]}:{d["class"] or "-"}' for d in sel if d["lrgasp_type"] == ty)) for ty in ("FSM", "NIC", "NNC")}}
            res[f"{p}_{setting}"] = summ
            print(f"== {p} {setting}: " + "; ".join(f"{v}: {s['called_exact']}/{s['n']} called, confident {s['confident']}, artifact {s['artifact']}, {s['class_counts']}" for v, s in summ.items()))
    with open(f"{W}/targets/per_target_verdicts.tsv", "w") as o:
        w = csv.DictWriter(o, delimiter="\t", fieldnames=list(rows[0].keys())); w.writeheader()
        for d in rows: d["isoforms"] = ",".join(d["isoforms"]); w.writerow(d)
    if a.emit_metrics:
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{W}/pb/adj_sr.provenance.log") if "\t" in l and not l.startswith("#"))
        env = {"protocol": "lrgasp_wtc11", "scope": "lrgasp_wtc11_validated_isoforms_pacbio_ont_illumina", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "REAL public LRGASP WTC11 (ENCODE): PacBio cDNA ENCSR507JOF (3 reps), ONT cDNA ENCSR539ZXJ (rep 1), "
                                  "Illumina ENCSR673UKZ (rep 1) -> minimap2 / STAR 2-pass, IsoQuant 3.13, SQANTI3 6.0.1 vs GENCODE v49; "
                                  "truth = LRGASP Supplementary Data 15 (amplicon-sequencing validation of 89 WTC11 targets) joined to the "
                                  "WTC11 consolidated models BED (cgl.gi.ucsc.edu/data/LRGASP). No private data.",
               "command": "bash benchmark/lrgasp_wtc11/run.sh all; benchmark/lrgasp_wtc11/score_lrgasp.py --work W --emit-metrics ...",
               "metrics": res,
               "notes": "Per platform and short-read setting: of the LRGASP WTC11 validation targets whose intron chain our IsoQuant "
                        "run reported exactly, the PanIsoGuard class distribution for supported / likely / unsupported targets. "
                        "The validation set is mostly positives (79 supported, 3 likely, 7 unsupported), so this is primarily a "
                        "check that experimentally supported isoforms are not called ARTIFACT and are kept, on real data with "
                        "both platforms; it is not a precision estimate."}
        os.makedirs(os.path.dirname(a.emit_metrics), exist_ok=True)
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
