#!/usr/bin/env python3
"""Score PanIsoGuard's verdicts on the LRGASP WTC11 experimentally validated isoforms.

Two target sets (Supplementary Data 15): `wtc11` (NIC/NNC/FSM models of the WTC11
consolidated set, 89 targets) and `gencode` (GENCODE-reviewed known/novel/rejected models,
79 targets). For every target (an intron chain), per platform (pb = PacBio cDNA, ont = ONT
cDNA) and short-read setting (sr = with the Illumina SJ.tab, nosr = without): did our caller
(IsoQuant) report the same intron chain, and if so what class / verdict / reason did
PanIsoGuard give? Targets are 'supported', 'likely' or 'unsupported' by amplicon sequencing.
Also counted: longer models containing the target chain, and SQANTI3's RTS_stage sign on
validated isoforms (a false sign when the isoform is experimentally supported).

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


def load_targets(path, setname):
    out = []
    for t in csv.DictReader(open(path), delimiter="\t"):
        if t["chain_key"] in ("NA", "monoexonic"): continue
        chrom, strand, rest = t["chain_key"].split("|")
        t["set"] = setname
        t["_key"] = (chrom, strand, tuple(tuple(int(x) for x in p.split("-")) for p in rest.split(",")))
        out.append(t)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--emit-metrics")
    a = ap.parse_args(); W = a.work
    sets = [("wtc11", f"{W}/targets/targets.tsv")]
    if os.path.exists(f"{W}/targets/targets_gencode.tsv"): sets.append(("gencode", f"{W}/targets/targets_gencode.tsv"))
    ours_cache = {}; adj_cache = {}
    res = {}; rows = []
    for setname, tpath in sets:
        targets = load_targets(tpath, setname)
        for p in ("pb", "ont"):
            if p not in ours_cache:
                ours = gtf_chains(f"{W}/{p}/isoquant/iq/iq.transcript_models.gtf")
                by_cs = collections.defaultdict(list)
                for iso, k in ours.items(): by_cs[(k[0], k[1])].append((iso, k[2]))
                ours_cache[p] = by_cs
            by_cs = ours_cache[p]
            exact = {}; contain = {}
            for t in targets:
                chrom, strand, chain = t["_key"]; n = len(chain); ex_l = []; co_l = []
                for iso, c in by_cs[(chrom, strand)]:
                    if c == chain: ex_l.append(iso)
                    elif len(c) > n and any(c[i:i + n] == chain for i in range(len(c) - n + 1)): co_l.append(iso)
                exact[t["target"]] = ex_l; contain[t["target"]] = co_l
            for setting in ("sr", "nosr"):
                key = (p, setting)
                if key not in adj_cache:
                    adj = {r["isoform_id"]: r for r in csv.DictReader(open(f"{W}/{p}/adj_{setting}.adjudicated.tsv"), delimiter="\t")}
                    rts = {}
                    for l in open(f"{W}/{p}/adj_{setting}.attribution.jsonl"):
                        d = json.loads(l); rts[d["isoform"]] = d["evidence"].get("rts_stage")
                    adj_cache[key] = (adj, rts)
                adj, rts = adj_cache[key]
                detail = []
                for t in targets:
                    tg = t["target"]; isos = exact[tg]
                    if not isos:
                        status, cls, reason, cat, rs = "not_called", "", "", "", ""
                    else:
                        r = adj[isos[0]]; status = "called"; cls = r["confidence_class"]; reason = r["verdict_reason"]
                        cat = r["structural_category"]; rs = rts.get(isos[0])
                    detail.append({"set": setname, "platform": p, "short_reads": setting, "target": tg, "transcript": t["transcript"],
                                   "lrgasp_type": t["type"], "validation": t["category"], "detail": t.get("detail", ""), "status": status,
                                   "isoforms": isos[:3], "n_containing": len(contain[tg]), "sqanti3_category": cat, "class": cls,
                                   "reason": reason, "rts_stage": rs, "pbPlatformCount": t.get("pbPlatformCount", ""),
                                   "ontPlatformCount": t.get("ontPlatformCount", ""), "mediaCpm": t.get("mediaCpm", "")})
                rows += detail
                summ = {}
                for val in ("supported", "likely", "unsupported"):
                    sel = [d for d in detail if d["validation"] == val]
                    called = [d for d in sel if d["status"] == "called"]
                    summ[val] = {"n": len(sel), "called_exact": len(called), "contained_only": sum(1 for d in sel if d["status"] != "called" and d["n_containing"] > 0),
                                 "class_counts": dict(collections.Counter(d["class"] for d in called)),
                                 "confident": sum(d["class"] in ("NOVEL_CONFIRMED", "NOVEL_PROBABLE", "KNOWN") for d in called),
                                 "artifact": sum(d["class"] == "ARTIFACT" for d in called),
                                 "by_type": {ty: dict(collections.Counter(f'{d["status"]}:{d["class"] or "-"}' for d in sel if d["lrgasp_type"] == ty))
                                             for ty in sorted({d["lrgasp_type"] for d in sel})}}
                called_sup = [d for d in detail if d["validation"] == "supported" and d["status"] == "called"]
                summ["rts_sign_on_validated_called"] = {"n": len(called_sup), "rts_stage_true": sum(d["rts_stage"] is True for d in called_sup)}
                res[f"{setname}_{p}_{setting}"] = summ
                print(f"== {setname} {p} {setting}: " + "; ".join(
                    f"{v}: {s['called_exact']}/{s['n']} called, confident {s['confident']}, artifact {s['artifact']}, {s['class_counts']}"
                    for v, s in summ.items() if v in ("supported", "likely", "unsupported")) +
                    f"; RTS on validated called {summ['rts_sign_on_validated_called']['rts_stage_true']}/{summ['rts_sign_on_validated_called']['n']}")
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
                                  "truth = LRGASP Supplementary Data 15 (amplicon-sequencing validation): sheet wtc11-support (89 WTC11 "
                                  "consolidated-model targets, chains from WTC11_consolidated.bed.gz) and sheet gencode-support (79 "
                                  "GENCODE-reviewed known/novel/rejected models, chains from LRGASP_manual_annotation.human.all.bed.gz + "
                                  "human_GENCODE_tmerge_transcripts.bed.gz; 76 resolvable). No private data.",
               "command": "bash benchmark/lrgasp_wtc11/run.sh all; benchmark/lrgasp_wtc11/score_lrgasp.py --work W --emit-metrics ...",
               "metrics": res,
               "notes": "Per target set, platform and short-read setting: of the targets whose intron chain our IsoQuant run reported "
                        "exactly, the PanIsoGuard class distribution for supported / likely / unsupported targets, plus how many "
                        "experimentally supported isoforms carry SQANTI3's RTS_stage sign. Mostly positives (wtc11: 79 supported / "
                        "7 unsupported; gencode: 48 supported / 23 unsupported), so a recall-side check first; the negatives the caller "
                        "did reconstruct are the only real-data precision evidence."}
        os.makedirs(os.path.dirname(a.emit_metrics), exist_ok=True)
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
