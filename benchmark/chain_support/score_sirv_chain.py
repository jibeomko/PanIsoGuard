#!/usr/bin/env python3
"""Whole-chain read support on the SIRV-Set4 spike-in truth (second truth set, PBSIM3 reads on
the Lexogen SIRV-Set4 reference, FLAIR + IsoQuant from ../multicaller's SIRV run).

(1) Per caller, full-chain read counts for genuine (hidden-chain truth) vs false novel chains;
(2) the engine ON vs OFF on FLAIR, with and without the simulated short-read SJ.tab, scored
as in score_chain.py. Truth matching is strand-agnostic on the exact intron chain (the SIRV
truth carries no strand), introns as 0-based half-open [exon_end, next_exon_start).

Usage: score_sirv_chain.py --work /path/to/pig_sirv/work --mc /path/to/pig_sirv_mc --out OUTDIR --pig panisoguard --cfg rules.toml [--emit-metrics m.json]
"""
import argparse, collections, json, os, statistics, subprocess, sys
import pysam
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "sqanti3_filter_h2h"))
from score_h2h import load, metrics, table  # noqa: E402
from score_chain import risk_coverage, DECISIVE, POSITIVE  # noqa: E402


def gtf_chains(path):
    ex = collections.defaultdict(list); chrom = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or 'transcript_id "' not in f[8]: continue
        t = f[8].split('transcript_id "')[1].split('"')[0]; ex[t].append((int(f[3]), int(f[4]))); chrom[t] = f[0]
    return {t: (chrom[t], tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)))
            for t, e in ((t, sorted(e)) for t, e in ex.items()) if len(e) > 1}


def read_introns(r):
    out = []; pos = r.reference_start
    for op, n in r.cigartuples:
        if op == 3: out.append((pos, pos + n)); pos += n
        elif op in (0, 2, 7, 8): pos += n
    return tuple(out)


def contiguous(sub, full):
    k = len(sub); return any(full[i:i + k] == sub for i in range(len(full) - k + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--mc", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--pig", required=True); ap.add_argument("--cfg", required=True); ap.add_argument("--emit-metrics")
    a = ap.parse_args(); W = a.work; os.makedirs(a.out, exist_ok=True)
    truth = {}
    for l in open(f"{W}/truth_hidden_chains.tsv"):
        tid, contig, chain = l.rstrip("\n").split("\t")
        truth[(contig, tuple(tuple(int(x) for x in p.split("-")) for p in chain.split(";")))] = tid
    known = set()
    for t, (c, ch) in gtf_chains(f"{W}/reduced.gtf").items(): known.add((c, ch))
    bam = pysam.AlignmentFile(f"{W}/aln.bam")
    res = {"n_hidden_truth": len(truth), "per_caller": {}}
    print("caller    novel genuine false | genuine full>=3  false full>=3 | medians")
    callers = {"flair": f"{W}/flair.isoforms.gtf", "isoquant": f"{a.mc}/isoquant_out/iq/iq.transcript_models.gtf"}
    labels = {}
    for name, gtf in callers.items():
        rows = []; lab = {}
        for t, (c, ch) in gtf_chains(gtf).items():
            if (c, ch) in known: continue
            l = "genuine_novel" if (c, ch) in truth else "false_novel"; lab[t] = l
            n_full = sum(1 for r in bam.fetch(c, max(0, ch[0][0] - 1), ch[0][0] + 1)
                         if not (r.is_secondary or r.is_supplementary or r.is_unmapped) and contiguous(ch, read_introns(r)))
            rows.append((l, n_full))
        labels[name] = lab
        g = [n for l, n in rows if l == "genuine_novel"]; f = [n for l, n in rows if l == "false_novel"]
        res["per_caller"][name] = {"novel": len(rows), "genuine": len(g), "false": len(f),
                                   "genuine_full_ge3": sum(n >= 3 for n in g), "false_full_ge3": sum(n >= 3 for n in f),
                                   "genuine_median": statistics.median(g) if g else None, "false_median": statistics.median(f) if f else None}
        print(f"{name:9s} {len(rows):5d} {len(g):7d} {len(f):5d} | {sum(n>=3 for n in g):7d}/{len(g):<6d} {sum(n>=3 for n in f):6d}/{len(f):<5d} | {statistics.median(g) if g else '-'} / {statistics.median(f) if f else '-'}")
    # engine ON / OFF on FLAIR (the caller with enough novel calls), +- short reads
    cfg_off = f"{a.out}/cfg_off.toml"
    open(cfg_off, "w").write(open(a.cfg).read().replace("min_full_reads = 3", "min_full_reads = 1000000"))
    common = ["--classification", f"{W}/sq/s_classification.txt", "--isoforms-gtf", f"{W}/flair.isoforms.gtf", "--ref-gtf", f"{W}/reduced.gtf",
              "--bam", f"{W}/aln.bam", "--reference", f"{W}/sirv.fa"]
    lab = labels["flair"]; ids = sorted(lab)
    res["engine_flair"] = {}
    for setting, extra in (("sr", ["--sj-tab", f"{W}/real.SJ.tab"]), ("nosr", [])):
        out = {}
        cls = {}
        for state, cfg in (("off", cfg_off), ("on", a.cfg)):
            pfx = f"{a.out}/flair.{setting}.{state}"
            subprocess.run([a.pig, "adjudicate", "--config", cfg, *common, *extra, "--out-prefix", pfx], check=True, capture_output=True)
            p = pfx + ".adjudicated.tsv"; pred = load("pig", p)
            cls[state] = {r["isoform_id"]: r["confidence_class"] for r in table(p)}
            m = metrics(pred, ids, lab)
            out[state] = {k: (None if m[k] != m[k] else round(m[k], 3)) for k in ("prec", "rec", "f1", "auprc")}
            out[state].update(TP=m["tp"], FP=m["fp"], **risk_coverage(cls[state], ids, lab))
        trans = collections.Counter((lab[t], cls["off"][t], cls["on"][t]) for t in ids if cls["off"][t] != cls["on"][t])
        out["transitions_off_to_on"] = {f"{l}:{x}->{y}": n for (l, x, y), n in sorted(trans.items())}
        out["lifted_to_confident"] = {"genuine": sum(n for (l, x, y), n in trans.items() if l == "genuine_novel" and y in POSITIVE),
                                      "false": sum(n for (l, x, y), n in trans.items() if l == "false_novel" and y in POSITIVE)}
        res["engine_flair"][setting] = out
        print(f"flair {setting}: prec {out['off']['prec']}->{out['on']['prec']} rec {out['off']['rec']}->{out['on']['rec']} lifted {out['lifted_to_confident']} risk@cov {out['off']['selective_risk']}@{out['off']['coverage']} -> {out['on']['selective_risk']}@{out['on']['coverage']}")
    # SQANTI3 categories of the genuine FLAIR chains: on SIRV most land outside NIC/NNC (antisense:
    # the SIRV annotation's strand vs FLAIR's), so the engine-level ON/OFF cannot lift them.
    cat = {r["isoform"]: r["structural_category"] for r in table(f"{W}/sq/s_classification.txt")}
    res["flair_genuine_sqanti3_categories"] = dict(collections.Counter(cat.get(t, "NA") for t, l in labels["flair"].items() if l == "genuine_novel"))
    if a.emit_metrics:
        env = json.load(open(a.emit_metrics)); env["metrics"]["sirv_set4"] = res
        marker = " SIRV-Set4 (second truth set"
        env["notes"] = env["notes"].split(marker)[0]
        env["notes"] += (marker + ", PBSIM3 reads): read level replicates, " +
                         "; ".join(f"{c}: genuine {v['genuine_full_ge3']}/{v['genuine']} vs false {v['false_full_ge3']}/{v['false']} with >= 3 full-chain reads"
                                   for c, v in res["per_caller"].items()) +
                         ". Engine level not evaluable there: SQANTI3 classifies most FLAIR SIRV models antisense/intergenic "
                         f"(genuine chains by category: {res['flair_genuine_sqanti3_categories']}), which PanIsoGuard does not adjudicate, "
                         "and the adjudicated NNC carry mapping/noncanonical signs from the simulated reads, so ON == OFF.")
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("updated", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
