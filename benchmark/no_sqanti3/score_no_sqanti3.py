#!/usr/bin/env python3
"""Compare `adjudicate --no-sqanti3` with the SQANTI3-based run on the chr22 five-caller set.

For each caller and short-read setting: the structural categories of the computed rows against
SQANTI3's (confusion over the isoforms both classify), the agreement of all_canonical and of
perc_A_downstream_TTS (within 1 point), the agreement of the confidence class, and the
precision / recall of genuine-novel detection of both modes against the SQANTI-SIM truth
(positive = NOVEL_CONFIRMED / NOVEL_PROBABLE, over multi-exon isoforms the truth labels
genuine_novel or false_novel).

Usage: score_no_sqanti3.py --out DIR --truth truth.truth.tsv --gtf NAME=PATH ... --cls NAME=PATH ...
                           [--tool-version V] [--emit-metrics m.json]
"""
import argparse, collections, csv, datetime, json, sys

CONF = {"NOVEL_CONFIRMED", "NOVEL_PROBABLE"}


def chain_keys(gtf):
    ex, strand, chrom = collections.defaultdict(list), {}, {}
    for line in open(gtf):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        tid = f[8].split('transcript_id "')[1].split('"')[0]
        ex[tid].append((int(f[3]), int(f[4]))); strand[tid] = f[6]; chrom[tid] = f[0]
    keys = {}
    for tid, e in ex.items():
        e = sorted(e)
        intr = [(e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)]
        keys[tid] = f"{chrom[tid]}|{strand[tid]}|" + ",".join(f"{s}-{t}" for s, t in intr) if intr else None
    return keys


def table(path, key):
    with open(path) as f:
        return {r[key]: r for r in csv.DictReader(f, delimiter='\t')}


def na(v):
    return v in ("", "NA", "NaN", "nan", ".")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--truth", required=True)
    ap.add_argument("--gtf", action="append", required=True); ap.add_argument("--cls", action="append", required=True)
    ap.add_argument("--tool-version", default=None); ap.add_argument("--emit-metrics")
    a = ap.parse_args()
    gtfs = dict(x.split("=", 1) for x in a.gtf); clss = dict(x.split("=", 1) for x in a.cls)
    truth = dict(l.rstrip('\n').split('\t') for l in open(a.truth))
    res = {}
    for caller in gtfs:
        keys = chain_keys(gtfs[caller]); sq = table(clss[caller], "isoform")
        for setting in ("nosr", "sr1M"):
            s3 = table(f"{a.out}/{caller}.{setting}.sqanti3.adjudicated.tsv", "isoform_id")
            cp = table(f"{a.out}/{caller}.{setting}.computed.adjudicated.tsv", "isoform_id")
            ids = [i for i in cp if i in s3]
            multi = [i for i in ids if keys.get(i)]          # multi-exon: the only isoforms PanIsoGuard adjudicates
            confusion = collections.Counter((sq[i]["structural_category"], cp[i]["structural_category"]) for i in multi if i in sq)
            cat_agree = sum(v for (x, y), v in confusion.items() if x == y)
            mono = collections.Counter(sq[i]["structural_category"] for i in ids if i in sq and not keys.get(i))
            canon = collections.Counter()
            jl3 = {json.loads(l)["isoform"]: json.loads(l) for l in open(f"{a.out}/{caller}.{setting}.sqanti3.attribution.jsonl")}
            jlc = {json.loads(l)["isoform"]: json.loads(l) for l in open(f"{a.out}/{caller}.{setting}.computed.attribution.jsonl")}
            for i in ids:
                e3, ec = jl3[i]["evidence"], jlc[i]["evidence"]
                if e3["noncanonical"] is not None and ec["noncanonical"] is not None:
                    canon["both_evaluable"] += 1
                    canon["agree"] += e3["noncanonical"] == ec["noncanonical"]
                elif ec["noncanonical"] is None and e3["noncanonical"] is not None:
                    canon["computed_na"] += 1
            class_agree = sum(s3[i]["confidence_class"] == cp[i]["confidence_class"] for i in ids)
            moves = collections.Counter((s3[i]["confidence_class"], cp[i]["confidence_class"]) for i in ids if s3[i]["confidence_class"] != cp[i]["confidence_class"])
            acc = {}
            for mode, tab in (("sqanti3", s3), ("computed", cp)):
                tp = fp = fn = tn = 0
                for i, r in tab.items():
                    k = keys.get(i)
                    if not k: continue
                    lab = truth.get(k, "false_novel")        # a multi-exon chain matching no simulated transcript is false
                    if lab not in ("genuine_novel", "false_novel"): continue
                    pos = r["confidence_class"] in CONF
                    if lab == "genuine_novel": tp += pos; fn += not pos
                    else: fp += pos; tn += not pos
                acc[mode] = dict(TP=tp, FP=fp, FN=fn, TN=tn,
                                 precision=round(tp / (tp + fp), 3) if tp + fp else None,
                                 recall=round(tp / (tp + fn), 3) if tp + fn else None)
            res[f"{caller}.{setting}"] = dict(
                isoforms=len(ids), multi_exon=len(multi), mono_exon_sqanti3_categories=dict(mono),
                category_agreement=round(cat_agree / max(1, sum(confusion.values())), 4),
                category_confusion={f"{x}->{y}": v for (x, y), v in sorted(confusion.items()) if x != y},
                all_canonical=dict(canon), class_agreement=round(class_agree / max(1, len(ids)), 4),
                class_moves={f"{x}->{y}": v for (x, y), v in moves.most_common()}, accuracy=acc)
            print(f"{caller} {setting}: n {len(ids)} ({len(multi)} multi-exon), categories agree {res[f'{caller}.{setting}']['category_agreement']}, "
                  f"classes agree {res[f'{caller}.{setting}']['class_agreement']}, "
                  f"P/R sqanti3 {acc['sqanti3']['precision']}/{acc['sqanti3']['recall']} computed {acc['computed']['precision']}/{acc['computed']['recall']}")
    if a.emit_metrics:
        env = {"protocol": "no_sqanti3", "scope": "gencode_v49_chr22_5caller", "status": "tracked",
               "tool_version": a.tool_version, "ruleset_version": None,
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
               "data_provenance": "SQANTI-SIM GENCODE v49 chr22 (PBSIM3 HiFi, one shared minimap2 alignment) called by FLAIR, "
                                  "IsoQuant, Bambu, ESPRESSO and TALON (the benchmark/multicaller set), each with its SQANTI3 6.0.1 "
                                  "classification; short reads = 1 M wgsim pairs -> STAR (sqanti3_filter_h2h) or none. Public/simulated.",
               "command": "bash benchmark/no_sqanti3/run.sh; benchmark/no_sqanti3/score_no_sqanti3.py --out ... --truth ... --gtf ... --cls ... --emit-metrics ...",
               "source": None,
               "notes": "adjudicate --no-sqanti3 (categories from the chains and the reference GTF, splice motifs and the 3' A content "
                        "from the genome, RT-switch not evaluable) against the SQANTI3-based run on the same isoforms: category "
                        "agreement over the isoforms SQANTI3 classifies, all_canonical agreement where both are evaluable, confidence "
                        "class agreement, and genuine-novel precision / recall of both modes against the truth.",
               "metrics": res}
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics, file=sys.stderr)


if __name__ == "__main__":
    main()
