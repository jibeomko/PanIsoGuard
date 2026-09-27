#!/usr/bin/env python3
"""Head-to-head: PanIsoGuard vs the SQANTI3 rules / ML filter on SQANTI-SIM truth.

Same truth labels as ../sqanti_sim/score.py (exact intron-chain match):
genuine_novel = chain of a simulated transcript deleted from the reduced annotation,
false_novel = multi-exon chain matching no simulated transcript. Task: genuine vs
false among multi-exon isoforms, reported overall (non-FSM) and per SQANTI category.

Methods are NAME=KIND:PATH, KIND one of
  pig     PanIsoGuard adjudicated.tsv   positive = HIGH/MEDIUM_CONF_NOVEL, score = class rank
  piglen  PanIsoGuard adjudicated.tsv   positive = anything but ARTIFACT (keep-filter semantics)
  rules   SQANTI3 *_RulesFilter_classification.txt   positive = filter_result == Isoform
  ml      SQANTI3 *_ML_classification.txt            positive = filter_result == Isoform, score = POS_MLprob
  mincov  SQANTI3 classification run with -c         positive = min_cov >= 3 (every junction has >= 3 unique
                                                     SR reads, the same gate as sj_min_uniq_reads)

Usage: score_h2h.py <flair.isoforms.gtf> <truth.tsv> <sqanti_classification.txt> NAME=KIND:PATH ...
"""
import sys
from collections import defaultdict
from sklearn.metrics import average_precision_score

flair_gtf, truth_tsv, cls_tsv = sys.argv[1:4]
methods = [a.split('=', 1) for a in sys.argv[4:]]

def tid_of(a):
    return a.split('transcript_id "')[1].split('"')[0]

ex = defaultdict(list); strand = {}; chrom = {}
for line in open(flair_gtf):
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    t = tid_of(f[8]); ex[t].append((int(f[3]), int(f[4]))); strand[t] = f[6]; chrom[t] = f[0]

def chain_key(t):
    e = sorted(ex[t])
    intr = [(e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)]
    return f"{chrom[t]}|{strand[t]}|" + ",".join(f"{s}-{ee}" for s, ee in intr) if intr else None

truth = dict(l.rstrip('\n').split('\t') for l in open(truth_tsv))

def table(path):
    it = (l.rstrip('\n').split('\t') for l in open(path))
    hdr = next(it)
    return [dict(zip(hdr, f)) for f in it]

cat = {r['isoform']: r['structural_category'] for r in table(cls_tsv)}
label = {}
for t in cat:
    ck = chain_key(t) if t in ex else None
    if ck is not None:
        label[t] = truth.get(ck, 'false_novel')
label = {t: l for t, l in label.items() if l in ('genuine_novel', 'false_novel')}

RANK = {"HIGH_CONF_NOVEL": 1.0, "MEDIUM_CONF_NOVEL": 0.75, "AMBIGUOUS": 0.5, "LOW_CONF_PARTIAL": 0.25,
        "HIGH_CONF_KNOWN": 0.1, "ARTIFACT": 0.0, "PAN_REF_RESCUED_FALSE_NOVEL": 0.0}

def load(kind, path):
    """-> {isoform: (predicted_positive, score)}"""
    out = {}
    if kind in ('pig', 'piglen'):
        for r in table(path):
            c = r['confidence_class']
            pos = c in ("HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL") if kind == 'pig' else c != "ARTIFACT"
            out[r['isoform_id']] = (pos, RANK.get(c, 0.0))
    elif kind in ('rules', 'ml'):
        for r in table(path):
            pos = r['filter_result'] == 'Isoform'
            ml = kind == 'ml' and r['POS_MLprob'] != 'NA'  # NA only on mono-exonic (unlabelled)
            out[r['isoform']] = (pos, float(r['POS_MLprob']) if ml else float(pos))
    elif kind == 'mincov':
        for r in table(path):
            pos = r['min_cov'] not in ('NA', '') and float(r['min_cov']) >= 3
            out[r['isoform']] = (pos, float(pos))
    else:
        sys.exit(f"unknown kind {kind}")
    return out

def metrics(pred, ids):
    y = [label[t] == 'genuine_novel' for t in ids]
    p = [pred[t][0] for t in ids]; s = [pred[t][1] for t in ids]
    tp = sum(a and b for a, b in zip(y, p)); fp = sum(b and not a for a, b in zip(y, p))
    fn = sum(a and not b for a, b in zip(y, p)); tn = sum(not a and not b for a, b in zip(y, p))
    div = lambda a, b: a / b if b else float('nan')
    prec, rec = div(tp, tp + fp), div(tp, tp + fn)
    auprc = average_precision_score(y, s) if 0 < sum(y) < len(y) else float('nan')
    return dict(n=len(ids), pos=sum(y), tp=tp, fp=fp, prec=prec, rec=rec, spec=div(tn, tn + fp),
                f1=div(2 * prec * rec, prec + rec), auprc=auprc)

preds = {name: load(*spec.split(':', 1)) for name, spec in methods}
for name, pred in preds.items():
    missing = set(label) - set(pred)
    if missing:
        sys.exit(f"{name}: {len(missing)} labelled isoforms missing from output")

strata = [("non-FSM (novel candidates)", [t for t in label if cat[t] != 'full-splice_match'])]
for c in sorted({cat[t] for t in label}):
    strata.append((c, [t for t in label if cat[t] == c]))

cols = ["prec", "rec", "spec", "f1", "auprc"]
for title, ids in strata:
    m0 = metrics(next(iter(preds.values())), ids)
    print(f"\n=== {title}: n={m0['n']} genuine={m0['pos']} false={m0['n'] - m0['pos']} ===")
    print("method".ljust(16) + "".join(c.rjust(8) for c in ["TP", "FP"] + cols))
    for name, pred in preds.items():
        m = metrics(pred, ids)
        print(name.ljust(16) + f"{m['tp']:8d}{m['fp']:8d}" + "".join(f"{m[c]:8.3f}" for c in cols))
