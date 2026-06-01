#!/usr/bin/env python3
"""Score PanIsoGuard verdicts against SQANTI-SIM truth.

Truth per FLAIR isoform (matched by exact intron chain):
  genuine_novel  - chain == a simulated transcript that was DELETED from the
                   reduced annotation (a real isoform the caller must rediscover)
  known          - chain == a simulated transcript still in the reduced annotation
  false_novel    - multi-exon chain matching NO simulated transcript
                   (a FLAIR/alignment artifact)
  monoexonic     - no intron chain (excluded from the novel-detection task)

Reports: confusion (truth x confidence class), genuine-novel detection
precision/recall at the decisive threshold (HIGH/MEDIUM_CONF_NOVEL = predicted
genuine), per-SQANTI-category recall, and AUPRC for genuine-vs-false ranking.

Usage: score.py <flair.isoforms.gtf> <truth.tsv> <adjudicated.tsv>
"""
import sys
from collections import defaultdict, Counter

flair_gtf, truth_tsv, adj_tsv = sys.argv[1:4]

def tid_of(a):
    return a.split('transcript_id "')[1].split('"')[0]

def chain_key(exlist, chrom, strand):
    e = sorted(exlist)
    intr = [(e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)]
    if not intr:
        return None
    return f"{chrom}|{strand}|" + ",".join(f"{s}-{ee}" for s, ee in intr)

# FLAIR isoform -> chain key
ex = defaultdict(list); strand = {}; chrom = {}
for line in open(flair_gtf):
    if line.startswith('#'):
        continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon':
        continue
    t = tid_of(f[8]); ex[t].append((int(f[3]), int(f[4]))); strand[t] = f[6]; chrom[t] = f[0]
iso_chain = {t: chain_key(ex[t], chrom[t], strand[t]) for t in ex}

# truth map
truth = {}
for line in open(truth_tsv):
    k, v = line.rstrip('\n').split('\t'); truth[k] = v

# adjudicated verdicts
hdr = None; verdict = {}; cat = {}
for line in open(adj_tsv):
    f = line.rstrip('\n').split('\t')
    if hdr is None:
        hdr = {c: i for i, c in enumerate(f)}; continue
    iid = f[hdr['isoform_id']]
    verdict[iid] = f[hdr['confidence_class']]
    cat[iid] = f[hdr['structural_category']]

# assign truth label per isoform
def truth_of(iid):
    ck = iso_chain.get(iid)
    if ck is None:
        return 'monoexonic'
    return truth.get(ck, 'false_novel')

rows = [(iid, truth_of(iid), verdict.get(iid, 'NA'), cat.get(iid, 'NA')) for iid in verdict]

# ---- confusion: truth x confidence class ----
conf = defaultdict(Counter)
for _, tl, vc, _ in rows:
    conf[tl][vc] += 1
classes = ["HIGH_CONF_KNOWN","HIGH_CONF_NOVEL","MEDIUM_CONF_NOVEL","LOW_CONF_PARTIAL","AMBIGUOUS","ARTIFACT"]
print("=== confusion (truth rows x PanIsoGuard class) ===")
print("truth".ljust(14) + "".join(c[:9].rjust(10) for c in classes) + "   total")
for tl in ["genuine_novel","known","false_novel","monoexonic"]:
    c = conf[tl]; print(tl.ljust(14) + "".join(str(c[k]).rjust(10) for k in classes) + str(sum(c.values())).rjust(8))

# ---- genuine-novel detection (decisive): predicted genuine = HIGH/MEDIUM_CONF_NOVEL ----
# restrict to the discrimination task: truth in {genuine_novel, false_novel} (multi-exon)
GEN = {"HIGH_CONF_NOVEL","MEDIUM_CONF_NOVEL"}
tp = sum(1 for _,tl,vc,_ in rows if tl=='genuine_novel' and vc in GEN)
fn = sum(1 for _,tl,vc,_ in rows if tl=='genuine_novel' and vc not in GEN)
fp = sum(1 for _,tl,vc,_ in rows if tl=='false_novel'  and vc in GEN)
tn = sum(1 for _,tl,vc,_ in rows if tl=='false_novel'  and vc not in GEN)
prec = tp/(tp+fp) if tp+fp else float('nan')
rec  = tp/(tp+fn) if tp+fn else float('nan')
spec = tn/(tn+fp) if tn+fp else float('nan')
print("\n=== genuine-novel detection (predicted genuine = HIGH/MEDIUM_CONF_NOVEL) ===")
print(f"genuine_novel={tp+fn}  false_novel={fp+tn}")
print(f"precision={prec:.3f}  recall={rec:.3f}  specificity(false rejected)={spec:.3f}")

# ---- per-SQANTI-category recall of genuine novels ----
catrec = defaultdict(lambda:[0,0])
for _,tl,vc,ct in rows:
    if tl=='genuine_novel':
        catrec[ct][1]+=1
        if vc in GEN: catrec[ct][0]+=1
print("\n=== genuine-novel recall by SQANTI category ===")
for ct,(hit,tot) in sorted(catrec.items()):
    print(f"  {ct:28s} {hit}/{tot}  recall={hit/tot:.3f}")

# ---- AUPRC for genuine-vs-false ranking ----
score = {"HIGH_CONF_NOVEL":1.0,"MEDIUM_CONF_NOVEL":0.75,"AMBIGUOUS":0.5,
         "LOW_CONF_PARTIAL":0.25,"ARTIFACT":0.0,"HIGH_CONF_KNOWN":0.1,"PAN_REF_RESCUED_FALSE_NOVEL":0.0,"NA":0.0}
y=[]; s=[]
for _,tl,vc,_ in rows:
    if tl in ('genuine_novel','false_novel'):
        y.append(1 if tl=='genuine_novel' else 0); s.append(score.get(vc,0.0))
try:
    from sklearn.metrics import average_precision_score
    print(f"\n=== AUPRC (genuine-novel vs false, n={len(y)}, positives={sum(y)}) ===")
    print(f"AUPRC = {average_precision_score(y, s):.4f}   (baseline = {sum(y)/len(y):.4f})")
except Exception as e:
    print("AUPRC skipped:", e)
