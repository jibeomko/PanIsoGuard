import re
from collections import Counter
tid_re=re.compile(r'transcript_id "([^"]+)"')
# flair isoform chains (intron tuples, same convention as prep: (exon_i_end_1based, exon_{i+1}_start_1based-1))
ex={}
for ln in open("flair.isoforms.gtf"):
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="exon": continue
    m=tid_re.search(f[8])
    if not m: continue
    ex.setdefault(m.group(1),[]).append((int(f[3]),int(f[4])))
def chain(e):
    e=sorted(e); return tuple((e[i][1], e[i+1][0]-1) for i in range(len(e)-1))
iso_chain={t:chain(e) for t,e in ex.items()}
# hidden (genuine) chains
hidden=set()
for ln in open("truth_hidden_chains.tsv"):
    t,ch=ln.rstrip("\n").split("\t")
    if ch: hidden.add(tuple(tuple(int(x) for x in p.split("-")) for p in ch.split(";")))
# adjudication
cls={}; cat={}; h=None
for ln in open("e2e_adj.adjudicated.tsv"):
    f=ln.rstrip("\n").split("\t")
    if h is None: h={n:k for k,n in enumerate(f)}; continue
    iso=f[h["isoform_id"]]; cls[iso]=f[h["confidence_class"]]; cat[iso]=f[h["structural_category"]]
# recovery: how many hidden chains exactly reconstructed by FLAIR
recovered=sum(1 for c in iso_chain.values() if c in hidden)
print(f"hidden(genuine) transcripts={len(hidden)}  FLAIR isoforms reconstructing a hidden chain (exact)={recovered}")
# score the SQANTI-novel set (NIC/NNC) by truth
def pred(c):
    if c in ("HIGH_CONF_NOVEL","MEDIUM_CONF_NOVEL"): return "genuine"
    if c in ("LOW_CONF_PARTIAL","ARTIFACT"): return "false"
    return "abstain"
conf=Counter(); bytruth_class=Counter()
for iso,c in cat.items():
    if c not in ("novel_in_catalog","novel_not_in_catalog"): continue
    t = "genuine" if iso_chain.get(iso) in hidden else "false"
    conf[(t,pred(cls[iso]))]+=1
    bytruth_class[(t,cls[iso])]+=1
print("\n=== SQANTI-novel (NIC/NNC) scored by intron-chain truth ===")
print(f"{'':9}pred_genuine pred_false abstain")
for t in ("genuine","false"):
    print(f"{t:9}{conf[(t,'genuine')]:11}{conf[(t,'false')]:11}{conf[(t,'abstain')]:8}")
TP,FP=conf[("genuine","genuine")],conf[("false","genuine")]
FN,TN=conf[("genuine","false")],conf[("false","false")]
prec=TP/(TP+FP) if TP+FP else 0; spec=TN/(TN+FP) if TN+FP else 0; rec=TP/(TP+FN) if TP+FN else 0
print(f"\nprecision(genuine)={prec:.3f}  specificity(false)={spec:.3f}  genuine-recall(decisive)={rec:.3f}")
print(f"TP={TP} FP={FP} TN={TN} FN={FN}  genuine_total={TP+FN+conf[('genuine','abstain')]} false_total={TN+FP+conf[('false','abstain')]}")
print("\n=== class breakdown by truth ===")
for t in ("genuine","false"):
    cc=", ".join(f"{k}={v}" for (tt,k),v in sorted(bytruth_class.items()) if tt==t)
    print(f"  {t}: {cc}")
