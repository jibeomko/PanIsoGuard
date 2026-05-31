import re,random
random.seed(7)
import sys
gtf=sys.argv[1]
# read chr22 sequence
seq=[]
for ln in open(sys.argv[2]):
    if ln.startswith(">"): continue
    seq.append(ln.strip())
S="".join(seq).upper()
comp=str.maketrans("ACGTN","TGCAN")
def rc(s): return s.translate(comp)[::-1]
tid_re=re.compile(r'transcript_id "([^"]+)"')
tx={}
for ln in open(gtf):
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or f[2]!="exon": continue
    m=tid_re.search(f[8]);
    if not m: continue
    d=tx.setdefault(m.group(1),{"strand":f[6],"exons":[],"lines":[]})
    d["exons"].append((int(f[3]),int(f[4])))
multi={t:d for t,d in tx.items() if len(d["exons"])>=3}   # >=3 exons -> richer chains
for d in multi.values(): d["exons"].sort()
chosen=sorted(multi)[:150]                                 # 150 transcripts
hidden=set(chosen[::5][:30])                               # 30 hidden = genuine-novel truth
def tseq(d):
    s="".join(S[a-1:b] for (a,b) in d["exons"])            # 1-based inclusive exons
    return rc(s) if d["strand"]=="-" else s
# pbsim transcript input (all chosen): id  count_sense  count_antisense  sequence
with open(sys.argv[3]+"/sim.transcript","w") as o:
    n=0
    for t in chosen:
        s=tseq(multi[t])
        if len(s)<300: continue
        o.write(f"{t}\t40\t0\t{s}\n"); n+=1
print("sim transcripts:",n,"hidden(genuine):",len(hidden))
# reduced annotation = full chr22 GTF minus hidden transcripts
with open(sys.argv[3]+"/reduced.gtf","w") as o:
    for ln in open(gtf):
        m=tid_re.search(ln)
        if m and m.group(1) in hidden: continue
        o.write(ln)
# truth chains (hidden transcripts' intron chains, 0-based half-open) for matching
def introns(d):
    e=d["exons"]; return [(e[i][1], e[i+1][0]-1) for i in range(len(e)-1)]
with open(sys.argv[3]+"/truth_hidden_chains.tsv","w") as o:
    for t in hidden:
        ch=";".join(f"{a}-{b}" for (a,b) in introns(multi[t]))
        o.write(f"{t}\t{ch}\n")
# SJ.tab: all real introns of chosen transcripts (short-read support)
seen=set()
with open(sys.argv[3]+"/real.SJ.tab","w") as o:
    sc={"+":1,"-":2}
    for t in chosen:
        d=multi[t]
        for (a,b) in introns(d):
            k=(a,b)
            if b<a or k in seen: continue
            seen.add(k); o.write(f"chr22\t{a+1}\t{b}\t{sc[d['strand']]}\t1\t1\t20\t0\t30\n")
print("reduced.gtf + truth + SJ written; real junctions:",len(seen))
