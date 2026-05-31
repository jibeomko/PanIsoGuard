import subprocess, sys
bam=sys.argv[1]
out=open(sys.argv[2],"w")
import re
cig_re=re.compile(r'(\d+)([MIDNSHP=X])')
p=subprocess.Popen(["samtools","view","-F","0x904",bam],stdout=subprocess.PIPE,text=True)  # primary only
for line in p.stdout:
    f=line.rstrip("\n").split("\t")
    flag=int(f[1]); chrom=f[2]; pos0=int(f[3])-1; name=f[0]; cigar=f[5]
    if chrom=="*" or cigar=="*": continue
    strand="-" if (flag&16) else "+"
    blocks=[]; bs=pos0; cur=0; ref=pos0
    # build exon blocks split by N
    for n,op in cig_re.findall(cigar):
        n=int(n)
        if op in "M=X":
            cur+=n; ref+=n
        elif op=="D":
            cur+=n; ref+=n
        elif op=="N":
            blocks.append((bs, bs+cur)); ref+=n; bs=ref; cur=0
        # I,S,H,P: no ref consume
    blocks.append((bs, bs+cur))
    chromStart=blocks[0][0]; chromEnd=blocks[-1][1]
    sizes=",".join(str(e-s) for s,e in blocks)+","
    starts=",".join(str(s-chromStart) for s,e in blocks)+","
    out.write(f"{chrom}\t{chromStart}\t{chromEnd}\t{name}\t0\t{strand}\t{chromStart}\t{chromEnd}\t0\t{len(blocks)}\t{sizes}\t{starts}\n")
out.close()
