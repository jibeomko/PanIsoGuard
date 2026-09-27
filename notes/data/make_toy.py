#!/usr/bin/env python3
"""Write the toy locus used by every study note (notes/*.md).

One gene (GENE_T, + strand) on a 1,700 bp contig `chrT`. All coordinates below are
1-based inclusive, like GTF. The script writes the hand-designed inputs; build_toy.sh
then runs STAR, SQANTI3 and samtools on them to produce the rest.

    python3 make_toy.py          # run from notes/data/
"""
import random

random.seed(1)
L = 1700
seq = [random.choice("ACGT") for _ in range(L)]


def put(pos, s):
    """Write string s starting at 1-based position pos."""
    for i, c in enumerate(s):
        seq[pos - 1 + i] = c


# Exons of the reference gene (100 bp each, introns 200 bp).
E1, E2, E3, E4, E5 = (101, 200), (401, 500), (701, 800), (1001, 1100), (1301, 1400)

# Splice-site dinucleotides. Donor = first two intron bases, acceptor = last two.
for donor in (201, 501, 801, 1101, 941):        # 941: donor of the extra exon in iso_C
    put(donor, "GT")
for acc_end in (400, 700, 1000, 1300, 860, 1020, 1004):
    put(acc_end - 1, "AG")
put(729, "AC")   # iso_A acceptor: non-canonical on the reference ...
put(449, "TT")   # iso_F acceptor: non-canonical, and stays that way

ref = "".join(seq)
hap1 = ref[:729] + "G" + ref[730:]   # ... one SNV (C->G at 730) makes it AG on haplotype 1
hap2 = ref                           # haplotype 2 carries the reference allele


def fasta(path, s):
    with open(path, "w") as f:
        f.write(">chrT\n")
        for i in range(0, len(s), 60):
            f.write(s[i:i + 60] + "\n")


fasta("toy.fa", ref)
fasta("hap1.fa", hap1)
fasta("hap2.fa", hap2)

# Reference annotation: TX1 uses all five exons, TX2 skips E2, TX3 skips E4.
REF = {
    "TX1": [E1, E2, E3, E4, E5],
    "TX2": [E1, E3, E4, E5],
    "TX3": [E1, E2, E3, E5],
}
# Caller 1 (the caller we adjudicate).
CALLER1 = {
    "iso_known": [E1, E2, E3, E4, E5],               # = TX1
    "iso_ism":   [E2, E3, E4, E5],                   # truncated TX1
    "iso_A":     [E1, E2, (731, 800), E4, E5],       # novel acceptor 730 (the protagonist)
    "iso_B":     [E1, E2, E3, (1021, 1100), E5],     # novel acceptor 1020
    "iso_C":     [E1, E2, E3, (861, 940), E4, E5],   # extra exon: two novel junctions
    "iso_D":     [E1, E3, E5],                       # known junctions, new combination
    "iso_E":     [E1, E2, E4, E5],                   # skips E3: known sites, new intron
    "iso_F":     [E1, (451, 500), E3, E4, E5],       # novel acceptor 450 (TT)
    "iso_G":     [E1, E2, E3, (1005, 1100), E5],     # acceptor 4 bp off the known 1000
}
# Caller 2: a second tool run on the same sample (different read ends, some shared chains).
CALLER2 = {
    "c2_1": [(121, 200), E2, E3, E4, (1301, 1380)],  # same chain as iso_known
    "c2_2": [(131, 200), E2, E3, (1021, 1100), E5],  # same chain as iso_B
    "c2_3": [E1, E2, E3, (1005, 1100), E5],          # same chain as iso_G
    "c2_4": [E1, E2, E3, (1051, 1100), E5],          # caller-2-only novel acceptor
}


def gtf(path, source, models):
    with open(path, "w") as f:
        for tid, exons in models.items():
            s, e = exons[0][0], exons[-1][1]
            attr = f'gene_id "GENE_T"; transcript_id "{tid}";'
            f.write(f"chrT\t{source}\ttranscript\t{s}\t{e}\t.\t+\t.\t{attr}\n")
            for a, b in exons:
                f.write(f"chrT\t{source}\texon\t{a}\t{b}\t.\t+\t.\t{attr}\n")


gtf("reference.gtf", "toyref", REF)
gtf("caller1.gtf", "caller1", CALLER1)
gtf("caller2.gtf", "caller2", CALLER2)


def cigar(exons, dels=None):
    """M for exons, N for introns; dels = {exon_index: (offset, length)} inserts a deletion."""
    ops = []
    for i, (a, b) in enumerate(exons):
        if i:
            ops.append(f"{a - exons[i - 1][1] - 1}N")
        n = b - a + 1
        if dels and i in dels:
            off, dl = dels[i]
            ops += [f"{off}M", f"{dl}D", f"{n - off - dl}M"]
        else:
            ops.append(f"{n}M")
    return "".join(ops)


# Long reads (the BAM for the mapping axis): full-length alignments, sequence omitted ('*').
# PanIsoGuard reads only FLAG, MAPQ and CIGAR.
reads = []
for iso in ("iso_known", "iso_A", "iso_B", "iso_C", "iso_E", "iso_F"):
    for k in range(5):
        reads.append((f"{iso}.r{k + 1}", 0, CALLER1[iso], 60, None))
for k in range(4):
    # iso_G: 3 of 4 reads carry a 3 bp deletion 5 bp after the junction (inside exon 4').
    reads.append((f"iso_G.r{k + 1}", 0, CALLER1["iso_G"], 60, {3: (5, 3)} if k < 3 else None))

with open("long_reads.sam", "w") as f:
    f.write("@HD\tVN:1.6\tSO:coordinate\n")
    f.write(f"@SQ\tSN:chrT\tLN:{L}\n")
    for name, flag, exons, mapq, dels in sorted(reads, key=lambda r: r[2][0][0]):
        f.write(f"{name}\t{flag}\tchrT\t{exons[0][0]}\t{mapq}\t{cigar(exons, dels)}\t*\t0\t0\t*\t*\n")


def spliced(exons, s):
    return "".join(s[a - 1:b] for a, b in exons)


# Short reads (100 bp, single-end) that each span exactly one junction of one isoform.
# junction k of an isoform lies between exons k and k+1; n reads per junction.
SHORT = [
    ("iso_known", 0, 30), ("iso_known", 1, 25), ("iso_known", 2, 20), ("iso_known", 3, 30),
    ("iso_A", 1, 9),     # the novel junction of iso_A (read from haplotype 1)
    ("iso_B", 2, 12),
    ("iso_C", 2, 7), ("iso_C", 3, 1),
    ("iso_E", 1, 6),
    ("iso_D", 0, 10), ("iso_D", 1, 10),
]
rng = random.Random(11)
with open("short_reads.fq", "w") as f:
    n = 0
    for iso, k, count in SHORT:
        exons = CALLER1[iso]
        tx = spliced(exons, hap1)
        cut = sum(b - a + 1 for a, b in exons[:k + 1])      # transcript offset of junction k
        for _ in range(count):
            up = rng.randint(35, 65)                        # bases upstream of the junction
            r = tx[cut - up:cut - up + 100]
            n += 1
            f.write(f"@{iso}.j{k + 1}.{n}\n{r}\n+\n{'I' * len(r)}\n")

print(f"wrote toy.fa hap1.fa hap2.fa reference.gtf caller1.gtf caller2.gtf "
      f"long_reads.sam ({len(reads)} reads) short_reads.fq ({n} reads)")
print("ref 729-730:", ref[728:730], "| hap1 729-730:", hap1[728:730])
