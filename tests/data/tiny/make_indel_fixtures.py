#!/usr/bin/env python3
"""Builds the indel_* fixtures next to this file (needs bgzip, tabix and samtools on PATH).

One contig, chrI, with seven 30-bp introns between 40-bp exons. Each intron carries one case of
the variant axis; the unit tests in tests/unit/vcf_haplotypes_test.cpp say what each must give.

  intron            reference   the sample                                   reads in indel.bam
  A [40,70)         GA..AG      1-bp deletion makes the donor GT (1|1)        3
  B [110,140)       GT..TG      1-bp insertion makes the acceptor AG (1|1)    3
  C [180,210)       GT..TC      TTAG inserted exactly at the acceptor (1|1)   4, without the insertion
  G [250,280)       GT..TC      the same                                      4, carrying the insertion
  D [320,350)       CC..TT      a deletion equal to the intron (0|1)          4, haplotype 2
  E [390,420)       GG..AG      G>T makes the donor GT (0|1)                  4, haplotype 2
  F [460,490)       GG..AG      the same                                      4, haplotype 1

Heterozygous single-base sites in the exons (0|1) tell the two haplotypes apart.

Three more files on the same reference:
  indel_clash.vcf.gz      E: the G>T and two different insertions at that exon boundary, one haplotype;
                          F: the G>T and one insertion at the boundary written twice (once with a flanking base)
  indel_ps_string.vcf.gz  FORMAT/PS declared as a String: the two G>T in the sets PATMAT and 9, a third
                          heterozygous site in PATMAT
  indel_ps_float.vcf.gz   the same with FORMAT/PS declared as a Float
  indel_bad_row.vcf.gz    the two G>T after a row whose GT is not a genotype
  indel_set10.vcf.gz      the two G>T on haplotype 1, in phase set 10
  indel_set20.vcf.gz      the G>T of E alone, on haplotype 1, in phase set 20
  indel_unphased2.vcf.gz  the two G>T as 0/1
  chain_shallow.bam       three reads that carry the chain of E and F end to end
  chain_deep.bam          three such reads and nine that span the chain without splicing

long_ref.fa is another contig, chrL, with one 400-bp intron [100, 500). long_del.vcf.gz deletes it
(1|0); long_del_bad.vcf.gz is the same record with one base of its REF, 200 bases in, changed.
"""
import os
import random
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
rng = random.Random(11)
INTRONS = {"A": 40, "B": 110, "C": 180, "G": 250, "D": 320, "E": 390, "F": 460}
N = 530

seq = []
while len(seq) < N:   # no base twice in a row: nothing can slide by accident
    b = rng.choice("ACGT")
    if not seq or seq[-1] != b: seq.append(b)


def plant(pos, s):
    seq[pos:pos + len(s)] = s


plant(39, "CGATC");  plant(66, "TCAGC")         # A: donor GA, the deletion of A leaves GT
plant(109, "CGTAC"); plant(136, "CATGC")        # B: acceptor TG, an A between them gives AG
for s in (180, 250):                            # C and G: acceptor TC, exon starts with G
    plant(s - 1, "CGTAC"); plant(s + 26, "CATCGA")
plant(319, "ACCTA"); plant(346, "GATTCA")       # D: first intron base C = first base after it
plant(389, "CGGAC"); plant(416, "TCAGC")        # E
plant(459, "CGGAC"); plant(486, "TCAGC")        # F
ref = "".join(seq)
assert ref[320] == ref[350] and ref[321] != ref[351]

# (0-based position, REF, ALT, GT)
variants = [
    (40, "GA", "G", "1|1"),
    (138, "T", "TA", "1|1"),
    (209, "C", "CTTAG", "1|1"),
    (279, "C", "CTTAG", "1|1"),
    (320, ref[320:351], ref[320], "0|1"),       # the deletion written one base to the right of the junction
    (391, "G", "T", "0|1"),
    (461, "G", "T", "0|1"),
]
SITES = [300, 305, 360, 365, 380, 430, 435, 450, 500, 505]
for p in SITES:
    alt = {"A": "C", "C": "A", "G": "T", "T": "G"}[ref[p]]
    variants.append((p, ref[p], alt, "0|1"))
variants.sort()
for p, r, a, g in variants: assert ref[p:p + len(r)] == r, (p, r, ref[p:p + len(r)])


def write_fasta(name, contig, s):
    with open(os.path.join(HERE, name), "w") as f:
        f.write(f">{contig}\n")
        for i in range(0, len(s), 60): f.write(s[i:i + 60] + "\n")
    with open(os.path.join(HERE, name + ".fai"), "w") as f:
        f.write(f"{contig}\t{len(s)}\t{len(contig) + 2}\t60\t61\n")


def write_vcf(name, rows, ps_type=None, contig=("chrI", N)):
    path = os.path.join(HERE, name)
    with open(path[:-3], "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=%s,length=%d>\n" % contig)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        if ps_type: f.write(f'##FORMAT=<ID=PS,Number=1,Type={ps_type},Description="Phase set">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        fmt = "GT:PS" if ps_type else "GT"
        for p, r, a, g in rows: f.write(f"{contig[0]}\t{p + 1}\t.\t{r}\t{a}\t.\tPASS\t.\t{fmt}\t{g}\n")
    subprocess.check_call(["bgzip", "-f", path[:-3]])
    subprocess.check_call(["tabix", "-f", "-p", "vcf", path])


write_fasta("indel_ref.fa", "chrI", ref)
write_fasta("indel_hap_short.fa", "chrI", ref[:100] + ref[101:])   # one base shorter: other coordinates
write_vcf("indel.vcf.gz", variants)
write_vcf("indel_clash.vcf.gz", sorted([
    (389, ref[389], ref[389] + "C", "0|1"), (389, ref[389], ref[389] + "AA", "0|1"), (391, "G", "T", "0|1"),
    (458, ref[458:460], ref[458:460] + "C", "0|1"), (459, ref[459], ref[459] + "C", "0|1"), (461, "G", "T", "0|1")]))
write_vcf("indel_ps_string.vcf.gz", [(391, "G", "T", "0|1:PATMAT"), (430, ref[430], {"A": "C", "C": "A", "G": "T", "T": "G"}[ref[430]], "0|1:PATMAT"),
                                     (461, "G", "T", "0|1:9")], ps_type="String")
write_vcf("indel_ps_float.vcf.gz", [(391, "G", "T", "0|1:5"), (461, "G", "T", "0|1:9")], ps_type="Float")
write_vcf("indel_bad_row.vcf.gz", [(380, ref[380], "N", "x|y"), (391, "G", "T", "0|1"), (461, "G", "T", "0|1")])
write_vcf("indel_set10.vcf.gz", [(391, "G", "T", "1|0:10"), (461, "G", "T", "1|0:10")], ps_type="Integer")
write_vcf("indel_set20.vcf.gz", [(391, "G", "T", "1|0:20")], ps_type="Integer")
write_vcf("indel_unphased2.vcf.gz", [(391, "G", "T", "0/1"), (461, "G", "T", "0/1")])

lrng = random.Random(23)
lseq = []
while len(lseq) < 600:
    b = lrng.choice("ACGT")
    if not lseq or lseq[-1] != b: lseq.append(b)
lseq[100:102] = "CC"; lseq[498:500] = "TT"       # not a splice motif, and no base to slide the deletion over
lref = "".join(lseq)
assert lref[99] != lref[499] and lref[100] != lref[500]
write_fasta("long_ref.fa", "chrL", lref)
whole = lref[99:500]
other = whole[:200] + {"A": "C", "C": "A", "G": "T", "T": "G"}[whole[200]] + whole[201:]
for name, r in (("long_del.vcf.gz", whole), ("long_del_bad.vcf.gz", other)):
    write_vcf(name, [(99, r, r[0], "1|0")], contig=("chrL", len(lref)))

hap2 = list(ref)                                # haplotype 2 in reference coordinates (exons only matter)
for p in SITES: hap2[p] = {"A": "C", "C": "A", "G": "T", "T": "G"}[ref[p]]
hap2 = "".join(hap2)


def read(name, start, hap, inserted=""):
    s = hap[start - 30:start] + inserted + hap[start + 30:start + 60]
    cigar = "30M30N" + (f"{len(inserted)}I" if inserted else "") + "30M"
    return f"{name}\t0\tchrI\t{start - 30 + 1}\t60\t{cigar}\t*\t0\t0\t{s}\t*"


reads = []
for tag, hap, n, ins in (("A", ref, 3, ""), ("B", ref, 3, ""), ("C", ref, 4, ""), ("G", ref, 4, "TTAG"),
                         ("D", hap2, 4, ""), ("E", hap2, 4, ""), ("F", ref, 4, "")):
    reads += [read(f"{tag}{i}", INTRONS[tag], hap, ins) for i in range(n)]
sam = os.path.join(HERE, "indel.sam")
with open(sam, "w") as f:
    f.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrI\tLN:%d\n" % N)
    f.write("\n".join(reads) + "\n")
subprocess.check_call(["samtools", "sort", "--no-PG", "-O", "bam", "-o", os.path.join(HERE, "indel.bam"), sam])
subprocess.check_call(["samtools", "index", os.path.join(HERE, "indel.bam")])

# two libraries with the same three whole-chain reads of E-F and another depth
whole_chain = lambda i: f"w{i}\t0\tchrI\t361\t60\t30M30N40M30N30M\t*\t0\t0\t{ref[360:390] + ref[420:460] + ref[490:520]}\t*"
unspliced = lambda i: f"u{i}\t0\tchrI\t361\t60\t160M\t*\t0\t0\t{ref[360:520]}\t*"
for name, lines in (("chain_shallow", [whole_chain(i) for i in range(3)]),
                    ("chain_deep", [whole_chain(i) for i in range(3)] + [unspliced(i) for i in range(9)])):
    sam = os.path.join(HERE, name + ".sam")
    with open(sam, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chrI\tLN:%d\n" % N)
        f.write("\n".join(lines) + "\n")
    subprocess.check_call(["samtools", "sort", "--no-PG", "-O", "bam", "-o", os.path.join(HERE, name + ".bam"), sam])
    subprocess.check_call(["samtools", "index", os.path.join(HERE, name + ".bam")])
    os.remove(sam)
