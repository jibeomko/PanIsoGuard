#!/usr/bin/env python3
"""Label arm A's novel junctions by what the same reads do on the personal genome (arm B).

For every novel junction j of the GRCh38 pipeline (SQANTI3 junctions file of arm A), take
the primary alignments that splice exactly at j on GRCh38, and look at the same reads'
primary alignments on the personal genome. A read's candidate intron there is the one
whose *read offset* (the position within the read where the intron starts) is closest to
j's; its genomic identity is then verified in GRCh38 coordinates (directly for an arm that
shares them, --same-coordinates; through `paftools.js liftover` on the assembly-to-GRCh38
PAF otherwise, --paf). Each read gets one class:

  same-canonical     same read offset and same GRCh38 coordinates; the intron is canonical
                     (GT-AG / GC-AG / AT-AC) on the personal genome
  same-noncanonical  same intron, still non-canonical on the personal genome
  shifted-known      the read is split within --window bp of the same read position but the
                     intron differs, and the personal-genome intron is annotated (lifted GENCODE)
  shifted-novel      as above, not annotated
  unspliced          the read aligns through that position at j's locus with no intron within
                     --window bp (the GRCh38 "intron" is not a splice in the person)
  unassessable       read unmapped / not primary on the personal genome, aligned at another
                     locus, position outside its aligned span, or coordinates that cannot be
                     verified

The junction's label is the majority class over its assessable reads (>= --min-reads).
Reference-bias truth for the scorer: same-canonical (when j was non-canonical on GRCh38),
shifted-known, and unspliced; same-noncanonical with a canonical GRCh38 junction is the
blind spot in the other direction (a motif *disrupted* in the person).

The motif on the personal genome is read on the transcript strand of arm A's junction,
flipped when the read maps in the opposite orientation there (an inverted assembly contig).
MAPQ is not used: on a diploid reference nearly every read is ambiguous between the two
haplotypes, which is expected, not a placement problem.

Usage:
  classify_junctions.py --junctions armA/sq/iq_junctions.txt --bam-a rna.aln.bam \
      --bam-b armB1/rna.aln.bam --fasta-b armB1/diploid.fa --gtf-b armB1/diploid.gtf \
      --same-coordinates --out armB1/junction_classes.tsv
  classify_junctions.py ... --bam-b armB2/rna.aln.bam --fasta-b armB2/diploid.fa \
      --gtf-b armB2/diploid.gtf --paf paf/mat.asm5.paf paf/pat.asm5.paf --out armB2/...
"""
import argparse, collections, json, os, re, subprocess, sys, tempfile
import pysam

CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
COMP = str.maketrans("ACGTacgt", "TGCAtgca")
rc = lambda s: s.translate(COMP)[::-1]
M, I, D, N, S, H, EQ, X = 0, 1, 2, 3, 4, 5, 7, 8
CLASSES = ["same-canonical", "same-noncanonical", "shifted-known", "shifted-novel", "unspliced", "unassessable"]


def introns_with_offsets(aln):
    """[(start0, end, q)] per N op (reference coords + read offset of the splice, in the ORIGINAL
    read orientation), the aligned query span in that orientation, and the read length."""
    out, q, r = [], 0, aln.reference_start
    qa0 = qa1 = None
    for op, n in aln.cigartuples:
        if op in (M, EQ, X):
            if qa0 is None: qa0 = q
            q += n; r += n; qa1 = q
        elif op == I:
            q += n
        elif op == D:
            r += n
        elif op == N:
            out.append((r, r + n, q)); r += n
        elif op in (S, H):
            q += n
    L = q
    if aln.is_reverse:
        out = [(a, b, L - qq) for a, b, qq in out]
        qa0, qa1 = L - qa1, L - qa0
    return out, (qa0, qa1), L


def motif(fa, contig, start0, end, strand):
    try:
        d, a = fa.fetch(contig, start0, start0 + 2).upper(), fa.fetch(contig, end - 2, end).upper()
    except (KeyError, ValueError):
        return "NA", False
    if strand == "-":
        d, a = rc(a), rc(d)
    return d + a, (d, a) in CANON


def gtf_introns(path):
    """Set of (contig, start0, end) introns of every transcript in a GTF."""
    ex = collections.defaultdict(list)
    tid = re.compile(r'transcript_id "([^"]+)"')
    with open(path) as f:
        for ln in f:
            if ln.startswith("#"): continue
            p = ln.split("\t", 9)
            if len(p) < 9 or p[2] != "exon": continue
            m = tid.search(p[8])
            if m: ex[(p[0], m.group(1))].append((int(p[3]), int(p[4])))
    introns = set()
    for (contig, _), exons in ex.items():
        exons.sort()
        for (s1, e1), (s2, e2) in zip(exons, exons[1:]):
            if s2 - 1 > e1: introns.add((contig, e1, s2 - 1))
    return introns


def liftover(intervals, pafs, paftools):
    """(full, partial): {(contig, start0, end): (chrom, start0, end)} through the first PAF that
    lifts it. paftools marks an interval whose 5' or 3' end falls outside the alignment block
    with _t5 / _t3 and clamps that end; those go to `partial` (good for "which locus", not for
    exact coordinates)."""
    out, partial = {}, {}
    todo = list(intervals)
    for paf in pafs:
        if not todo: break
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as bed:
            for c, s, e in todo: bed.write(f"{c}\t{s}\t{e}\t{c}|{s}|{e}\n")
        try:
            res = subprocess.run([paftools, "liftover", paf, bed.name], capture_output=True, text=True, check=True)
        finally:
            os.unlink(bed.name)
        for ln in res.stdout.splitlines():
            p = ln.split("\t")
            if len(p) < 4: continue
            # paftools ignores the BED name and writes <contig>_<start>_<end>[_t5][_t3]
            # (assembly contig names carry no underscore)
            name, trunc = p[3], False
            while name.endswith(("_t5", "_t3")):
                name, trunc = name[:-3], True
            c, s, e = name.rsplit("_", 2)
            k = (c, int(s), int(e))
            if trunc:
                partial.setdefault(k, (p[0], int(p[1]), int(p[2])))
            else:
                out[k] = (p[0], int(p[1]), int(p[2]))
        todo = [k for k in todo if k not in out]
    return out, partial


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--junctions", required=True, help="SQANTI3 *_junctions.txt of arm A (GRCh38)")
    ap.add_argument("--bam-a", required=True, help="arm A BAM (reads on GRCh38)")
    ap.add_argument("--bam-b", required=True, help="arm B BAM (same reads on the personal genome)")
    ap.add_argument("--fasta-b", required=True, help="arm B genome (faidx)")
    ap.add_argument("--gtf-b", required=True, help="arm B annotation (introns = known)")
    ap.add_argument("--out", required=True, help="per-junction TSV; summary JSON goes to <out>.summary.json")
    ap.add_argument("--min-reads", type=int, default=2, help="assessable reads needed to label a junction")
    ap.add_argument("--window", type=int, default=100, help="read-offset window for 'shifted' (bp)")
    ap.add_argument("--same-coordinates", action="store_true",
                    help="arm B shares GRCh38 coordinates; contigs are <chrom>_<hap>")
    ap.add_argument("--paf", nargs="*", default=[], help="assembly->GRCh38 PAFs (minimap2 -cx asm5 --cs) for liftover")
    ap.add_argument("--paftools", default=os.path.expanduser("~/miniconda3/envs/flair300/bin/paftools.js"))
    ap.add_argument("--locus-window", type=int, default=100000, help="how far from j a read may lie and still count as 'at the locus'")
    ap.add_argument("--lift-slop", type=int, default=2, help="bp tolerance when comparing lifted coordinates")
    a = ap.parse_args()
    if not a.same_coordinates and not a.paf:
        sys.exit("need --same-coordinates or --paf")

    # 1. arm A novel junctions from SQANTI3
    J = {}   # (chrom, start0, end) -> dict
    with open(a.junctions) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if p[ix["junction_category"]] != "novel": continue
            k = (p[ix["chrom"]], int(p[ix["genomic_start_coord"]]) - 1, int(p[ix["genomic_end_coord"]]))
            j = J.setdefault(k, {"strand": p[ix["strand"]], "canonical_A": p[ix["canonical"]] == "canonical",
                                 "motif_A": p[ix["splice_site"]], "isoforms": set()})
            j["isoforms"].add(p[ix["isoform"]])
    print(f"[1] {len(J)} novel junctions in arm A", file=sys.stderr)

    # 2. reads spanning them on GRCh38: read -> [(jkey, q_offset, is_reverse)]
    want = collections.defaultdict(list)
    n_reads = 0
    with pysam.AlignmentFile(a.bam_a) as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: continue
            n_reads += 1
            chrom = aln.reference_name
            for s0, e, q in introns_with_offsets(aln)[0]:
                k = (chrom, s0, e)
                if k in J: want[aln.query_name].append((k, q, aln.is_reverse))
    print(f"[2] {n_reads} primary alignments on GRCh38; {len(want)} reads span a novel junction", file=sys.stderr)

    # 3. the same reads on the personal genome: collect each (read, j) observation
    known_B = gtf_introns(a.gtf_b)
    print(f"[3] {len(known_B)} annotated introns on the personal genome", file=sys.stderr)
    obs = []      # (jkey, contig, strand_B, cand or None, span (s, e) or None, q_inside)
    seen = set()
    with pysam.AlignmentFile(a.bam_b) as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: continue
            name = aln.query_name
            if name not in want: continue
            seen.add(name)
            contig = aln.reference_name
            introns, (qa0, qa1), L = introns_with_offsets(aln)
            for k, qA, revA in want[name]:
                strand = J[k]["strand"] if aln.is_reverse == revA else ("-" if J[k]["strand"] == "+" else "+")
                best = None
                for s0, e, qB in introns:
                    d = abs(qB - qA)
                    if d <= a.window and (best is None or d < best[0]): best = (d, s0, e)
                inside = qa0 is not None and qa0 + 20 <= qA <= qa1 - 20
                obs.append((k, contig, strand, best, (aln.reference_start, aln.reference_end), inside))
    print(f"[4] {len(seen)} of {len(want)} spanning reads found on the personal genome; {len(obs)} observations", file=sys.stderr)

    # 4. verify genomic identity in GRCh38 coordinates
    def to_grch38_direct(contig, s, e):
        return (contig.rsplit("_", 1)[0], s, e)
    lifted, lifted_partial = {}, {}
    if not a.same_coordinates:
        ivs = set()
        for k, contig, strand, best, span, inside in obs:
            if best is not None: ivs.add((contig, best[1], best[2]))
            ivs.add((contig, span[0], span[1]))
        lifted, lifted_partial = liftover(ivs, a.paf, a.paftools)
        print(f"[5] lifted {len(lifted)} of {len(ivs)} intervals to GRCh38 ({len(lifted_partial)} partial)", file=sys.stderr)
    def to_grch38(contig, s, e, approximate=False):
        if a.same_coordinates: return to_grch38_direct(contig, s, e)
        g = lifted.get((contig, s, e))
        return g if g is not None or not approximate else lifted_partial.get((contig, s, e))

    # 5. classify observations
    fa = pysam.FastaFile(a.fasta_b)
    per_j = collections.defaultdict(collections.Counter)
    detail = collections.defaultdict(list)
    for k, contig, strand, best, span, inside in obs:
        if best is not None:
            g = to_grch38(contig, best[1], best[2])
            if g is None:
                cls = "unassessable"
            elif g[0] == k[0] and abs(g[1] - k[1]) <= a.lift_slop and abs(g[2] - k[2]) <= a.lift_slop:
                mo, can = motif(fa, contig, best[1], best[2], strand)
                cls = "same-canonical" if can else "same-noncanonical"
                detail[k].append((cls, contig, best[1], best[2], mo, (contig, best[1], best[2]) in known_B))
            elif g[0] == k[0] and abs(g[1] - k[1]) <= a.locus_window:
                kn = (contig, best[1], best[2]) in known_B
                mo, _ = motif(fa, contig, best[1], best[2], strand)
                cls = "shifted-known" if kn else "shifted-novel"
                detail[k].append((cls, contig, best[1], best[2], mo, kn))
            else:
                cls = "unassessable"                      # the read is somewhere else
        else:
            g = to_grch38(contig, span[0], span[1], approximate=True)
            at_locus = g is not None and g[0] == k[0] and g[1] - a.locus_window <= k[1] <= g[2] + a.locus_window
            if at_locus and inside:
                cls = "unspliced"
                detail[k].append((cls, contig, -1, -1, "NA", False))
            else:
                cls = "unassessable"
        per_j[k][cls] += 1
    for name, hits in want.items():          # reads with no primary alignment on the personal genome
        if name not in seen:
            for k, _, _ in hits: per_j[k]["unassessable"] += 1

    # 6. per-junction label
    summary = collections.Counter()
    with open(a.out, "w") as out:
        out.write("\t".join(["chrom", "start0", "end", "strand", "canonical_A", "motif_A", "n_isoforms_A", "n_reads_A",
                             "label"] + ["n_" + c for c in CLASSES] + ["B_contig", "B_start0", "B_end", "B_motif", "B_known"]) + "\n")
        for k, j in sorted(J.items()):
            c = per_j.get(k, collections.Counter())
            n_A = sum(c.values())
            assessable = {cl: n for cl, n in c.items() if cl != "unassessable"}
            if n_A == 0:
                label = "no-reads"
            elif sum(assessable.values()) < a.min_reads:
                label = "unassessable"
            else:
                label = max(assessable.items(), key=lambda t: (t[1], t[0]))[0]
            rep = next((d for d in detail.get(k, []) if d[0] == label), None)
            bcol = [rep[1], str(rep[2]), str(rep[3]), rep[4], str(rep[5])] if rep else ["NA"] * 5
            summary[(label, "canonical_A" if j["canonical_A"] else "noncanonical_A")] += 1
            out.write("\t".join([k[0], str(k[1]), str(k[2]), j["strand"], str(j["canonical_A"]), j["motif_A"],
                                 str(len(j["isoforms"])), str(n_A), label] + [str(c.get(cl, 0)) for cl in CLASSES] + bcol) + "\n")
    js = {"n_novel_junctions_A": len(J), "n_spanning_reads_A": len(want), "n_reads_found_B": len(seen),
          "labels": {f"{l}|{ca}": n for (l, ca), n in sorted(summary.items())}}
    json.dump(js, open(a.out + ".summary.json", "w"), indent=1)
    print(json.dumps(js, indent=1), file=sys.stderr)


if __name__ == "__main__":
    main()
