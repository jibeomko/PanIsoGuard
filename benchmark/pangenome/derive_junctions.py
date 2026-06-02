#!/usr/bin/env python3
"""Derive graph-supported splice junctions from a `vg deconstruct` VCF (GATE-1).

A population **deletion** carried on a pangenome haplotype is the genomic form of a
reference-bias "novel intron": an RNA read from a haplotype that carries the deletion,
aligned to the linear reference, spans the deleted span as an apparent novel intron
(an N/skip), even though no new splice site exists. So each deletion allele in the
graph's deconstruct VCF (relative to the linear reference path) becomes a
graph-supported junction at the deleted span.

  VCF deletion (left-anchored):  POS  REF=anchor+deleted  ALT=anchor
  deleted span (1-based incl.):  [POS + len(ALT) , POS + len(REF) - 1]
  -> PanIsoGuard junction at that intron, for BOTH strands (a genomic deletion is
     strand-agnostic; the isoform's novel intron may be on either strand).
  n_haplotypes = AC (population allele count) -> feeds the min_haplotypes gate.

Independence (non-circular) holds only because the adjudicated sample is OUT of the
graph — verify e.g. `vg paths --list-samples` does not include it. The output is the
`--pangenome-junctions` file consumed by PanIsoGuard's pangenome reference-bias axis.

Usage: derive_junctions.py <deconstruct.vcf> <out.tsv> [min_del_len=50]
"""
import sys


def parse_info_ac(info: str, n_alts: int):
    """Per-allele AC list (population allele count); fall back to NS, then 1."""
    ac = ns = None
    for kv in info.split(";"):
        if kv.startswith("AC="):
            ac = kv[3:]
        elif kv.startswith("NS="):
            ns = kv[3:]
    if ac is not None:
        vals = ac.split(",")
        out = []
        for i in range(n_alts):
            try:
                out.append(max(1, int(vals[i])))
            except (ValueError, IndexError):
                out.append(1)
        return out
    base = 1
    if ns is not None:
        try:
            base = max(1, int(ns))
        except ValueError:
            base = 1
    return [base] * n_alts


def main() -> int:
    if len(sys.argv) < 3:
        sys.stderr.write("usage: derive_junctions.py <deconstruct.vcf> <out.tsv> [min_del_len=50]\n")
        return 2
    vcf, out_path = sys.argv[1], sys.argv[2]
    min_del = int(sys.argv[3]) if len(sys.argv) > 3 else 50

    seen = {}  # (chrom, start, end) -> n_haplotypes (keep the max support seen)
    n_records = n_del = 0
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        n_records += 1
        # contig is e.g. "GRCh38#0#chr22" -> use the assembly contig name "chr22"
        chrom = f[0].split("#")[-1]
        try:
            pos = int(f[1])
        except ValueError:
            continue
        ref, alts = f[3], f[4].split(",")
        ac = parse_info_ac(f[7], len(alts))
        for i, alt in enumerate(alts):
            del_len = len(ref) - len(alt)
            if del_len < min_del:
                continue
            # left-anchored deletion: deleted span is REF beyond the shared ALT prefix
            start = pos + len(alt)          # 1-based first deleted base
            end = pos + len(ref) - 1        # 1-based last deleted base
            if start > end:
                continue
            n_del += 1
            key = (chrom, start, end)
            if ac[i] > seen.get(key, 0):
                seen[key] = ac[i]

    with open(out_path, "w") as o:
        o.write("# graph-supported junctions from deletions >= %dbp (derive_junctions.py)\n" % min_del)
        o.write("# chrom\tstart\tend\tstrand\tn_haplotypes  (1-based inclusive intron; both strands emitted)\n")
        for (chrom, start, end), nhap in sorted(seen.items()):
            o.write(f"{chrom}\t{start}\t{end}\t+\t{nhap}\n")
            o.write(f"{chrom}\t{start}\t{end}\t-\t{nhap}\n")

    sys.stderr.write(f"records={n_records}  deletions>= {min_del}bp={n_del}  "
                     f"unique junctions={len(seen)} (x2 strands = {2*len(seen)} rows)\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
