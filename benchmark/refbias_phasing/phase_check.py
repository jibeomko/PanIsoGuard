#!/usr/bin/env python3
"""Non-circular check of the reference-bias rescue: are junction reads from the right haplotype?

The rescue labels a novel junction reference bias when it is non-canonical on GRCh38 but
canonical on one of the individual's haplotypes. The benchmarks' truth label is that same
criterion, so they cannot test it (docs/validation.md). The motif-creating base is intronic,
so reads that use the junction never contain it. This script tests the claim another way:

  For a HETEROZYGOUS reference-bias junction (canonical on haplotype H only), long reads that
  use the junction should carry H's alleles at the gene's other heterozygous SNVs. Reads from
  the other haplotype O should not use it, while O is still expressed at the locus (non-junction
  reads carry O alleles).

Heterozygous sites: SNVs where the two assembly haplotypes (paftools call VCFs) differ; a site
absent from one VCF is taken as the reference allele there (~91% of GRCh38 is covered by exactly
one contig per haplotype). Phase comes from the DNA assembly only. A read is assigned to a
haplotype when its alleles favour it by >= MARGIN sites, so one sequencing error or one bad site
cannot flip a read.

Usage: phase_check.py LABEL:WORKDIR [LABEL:WORKDIR ...] [--emit-metrics metrics.json]
  WORKDIR holds refbias_unique.tsv, rna.aln.bam(.bai), pat.snv.vcf.gz, mat.snv.vcf.gz
  (the benchmark/hg03516_refbias pipeline output)
"""
import collections
import json
import sys

import pysam
from scipy.stats import fisher_exact

MARGIN = 2                # a read's alleles must favour one haplotype by >= this many sites
MIN_PHASED_J = 3          # junction reads that must be phased for a junction to be testable
MIN_O_NONJ = 2            # other-haplotype reads at the locus, i.e. O is expressed there


def load_junctions(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)][1:]
    return [(c, int(s), int(e), st, h.split(",")) for c, s, e, st, h in rows]


def uses(read, s, e, window=10):
    """(read has an N op exactly == [s, e), read has an I/D within `window` bp of s or e)"""
    pos, hit, indel = read.reference_start, False, False
    for op, n in read.cigartuples:
        if op in (1, 2) and (abs(pos - s) <= window or abs(pos - e) <= window):
            indel = True
        if op == 3 and pos == s and pos + n == e:
            hit = True
        if op in (0, 2, 3, 7, 8):
            pos += n
    return hit, indel


def check(wd, chrom, s, e, hap):
    bam = pysam.AlignmentFile(f"{wd}/rna.aln.bam")
    vcf = {h: pysam.VariantFile(f"{wd}/{h}.snv.vcf.gz") for h in ("pat", "mat")}
    reads = []
    for r in bam.fetch(chrom, s, e):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        if r.reference_start < s and r.reference_end > e:        # spans the whole intron
            reads.append(uses(r, s, e) + (r,))
    if not reads:
        return None
    lo, hi = min(r.reference_start for _, _, r in reads), max(r.reference_end for _, _, r in reads)
    alt = {h: {v.pos - 1: (v.ref, v.alts[0]) for v in vcf[h].fetch(chrom, lo, hi)} for h in vcf}
    sites = {}
    for p in set(alt["pat"]) | set(alt["mat"]):
        ref = (alt["pat"].get(p) or alt["mat"].get(p))[0]
        a_pat = alt["pat"][p][1] if p in alt["pat"] else ref
        a_mat = alt["mat"][p][1] if p in alt["mat"] else ref
        if a_pat != a_mat and len(ref) == 1:
            sites[p] = {"pat": a_pat, "mat": a_mat}
    other = "mat" if hap == "pat" else "pat"
    tab, covered = collections.Counter(), set()                  # (is_j, "H"/"O"/"unphased")
    for is_j, _, r in reads:
        q = r.query_sequence
        b = {rp: q[qp] for qp, rp in r.get_aligned_pairs(matches_only=True) if rp in sites}
        covered |= b.keys()
        vh = sum(base == sites[p][hap] for p, base in b.items())
        vo = sum(base == sites[p][other] for p, base in b.items())
        tab[(is_j, "H" if vh - vo >= MARGIN else "O" if vo - vh >= MARGIN else "unphased")] += 1
    n_j = tab[(True, "H")] + tab[(True, "O")] + tab[(True, "unphased")]
    return {"het_sites": len(covered), "J_reads": n_j,
            "J_indel_near": round(sum(j and d for j, d, _ in reads) / n_j, 2) if n_j else None,
            "J_H": tab[(True, "H")], "J_O": tab[(True, "O")],
            "nonJ_H": tab[(False, "H")], "nonJ_O": tab[(False, "O")]}


def judge(wd, chrom, s, e, strand, haps):
    row = {"junction": f"{chrom}:{s}-{e}:{strand}", "haplotypes": ",".join(haps)}
    r = check(wd, chrom, s, e, haps[0]) or {}
    row.update(r)
    if len(haps) != 1:
        row["verdict"] = "untestable (homozygous)"
        for k in ("J_H", "J_O", "nonJ_H", "nonJ_O", "het_sites"):
            row.pop(k, None)
        return row
    phased_j = r.get("J_H", 0) + r.get("J_O", 0)
    if phased_j < MIN_PHASED_J:
        row["verdict"] = "uninformative (too few phased junction reads)"
    elif r["J_O"] > 0:
        row["verdict"] = "CONTRADICTED (junction reads from the non-canonical haplotype)"
    elif r["nonJ_O"] < MIN_O_NONJ:
        row["verdict"] = "uninformative (other haplotype not seen at the locus)"
    else:
        row["verdict"] = "consistent"
    if phased_j and r["nonJ_H"] + r["nonJ_O"]:
        row["fisher_p"] = fisher_exact([[r["J_H"], r["J_O"]], [r["nonJ_H"], r["nonJ_O"]]],
                                       alternative="greater")[1]
    return row


def main():
    pairs = [a.split(":", 1) for a in sys.argv[1:] if ":" in a and not a.startswith("--")]
    per = {}
    for label, wd in pairs:
        rows = [judge(wd, *j) for j in load_junctions(f"{wd}/refbias_unique.tsv")]
        for row in rows:
            keys = ("het_sites", "J_reads", "J_indel_near", "J_H", "J_O", "nonJ_H", "nonJ_O")
            nums = " ".join(f"{k}={row[k]}" for k in keys if k in row)
            p = f" p={row['fisher_p']:.2g}" if "fisher_p" in row else ""
            print(f"{label} {row['junction']:32s} {row['haplotypes']:7s} {nums}{p} -> {row['verdict']}")
        tested = [r for r in rows if r["verdict"].startswith(("consistent", "CONTRADICTED"))]
        per[label] = {
            "junctions": len(rows),
            "heterozygous": sum(1 for r in rows if "," not in r["haplotypes"]),
            "verdicts": dict(collections.Counter(r["verdict"].split(" (")[0] for r in rows)),
            "tested_junction_reads_motif_haplotype": sum(r["J_H"] for r in tested),
            "tested_junction_reads_other_haplotype": sum(r["J_O"] for r in tested),
            "tested_locus_reads_motif_haplotype": sum(r["nonJ_H"] for r in tested),
            "tested_locus_reads_other_haplotype": sum(r["nonJ_O"] for r in tested),
            "rows": rows,
        }
        print(f"{label} summary: {per[label]['verdicts']} | junction reads motif/other haplotype "
              f"{per[label]['tested_junction_reads_motif_haplotype']}/{per[label]['tested_junction_reads_other_haplotype']}")
    if "--emit-metrics" in sys.argv:
        env = {
            "protocol": "refbias_phasing", "scope": "hprc_r2_" + "_".join(l.lower() for l, _ in pairs),
            "status": "tracked", "tool_version": None, "ruleset_version": None,
            "generated_utc": None, "source": None,
            "data_provenance": "the benchmark/hg03516_refbias work dirs: public HPRC Release 2 PacBio Kinnex "
                               "Iso-Seq (minimap2 splice:hq to GRCh38) and trio-phased HiFi assembly SNVs "
                               "(pat/mat, paftools call). No PanIsoGuard run: tests the rescue's claim directly.",
            "command": "benchmark/refbias_phasing/phase_check.py " + " ".join(f"{l}:<workdir>" for l, _ in pairs)
                       + " --emit-metrics benchmark/results/refbias_phasing/metrics.json",
            "metrics": per,
            "notes": f"MARGIN={MARGIN}, MIN_PHASED_J={MIN_PHASED_J}, MIN_O_NONJ={MIN_O_NONJ}. 'motif haplotype' = "
                     "the haplotype on which the junction is canonical. Homozygous junctions cannot be phased.",
        }
        json.dump(env, open(sys.argv[sys.argv.index("--emit-metrics") + 1], "w"), indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
