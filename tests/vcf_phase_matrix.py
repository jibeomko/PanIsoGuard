#!/usr/bin/env python3
"""Property check of the VCF haplotype route on one junction whose canonical motif needs two
variants: a donor G>T and an acceptor A>G on a GG-AA reference intron. Every pair of diploid
genotypes (0|0, 1|1, 1|0, 0|1, 1|., .|1, 0|., .|0, .|.) is run with the two records in one phase
set, in different phase sets, and unphased: 243 inputs. For each, the verdict must be REFERENCE_BIAS
exactly when a canonical haplotype is guaranteed under every completion of the missing alleles and,
without a shared phase set, under either relative orientation of the two records. In the two modes
without linkage the verdict must also be the same when the second record's alleles are written in
the other order: the file carries the same information either way.

Usage: vcf_phase_matrix.py <panisoguard> [--bgzip PATH] [--tabix PATH] [--work DIR]
Exit 0 when every case holds, 1 otherwise; prints the failures.
"""
import argparse, csv, itertools, json, os, shutil, subprocess, sys, tempfile

GENOTYPES = ["0|0", "1|1", "1|0", "0|1", "1|.", ".|1", "0|.", ".|0", ".|."]
HEADER = ('##fileformat=VCFv4.2\n##contig=<ID=chrT2,length=100>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
          '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')


def completions(gt):
    return list(itertools.product(*[([0, 1] if x == "." else [int(x)]) for x in gt.replace("/", "|").split("|")]))


def canonical_guaranteed(a, b, linked):
    """Some haplotype carries both ALTs under every completion (and, without linkage, either orientation)."""
    for x, y in itertools.product(completions(a), completions(b)):
        for z in ([y] if linked else [y, y[::-1]]):
            if not any(x[h] == 1 and z[h] == 1 for h in range(2)):
                return False
    return True


def swapped(gt):
    sep = "/" if "/" in gt else "|"
    return sep.join(reversed(gt.split(sep)))


def write_inputs(work):
    seq = ["A"] * 100
    seq[20] = seq[21] = "G"          # donor GG at 0-based 20-21 (the intron is {20,40})
    seq[38] = seq[39] = "A"          # acceptor AA at 38-39
    open(f"{work}/ref.fa", "w").write(">chrT2\n" + "".join(seq) + "\n")
    # an annotated transcript on the locus with another intron, so the caller's intron is novel in an annotated locus
    open(f"{work}/ref.gtf", "w").write('chrT2\tx\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                                       'chrT2\tx\texon\t50\t60\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    open(f"{work}/caller.gtf", "w").write('chrT2\tx\texon\t1\t20\t.\t+\t.\tgene_id "G1"; transcript_id "iso1";\n'
                                          'chrT2\tx\texon\t41\t60\t.\t+\t.\tgene_id "G1"; transcript_id "iso1";\n')
    open(f"{work}/sq.tsv", "w").write("isoform\tchrom\tstrand\tstructural_category\tall_canonical\tRTS_stage\tperc_A_downstream_TTS\n"
                                      "iso1\tchrT2\t+\tnovel_not_in_catalog\tnon_canonical\tFALSE\t0\n")


def run_case(pig, bgzip, tabix, work, name, gt_donor, gt_acceptor, mode):
    ps_a, ps_b = ("1", "1") if mode == "same_ps" else ("1", "2")
    if mode == "unphased":
        rows = f"chrT2\t22\t.\tG\tT\t.\tPASS\t.\tGT\t{gt_donor}\nchrT2\t40\t.\tA\tG\t.\tPASS\t.\tGT\t{gt_acceptor}\n"
    else:
        rows = f"chrT2\t22\t.\tG\tT\t.\tPASS\t.\tGT:PS\t{gt_donor}:{ps_a}\nchrT2\t40\t.\tA\tG\t.\tPASS\t.\tGT:PS\t{gt_acceptor}:{ps_b}\n"
    vcf = f"{work}/{name}.vcf"
    open(vcf, "w").write(HEADER + rows)
    with open(vcf + ".gz", "wb") as f:
        subprocess.run([bgzip, "-c", vcf], stdout=f, check=True)
    subprocess.run([tabix, "-f", "-p", "vcf", vcf + ".gz"], check=True)
    prefix = f"{work}/{name}"
    subprocess.run([pig, "adjudicate", "--classification", f"{work}/sq.tsv", "--isoforms-gtf", f"{work}/caller.gtf",
                    "--ref-gtf", f"{work}/ref.gtf", "--reference", f"{work}/ref.fa", "--haplotype-vcf", vcf + ".gz",
                    "--haplotype-provenance", "wgs", "--out-prefix", prefix], check=True, capture_output=True)
    r = json.loads(open(prefix + ".attribution.jsonl").readline())
    with open(prefix + ".junctions.tsv") as f:
        j = next(csv.DictReader(f, delimiter="\t"))
    return r["confidence_class"], r["verdict_reason"], j["hap_motifs"], j["variant_phase"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pig"); ap.add_argument("--bgzip", default=shutil.which("bgzip")); ap.add_argument("--tabix", default=shutil.which("tabix"))
    ap.add_argument("--work")
    a = ap.parse_args()
    if not a.bgzip or not a.tabix:
        sys.exit("bgzip and tabix are needed (htslib)")
    work = a.work or tempfile.mkdtemp(prefix="phase_matrix_")
    os.makedirs(work, exist_ok=True)
    write_inputs(work)
    failures, n = [], 0
    for mode in ("same_ps", "different_ps", "unphased"):
        linked = mode == "same_ps"
        for i, ga in enumerate(GENOTYPES):
            for j, gb in enumerate(GENOTYPES):
                gd = ga.replace("|", "/") if mode == "unphased" else ga
                gc = gb.replace("|", "/") if mode == "unphased" else gb
                cls, reason, motifs, phase = run_case(a.pig, a.bgzip, a.tabix, work, f"{mode}_{i}_{j}", gd, gc, mode)
                n += 1
                guaranteed = canonical_guaranteed(ga, gb, linked)
                if (cls == "REFERENCE_BIAS") != guaranteed:
                    failures.append(f"{mode} {gd} {gc}: {cls} ({reason}; motifs {motifs}, phase {phase}) but canonical guaranteed = {guaranteed}")
                if not linked and gc != swapped(gc):
                    cls2, reason2, motifs2, phase2 = run_case(a.pig, a.bgzip, a.tabix, work, f"{mode}_{i}_{j}_swapped", gd, swapped(gc), mode)
                    n += 1
                    if (cls2, reason2) != (cls, reason):
                        failures.append(f"{mode} {gd} {gc} vs {swapped(gc)}: {cls}/{reason} vs {cls2}/{reason2} (motifs {motifs} vs {motifs2})")
    print(f"{n} runs, {len(failures)} failures")
    for f in failures:
        print("  " + f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
