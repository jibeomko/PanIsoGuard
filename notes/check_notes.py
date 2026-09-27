#!/usr/bin/env python3
"""Recompute the study notes' toy verdicts without PanIsoGuard and compare with the binary.

The explanations are in notes/*.md (Korean). This script rebuilds every verdict of the toy
locus (notes/data/) from the raw input files with the Python standard library only, runs
`panisoguard adjudicate` / `combine` on the same inputs, and stops with an error if any
isoform differs in novelty support, mechanism, confidence class or junction counts. It
also checks the headline numbers the notes quote.

    python3 notes/check_notes.py [path/to/panisoguard]      # default: build/panisoguard
"""
import os
import re
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
PIG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "build", "panisoguard")
D = lambda name: os.path.join(DATA, name)

NOVEL = ("novel_in_catalog", "novel_not_in_catalog")
CFG = dict(min_uniq=3, need_canonical=True, min_callers=2, min_mapq=20, clip_bp=20, window=10,
           gates={"low_mapq": 0.5, "supplementary": 0.5, "indel_near": 0.5, "softclip": 1.01})


# ---------------------------------------------------------------- readers (all -> 0-based half-open)
def gtf_chains(path):
    tx = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            tx.setdefault(t, [f[0], f[6], []])[2].append((int(f[3]), int(f[4])))
    out = {}
    for t, (chrom, strand, ex) in tx.items():
        ex.sort()
        out[t] = (chrom, strand, [(ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1)])
    return out


def read_sj(path):
    sj = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        key = (f[0], {"1": "+", "2": "-"}.get(f[3], "."), (int(f[1]) - 1, int(f[2])))
        r = sj.setdefault(key, {"motif": 0, "n_uniq": 0})
        r["n_uniq"] += int(f[6])                                    # duplicate rows are summed
        r["motif"] = max(r["motif"], int(f[4]))
    return sj


def read_fasta(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def read_sam(path):
    reads = []
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.split("\t")
        reads.append((int(f[1]), int(f[3]) - 1, int(f[4]), re.findall(r"(\d+)([MIDNSHP=X])", f[5])))
    return reads


def read_classification(path):
    rows = [l.rstrip("\n").split("\t") for l in open(path)]
    return {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}


# ---------------------------------------------------------------- evidence
def bam_features(reads, jx):
    n = low = supp = clip = indel = 0
    near = lambda p: abs(p - jx[0]) <= CFG["window"] or abs(p - jx[1]) <= CFG["window"]
    for flag, pos, mapq, cig in reads:
        if flag & 4:
            continue
        span = c = d = False
        for i, (ln, op) in enumerate(cig):
            ln = int(ln)
            if op == "S" and ln >= CFG["clip_bp"] and i in (0, len(cig) - 1):
                c = True
            elif op == "N":
                span |= (pos, pos + ln) == jx
                pos += ln
            elif op == "D":
                d |= near(pos)
                pos += ln
            elif op == "I":
                d |= near(pos)
            elif op in "M=X":
                pos += ln
        if span:
            n += 1
            low += mapq < CFG["min_mapq"]
            supp += bool(flag & (256 | 2048))
            clip += c
            indel += d
    return n, {"low_mapq": low, "supplementary": supp, "softclip": clip, "indel_near": indel}


CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}


def canonical(seq, s, e, strand):
    left, right = seq[s:s + 2], seq[e - 2:e]
    if strand == "-":
        comp = str.maketrans("ACGT", "TGCA")
        left, right = right[::-1].translate(comp), left[::-1].translate(comp)
    return (left, right) in CANON


def fingerprint(chrom, strand, introns):
    h = 1469598103934665603
    for b in chrom.encode() + strand.encode() + b"".join(struct.pack("<qq", s, e) for s, e in sorted(introns)):
        h = ((h ^ b) * 1099511628211) % 2 ** 64
    return h


def verdict(iso, row, chain, inp):
    """Return (novelty_support, primary_mechanism, confidence_class, n_novel, n_sr, circular)."""
    catalog, sj, reads = inp.get("catalog"), inp.get("sj"), inp.get("reads")
    n_novel = n_sr = n_created = n_pan = 0
    spanning, worst = 0, dict.fromkeys(CFG["gates"], 0.0)
    if chain and catalog is not None:
        chrom, strand, introns = chain
        for jx in introns:
            if (chrom, strand, jx) in catalog:
                continue
            n_novel += 1
            if "ref_seq" in inp:
                hap_ok = [canonical(h, *jx, strand) for h in inp["haps"]]
                n_created += not canonical(inp["ref_seq"], *jx, strand) and any(hap_ok)
            if "pangenome" in inp:
                n_pan += (chrom, strand, jx) in inp["pangenome"]
            if sj is not None:
                r = sj.get((chrom, strand, jx))
                n_sr += bool(r and r["n_uniq"] >= CFG["min_uniq"] and (r["motif"] >= 1 or not CFG["need_canonical"]))
            if reads is not None:
                n, hits = bam_features(reads, jx)
                spanning += n
                if n:
                    for k in worst:
                        worst[k] = max(worst[k], hits[k] / n)
    cat = row["structural_category"]
    if cat == "full-splice_match":
        return "SUPPORTED", "none", "HIGH_CONF_KNOWN", n_novel, n_sr, False
    if cat == "incomplete-splice_match":
        return "UNKNOWN", "none", "LOW_CONF_PARTIAL", n_novel, n_sr, False
    if cat not in NOVEL:
        return "UNKNOWN", "none", "AMBIGUOUS", n_novel, n_sr, False
    for axis, n_ok, mech in (("pangenome", n_pan, "population_known"), ("ref_seq", n_created, "variant_created")):
        if axis in inp and n_novel and n_ok == n_novel:
            circular = inp[axis + "_circular"]
            return ("UNKNOWN", mech, "AMBIGUOUS" if circular else "PAN_REF_RESCUED_FALSE_NOVEL",
                    n_novel, n_sr, circular)
    if sj is None or not chain or n_novel == 0:
        sup = "UNKNOWN"
    else:
        sup = "SUPPORTED" if n_sr == n_novel else "UNSUPPORTED" if n_sr == 0 else "PARTIAL"
    perc_a = row.get("perc_A_downstream_TTS", "NA")
    if reads is not None and spanning and any(worst[k] > g for k, g in CFG["gates"].items()):
        mech = "mapping_or_repeat"
    elif row.get("all_canonical") == "non_canonical":
        mech = "noncanonical"
    elif row.get("RTS_stage", "").upper() == "TRUE":
        mech = "rt_switch"
    elif perc_a not in ("NA", "") and float(perc_a) >= 60:
        mech = "degradation"
    else:
        mech = "none"
    if sup == "UNKNOWN":
        n_callers = inp.get("callers", {}).get(iso, 0)
        if mech == "mapping_or_repeat":
            cls = "ARTIFACT"
        elif n_callers >= CFG["min_callers"]:
            cls = "MEDIUM_CONF_NOVEL" if mech == "none" else "LOW_CONF_PARTIAL"
        else:
            cls = "AMBIGUOUS"
    elif sup == "SUPPORTED":
        cls = "HIGH_CONF_NOVEL" if mech == "none" else "MEDIUM_CONF_NOVEL"
    elif sup == "PARTIAL":
        cls = "LOW_CONF_PARTIAL"
    else:
        cls = "LOW_CONF_PARTIAL" if mech == "none" else "ARTIFACT"
    return sup, mech, cls, n_novel, n_sr, False


# ---------------------------------------------------------------- runs
def run(args, out):
    subprocess.run([PIG, "adjudicate", *args, "--out-prefix", out], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rows = [l.rstrip("\n").split("\t") for l in open(out + ".adjudicated.tsv")]
    col = {c: i for i, c in enumerate(rows[0])}
    keys = ("novelty_support", "primary_mechanism", "confidence_class", "n_novel_junctions",
            "n_novel_jx_sr_supported", "circularity_flag")
    return {r[0]: (r[col[keys[0]]], r[col[keys[1]]], r[col[keys[2]]], int(r[col[keys[3]]]),
                   int(r[col[keys[4]]]), r[col[keys[5]]] == "true") for r in rows[1:]}


def main():
    cls = read_classification(D("caller1_classification.txt"))
    chains = gtf_chains(D("caller1.gtf"))
    catalog = {(c, s, j) for c, s, js in gtf_chains(D("reference.gtf")).values() for j in js}
    sj, reads = read_sj(D("short_reads.SJ.out.tab")), read_sam(D("long_reads.sam"))
    ref_seq, haps = read_fasta(D("toy.fa")), [read_fasta(D("hap1.fa")), read_fasta(D("hap2.fa"))]

    groups = {}
    for caller in ("caller1", "caller2"):
        for t, (c, s, js) in gtf_chains(D(caller + ".gtf")).items():
            groups.setdefault(fingerprint(c, s, js), set()).add((caller, t))
    callers = {t: len({c for c, _ in g}) for g in groups.values() for c, t in g if c == "caller1"}

    base = ["--classification", D("caller1_classification.txt"), "--isoforms-gtf", D("caller1.gtf")]
    R, S, M = ["--ref-gtf", D("reference.gtf")], ["--sj-tab", D("short_reads.SJ.out.tab")], ["--bam", D("long_reads.bam")]
    H = ["--reference", D("toy.fa"), "--reference-haplotype", D("hap1.fa"), "--reference-haplotype", D("hap2.fa")]
    hap = dict(ref_seq=ref_seq, haps=haps)
    configs = [  # name, CLI args, what the Python side sees
        ("step0", [], {}),
        ("step1", R, dict(catalog=catalog)),
        ("step2", R + S, dict(catalog=catalog, sj=sj)),
        ("step3", R + S + M, dict(catalog=catalog, sj=sj, reads=reads)),
        ("step4", R + S + M + H, dict(catalog=catalog, sj=sj, reads=reads, ref_seq_circular=True, **hap)),
        ("step5", R + S + M + H + ["--haplotype-provenance", "wgs"],
         dict(catalog=catalog, sj=sj, reads=reads, ref_seq_circular=False, **hap)),
        ("pangenome", R + S + ["--pangenome-junctions", None, "--pangenome-provenance", "population"],
         dict(catalog=catalog, sj=sj, pangenome={("chrT", "+", (500, 730))}, pangenome_circular=False)),
        ("consensus", R + ["--caller-support", None], dict(catalog=catalog, callers=callers)),
        ("consensus_bam", R + M + ["--caller-support", None], dict(catalog=catalog, reads=reads, callers=callers)),
    ]
    failures = 0
    with tempfile.TemporaryDirectory() as tmp:
        pan = os.path.join(tmp, "pangenome.tsv")
        open(pan, "w").write("chrT\t501\t730\t+\t12\n")
        matrix = os.path.join(tmp, "matrix.tsv")
        subprocess.run([PIG, "combine", "--gtf", "caller1:" + D("caller1.gtf"), "--gtf", "caller2:" + D("caller2.gtf"),
                        "--ref-gtf", D("reference.gtf"), "--out", matrix], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # An isoform id that does not join must be reported, not silently held AMBIGUOUS (01, 06).
        renamed = os.path.join(tmp, "renamed.gtf")
        open(renamed, "w").write(open(D("caller1.gtf")).read().replace('"iso_A"', '"FLAIR_iso_A"'))
        err = subprocess.run([PIG, "adjudicate", *base[:2], "--isoforms-gtf", renamed,
                              "--out-prefix", os.path.join(tmp, "renamed")], capture_output=True, text=True).stderr
        warned = 'WARNING: 1 of 9 SQANTI3 isoform id(s) (e.g. "iso_A")' in err
        results = {}
        for name, args, inp in configs:
            args = [pan if a is None and name == "pangenome" else matrix if a is None else a for a in args]
            got = run(base + args, os.path.join(tmp, name))
            results[name] = got
            for iso, row in cls.items():
                want = verdict(iso, row, chains.get(iso), inp)
                if got[iso] != want:
                    failures += 1
                    print(f"MISMATCH {name} {iso}: PanIsoGuard {got[iso]} vs by hand {want}")
            print(f"{name:14s} {len(cls)} isoforms, iso_A -> {got['iso_A'][2]}")

    # The headline numbers quoted in the notes.
    step = [results[f"step{i}"]["iso_A"][2] for i in range(6)]
    expect = ["AMBIGUOUS", "AMBIGUOUS", "ARTIFACT", "ARTIFACT", "AMBIGUOUS", "PAN_REF_RESCUED_FALSE_NOVEL"]
    checks = [
        ("iso_A across the six steps (00)", step == expect),
        ("iso_C is PARTIAL -> LOW_CONF_PARTIAL (02)", results["step2"]["iso_C"][:3] == ("PARTIAL", "none", "LOW_CONF_PARTIAL")),
        ("iso_G is a mapping artifact (03)", results["step3"]["iso_G"][1:3] == ("mapping_or_repeat", "ARTIFACT")),
        ("iso_A novel junction absent from SJ.out.tab (02)", ("chrT", "+", (500, 730)) not in sj),
        ("fingerprint of iso_known (06)", fingerprint("chrT", "+", chains["iso_known"][2]) == 4486241254994410266),
        ("unjoined isoform id is warned about (01)", warned),
        ("consensus lifts iso_G, BAM drops it (06)",
         results["consensus"]["iso_G"][2] == "MEDIUM_CONF_NOVEL" and results["consensus_bam"]["iso_G"][2] == "ARTIFACT"),
    ]
    counts = {1.0: (385, 4), 0.75: (45, 1), 0.5: (146, 6), 0.25: (154, 61), 0.1: (0, 8), 0.0: (0, 68)}
    tp = fp = 0
    ap, pos = 0.0, sum(g for g, _ in counts.values())
    for score in sorted(counts, reverse=True):                     # average precision over tied class ranks
        g, f = counts[score]
        tp, fp = tp + g, fp + f
        ap += g / pos * tp / (tp + fp)
    checks.append(("AUPRC 0.9712 from the SQANTI-SIM confusion table (07)", round(ap, 4) == 0.9712))
    for name, ok in checks:
        print(f"{'ok  ' if ok else 'FAIL'} {name}")
        failures += not ok
    if failures:
        sys.exit(f"{failures} check(s) failed")
    print("All checks passed.")


if __name__ == "__main__":
    main()
