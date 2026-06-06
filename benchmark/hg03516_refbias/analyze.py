#!/usr/bin/env python3
"""Consolidate the HG03516 reference-bias scan + adjudication into the metrics envelope.

Audited result (an independent 5-angle adversarial audit corrected the raw numbers):
drops one displaced-coordinate hit, uses exact-coordinate read support, and frames the
GM12878 contrast as directional (not a controlled rate). See README.md for caveats.

Usage: analyze.py <work_dir> <gencode.gtf> <rna.aln.bam> [--emit-metrics out.json]
  work_dir holds refbias_unique.tsv + scan/<chrom>.<hap>/{vt.truth,vt.gtf,adj_wgs.*,adj_unk.*}
"""
import sys, os, re, glob, subprocess, statistics, json
from collections import defaultdict

WORK, GENCODE, BAM = sys.argv[1], sys.argv[2], sys.argv[3]
EMIT = sys.argv[5] if len(sys.argv) > 5 and sys.argv[4] == "--emit-metrics" else None
SAM = os.environ.get("SAMTOOLS", "samtools")
# One hit excluded by the audit: a 4bp-displaced IsoQuant model coordinate with 0 exact
# spanning reads (its reads splice 4bp away to an already-canonical GRCh38 site).
DROP = ("chr17", 3681978, 3685514)

rows = [l.strip().split('\t') for l in open(f"{WORK}/refbias_unique.tsv")][1:]
rows = [r for r in rows if not (r[0] == DROP[0] and int(r[1]) == DROP[1] and int(r[2]) == DROP[2])]
hom = sum(1 for r in rows if r[4] == "mat,pat")

# GENCODE v49 intron set (exact)
ex = defaultdict(list); chrm = {}; tid = re.compile(r'transcript_id "([^"]+)"')
for ln in open(GENCODE):
    if ln.startswith('#'): continue
    f = ln.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] != 'exon': continue
    m = tid.search(f[8])
    if m: ex[m.group(1)].append((int(f[3]), int(f[4]))); chrm[m.group(1)] = f[0]
gintr = defaultdict(set)
for t, exs in ex.items():
    e = sorted(exs)
    for i in range(len(e) - 1): gintr[chrm[t]].add((e[i][1], e[i + 1][0] - 1))
novel = sum(1 for r in rows if (int(r[1]), int(r[2])) not in gintr.get(r[0], set()))

# exact-coordinate spanning Iso-Seq reads (primary only)
cig = re.compile(r'(\d+)([MIDNSHP=X])')
def span(chrom, a, b):
    out = subprocess.run([SAM, "view", BAM, f"{chrom}:{max(1,a-1)}-{b+1}"], capture_output=True, text=True).stdout
    n = 0
    for ln in out.splitlines():
        f = ln.split('\t')
        if len(f) < 6 or (int(f[1]) & 0x900): continue
        rp = int(f[3]) - 1
        for L, op in cig.findall(f[5]):
            L = int(L)
            if op == 'N' and rp == a and rp + L == b: n += 1; break
            if op in 'MDN=X': rp += L
    return n
sup = [span(r[0], int(r[1]), int(r[2])) for r in rows]

# denominator: unique novel junctions scanned
uniqj = set()
for g in glob.glob(f"{WORK}/scan/chr*.novel.gtf"):
    chrom = g.split('/')[-1].split('.')[0]; e2 = defaultdict(list)
    for ln in open(g):
        f = ln.rstrip('\n').split('\t'); iso = f[8].split('transcript_id "')[1].split('"')[0]
        e2[iso].append((int(f[3]), int(f[4])))
    for iso, e in e2.items():
        e = sorted(e)
        for i in range(len(e) - 1): uniqj.add((chrom, e[i][1], e[i + 1][0] - 1))
D = len(uniqj)

# rescue (wgs) + firewall (unknown) from the axis-engaged per-chrom adjudications
def tally(tag):
    resc = fr = held = 0
    for d in glob.glob(f"{WORK}/scan/chr*.pat") + glob.glob(f"{WORK}/scan/chr*.mat"):
        adj = f"{d}/adj_{tag}.adjudicated.tsv"
        if not os.path.exists(adj): continue
        truth = {l.split('\t')[0]: l.strip().split('\t')[1] for l in open(f"{d}/vt.truth")}
        exj = defaultdict(list); chrom = d.split('/')[-1].split('.')[0]
        for ln in open(f"{d}/vt.gtf"):
            f = ln.rstrip('\n').split('\t'); iso = f[8].split('transcript_id "')[1].split('"')[0]
            exj[iso].append((int(f[3]), int(f[4])))
        jc = {i: (chrom, sorted(v)[0][1], sorted(v)[1][0] - 1) for i, v in exj.items() if len(v) >= 2}
        hdr = None
        for ln in open(adj):
            f = ln.rstrip('\n').split('\t')
            if hdr is None: hdr = {c: i for i, c in enumerate(f)}; continue
            iso = f[hdr['isoform_id']]; cls = f[hdr['confidence_class']]; lab = truth.get(iso)
            if jc.get(iso) == DROP: continue
            if tag == "wgs":
                if lab == "CREATED" and cls == "PAN_REF_RESCUED_FALSE_NOVEL": resc += 1
                if lab in ("DISRUPTED", "CONTROL_canon") and cls == "PAN_REF_RESCUED_FALSE_NOVEL": fr += 1
            elif lab == "CREATED" and cls == "AMBIGUOUS": held += 1
    return resc, fr, held
rescued_rows, false_resc, _ = tally("wgs")
_, _, held_rows = tally("unk")
lam = 1364 * len(rows) / D
import math
m = {
    "individual": "HG03516 (ESN, Esan in Nigeria; AFR superpopulation; HPRC Release 2 LCL)",
    "rna": "PacBio Kinnex Iso-Seq (flnc), 11.34M full-length reads aligned (minimap2 splice:hq)",
    "novel_transcript_models_isoquant": 37568,
    "unique_novel_junctions_scanned": D,
    "reference_bias_junctions": len(rows),
    "novel_vs_gencode_v49": novel,
    "annotated_noncanonical": len(rows) - novel,
    "homozygous_both_haplotypes": hom,
    "heterozygous_one_haplotype": len(rows) - hom,
    "all_trace_to_assembly_snv": True,
    "read_support_exact": {"ge2_spanning_reads": sum(1 for x in sup if x >= 2), "of": len(sup),
                            "median": int(statistics.median(sup)), "max": max(sup), "min": min(sup)},
    "rescue_independent_dna_provenance": {"created_rows_rescued": rescued_rows,
                                           "unique_junctions": len(rows), "false_rescues": false_resc},
    "firewall_circular_risk_provenance": {"created_rows_held_ambiguous": held_rows},
    "per_junction_rate": round(len(rows) / D, 6),
    "contrast_gm12878_directional": {"reference_bias_junctions": 0, "junctions_scanned": 1364,
                                      "expected_at_hg03516_rate": round(lam, 2), "p_observe_zero": round(math.exp(-lam), 3)},
}
print(json.dumps(m, indent=2))
if EMIT:
    out = {
        "protocol": "hg03516_refbias", "scope": "hg03516_esn_afr_real_longread",
        "status": "tracked", "tool_version": "0.0.3", "ruleset_version": "builtin-0.0.1",
        "data_provenance": "REAL public HPRC Release 2 HG03516 (ESN/Nigeria, AFR LCL): PacBio Kinnex Iso-Seq "
            "(s3://human-pangenomics .../HG03516/raw_data/PacBio_Kinnex/*.flnc.bam) + HiFi de-novo phased "
            "assembly (pat/mat, GCA_018469415.2/GCA_018469425.2). Personal SNV-consensus haplotype = assembly "
            "minimap2 asm5 + paftools.js call vs GRCh38, SNVs only (length-matched). GENCODE v49. No private data.",
        "command": "see run.sh: download flnc+assemblies; asm5+paftools->SNV-consensus haplotypes; minimap2 "
            "splice:hq + IsoQuant -> novel junctions; scan_variant_axis.py per chrom vs pat+mat haplotypes; "
            "adjudicate --reference --reference-haplotype --haplotype-provenance {wgs,unknown}; analyze.py",
        "metrics": m,
        "notes": "POSITIVE real-data application: PanIsoGuard's reference-bias rescue FIRES on real long-read RNA "
            "from a divergent (African, ESN) individual -- 40 caller-reported novel splice junctions are reference "
            "bias (non-canonical on GRCh38, canonical on HG03516's own HiFi-assembly haplotype via a real personal "
            "SNV), 35/40 novel vs GENCODE (a reference-only pipeline would report them as discoveries). All 40 "
            "rescued (PAN_REF_RESCUED_FALSE_NOVEL) under independent-DNA provenance with 0 false rescues; the "
            "circularity firewall holds all under circular-risk provenance. Independently adversarially audited (5 "
            "angles); REQUIRED CAVEATS, all stated honestly: (1) one hit (chr17:3681978-3685514, P2RX5) was DROPPED "
            "-- a 4bp-displaced IsoQuant model coordinate with 0 exact spanning reads whose reads splice to an "
            "already-canonical GRCh38 site. (2) The HG03516-vs-GM12878(0) contrast is DIRECTIONAL, NOT a controlled "
            "rate: ~72x denominator asymmetry (98,672 vs 1,364 junctions), different chemistry (PacBio vs ONT), "
            "caller set (IsoQuant vs multi-caller), and variant source (HiFi assembly genome-wide vs GIAB v4.2.1 "
            "high-confidence SNPs excluding segdups); at HG03516's rate GM12878's 0 is expected (E=0.55, P(0)=0.58). "
            "The apples-to-apples 'scales with divergence' evidence is giab_cohort_rescue (137 across 4 genomes). "
            "(3) 'Independent-DNA provenance' = HiFi de-novo ASSEMBLY variant calls (paftools.js, GT 1/1 + QUAL 60 "
            "placeholders, SNV-only so indel-mediated reference bias is invisible), not short-read WGS genotyping; "
            "the firewall's independence is MOLECULE independence (DNA vs RNA), which holds, but both arms share the "
            "GRCh38 + minimap2 backbone. Mitigations: 40/40 independent exact read support (median 5); audit found "
            "0/40 hits in segdup/paralog/low-complexity regions and the restoring SNVs collinear on one assembly contig.",
    }
    json.dump(out, open(EMIT, "w"), indent=2, sort_keys=True)
    print(f"\nwrote {EMIT}", file=sys.stderr)
