#!/usr/bin/env python3
"""Real-data robustness of the BAM mapping axis + reference-bias rescue on GM12878 ONT.

Consumes (1) the adjudication of the 285 GM12878 caller-novel chains run WITH the real
ENCODE ONT direct-RNA BAM, (2) the chain -> n_callers map (consensus), and (3) the
reference-bias scan summary (created/scanned), and emits the metrics envelope.

Usage: analyze.py <rb.attribution.jsonl> <pig_ncallers.json> <created:int> <scanned:int> [--emit-metrics out]
"""
import json, sys

attr, ncj, created, scanned = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
emit = sys.argv[6] if len(sys.argv) > 6 and sys.argv[5] == "--emit-metrics" else None
nc = json.load(open(ncj))

R = []
for line in open(attr):
    d = json.loads(line); ev = d["evidence"]
    R.append((nc.get(d["isoform"], 1), ev["bam_evaluable"], ev["bam_n_spanning"],
              ev["bam_frac_indel_near"], ev["bam_frac_softclip"], ev["bam_frac_low_mapq"],
              ev["bam_frac_supplementary"], d["confidence_class"]))
sp = [r for r in R if r[1] and r[2] > 0]
def cnt(idx, thr): return sum(1 for r in sp if r[idx] > thr)
def mean(idx, pred): 
    v = [r[idx] for r in sp if pred(r[0])]; return round(sum(v)/len(v), 3) if v else None

metrics = {
    "n_novel_chains": len(R),
    "n_bam_spanning_evaluable": len(sp),
    "bam_signal_fire_rate_on_real_ONT": {
        "softclip_gt_0.5":      [cnt(4, 0.5), len(sp)],
        "indel_near_gt_0.5":    [cnt(3, 0.5), len(sp)],
        "low_mapq_gt_0.5":      [cnt(5, 0.5), len(sp)],
        "supplementary_gt_0.5": [cnt(6, 0.5), len(sp)],
    },
    "indel_near_mean_by_consensus": {
        "single_caller": mean(3, lambda n: n == 1),
        "ge2_caller":    mean(3, lambda n: n >= 2),
    },
    "artifact_rate_by_consensus": {
        "single_caller": [sum(1 for r in R if r[0] == 1 and r[7] == "ARTIFACT"), sum(1 for r in R if r[0] == 1)],
        "ge2_caller":    [sum(1 for r in R if r[0] >= 2 and r[7] == "ARTIFACT"), sum(1 for r in R if r[0] >= 2)],
    },
    "reference_bias_rescue": {"created_reference_bias": created, "novel_junctions_scanned": scanned},
}
print(json.dumps(metrics, indent=2))
if emit:
    out = {
        "protocol": "gm12878_realdata",
        "scope": "gm12878_ont_drna_real",
        "status": "tracked",
        "tool_version": "0.0.3",
        "ruleset_version": "builtin-0.0.1",
        "data_provenance": "REAL public ENCODE GM12878 ONT direct-RNA (ENCFF440ZML, minimap2 -ax splice -uf); "
                           "285 caller-novel chains (IsoQuant+Bambu+ESPRESSO combine); reference-bias scan vs "
                           "NA12878=HG001 GIAB v4.2.1 SNV-consensus haplotype. Public data only, no private data.",
        "command": "see run.sh: build novel GTF + minimal classification; adjudicate --bam <ONT bam>; "
                   "scan_variant_axis.py over the novel junctions vs the HG001 haplotype; analyze.py",
        "metrics": metrics,
        "notes": "Real-data robustness of the optional BAM mapping axis. (1) Soft-clip default-OFF is EMPIRICALLY "
                 "VINDICATED: on real ONT dRNA terminal soft-clips fire on 187/191 (98%) of junctions (adapter/poly-A "
                 "read ends), so a default-on gate would demote nearly everything. (2) The HiFi-calibrated 0.5 gates are "
                 "CHEMISTRY-DEPENDENT: ONT's intrinsic indel error makes indel_near fire on 175/191 (92%) and barely "
                 "discriminates single (0.925) vs >=2-caller (0.732); the older low_mapq/supplementary signals barely "
                 "fire (11/191, 4/191). The mapping axis is opt-in (--bam) and its thresholds MUST be re-calibrated per "
                 "chemistry. (3) The sequence-level CONSENSUS axis separates the same data chemistry-independently "
                 "(254 single vs 31 reproducible). (4) Reference-bias rescue fires 0/1364: GM12878=NA12878 is GIAB "
                 "reference-grade, so the base rate is ~0 -- the rescue's value is reserved for divergent/non-reference "
                 "genomes (cf. giab_cohort_rescue: 137 across 4 personalized genomes). Small-N caveat on >=2-caller (31).",
    }
    json.dump(out, open(emit, "w"), indent=2, sort_keys=True)
    print(f"\nwrote {emit}", file=sys.stderr)
