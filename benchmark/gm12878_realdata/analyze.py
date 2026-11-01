#!/usr/bin/env python3
"""Real-data robustness of the BAM mapping axis + reference-bias rescue on GM12878 ONT.

Consumes (1) the adjudication of the 285 GM12878 caller-novel chains run WITH the real
ENCODE ONT direct-RNA BAM, (2) the chain -> n_callers map (consensus), and (3) the
reference-bias scan summary (created/scanned), and emits the metrics envelope.

Usage: analyze.py <rb.attribution.jsonl> <pig_ncallers.json> <created:int> <scanned:int>
                  [--ont rb_ont.attribution.jsonl] [--wobble chain_wobble.json] [--emit-metrics out]
"""
import collections, json, statistics, sys

attr, ncj, created, scanned = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
def opt(flag): return sys.argv[sys.argv.index(flag) + 1] if flag in sys.argv else None
emit, ont_attr, wobble = opt("--emit-metrics"), opt("--ont"), opt("--wobble")
nc = json.load(open(ncj))


def class_counts(path):
    return dict(collections.Counter(json.loads(l)["confidence_class"] for l in open(path)))


def chain_stats(path):
    rows = [json.loads(l) for l in open(path)]
    full = [r["evidence"].get("chain_n_full_reads") or 0 for r in rows]
    by = {}
    for name, pred in (("single_caller", lambda n: n < 2), ("ge2_caller", lambda n: n >= 2)):
        sel = [r for r in rows if pred(nc.get(r["isoform"], 1))]
        by[name] = {"n": len(sel), "chain_full_ge3": sum((r["evidence"].get("chain_n_full_reads") or 0) >= 3 for r in sel),
                    "NOVEL_PROBABLE": sum(r["confidence_class"] == "NOVEL_PROBABLE" for r in sel),
                    "mapping_sign_held": sum(r["primary_mechanism"] == "mapping_or_repeat" for r in sel)}
    return {"chain_full_reads_median": statistics.median(full), "chains_full_ge1": sum(x >= 1 for x in full),
            "chains_full_ge3": sum(x >= 3 for x in full), "by_consensus": by}

R = []
for line in open(attr):
    d = json.loads(line); ev = d["evidence"]
    R.append((nc.get(d["isoform"], 1), ev["bam_evaluable"], ev["bam_n_spanning"],
              ev["bam_frac_indel_near"], ev["bam_frac_low_mapq"], ev["bam_frac_supplementary"],
              d["primary_mechanism"]))
sp = [r for r in R if r[1] and r[2] > 0]
def cnt(idx, thr): return sum(1 for r in sp if r[idx] > thr)
def mean(idx, pred): 
    v = [r[idx] for r in sp if pred(r[0])]; return round(sum(v)/len(v), 3) if v else None

metrics = {
    "n_novel_chains": len(R),
    "n_bam_spanning_evaluable": len(sp),
    "bam_signal_fire_rate_on_real_ONT": {
        "indel_near_gt_0.5":    [cnt(3, 0.5), len(sp)],
        "low_mapq_gt_0.5":      [cnt(4, 0.5), len(sp)],
        "supplementary_gt_0.5": [cnt(5, 0.5), len(sp)],
    },
    "indel_near_mean_by_consensus": {
        "single_caller": mean(3, lambda n: n == 1),
        "ge2_caller":    mean(3, lambda n: n >= 2),
    },
    # the mapping sign under the HiFi defaults (UNKNOWN support here, so the sign alone decides the
    # class: ARTIFACT before ruleset 0.0.5, UNCONFIRMED / mapping_without_short_reads since)
    "mapping_sign_rate_by_consensus": {
        "single_caller": [sum(1 for r in R if r[0] == 1 and r[6] == "mapping_or_repeat"), sum(1 for r in R if r[0] == 1)],
        "ge2_caller":    [sum(1 for r in R if r[0] >= 2 and r[6] == "mapping_or_repeat"), sum(1 for r in R if r[0] >= 2)],
    },
    "reference_bias_rescue": {"created_reference_bias": created, "novel_junctions_scanned": scanned},
    "classes_hifi_preset": class_counts(attr),
}
if ont_attr:
    metrics["classes_ont_preset"] = class_counts(ont_attr)
    metrics["chain_support_on_real_ONT"] = chain_stats(ont_attr)
if wobble:
    metrics["chain_support_on_real_ONT"]["splice_site_tolerance"] = json.load(open(wobble))
print(json.dumps(metrics, indent=2))
if emit:
    prov = dict(l.rstrip("\n").split("\t", 1) for l in open(attr.replace(".attribution.jsonl", ".provenance.log"))
                if "\t" in l and not l.startswith("#"))
    fr = metrics["bam_signal_fire_rate_on_real_ONT"]
    cs = metrics["chain_support_on_real_ONT"]
    hifi = metrics.get("classes_hifi_preset", {})
    ont = metrics.get("classes_ont_preset", {})
    wob_ge3 = (cs.get("splice_site_tolerance") or {}).get("by_tolerance_bp", {}).get("0", {}).get("chains_ge3", cs["chains_full_ge3"])
    out = {
        "protocol": "gm12878_realdata",
        "scope": "gm12878_ont_drna_real",
        "status": "tracked",
        "tool_version": prov.get("tool_version"),
        "ruleset_version": prov.get("ruleset_version"),
        "data_provenance": "real public ENCODE GM12878 ONT direct-RNA (ENCFF440ZML, minimap2 -ax splice -uf); "
                           "285 caller-novel chains (IsoQuant+Bambu+ESPRESSO combine); reference-bias scan vs "
                           "NA12878=HG001 GIAB v4.2.1 SNV-consensus haplotype. Public data only, no private data.",
        "command": "see run.sh: build novel GTF + minimal classification; adjudicate --bam <ONT bam>; "
                   "scan_variant_axis.py over the novel junctions vs the HG001 haplotype; analyze.py",
        "metrics": metrics,
        "notes": (
            "Real-data robustness of the optional BAM mapping axis. (1) The terminal soft-clip fraction, never a "
            "default gate and removed in 0.0.8, fired on 134/191 (70%) of junctions on real ONT dRNA in the last run "
            "that measured it (adapter/poly-A read ends). (2) The HiFi-calibrated 0.5 gates are chemistry-dependent: "
            f"ONT's intrinsic indel error makes indel_near fire on {fr['indel_near_gt_0.5'][0]}/{fr['indel_near_gt_0.5'][1]} "
            f"({round(100 * fr['indel_near_gt_0.5'][0] / fr['indel_near_gt_0.5'][1])}%) and it does not separate single "
            f"({metrics['indel_near_mean_by_consensus']['single_caller']}) from >=2-caller "
            f"({metrics['indel_near_mean_by_consensus']['ge2_caller']}) chains; the older low_mapq/supplementary signals "
            f"barely fire ({fr['low_mapq_gt_0.5'][0]}/{fr['low_mapq_gt_0.5'][1]}, {fr['supplementary_gt_0.5'][0]}/"
            f"{fr['supplementary_gt_0.5'][1]}). The mapping axis is opt-in (--bam) and its thresholds need re-calibration "
            "per chemistry. (3) The sequence-level consensus axis separates the same data independently of chemistry "
            f"({cs['by_consensus']['single_caller']['n']} single vs {cs['by_consensus']['ge2_caller']['n']} reproducible). "
            f"(4) Reference-bias rescue fires {created}/{scanned}: GM12878=NA12878 is GIAB reference-grade, so the base "
            "rate is ~0; the rescue's value lies with divergent / non-reference genomes (cf. giab_cohort_rescue: 137 "
            "across 4 personalized genomes). Small-N caveat on >=2-caller (31). (5) 0.0.7 whole-chain read support works "
            "on real ONT dRNA with exact coordinates: 94.5% of splice sites among near-matching reads are exact, "
            f"{cs['chains_full_ge3']}/{metrics['n_novel_chains']} chains have >= 3 exact full-chain reads on the "
            f"isoform's strand (chain_wobble.py counts {wob_ge3} without the strand check; a +-5 bp tolerance would add "
            f"about 20), {round(100 * cs['by_consensus']['ge2_caller']['chain_full_ge3'] / cs['by_consensus']['ge2_caller']['n'])}% of "
            f">=2-caller chains vs {round(100 * cs['by_consensus']['single_caller']['chain_full_ge3'] / cs['by_consensus']['single_caller']['n'])}% "
            "of single-caller ones; no tolerance is needed. With --platform ont the indel-near gate is off and these "
            f"chains become NOVEL_PROBABLE ({ont.get('NOVEL_PROBABLE', 0)}) instead of being held UNCONFIRMED by the "
            f"mapping sign ({hifi.get('UNCONFIRMED', 0)} under the HiFi defaults; ARTIFACT before ruleset 0.0.5). No truth "
            "here: this is robustness, not accuracy."),
    }
    json.dump(out, open(emit, "w"), indent=2, sort_keys=True)
    print(f"\nwrote {emit}", file=sys.stderr)
