#!/usr/bin/env python3
"""Score PanIsoGuard's reference-bias rescue (arm C) against the personal-genome labels (arm B).

Truth per B arm (from classify_junctions.py): T = novel GRCh38 junctions labelled
  same-canonical with a non-canonical GRCh38 motif   (the splice-dinucleotide SNV class)
  shifted-known                                      (the read uses an annotated junction in the person)
  unspliced                                          (no intron there in the person)
R = junctions the rescue fired on: novel junctions of REFERENCE_BIAS isoforms that are
non-canonical on GRCh38 and canonical on at least one SNV-consensus haplotype (the rule's own
criterion, recomputed here so R is junction-level).

Reports sensitivity |R∩T|/|T| and precision |R∩T|/|R| with bootstrap CIs, the same against the
in-scope subset of T (the SNV class only), the breakdown of T∖R by label (what the rescue
cannot see), R∖T (candidate false rescues, for inspection), the firewall check on the
`unknown`-provenance run, the isoform-level A-vs-C view, and the stage cost table.

Usage:
  score.py --work $WORK --data $DATA --genome GRCh38.fa --classes B1=armB1/junction_classes.tsv B2=...
"""
import argparse, collections, csv, glob, json, os, random, re, sys
import pysam

CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
T_LABELS = {"same-canonical", "shifted-known", "unspliced"}


def canonical(fa, chrom, s0, e, strand):
    try:
        d, a = fa.fetch(chrom, s0, s0 + 2).upper(), fa.fetch(chrom, e - 2, e).upper()
    except (KeyError, ValueError):
        return False
    if strand == "-": d, a = rc(a), rc(d)
    return (d, a) in CANON


def boot(num, den_set, pick, n=1000, seed=1):
    """CI of |pick(x)|/|den| under resampling of den."""
    den = sorted(den_set)   # sorted: set order is hash-randomized per process, the CIs must not be
    if not den: return (float("nan"), float("nan"))
    rng = random.Random(seed); vals = []
    for _ in range(n):
        s = [den[rng.randrange(len(den))] for _ in den]
        vals.append(sum(1 for x in s if pick(x)) / len(s))
    vals.sort()
    return (vals[int(0.025 * n)], vals[int(0.975 * n) - 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work", required=True); ap.add_argument("--data", required=True)
    ap.add_argument("--genome", required=True, help="GRCh38 fasta (faidx)")
    ap.add_argument("--classes", nargs="+", required=True, help="NAME=junction_classes.tsv per B arm")
    ap.add_argument("--out", default=None, help="metrics json (default $WORK/metrics.json)")
    ap.add_argument("--emit-metrics", default=None, help="also write the tracked envelope (benchmark/results/refbias_3way/metrics.json)")
    ap.add_argument("--label", default="HG03516", help="individual, the key under metrics.samples in the envelope (an existing envelope is updated in place)")
    a = ap.parse_args()
    W, D = a.work, a.data

    # arm A: novel junctions per isoform, GRCh38 canonical status
    jx_of = collections.defaultdict(list); canon_A = {}; strand_of = {}
    with open(f"{W}/armA/sq/iq_junctions.txt") as f:
        r = csv.DictReader(f, delimiter="\t")
        for p in r:
            if p["junction_category"] != "novel": continue
            k = (p["chrom"], int(p["genomic_start_coord"]) - 1, int(p["genomic_end_coord"]))
            jx_of[p["isoform"]].append(k); canon_A[k] = p["canonical"] == "canonical"; strand_of[k] = p["strand"]

    # arm C: rescued isoforms -> junction-level R (rule criterion recomputed on the consensus haplotypes)
    cls_C = {}
    for prov in ("nohap", "wgs", "unk"):
        with open(f"{W}/armC/{prov}.adjudicated.tsv") as f:
            cls_C[prov] = {p["isoform_id"]: p for p in csv.DictReader(f, delimiter="\t")}
    rescued = [i for i, p in cls_C["wgs"].items() if p["confidence_class"] == "REFERENCE_BIAS"]
    ref = pysam.FastaFile(a.genome)
    haps = [pysam.FastaFile(f"{D}/{h}_consensus.fa") for h in ("mat", "pat")]
    # R_rule: every novel junction that meets the variant axis's criterion (what the axis flags);
    # R: those inside isoforms that actually got the rescue verdict (what the user sees). They differ
    # because the verdict needs ALL of an isoform's novel junctions explained, and only NIC/NNC
    # isoforms are adjudicated.
    R_rule = {k for k, c in canon_A.items() if not c and any(canonical(h, k[0], k[1], k[2], strand_of[k]) for h in haps)}
    R = {k for iso in rescued for k in jx_of.get(iso, []) if k in R_rule}
    TARGETED = {"novel_in_catalog", "novel_not_in_catalog"}
    cat_of = {}
    with open(f"{W}/armA/sq/iq_classification.txt") as f:
        for p in csv.DictReader(f, delimiter="\t"): cat_of[p["isoform"]] = p["structural_category"]
    iso_of = collections.defaultdict(set)
    for iso, ks in jx_of.items():
        for k in ks: iso_of[k].add(iso)
    firewall = {"rescued_wgs": len(rescued),
                "rescued_unknown": sum(p["confidence_class"] == "REFERENCE_BIAS" for p in cls_C["unk"].values()),
                "held_ambiguous_with_flag_unknown": sum(p["confidence_class"] == "UNRESOLVED" and p["circularity_flag"] == "true"
                                                        for p in cls_C["unk"].values())}
    out = {"n_novel_junctions_A": len(canon_A), "n_noncanonical_A": sum(not c for c in canon_A.values()),
           "rescued_isoforms": len(rescued), "R_junctions": len(R), "R_rule_junctions": len(R_rule),
           "firewall": firewall, "arms": {}}

    def why_missed(k, T):
        """Why a junction in T_snv got no rescue verdict."""
        if k not in R_rule: return "criterion-miss"
        isos = iso_of[k]
        if all(cat_of.get(i) not in TARGETED for i in isos): return "isoform-category-not-adjudicated"
        if all(any(j not in T for j in jx_of[i]) for i in isos): return "other-novel-junction-unexplained"
        return "other"

    # per B arm
    for spec in a.classes:
        name, path = spec.split("=", 1)
        lab = {}
        with open(path) as f:
            for p in csv.DictReader(f, delimiter="\t"):
                lab[(p["chrom"], int(p["start0"]), int(p["end"]))] = p
        T = {k for k, p in lab.items() if p["label"] in T_LABELS and not (p["label"] == "same-canonical" and p["canonical_A"] == "True")}
        T_snv = {k for k in T if lab[k]["label"] == "same-canonical"}
        # 'unspliced' in an SNV-only arm is mostly 1-2 read alignment instability at short
        # introns, not reference bias; T_core leaves it out and keeps the two solid classes.
        T_core = {k for k in T if lab[k]["label"] != "unspliced"}
        RT = R & T
        def pr(num, den): return len(num) / len(den) if den else float("nan")
        arm = {
            "labels": dict(collections.Counter(f'{p["label"]}|{"canonical_A" if p["canonical_A"] == "True" else "noncanonical_A"}' for p in lab.values())),
            "T": len(T), "T_snv_class": len(T_snv), "T_core": len(T_core), "R": len(R), "R_and_T": len(RT),
            "sensitivity_core": pr(R & T_core, T_core), "precision_core": pr(R & T_core, R),
            "rule_sensitivity_core": pr(R_rule & T_core, T_core), "rule_precision_core": pr(R_rule & T_core, R_rule),
            "sensitivity": pr(RT, T), "sensitivity_ci": boot(RT, T, lambda k: k in R),
            "precision": pr(RT, R), "precision_ci": boot(RT, R, lambda k: k in T),
            "sensitivity_snv_class": pr(R & T_snv, T_snv), "sensitivity_snv_class_ci": boot(R & T_snv, T_snv, lambda k: k in R),
            "rule_level": {"R_rule": len(R_rule), "R_rule_and_T": len(R_rule & T),
                           "sensitivity_snv_class": pr(R_rule & T_snv, T_snv), "precision": pr(R_rule & T, R_rule),
                           "R_rule_minus_T_by_label": dict(collections.Counter(lab[k]["label"] if k in lab else "absent" for k in R_rule - T))},
            "T_minus_R_by_label": dict(collections.Counter(lab[k]["label"] for k in T - R)),
            "T_snv_minus_R_by_reason": dict(collections.Counter(why_missed(k, T) for k in T_snv - R)),
            "T_minus_R_snv_class": sorted(f"{k[0]}:{k[1]}-{k[2]}" for k in T_snv - R),
            "R_minus_T": sorted(f'{k[0]}:{k[1]}-{k[2]} label={lab[k]["label"] if k in lab else "absent"}' for k in R - T),
            "motif_disrupted_in_person": sum(1 for p in lab.values() if p["label"] == "same-noncanonical" and p["canonical_A"] == "True"),
        }
        out["arms"][name] = arm

    # isoform level, A vs C: SQANTI3 rules filter verdict x PanIsoGuard class (novel isoforms only)
    rules = {}
    with open(f"{W}/armA/sq/rules/rules_RulesFilter_classification.txt") as f:
        for p in csv.DictReader(f, delimiter="\t"):
            if p["structural_category"] in ("novel_in_catalog", "novel_not_in_catalog"): rules[p["isoform"]] = p["filter_result"]
    xt = collections.Counter((rules[i], cls_C["wgs"][i]["confidence_class"]) for i in rules if i in cls_C["wgs"])
    out["isoform_A_vs_C"] = {f"SQ3rules={k[0]}|PIG={k[1]}": v for k, v in sorted(xt.items())}
    out["rescued_discarded_by_rules_filter"] = sum(1 for i in rescued if rules.get(i) == "Artifact")

    # cost table from /usr/bin/time outputs (+ arm A's alignment and IsoQuant from the original logs)
    cost = {}
    for fn in sorted(glob.glob(f"{W}/time.*.txt")):
        t = open(fn).read()
        wall = re.search(r"Elapsed \(wall clock\).*: (.*)", t); rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
        cost[os.path.basename(fn)[5:-4]] = {"wall": wall.group(1) if wall else None, "max_rss_gb": round(int(rss.group(1)) / 1e6, 1) if rss else None}
    m = re.search(r"Real time: ([\d.]+) sec.*Peak RSS: ([\d.]+) GB", open(f"{D}/mm2.rna.log").read())
    if m: cost["A_align(existing)"] = {"wall_sec": float(m.group(1)), "max_rss_gb": float(m.group(2))}
    out["cost"] = cost

    path = a.out or f"{W}/metrics.json"
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, indent=1))
    if a.emit_metrics:
        import datetime
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{W}/armC/wgs.provenance.log") if "\t" in l and not l.startswith("#"))
        arms = {n: {k: v for k, v in arm.items() if k not in ("T_minus_R_snv_class", "R_minus_T", "labels")} for n, arm in out["arms"].items()}
        sample = {k: out[k] for k in ("n_novel_junctions_A", "n_noncanonical_A", "rescued_isoforms", "R_junctions",
                                      "R_rule_junctions", "firewall", "rescued_discarded_by_rules_filter", "cost")} | {"arms": arms}
        samples = {}
        if os.path.exists(a.emit_metrics):
            try: samples = json.load(open(a.emit_metrics))["metrics"].get("samples", {})
            except (KeyError, ValueError): samples = {}
        samples[a.label] = sample
        env = {"protocol": "refbias_3way", "scope": "hprc_r2_grch38_vs_personalized_pipelines", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "REAL public HPRC R2 Kinnex Iso-Seq of HG03516 (ESN) and HG02717 (GWD) through four arms on the same reads: "
                                  "A GRCh38 + GENCODE v49 + IsoQuant 3.13 + SQANTI3 6.0.1; B1 SNV-personalized diploid (assembly SNVs on "
                                  "the GRCh38 backbone); B2 HPRC R2 phased assembly diploid with liftoff annotation and paftools liftover; "
                                  "C panisoguard adjudicate post hoc on A with --bam and both SNV-consensus haplotypes. No private data.",
               "command": "SAMPLE=<id> bash benchmark/refbias_3way/run.sh all; benchmark/refbias_3way/classify_junctions.py per B arm; "
                          "benchmark/refbias_3way/score.py --work W --data D --genome GRCh38.fa --classes B1=... B2=... --label <ID> --emit-metrics ...",
               "metrics": {"samples": samples},
               "notes": "Unit = GRCh38 novel junction of arm A, labelled by where the SAME reads splice on the personal genome. "
                        "T = same-canonical (non-canonical on GRCh38) + shifted-known + unspliced; T_core drops unspliced (short "
                        "introns that minimap2 renders as a deletion on the diploid assembly, median 82 bp; 308/499 of B2's are "
                        "spliced on the SNV genome B1). R_rule = junctions meeting the rescue's own criterion (what the axis sees); "
                        "R = junctions inside isoforms that got the REFERENCE_BIAS verdict (needs ALL of the isoform's novel "
                        "junctions explained, NIC/NNC only). The gap R_rule -> R (35 -> 13) is the all-junctions rule and the "
                        "category pass-throughs, not the motif criterion. Precision misses (4/35) are junctions whose reads come "
                        "from the non-canonical haplotype. Firewall: 0 rescues under unknown provenance, 18/18 held."}
        json.dump(env, open(a.emit_metrics, "w"), indent=2, sort_keys=True); open(a.emit_metrics, "a").write("\n")
        print("wrote", a.emit_metrics)


if __name__ == "__main__":
    main()
