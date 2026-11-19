#!/usr/bin/env python3
"""Score PanIsoGuard's reference-bias rescue (arm C) against the personal-genome labels (arm B).

Truth per B arm (from classify_junctions.py): T = novel GRCh38 junctions labelled
  same-canonical with a non-canonical GRCh38 motif   (the splice-dinucleotide SNV class)
  shifted-known                                      (the read uses an annotated junction in the person)
  unspliced                                          (no intron there in the person)
R = junctions the rescue fired on: novel junctions of REFERENCE_BIAS isoforms that are
non-canonical on GRCh38 and canonical on at least one SNV-consensus haplotype (the rule's own
criterion, recomputed here so R is junction-level).

With --vcf-arm (the run that reads every call of the assembly, insertions and deletions included,
from one haploid VCF per haplotype) the same T is scored against what the tool itself reports in
its junction table: V_rule = novel junctions it marks created or deleted and that the junction
reads do not contradict, V = those inside REFERENCE_BIAS isoforms. --call-set NAME=PREFIX scores
further runs the same way (a call set that does not come from the assembly: short-read WGS, with
and without its phasing) and says which of the assembly arm's confirmed junctions each one finds.

Reports sensitivity |R and T|/|T| and precision |R and T|/|R| with bootstrap CIs, the same against the
in-scope subset of T (the SNV class only), the breakdown of T minus R by label (what the rescue
cannot see), R minus T (candidate false rescues, for inspection), the firewall check on the
`unknown`-provenance run, the isoform-level A-vs-C view, and the stage cost table.

Usage:
  score.py --work $WORK --data $DATA --genome GRCh38.fa --classes B1=armB1/junction_classes.tsv B2=...
"""
import argparse, collections, csv, glob, gzip, json, os, random, re, sys
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
    ap.add_argument("--vcf-arm", default=None, help="output prefix of the --haplotype-vcf run (armC_vcf/wgs)")
    ap.add_argument("--call-set", action="append", default=[], help="NAME=PREFIX[,VCF] of a --haplotype-vcf run with another call set (repeatable); with the "
                                                                    "VCF it read, the deletion junctions it misses are looked up in that file")
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

    # a --haplotype-vcf run, read from the tool's junction table
    def vcf_run(prefix):
        jrow = {}
        with open(prefix + ".junctions.tsv") as f:
            for p in csv.DictReader(f, delimiter="\t"):
                if p["known"] == "false": jrow[(p["chrom"], int(p["start"]) - 1, int(p["end"]))] = p
        with open(prefix + ".adjudicated.tsv") as f: adj_v = {p["isoform_id"]: p for p in csv.DictReader(f, delimiter="\t")}
        prov_v = dict(l.rstrip("\n").split("\t", 1) for l in open(prefix + ".provenance.log") if "\t" in l and not l.startswith("#"))
        explained = {k: p["variant"] for k, p in jrow.items() if p["variant"] in ("created", "deleted")}
        V_rule = {k for k in explained if jrow[k]["read_phase"] != "contradicted"}
        rescued_v = [i for i, p in adj_v.items() if p["confidence_class"] == "REFERENCE_BIAS"]
        V = {k for iso in rescued_v for k in jx_of.get(iso, []) if k in V_rule}
        top = {"rescued_isoforms": len(rescued_v),
               "rescued_by_reason": dict(collections.Counter(adj_v[i]["verdict_reason"].split(";mapping_sign")[0] for i in rescued_v)),
               "held_phase_unknown": sum(p["verdict_reason"] == "haplotype_phase_unknown" for p in adj_v.values()),
               "explained_junctions": dict(collections.Counter(explained.values())),
               "junction_phase": dict(collections.Counter(jrow[k]["variant_phase"] for k in explained)),
               "junctions_unlinked": sum(p["variant_phase"] == "unlinked" for p in jrow.values()),
               "read_phase": {v: dict(collections.Counter(jrow[k]["read_phase"] for k in explained if explained[k] == v)) for v in ("created", "deleted")},
               "indels_applied_at_motifs": int(prov_v.get("haplotype_vcf_indels_applied_at_motifs", 0)),
               "unphased_het_applied_at_motifs": int(prov_v.get("haplotype_vcf_unphased_het_applied_at_motifs", 0)),
               "unapplied_alleles_at_motifs": int(prov_v.get("haplotype_vcf_unapplied_alleles_at_motifs", 0)),
               "V_rule": len(V_rule), "V": len(V)}
        return {"jrow": jrow, "adj": adj_v, "explained": explained, "V_rule": V_rule, "V": V, "rescued": rescued_v, "top": top}

    def vcf_against(run, lab, T, T_snv, reference=None):
        """One VCF run against one B arm. `reference`: the assembly-call run, to say what this call set finds of it."""
        by_label = lambda S: dict(collections.Counter(lab[k]["label"] if k in lab else "absent" for k in S))
        pr = lambda num, den: len(num) / len(den) if den else float("nan")
        explained, V_rule, V = run["explained"], run["V_rule"], run["V"]
        T_unspliced = {k for k in T if lab[k]["label"] == "unspliced"}
        created = {k for k in V_rule if explained[k] == "created"}; deleted = V_rule - created
        res = {"V_rule": len(V_rule), "V_rule_and_T": len(V_rule & T), "precision": pr(V_rule & T, V_rule), "sensitivity": pr(V_rule & T, T),
               "precision_ci": boot(V_rule & T, V_rule, lambda k: k in T), "sensitivity_ci": boot(V_rule & T, T, lambda k: k in V_rule),
               "created": {"n": len(created), "in_T_snv_class": len(created & T_snv), "precision": pr(created & T_snv, created),
                           "sensitivity_snv_class": pr(created & T_snv, T_snv), "other_by_label": by_label(created - T_snv)},
               "deleted": {"n": len(deleted), "unspliced": len(deleted & T_unspliced), "precision": pr(deleted & T_unspliced, deleted),
                           "sensitivity_unspliced": pr(deleted & T_unspliced, T_unspliced), "other_by_label": by_label(deleted - T_unspliced)},
               "contradicted_by_reads_by_label": by_label(set(explained) - V_rule),
               "V": len(V), "V_and_T": len(V & T), "precision_in_rescued_isoforms": pr(V & T, V),
               "T_snv_missed": sorted(f"{k[0]}:{k[1]}-{k[2]}" for k in T_snv - created)}
        if reference is not None:   # of what the assembly's own calls explain and the assembly confirms: found here, or why not
            ref_ok = reference["V_rule"] & T
            def state(k):
                p = run["jrow"].get(k)
                if p is None: return "junction not in this run"
                if p["read_phase"] == "contradicted": return "contradicted by the reads"
                if p["variant_phase"] == "unlinked": return "phase not given (unlinked)"
                return {"none": "no explaining variant in the call set", "NA": "a base at the motif unknown", "disrupted": "disrupted"}.get(p["variant"], p["variant"])
            res["of_the_assembly_arm_confirmed"] = {
                "n": len(ref_ok), "found": len(ref_ok & V_rule),
                "created": {"n": sum(reference["explained"][k] == "created" for k in ref_ok), "found": sum(reference["explained"][k] == "created" and k in V_rule for k in ref_ok)},
                "deleted": {"n": sum(reference["explained"][k] == "deleted" for k in ref_ok), "found": sum(reference["explained"][k] == "deleted" and k in V_rule for k in ref_ok)},
                "missed_by_state": dict(collections.Counter(state(k) for k in ref_ok - V_rule)),
                "explained_here_only": len((V_rule & T) - ref_ok)}
        return res

    main_run = vcf_run(a.vcf_arm) if a.vcf_arm else None
    vcf = main_run["top"] if main_run else None
    if vcf is not None: out["vcf_route"] = vcf
    def deletions_of(path):
        """Deleted reference intervals [a, b) of a single-sample VCF, per chromosome: sequence alleles and <DEL> with END."""
        out = collections.defaultdict(list)
        for line in gzip.open(path, "rt"):
            if line[0] == "#": continue
            f = line.rstrip("\n").split("\t"); pos = int(f[1])
            if not re.search(r"[1-9]", f[9].split(":")[0]): continue
            if f[4].startswith("<DEL"):
                m = re.search(r"(?:^|;)END=(\d+)", f[7])
                if m: out[f[0]].append((pos, int(m.group(1)), "symbolic"))
            elif not f[4].startswith("<") and len(f[3]) > len(f[4]):
                out[f[0]].append((pos + len(f[4]) - 1, pos - 1 + len(f[3]), "sequence"))
        return out

    SIZES = (("under_50_bp", 0, 50), ("50_to_99_bp", 50, 100), ("100_to_299_bp", 100, 300), ("300_bp_or_more", 300, 10 ** 9))
    def deletion_lookup(run, dels, ref_deleted):
        """The junctions that equal a deletion in the assembly's calls, by length: found by this call set, and for the
        others what the call set's file holds there (no deletion that overlaps; one of the same extent that the tool did
        not use, e.g. next to an allele it cannot apply; one with other coordinates, by reciprocal overlap)."""
        res = {}
        for name, lo, hi in SIZES:
            ks = [k for k in ref_deleted if lo <= k[2] - k[1] < hi]; why = collections.Counter(); found = 0
            for c, s0, e0 in ks:
                if (c, s0, e0) in run["V_rule"] and run["explained"][(c, s0, e0)] == "deleted": found += 1; continue
                best = None
                for d0, d1, kind in dels.get(c, []):
                    ov = min(e0, d1) - max(s0, d0)
                    if ov <= 0: continue
                    rec = ov / max(e0 - s0, d1 - d0)
                    if best is None or rec > best[0]: best = (rec, d0 == s0 and d1 == e0)
                why["no_deletion_overlaps" if best is None else "same_extent_not_used" if best[1]
                    else "other_coordinates_overlap_ge_80pct" if best[0] >= 0.8 else "other_coordinates_overlap_lt_80pct"] += 1
            res[name] = {"n": len(ks), "found": found, "missed": dict(why)}
        return res

    def same_person(path, chrom="chr22"):
        """Share of the call set's SNV alleles on one chromosome that the assembly's own calls also have: near 1 for the
        same person, about one half for another."""
        mine = set()
        for line in gzip.open(path, "rt"):
            if line.startswith(chrom + "\t"):
                f = line.split("\t")
                if len(f[3]) == 1 and len(f[4]) == 1: mine.add((int(f[1]), f[4]))
        asm = set()
        for h in ("mat", "pat"):
            for line in pysam.TabixFile(f"{D}/{h}.hap.vcf.gz").fetch(chrom):
                f = line.split("\t")
                if len(f[3]) == 1 and len(f[4]) == 1: asm.add((int(f[1]), f[4]))
        return {"chrom": chrom, "snv_alleles": len(mine), "also_in_the_assembly_calls": len(mine & asm), "share": round(len(mine & asm) / len(mine), 3) if mine else None}

    call_sets = {}; call_vcf = {}; identity = {}
    for spec in a.call_set:
        cname, rest = spec.split("=", 1); prefix, _, vcf_path = rest.partition(",")
        call_sets[cname] = vcf_run(prefix)
        if vcf_path: call_vcf[cname] = deletions_of(vcf_path); identity[cname] = same_person(vcf_path)
    for n, r in call_sets.items():
        if n in identity: r["top"]["same_person_check"] = identity[n]
    if call_sets: out["call_sets"] = {n: {"summary": r["top"], "arms": {}} for n, r in call_sets.items()}   # the same dict: filled further below

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
        if main_run is not None: arm["vcf_route"] = vcf_against(main_run, lab, T, T_snv)
        for cname, r in call_sets.items():
            out["call_sets"][cname]["arms"][name] = vcf_against(r, lab, T, T_snv, reference=main_run)
            if cname in call_vcf and main_run is not None:
                ref_deleted = {k for k in main_run["V_rule"] & T if main_run["explained"][k] == "deleted"}
                out["call_sets"][cname]["arms"][name]["deletion_junctions_by_length"] = deletion_lookup(r, call_vcf[cname], ref_deleted)
        out["arms"][name] = arm

    # isoform level, A vs C: SQANTI3 rules filter verdict x PanIsoGuard class (novel isoforms only)
    rules = {}
    with open(f"{W}/armA/sq/rules/rules_RulesFilter_classification.txt") as f:
        for p in csv.DictReader(f, delimiter="\t"):
            if p["structural_category"] in ("novel_in_catalog", "novel_not_in_catalog"): rules[p["isoform"]] = p["filter_result"]
    xt = collections.Counter((rules[i], cls_C["wgs"][i]["confidence_class"]) for i in rules if i in cls_C["wgs"])
    out["isoform_A_vs_C"] = {f"SQ3rules={k[0]}|PIG={k[1]}": v for k, v in sorted(xt.items())}
    out["rescued_discarded_by_rules_filter"] = sum(1 for i in rescued if rules.get(i) == "Artifact")
    for run in ([main_run] if main_run else []) + list(call_sets.values()):   # what SQANTI3's rules filter does with the isoforms a VCF run rescues, by reason
        run["top"]["rescued_sqanti3_rules"] = {}
        for i in run["rescued"]:
            d = run["top"]["rescued_sqanti3_rules"].setdefault(run["adj"][i]["verdict_reason"].split(";mapping_sign")[0], {})
            d[rules.get(i, "not_filtered")] = d.get(rules.get(i, "not_filtered"), 0) + 1

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
        if vcf is not None: sample["vcf_route"] = vcf
        if call_sets:
            sample["call_sets"] = {n: {"summary": c["summary"], "arms": {arm: {k: v for k, v in x.items() if k != "T_snv_missed"} for arm, x in c["arms"].items()}}
                                   for n, c in out["call_sets"].items()}
        samples = {}
        if os.path.exists(a.emit_metrics):
            try: samples = json.load(open(a.emit_metrics))["metrics"].get("samples", {})
            except (KeyError, ValueError): samples = {}
        samples[a.label] = sample
        env = {"protocol": "refbias_3way", "scope": "hprc_r2_grch38_vs_personalized_pipelines", "status": "tracked",
               "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
               "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
               "data_provenance": "real public HPRC R2 Kinnex Iso-Seq of HG03516 (ESN) and HG02717 (GWD) through four arms on the same reads: "
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
