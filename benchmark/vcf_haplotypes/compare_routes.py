#!/usr/bin/env python3
"""Compare the variant axis read from a VCF with the same axis read from consensus FASTAs, run by run.
For each LABEL:FASTA_PREFIX:VCF_PREFIX pair: isoforms whose confidence class and verdict reason agree,
junctions whose variant call (none / created / disrupted / NA) agrees and whose first-haplotype motif
agrees (the FASTA route carries one haplotype, the VCF route the sample's ploidy; haplotype 1 of the VCF is
the FASTA's: `bcftools consensus -H 1`), the VCF route's unphased-het and non-SNV counters from its
provenance log, and the junctions marked phased / unphased.

Usage: compare_routes.py LABEL:FASTA_PREFIX:VCF_PREFIX ... --emit-metrics out.json
"""
import csv, datetime, json, os, sys


def table(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def provenance(path):
    out = {}
    for l in open(path):
        if "\t" in l and not l.startswith("#"):
            k, _, v = l.rstrip("\n").partition("\t"); out.setdefault(k, v)
    return out


def main():
    args = sys.argv[1:]
    out = args[args.index("--emit-metrics") + 1]
    specs = [a for a in args if a.count(":") >= 2 and not a.startswith("--")]
    per = {}; prov = {}
    for spec in specs:
        label, fa, vc = spec.split(":", 2)
        p = per.setdefault(label, {"runs": 0, "isoforms": 0, "same_class": 0, "same_reason": 0, "reference_bias_fasta": 0, "reference_bias_vcf": 0,
                                   "junctions": 0, "same_variant_call": 0, "same_haplotype1_motif": 0, "variant_evaluable_vcf": 0,
                                   "isoforms_missing_in_vcf_run": 0, "isoforms_only_in_vcf_run": 0,
                                   "junctions_missing_in_vcf_run": 0, "junctions_only_in_vcf_run": 0,
                                   "phase": {"phased": 0, "unphased": 0, "NA": 0}, "unphased_het_applied": 0, "non_snv_skipped": 0, "missing_alleles": 0,
                                   "differences": []})
        fa_iso = {r["isoform_id"]: r for r in table(fa + ".adjudicated.tsv")}; vc_iso = {r["isoform_id"]: r for r in table(vc + ".adjudicated.tsv")}
        p["runs"] += 1; p["isoforms"] += len(fa_iso)   # the FASTA route's rows are the denominator
        p["reference_bias_fasta"] += sum(1 for r in fa_iso.values() if r["confidence_class"] == "REFERENCE_BIAS")   # over every row of each route
        p["reference_bias_vcf"] += sum(1 for r in vc_iso.values() if r["confidence_class"] == "REFERENCE_BIAS")
        extra = [i for i in vc_iso if i not in fa_iso]
        p["isoforms_only_in_vcf_run"] += len(extra)
        for i in extra[:5]: p["differences"].append(f"{label} {i}: only in the VCF run")
        for i, r in fa_iso.items():
            q = vc_iso.get(i)
            if q is None:
                p["isoforms_missing_in_vcf_run"] += 1   # counted against same_class / same_reason: a missing row is a disagreement
                if len(p["differences"]) < 20: p["differences"].append(f"{label} {i}: missing in the VCF run")
                continue
            p["same_class"] += r["confidence_class"] == q["confidence_class"]
            p["same_reason"] += r["verdict_reason"] == q["verdict_reason"]
            if r["confidence_class"] != q["confidence_class"] and len(p["differences"]) < 20:
                p["differences"].append(f"{label} {i}: {r['confidence_class']} ({r['verdict_reason']}) vs {q['confidence_class']} ({q['verdict_reason']})")
        key = lambda r: (r["chrom"], r["start"], r["end"], r["strand"])
        fa_jx = {key(r): r for r in table(fa + ".junctions.tsv")}; vc_jx = {key(r): r for r in table(vc + ".junctions.tsv")}
        p["junctions_only_in_vcf_run"] += sum(1 for k in vc_jx if k not in fa_jx)
        for k, r in fa_jx.items():
            p["junctions"] += 1   # every FASTA-route junction counts, present in the VCF run or not
            q = vc_jx.get(k)
            if q is None:
                p["junctions_missing_in_vcf_run"] += 1
                if len(p["differences"]) < 20: p["differences"].append(f"{label} junction {k}: missing in the VCF run")
                continue
            p["same_variant_call"] += r["variant"] == q["variant"]
            p["same_haplotype1_motif"] += r["hap_motifs"].split(",")[0] == q["hap_motifs"].split(",")[0]
            p["variant_evaluable_vcf"] += q["variant"] != "NA"
            p["phase"][q.get("variant_phase", "NA") if q.get("variant_phase", "NA") in p["phase"] else "NA"] += 1
        pv = provenance(vc + ".provenance.log"); prov = pv
        p["unphased_het_applied"] += int(pv.get("haplotype_vcf_unphased_het_applied_at_motifs", 0) or 0)
        p["non_snv_skipped"] += int(pv.get("haplotype_vcf_non_snv_skipped_at_motifs", 0) or 0)
        p["missing_alleles"] += int(pv.get("haplotype_vcf_missing_alleles_at_motifs", 0) or 0)
    total = {k: sum(p[k] for p in per.values()) for k in ("runs", "isoforms", "same_class", "same_reason", "junctions", "same_variant_call", "same_haplotype1_motif",
                                                           "isoforms_missing_in_vcf_run", "isoforms_only_in_vcf_run", "junctions_missing_in_vcf_run", "junctions_only_in_vcf_run")}
    for label, p in per.items():
        print(f"{label:<8} runs {p['runs']:>3}  isoforms {p['isoforms']:>5}  same class {p['same_class']:>5}  same reason {p['same_reason']:>5}  "
              f"REFERENCE_BIAS fasta/vcf {p['reference_bias_fasta']}/{p['reference_bias_vcf']}  junctions {p['junctions']:>5}  same variant call {p['same_variant_call']:>5}  "
              f"same hap-1 motif {p['same_haplotype1_motif']:>5}  missing/only-in-VCF isoforms {p['isoforms_missing_in_vcf_run']}/{p['isoforms_only_in_vcf_run']} "
              f"junctions {p['junctions_missing_in_vcf_run']}/{p['junctions_only_in_vcf_run']}  phase {p['phase']}  unphased het {p['unphased_het_applied']}  "
              f"non-SNV skipped {p['non_snv_skipped']}  missing alleles {p['missing_alleles']}")
        for d in p["differences"]: print("   ", d)
    print("total", total)
    if (total["isoforms"] != total["same_class"] or total["isoforms"] != total["same_reason"] or total["junctions"] != total["same_variant_call"]
            or total["junctions"] != total["same_haplotype1_motif"] or total["isoforms_only_in_vcf_run"] or total["junctions_only_in_vcf_run"]):
        print("ROUTES DIFFER", file=sys.stderr)
    env = {"protocol": "vcf_haplotypes", "scope": "hg002_chr22_giab_v4.2.1_and_hprc_r2_hg03516_hg02717_all_chromosomes", "status": "tracked",
           "tool_version": prov.get("tool_version"), "ruleset_version": prov.get("ruleset_version"),
           "generated_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), "source": None,
           "data_provenance": "The variant-axis runs already on disk, re-run with the haplotypes read from a VCF instead of consensus FASTAs: HG002 chr22 "
                              "(GIAB v4.2.1 biallelic SNVs; FASTA = bcftools consensus -H 1 of the same file; benchmark/hg002) and every per-chromosome scan "
                              "of HG03516 and HG02717 (paftools SNV calls per parental HPRC R2 assembly, consensus FASTA per haplotype; "
                              "benchmark/hg03516_refbias). Public GIAB / HPRC / GENCODE v49.",
           "command": "bash benchmark/vcf_haplotypes/run.sh; benchmark/vcf_haplotypes/compare_routes.py LABEL:FASTA_PREFIX:VCF_PREFIX ... --emit-metrics ...",
           "metrics": {"per_sample": per, "total": total},
           "notes": "Equivalence of the two haplotype inputs: a run from an indexed VCF must give the verdicts and the per-junction motifs a run from the "
                    "consensus FASTAs built from that VCF gives (haplotype 1 of the VCF = `bcftools consensus -H 1`). Only same-length substitutions are "
                    "applied from the VCF; indels and symbolic alleles at a motif are skipped and counted (non_snv_skipped), unphased heterozygous "
                    "alleles are applied by order and counted (unphased_het_applied) and the junction is marked unphased."}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(env, open(out, "w"), indent=2, sort_keys=True); open(out, "a").write("\n"); print("wrote", out)


if __name__ == "__main__":
    main()
