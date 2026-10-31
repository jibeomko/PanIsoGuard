#!/usr/bin/env python3
"""The rulebook of RuleEngine::evaluate() in 42 lines of Python (evaluate() below), replayed on the
worked locus and compared, isoform by isoform, with the verdicts the binary writes. The six runs add
one input at a time (reference GTF, short-read junctions, long-read BAM, haplotypes, provenance).
CTest runs this as `integration_rules_replay`.

    python3 rules_replay.py [path/to/panisoguard]      # default: ../../build/panisoguard
"""
import csv, json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
D = lambda name: os.path.join(HERE, "data", name)
PIG = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "build", "panisoguard")
NOVEL = ("novel_in_catalog", "novel_not_in_catalog")
MAP_GATES = {"bam_frac_low_mapq": 0.5, "bam_frac_supplementary": 0.5, "bam_frac_indel_near": 0.5}


def evaluate(cat, e, perc_a, min_callers=2, min_chain=3):
    """RuleEngine::evaluate() with the built-in defaults, on one attribution.jsonl record."""
    if cat == "full-splice_match":
        return "KNOWN"
    if cat == "incomplete-splice_match":                      # a known structure: only "did reads carry it whole?"
        if e["chain_evaluable"] and e["chain_n_full_reads"] >= min_chain and not (perc_a is not None and perc_a >= 60):
            return "KNOWN_FRAGMENT"
        return "UNCONFIRMED"
    if cat not in NOVEL:
        return "UNRESOLVED"
    rescues = [e[f"{a}_circular"] for a in ("pangenome", "variant") if e[f"{a}_evaluable"] and e[f"{a}_rescue"]]
    if rescues:                                                # reference-bias rescue first: independent evidence
        return "REFERENCE_BIAS" if not all(rescues) else "UNRESOLVED"   # on either axis promotes, else held
    n, k = e["n_novel_junctions"], e["n_novel_jx_sr_supported"]
    if not (e["sj_evaluable"] and e["chain_available"]) or n == 0:
        sup = "UNKNOWN"
    else:
        sup = "SUPPORTED" if k == n else "UNSUPPORTED" if k == 0 else "PARTIAL"
    if e["bam_evaluable"] and e["bam_n_spanning"] > 0 and any(e[f] > t for f, t in MAP_GATES.items()):
        mech = "mapping_or_repeat"
    elif e["noncanonical"]:
        mech = "noncanonical"
    elif e["rts_stage"]:
        mech = "rt_switch"
    elif perc_a is not None and perc_a >= 60:
        mech = "degradation"
    else:
        mech = "none"
    if sup == "UNKNOWN":
        weak = mech in ("none", "rt_switch")                     # the RT-switch sign does not demote a lift
        if mech == "mapping_or_repeat":
            return "UNCONFIRMED"                                 # held, not rejected, without short reads
        if e["consensus_evaluable"] and e["n_callers"] >= min_callers:
            return "NOVEL_PROBABLE" if weak else "UNCONFIRMED"
        if e["chain_evaluable"] and e["chain_n_full_reads"] >= min_chain:   # whole-chain long reads
            return "NOVEL_PROBABLE" if weak else "UNCONFIRMED"
        return "UNRESOLVED"
    if sup == "SUPPORTED":
        return "NOVEL_CONFIRMED" if mech == "none" else "NOVEL_PROBABLE"
    if sup == "PARTIAL":
        return "UNCONFIRMED"
    return "UNCONFIRMED" if mech == "none" else "ARTIFACT"


def main():
    base = ["--classification", D("caller1_classification.txt"), "--isoforms-gtf", D("caller1.gtf")]
    R, S, M = ["--ref-gtf", D("reference.gtf")], ["--sj-tab", D("short_reads.SJ.out.tab")], ["--bam", D("long_reads.bam")]
    H = ["--reference", D("toy.fa"), "--reference-haplotype", D("hap1.fa"), "--reference-haplotype", D("hap2.fa")]
    runs = [("step0", []), ("step1", R), ("step2", R + S), ("step3", R + S + M), ("step4", R + S + M + H),
            ("step5", R + S + M + H + ["--haplotype-provenance", "wgs"])]
    cls = {r["isoform"]: r for r in csv.DictReader(open(D("caller1_classification.txt")), delimiter="\t")}
    n = diff = 0
    with tempfile.TemporaryDirectory() as tmp:
        for name, args in runs:
            out = os.path.join(tmp, name)
            subprocess.run([PIG, "adjudicate", *base, *args, "--out-prefix", out], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            for line in open(out + ".attribution.jsonl"):
                r = json.loads(line)
                pa = cls[r["isoform"]].get("perc_A_downstream_TTS", "NA")
                mine = evaluate(r["structural_category"], r["evidence"], None if pa in ("NA", "") else float(pa))
                n += 1
                if mine != r["confidence_class"]:
                    diff += 1; print("DIFF", name, r["isoform"], mine, r["confidence_class"])
    print(f"{len(runs)} runs, {n} verdicts, {n - diff} identical to PanIsoGuard")
    sys.exit(1 if diff else 0)


if __name__ == "__main__":
    main()
