#!/usr/bin/env python3
"""Generate a synthetic-but-realistic adjudication output that exercises EVERY report
page, so `panisoguard-report` has a self-contained showcase fixture.

This does NOT run the engine; it writes a hand-curated <prefix>.attribution.jsonl +
.provenance.log whose records are internally consistent with the rule engine (the same
fields `adjudicate` emits), spanning all 7 confidence classes, every evidence axis
(short-read, BAM mapping, variant, pangenome, consensus) and the SQANTI3 bio priors.
Deterministic (no randomness). It is a DOCUMENTATION fixture, clearly labelled.

Usage: make_showcase.py [out_prefix]   (default: showcase)
"""
import json
import sys
from collections import Counter

prefix = sys.argv[1] if len(sys.argv) > 1 else "showcase"

# axis toggles for this fixture (all on, so every page renders)
AXES = dict(short_read="on", catalog="on", bam="on", variant="on", pangenome="on", consensus="on")


def rec(iso, cat, cls, support, mech, *, n_nov=0, n_sr=0,
        bam=None, variant=None, pangenome=None, n_callers=None,
        circular=False, bio=None, trace=None):
    """Build one attribution record with evidence consistent with (cls, support, mech)."""
    ev = dict(
        chain_available=True,
        sj_evaluable=(support != "UNKNOWN"),
        n_novel_junctions=n_nov,
        n_novel_jx_sr_supported=n_sr,
        rts_stage=(mech == "rt_switch"),
        noncanonical=(mech == "noncanonical"),
        bam_evaluable=bam is not None,
        bam_n_spanning=(bam or {}).get("n", 0),
        bam_frac_low_mapq=(bam or {}).get("low_mapq", 0.0),
        bam_frac_supplementary=(bam or {}).get("suppl", 0.0),
        bam_frac_softclip=(bam or {}).get("softclip", 0.0),
        bam_frac_indel_near=(bam or {}).get("indel", 0.0),
        variant_evaluable=variant is not None,
        variant_rescue=(variant or {}).get("rescue", False),
        variant_circular=(variant or {}).get("circular", False),
        pangenome_evaluable=pangenome is not None,
        pangenome_rescue=(pangenome or {}).get("rescue", False),
        pangenome_circular=(pangenome or {}).get("circular", False),
        n_novel_jx_pangenome=(pangenome or {}).get("n_pan", 0),
        consensus_evaluable=n_callers is not None,
        n_callers=(n_callers if n_callers is not None else 0),
    )
    b = bio or {}
    return dict(
        isoform=iso, structural_category=cat, confidence_class=cls,
        novelty_support=support, primary_mechanism=mech, circularity_flag=circular,
        evidence=ev,
        graph_trace=dict(novel_edges=n_nov, path_distance_to_reference=n_nov,
                         edges_sr_supported=n_sr, novel_site_type=("novel_splice_site"
                         if cat == "novel_not_in_catalog" else ("known_site_recombination"
                         if cat == "novel_in_catalog" else "none")),
                         pangenome_edge_support=(pangenome or {}).get("rescue", False),
                         variant_canonicalized_edge=(variant or {}).get("rescue", False)),
        bio_flags=dict(dist_to_CAGE_peak=b.get("cage"), dist_to_polyA_site=b.get("polya"),
                       predicted_NMD=b.get("nmd", "NA"), polyA_motif_found=b.get("motif", "NA")),
        rule_trace=trace or [f"project({support},{mech}) -> {cls}"],
    )


records = []
i = 0


def add(n, **kw):
    global i
    for _ in range(n):
        i += 1
        kw2 = dict(kw)
        kw2["iso"] = f"{kw['iso']}.{i}"
        records.append(rec(**kw2))


# --- known pass-through (FSM) ---
add(9, iso="FSM", cat="full-splice_match", cls="HIGH_CONF_KNOWN", support="SUPPORTED",
    mech="none", bio=dict(cage=5, polya=8, nmd="False", motif="AATAAA"),
    trace=["category=full-splice_match -> HIGH_CONF_KNOWN"])
# --- confident novel: SUPPORTED x none ---
add(7, iso="NNC_hi", cat="novel_not_in_catalog", cls="HIGH_CONF_NOVEL", support="SUPPORTED",
    mech="none", n_nov=3, n_sr=3, bam=dict(n=22, low_mapq=0.04, suppl=0.02, softclip=0.10),
    bio=dict(cage=12, polya=20, nmd="False", motif="AATAAA"))
# --- MEDIUM: SUPPORTED x noncanonical, and PARTIAL x none ---
add(5, iso="NNC_med1", cat="novel_not_in_catalog", cls="MEDIUM_CONF_NOVEL", support="SUPPORTED",
    mech="noncanonical", n_nov=2, n_sr=2, bio=dict(cage=40, polya=33))
add(4, iso="NIC_med2", cat="novel_in_catalog", cls="MEDIUM_CONF_NOVEL", support="PARTIAL",
    mech="none", n_nov=3, n_sr=2)
# --- LOW: PARTIAL x mechanism, UNSUPPORTED x none ---
add(4, iso="NNC_low1", cat="novel_not_in_catalog", cls="LOW_CONF_PARTIAL", support="PARTIAL",
    mech="rt_switch", n_nov=3, n_sr=1)
add(3, iso="NNC_low2", cat="novel_not_in_catalog", cls="LOW_CONF_PARTIAL", support="UNSUPPORTED",
    mech="none", n_nov=2, n_sr=0)
# --- ARTIFACT: UNSUPPORTED x mechanism (incl. a BAM mapping artifact) ---
add(3, iso="NNC_art", cat="novel_not_in_catalog", cls="ARTIFACT", support="UNSUPPORTED",
    mech="noncanonical", n_nov=2, n_sr=0, bio=dict(nmd="True"))
add(3, iso="NNC_map", cat="novel_not_in_catalog", cls="ARTIFACT", support="UNSUPPORTED",
    mech="mapping_or_repeat", n_nov=2, n_sr=0,
    bam=dict(n=14, low_mapq=0.86, suppl=0.71, softclip=0.40, indel=0.22),
    trace=["bam mapping artifact (low_mapq_frac=0.86) -> mechanism=mapping_or_repeat",
           "project(UNSUPPORTED,mapping_or_repeat) -> ARTIFACT"])
# --- reference-bias rescue: variant (promoted) + pangenome (promoted) ---
add(3, iso="RESC_var", cat="novel_not_in_catalog", cls="PAN_REF_RESCUED_FALSE_NOVEL",
    support="UNKNOWN", mech="variant_created", n_nov=2,
    variant=dict(rescue=True, circular=False),
    trace=["variant creates canonical splice motif on haplotype, non-canonical on "
           "reference -> PAN_REF_RESCUED_FALSE_NOVEL (reference bias)"])
add(2, iso="RESC_pan", cat="novel_not_in_catalog", cls="PAN_REF_RESCUED_FALSE_NOVEL",
    support="UNKNOWN", mech="population_known", n_nov=2,
    pangenome=dict(rescue=True, circular=False, n_pan=2),
    trace=["all 2/2 novel junctions realizable on a pangenome haplotype path -> "
           "PAN_REF_RESCUED_FALSE_NOVEL (reference bias, population_known)"])
# --- HELD by circularity firewall (the thesis): rescue fired but provenance circular-risk ---
add(2, iso="HELD_var", cat="novel_not_in_catalog", cls="AMBIGUOUS", support="UNKNOWN",
    mech="variant_created", n_nov=2, circular=True,
    variant=dict(rescue=True, circular=True),
    trace=["variant creates canonical motif on haplotype but provenance is circular-risk "
           "(RNA-derived/unknown) -> held AMBIGUOUS (not promoted)"])
add(1, iso="HELD_pan", cat="novel_not_in_catalog", cls="AMBIGUOUS", support="UNKNOWN",
    mech="population_known", n_nov=2, circular=True,
    pangenome=dict(rescue=True, circular=True, n_pan=2),
    trace=["all 2/2 novel junctions realizable on a pangenome path, but the junction-set "
           "provenance is sample-derived/unknown (circular-risk) -> held AMBIGUOUS"])
# --- consensus axis: UNKNOWN-support novels promoted by multi-caller agreement ---
add(5, iso="CONS_med", cat="novel_not_in_catalog", cls="MEDIUM_CONF_NOVEL", support="UNKNOWN",
    mech="none", n_callers=4,
    trace=["caller-consensus n_callers=4 >= 2 corroborates the novel chain (no short-read "
           "axis) -> consensus-supported", "project(UNKNOWN,none) -> MEDIUM_CONF_NOVEL"])
add(3, iso="CONS_low", cat="novel_not_in_catalog", cls="LOW_CONF_PARTIAL", support="UNKNOWN",
    mech="noncanonical", n_callers=2,
    trace=["caller-consensus n_callers=2 >= 2 corroborates the novel chain -> "
           "consensus-supported", "project(UNKNOWN,noncanonical) -> LOW_CONF_PARTIAL"])
add(4, iso="CONS_ambi", cat="novel_not_in_catalog", cls="AMBIGUOUS", support="UNKNOWN",
    mech="none", n_callers=1,
    trace=["novelty-support UNKNOWN reason=short_read/catalog_axis_absent",
           "project(UNKNOWN,none) -> AMBIGUOUS"])
# --- ISM pass-through ---
add(4, iso="ISM", cat="incomplete-splice_match", cls="LOW_CONF_PARTIAL", support="UNKNOWN",
    mech="none", trace=["category=incomplete-splice_match -> LOW_CONF_PARTIAL"])

# --- write outputs ---
with open(f"{prefix}.attribution.jsonl", "w") as fh:
    for r in records:
        fh.write(json.dumps(r) + "\n")

counts = Counter(r["confidence_class"] for r in records)
with open(f"{prefix}.provenance.log", "w") as fh:
    fh.write("# PanIsoGuard adjudicate provenance (SYNTHETIC SHOWCASE FIXTURE)\n")
    fh.write("tool_version\t0.0.3\nruleset_version\tbuiltin-0.0.1\nsqanti3_version_target\t6.0\n")
    fh.write("config\t<built-in defaults>\n")
    for k in ("classification", "isoforms", "ref_gtf", "sj_tab", "bam", "caller_support"):
        fh.write(f"{k}\tshowcase (synthetic demo)\n")
    for a, st in AXES.items():
        fh.write(f"axis.{a}\t{st}\n")
    fh.write("# confidence-class counts\n")
    for c in counts:
        fh.write(f"class.{c}\t{counts[c]}\n")
    fh.write(f"total\t{len(records)}\n")

sys.stderr.write(f"wrote {prefix}.attribution.jsonl + {prefix}.provenance.log "
                 f"({len(records)} records across {len(counts)} classes)\n")
