"""Data layer for the PanIsoGuard report: load + aggregate adjudication output.

Spec-independent. Reads the three files `panisoguard adjudicate` writes
(<prefix>.attribution.jsonl / .adjudicated.tsv / .provenance.log) into plain
structures the renderer aggregates. Pure stdlib (no third-party deps here).
"""
from __future__ import annotations

import json
import os
from collections import Counter, defaultdict


# Canonical ordering / tiering of the 7 confidence classes (engine order).
CONFIDENCE_ORDER = [
    "HIGH_CONF_KNOWN",
    "HIGH_CONF_NOVEL",
    "MEDIUM_CONF_NOVEL",
    "LOW_CONF_PARTIAL",
    "PAN_REF_RESCUED_FALSE_NOVEL",
    "AMBIGUOUS",
    "ARTIFACT",
]
NOVELTY_ORDER = ["SUPPORTED", "PARTIAL", "UNSUPPORTED", "UNKNOWN"]
MECHANISM_ORDER = [
    "none", "mapping_or_repeat", "noncanonical", "rt_switch",
    "degradation", "variant_created", "population_known",
]


def get_path(rec, dotted, default=None):
    """rec, 'evidence.n_callers' -> nested lookup (None-safe)."""
    cur = rec
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


def load_records(jsonl_path):
    """Parse <prefix>.attribution.jsonl into a list of dicts (skips blank lines)."""
    recs = []
    with open(jsonl_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            recs.append(json.loads(line))
    return recs


def load_provenance(log_path):
    """Parse <prefix>.provenance.log (tab-separated key/value, '#' comments) -> dict.

    Returns {meta: {...str...}, axes: {axis_name: state}, class_counts: {NAME: int},
    total: int|None}.
    """
    meta, axes, class_counts, total = {}, {}, {}, None
    if not os.path.exists(log_path):
        return dict(meta=meta, axes=axes, class_counts=class_counts, total=total)
    with open(log_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            key, val = parts[0], parts[1]
            if key == "total":
                try:
                    total = int(val)
                except ValueError:
                    pass
            elif key.startswith("class."):
                try:
                    class_counts[key[len("class."):]] = int(val)
                except ValueError:
                    pass
            elif key.startswith("axis."):
                axes[key[len("axis."):]] = val
            else:
                meta[key] = val
    return dict(meta=meta, axes=axes, class_counts=class_counts, total=total)


def count_field(records, dotted):
    """Counter of a (possibly nested) field's values across records."""
    c = Counter()
    for r in records:
        c[get_path(r, dotted)] += 1
    return c


def crosstab(records, row_field, col_field):
    """2D Counter[(row,col)] for a heatmap (e.g. novelty_support x primary_mechanism)."""
    ct = defaultdict(int)
    for r in records:
        ct[(get_path(r, row_field), get_path(r, col_field))] += 1
    return ct


def axis_evaluability(records):
    """How many records had each evidence axis evaluable (a waterfall of available evidence)."""
    axes = [
        ("short-read SJ", "evidence.sj_evaluable"),
        ("BAM mapping", "evidence.bam_evaluable"),
        ("variant", "evidence.variant_evaluable"),
        ("pangenome", "evidence.pangenome_evaluable"),
        ("consensus", "evidence.consensus_evaluable"),
    ]
    out = []
    n = len(records)
    for label, field in axes:
        k = sum(1 for r in records if get_path(r, field) is True)
        out.append((label, k, n))
    return out


def rescue_firewall(records):
    """Reference-bias rescue accounting: promoted vs held-by-firewall, per axis."""
    out = {}
    for axis in ("variant", "pangenome"):
        rescue = sum(1 for r in records if get_path(r, f"evidence.{axis}_rescue") is True)
        circ = sum(1 for r in records if get_path(r, f"evidence.{axis}_circular") is True)
        out[axis] = dict(rescue_fired=rescue, held_circular=circ)
    promoted = sum(1 for r in records
                   if get_path(r, "confidence_class") == "PAN_REF_RESCUED_FALSE_NOVEL")
    held = sum(1 for r in records if get_path(r, "circularity_flag") is True)
    out["promoted_total"] = promoted
    out["held_by_firewall_total"] = held
    return out


def numeric_values(records, dotted, predicate=None):
    """Collect non-null numeric values of a field (e.g. bio_flags.dist_to_CAGE_peak)."""
    vals = []
    for r in records:
        if predicate is not None and not predicate(r):
            continue
        v = get_path(r, dotted)
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            vals.append(v)
    return vals
