#!/usr/bin/env python3
"""PanIsoGuard adjudication report -- a SQANTI3-style multi-page PDF summary.

Companion tool (keeps the C++ core dependency-free, like SQANTI3's separate report
step). Consumes the files `panisoguard adjudicate` already writes -- it recomputes
nothing about the isoforms, it only visualizes the verdicts and the evidence behind
them:

    <prefix>.attribution.jsonl   per-isoform verdict + evidence + rule_trace (richest)
    <prefix>.provenance.log      run metadata, axis on/off, authoritative class tally

Output: <prefix>.report.pdf (landscape, one figure per page).

Pages: (1) executive summary + trust-at-a-glance, (2) the two-axis decision space,
(3) artifact-mechanism diagnostics, (4) reference-bias rescue & circularity firewall,
(5) multi-caller consensus [only if --caller-support was used], (6) SQANTI3 bio priors
pass-through [only if supplied], (7) provenance & defensibility appendix.

Requires matplotlib (optional dependency; install separately). Usage:
    python -m panisoguard.report --prefix run [--out run.report.pdf]
"""
from __future__ import annotations

import argparse
import os
import sys

from . import _report_data as D

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Rectangle
except ImportError:  # pragma: no cover
    sys.stderr.write(
        "panisoguard.report needs matplotlib (optional dependency).\n"
        "  pip install matplotlib   (or: conda install matplotlib)\n")
    raise

# --- locked semantic palette (learned once, reused on every page) -------------
CLASS_COLORS = {
    "HIGH_CONF_KNOWN": "#1A7F37",            # trusted catalog pass-through
    "HIGH_CONF_NOVEL": "#3FAE5A",            # confident novel
    "MEDIUM_CONF_NOVEL": "#8CC63F",          # confident novel
    "PAN_REF_RESCUED_FALSE_NOVEL": "#2C7FB8",  # reclassification (reference bias), NOT a grade
    "LOW_CONF_PARTIAL": "#F0A202",           # tentative
    "AMBIGUOUS": "#9E9E9E",                  # honest abstention (absence of evidence)
    "ARTIFACT": "#C1272D",                   # reject
}
CLASS_ORDER = [
    "HIGH_CONF_KNOWN", "HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL",
    "PAN_REF_RESCUED_FALSE_NOVEL", "LOW_CONF_PARTIAL", "AMBIGUOUS", "ARTIFACT",
]
MECH_COLORS = {
    "none": "#BDBDBD", "mapping_or_repeat": "#5D4037", "noncanonical": "#D08C3F",
    "rt_switch": "#E0B050", "degradation": "#EAD18A",
    "variant_created": "#2C7FB8", "population_known": "#4FA3D1",
}
AXIS_ON = "#2A9D8F"        # teal: evidence axis in play
AXIS_OFF = "#D9D9D9"       # grey: not evaluable
NOVEL_CATS = {"novel_in_catalog", "novel_not_in_catalog"}


# --- small drawing helpers ----------------------------------------------------
def _footer(fig, prov, prefix):
    m = prov["meta"]
    txt = (f"PanIsoGuard {m.get('tool_version', '?')}  ·  ruleset {m.get('ruleset_version', '?')}"
           f"  ·  SQANTI3 target {m.get('sqanti3_version_target', '?')}  ·  {os.path.basename(prefix)}")
    fig.text(0.5, 0.012, txt, ha="center", va="bottom", fontsize=6.5, color="#888888")


def _placeholder(ax, msg):
    """A skipped/conditional panel degrades to a visible grey-hatch callout, never empty axes."""
    ax.add_patch(Rectangle((0.04, 0.18), 0.92, 0.64, transform=ax.transAxes,
                           facecolor="#F2F2F2", edgecolor="#BBBBBB", hatch="//", lw=0.8))
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center",
            fontsize=8.5, color="#777777", wrap=True)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def _hbar(ax, labels, values, colors, title, hatches=None, annotate_pct_of=None):
    y = range(len(labels))
    bars = ax.barh(list(y), values, color=colors, edgecolor="white", height=0.74)
    if hatches:
        for b, h in zip(bars, hatches):
            if h:
                b.set_hatch(h)
    ax.set_yticks(list(y)); ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    total = annotate_pct_of if annotate_pct_of else 0
    vmax = max(values) if values and max(values) > 0 else 1
    for yi, v in zip(y, values):
        lab = f"{v:,}" + (f"  ({100 * v / total:.1f}%)" if total else "")
        ax.text(v + vmax * 0.01, yi, lab, va="center", fontsize=7.5, color="#333333")
    ax.set_title(title, fontsize=10, loc="left", weight="bold")
    ax.set_xlim(0, vmax * 1.18)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(length=0)


def _callout(ax, title, lines, accent="#333333"):
    ax.axis("off")
    ax.text(0.0, 1.0, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=10, weight="bold", color="#222222")
    y = 0.80
    for txt, color in lines:
        ax.text(0.02, y, txt, transform=ax.transAxes, ha="left", va="top",
                fontsize=9, color=color)
        y -= 0.135


def _stat_tiles(ax, tiles):
    """tiles: list of (value, label, color). Big-number tiles in a row."""
    ax.axis("off")
    n = len(tiles)
    for i, (val, label, color) in enumerate(tiles):
        x0 = i / n
        ax.add_patch(Rectangle((x0 + 0.01, 0.15), 1.0 / n - 0.02, 0.7,
                               transform=ax.transAxes, facecolor=color, alpha=0.16,
                               edgecolor=color, lw=1.2))
        ax.text(x0 + 0.5 / n, 0.62, str(val), transform=ax.transAxes, ha="center",
                va="center", fontsize=20, weight="bold", color=color)
        ax.text(x0 + 0.5 / n, 0.30, label, transform=ax.transAxes, ha="center",
                va="center", fontsize=7.8, color="#444444", wrap=True)


def _heatmap(ax, matrix, row_labels, col_labels, title, border_colors=None,
             cmap="Blues"):
    import numpy as np
    arr = np.array(matrix, dtype=float)
    ax.imshow(arr, cmap=cmap, aspect="auto")
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=7, rotation=35, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    vmax = arr.max() if arr.size and arr.max() > 0 else 1
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            v = int(arr[i, j])
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=7.5,
                        color="white" if arr[i, j] > vmax * 0.55 else "#222222")
            if border_colors is not None and border_colors[i][j] is not None:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=border_colors[i][j], lw=2.4))
    ax.set_title(title, fontsize=10, loc="left", weight="bold")
    ax.tick_params(length=0)


def _table(ax, col_labels, rows, title, col_widths=None, max_rows=14):
    ax.axis("off")
    ax.set_title(title, fontsize=10, loc="left", weight="bold")
    shown = rows[:max_rows]
    tbl = ax.table(cellText=shown, colLabels=col_labels, loc="center",
                   cellLoc="left", colWidths=col_widths)
    tbl.auto_set_font_size(False); tbl.set_fontsize(7)
    tbl.scale(1, 1.25)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        if r == 0:
            cell.set_facecolor("#EFEFEF"); cell.set_text_props(weight="bold")
    if len(rows) > max_rows:
        ax.text(0.5, -0.02, f"... {len(rows) - max_rows} more rows (see attribution.jsonl)",
                transform=ax.transAxes, ha="center", fontsize=6.5, color="#999999")


# --- projection rulebook (Page 2 heatmap borders) -----------------------------
def projected_class(support, mech):
    has_mech = mech not in (None, "none")
    if support == "SUPPORTED":
        return "HIGH_CONF_NOVEL" if not has_mech else "MEDIUM_CONF_NOVEL"
    if support == "PARTIAL":
        return "MEDIUM_CONF_NOVEL" if not has_mech else "LOW_CONF_PARTIAL"
    if support == "UNSUPPORTED":
        return "LOW_CONF_PARTIAL" if not has_mech else "ARTIFACT"
    # UNKNOWN: mapping artifact decisive; otherwise abstain (consensus exception annotated)
    if support == "UNKNOWN":
        return "ARTIFACT" if mech == "mapping_or_repeat" else "AMBIGUOUS"
    return "AMBIGUOUS"


# --- run context --------------------------------------------------------------
class Ctx:
    def __init__(self, records, prov):
        self.records = records
        self.prov = prov
        self.novel = [r for r in records
                      if D.get_path(r, "structural_category") in NOVEL_CATS]
        self.has_consensus = any(D.get_path(r, "evidence.consensus_evaluable") is True
                                 for r in records)
        self.has_bam = any(D.get_path(r, "evidence.bam_evaluable") is True for r in records)
        self.has_variant = any(D.get_path(r, "evidence.variant_evaluable") is True for r in records)
        self.has_pangenome = any(D.get_path(r, "evidence.pangenome_evaluable") is True for r in records)
        self.has_sr = any(D.get_path(r, "evidence.sj_evaluable") is True for r in records)
        cc = D.count_field(records, "confidence_class")
        self.jsonl_counts = {k: cc.get(k, 0) for k in CLASS_ORDER}
        # integrity: provenance authoritative tally vs JSONL recompute
        pc = prov["class_counts"]
        self.mismatch = bool(pc) and any(
            pc.get(k, 0) != self.jsonl_counts.get(k, 0) for k in set(pc) | set(self.jsonl_counts))


# --- PAGE 1: executive summary ------------------------------------------------
def page1(pdf, ctx, prefix):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("PanIsoGuard adjudication — executive summary",
                 fontsize=13, weight="bold", x=0.07, ha="left")
    fig.text(0.07, 0.935, "Orthogonal confidence layer — SQANTI3 features consumed as "
             "priors, never recomputed.", fontsize=8, style="italic", color="#666666")
    gs = fig.add_gridspec(2, 2, hspace=0.42, wspace=0.22,
                          left=0.07, right=0.96, top=0.88, bottom=0.07)

    # axis availability strip (custom: one lozenge per axis, fraction fill, k/n at right)
    ax = fig.add_subplot(gs[0, 0])
    axes = D.axis_evaluability(ctx.records)
    y = range(len(axes))
    for i, (lab, k, n) in enumerate(axes):
        frac = k / n if n else 0
        on = k > 0
        ax.barh(i, frac, color=AXIS_ON if on else AXIS_OFF, height=0.66,
                edgecolor="#BBBBBB", hatch=None if on else "//")
        ax.text(1.04, i, f"{k:,}/{n:,}", va="center", ha="left", fontsize=7.5, color="#333")
    ax.set_yticks(list(y)); ax.set_yticklabels([lab for lab, _, _ in axes], fontsize=8)
    ax.invert_yaxis(); ax.set_xlim(0, 1.25); ax.set_xticks([0, 0.5, 1.0])
    ax.set_title("Evidence axes in play (fraction of records evaluable)",
                 fontsize=10, loc="left", weight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(length=0)

    # headline verdict
    ax = fig.add_subplot(gs[0, 1])
    nov = ctx.novel
    nnov = len(nov)
    def cc_in(rs, names):
        return sum(1 for r in rs if D.get_path(r, "confidence_class") in names)
    trust = cc_in(nov, {"HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL"})
    look = cc_in(nov, {"LOW_CONF_PARTIAL", "AMBIGUOUS"})
    art = cc_in(nov, {"ARTIFACT"})
    resc = cc_in(nov, {"PAN_REF_RESCUED_FALSE_NOVEL"})
    _callout(ax, f"Headline — {nnov:,} novel isoforms", [
        (f"Trustworthy (HIGH/MEDIUM_CONF_NOVEL):  {trust:,}", CLASS_COLORS["HIGH_CONF_NOVEL"]),
        (f"Need a look (LOW_CONF_PARTIAL/AMBIGUOUS):  {look:,}", CLASS_COLORS["LOW_CONF_PARTIAL"]),
        (f"Look like artifacts (ARTIFACT):  {art:,}", CLASS_COLORS["ARTIFACT"]),
        (f"Reclassified as reference-bias (RESCUED):  {resc:,}", CLASS_COLORS["PAN_REF_RESCUED_FALSE_NOVEL"]),
    ])

    # confidence-class distribution (whole call set), authoritative tally
    ax = fig.add_subplot(gs[1, 0])
    counts = ctx.prov["class_counts"] or ctx.jsonl_counts
    present = [c for c in CLASS_ORDER if counts.get(c, 0) > 0]
    total = ctx.prov["total"] or sum(counts.values()) or 1
    _hbar(ax, present, [counts.get(c, 0) for c in present],
          [CLASS_COLORS[c] for c in present], "Confidence-class distribution (all isoforms)",
          annotate_pct_of=total)
    if ctx.mismatch:
        ax.text(0.5, 1.14, "⚠ provenance tally ≠ JSONL recompute", transform=ax.transAxes,
                ha="center", color="#C1272D", fontsize=8, weight="bold")

    # firewall scorecard
    ax = fig.add_subplot(gs[1, 1])
    fw = D.rescue_firewall(ctx.records)
    held, promoted = fw["held_by_firewall_total"], fw["promoted_total"]
    rate = held / (held + promoted) if (held + promoted) else float("nan")
    rate_s = "n/a" if rate != rate else f"{100 * rate:.0f}%"
    _stat_tiles(ax, [
        (held, "rescues HELD\n(circularity firewall)", "#2C7FB8"),
        (promoted, "rescues PROMOTED\n(PAN_REF_RESCUED)", "#1A7F37"),
        (rate_s, "firewall hold-rate\n(held / held+promoted)", "#5D4037"),
    ])
    ax.set_title("Circularity firewall — promotion conservatism", fontsize=10, loc="left", weight="bold")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# --- PAGE 2: two-axis decision space ------------------------------------------
def page2(pdf, ctx, prefix):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("The two-axis decision space — how each verdict was reached",
                 fontsize=13, weight="bold", x=0.07, ha="left")
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.28,
                          left=0.09, right=0.96, top=0.9, bottom=0.1)
    nov = ctx.novel

    # projection grid (centerpiece)
    ax = fig.add_subplot(gs[0, 0])
    rows, cols = D.NOVELTY_ORDER, D.MECHANISM_ORDER
    ct = D.crosstab(nov, "novelty_support", "primary_mechanism")
    mat = [[ct.get((rl, cl), 0) for cl in cols] for rl in rows]
    borders = [[CLASS_COLORS.get(projected_class(rl, cl)) for cl in cols] for rl in rows]
    _heatmap(ax, mat, rows, cols, "Axis A (novelty-support) × Axis B (mechanism)",
             border_colors=borders)
    ax.set_xlabel("artifact mechanism  ·  cell = isoform count, border = class the rulebook projects to",
                  fontsize=7.5)
    ax.set_ylabel("novelty support", fontsize=8)

    # Axis A by structural category
    ax = fig.add_subplot(gs[0, 1])
    cats = sorted({D.get_path(r, "structural_category") for r in nov})
    import numpy as np
    bottoms = np.zeros(len(cats))
    sup_colors = {"SUPPORTED": "#3FAE5A", "PARTIAL": "#8CC63F",
                  "UNSUPPORTED": "#C1272D", "UNKNOWN": "#9E9E9E"}
    for sup in D.NOVELTY_ORDER:
        vals = [sum(1 for r in nov if D.get_path(r, "structural_category") == c
                    and D.get_path(r, "novelty_support") == sup) for c in cats]
        ax.bar(range(len(cats)), vals, bottom=bottoms, label=sup,
               color=sup_colors[sup], edgecolor="white", width=0.7)
        bottoms += np.array(vals)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c.replace("novel_", "n_") for c in cats], fontsize=7, rotation=20, ha="right")
    ax.set_title("Novelty support by structural category", fontsize=10, loc="left", weight="bold")
    ax.legend(fontsize=6.5, ncol=2, loc="upper right", frameon=False)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    # Axis B mechanism distribution, priority order
    ax = fig.add_subplot(gs[1, 0])
    mc = D.count_field(nov, "primary_mechanism")
    present = [m for m in D.MECHANISM_ORDER if mc.get(m, 0) > 0]
    _hbar(ax, present, [mc.get(m, 0) for m in present],
          [MECH_COLORS[m] for m in present], "Primary mechanism (engine priority order)")

    # SR corroboration histogram (conditional)
    ax = fig.add_subplot(gs[1, 1])
    ratios = []
    for r in nov:
        if D.get_path(r, "evidence.sj_evaluable") is True:
            n = D.get_path(r, "evidence.n_novel_junctions", 0) or 0
            k = D.get_path(r, "evidence.n_novel_jx_sr_supported", 0) or 0
            if n > 0:
                ratios.append(k / n)
    if ctx.has_sr and ratios:
        ax.hist(ratios, bins=11, range=(0, 1), color="#3FAE5A", edgecolor="white")
        ax.axvline(0.5, ls="--", color="#C1272D", lw=1)
        ax.set_xlabel("per-isoform novel-junction SR-support ratio", fontsize=8)
        ax.set_title("Short-read corroboration of novel junctions", fontsize=10, loc="left", weight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    else:
        _placeholder(ax, "short-read axis not evaluable\n(no --sj-tab / --ref-gtf)")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# --- PAGE 3: mechanism diagnostics --------------------------------------------
def page3(pdf, ctx, prefix):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Artifact-mechanism diagnostics — why calls were down-weighted",
                 fontsize=13, weight="bold", x=0.07, ha="left")
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.28,
                          left=0.09, right=0.96, top=0.9, bottom=0.12)
    nov = ctx.novel

    # mechanism x confidence class
    ax = fig.add_subplot(gs[0, 0])
    rows = [m for m in D.MECHANISM_ORDER]
    cols = [c for c in CLASS_ORDER]
    ct = D.crosstab(nov, "primary_mechanism", "confidence_class")
    mat = [[ct.get((rl, cl), 0) for cl in cols] for rl in rows]
    keep_r = [i for i, rl in enumerate(rows) if sum(mat[i]) > 0]
    keep_c = [j for j, cl in enumerate(cols) if sum(mat[i][j] for i in range(len(rows))) > 0]
    mat2 = [[mat[i][j] for j in keep_c] for i in keep_r]
    _heatmap(ax, mat2, [rows[i] for i in keep_r], [cols[j] for j in keep_c],
             "Mechanism → resulting confidence class")

    # BAM mapping-evidence quality (conditional)
    ax = fig.add_subplot(gs[0, 1])
    if ctx.has_bam:
        gating, reported = [], []
        for r in ctx.records:
            if D.get_path(r, "evidence.bam_evaluable") is True:
                gating.append(D.get_path(r, "evidence.bam_frac_low_mapq", 0) or 0)
                reported.append(D.get_path(r, "evidence.bam_frac_softclip", 0) or 0)
        ax.hist(gating, bins=12, range=(0, 1), color="#5D4037", alpha=0.85,
                label="low_mapq frac (GATES verdict)", edgecolor="white")
        ax.hist(reported, bins=12, range=(0, 1), color="#D08C3F", alpha=0.5,
                label="softclip frac (reported, NOT verdict-driving)", hatch="//", edgecolor="white")
        ax.set_title("BAM mapping-evidence quality", fontsize=10, loc="left", weight="bold")
        ax.legend(fontsize=6.3, frameon=False)
        ax.set_xlabel("fraction of spanning reads", fontsize=8)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    else:
        _placeholder(ax, "BAM axis not evaluable\nmapping-artifact mechanism could not be assessed")

    # structural-signal evaluability mini-table (three-valued)
    ax = fig.add_subplot(gs[1, 0])
    def tri(field):
        t = sum(1 for r in ctx.records if D.get_path(r, field) is True)
        f = sum(1 for r in ctx.records if D.get_path(r, field) is False)
        n = len(ctx.records) - t - f
        return [str(t), str(f), str(n)]
    rows_t = [["noncanonical motif"] + tri("evidence.noncanonical"),
              ["RT-switch stage"] + tri("evidence.rts_stage"),
              ["caller chain available"] + tri("evidence.chain_available")]
    _table(ax, ["structural prior", "true", "false", "n/a"], rows_t,
           "Structural-signal availability (null = grey, never counted false)",
           col_widths=[0.4, 0.2, 0.2, 0.2])

    # BAM down-weighted register (conditional)
    ax = fig.add_subplot(gs[1, 1])
    mapping_recs = [r for r in ctx.records
                    if D.get_path(r, "primary_mechanism") == "mapping_or_repeat"]
    if ctx.has_bam and mapping_recs:
        rows_r = []
        for r in mapping_recs:
            rows_r.append([
                str(D.get_path(r, "isoform"))[:18],
                D.get_path(r, "confidence_class", "")[:14],
                str(D.get_path(r, "evidence.bam_n_spanning", "")),
                f"{D.get_path(r, 'evidence.bam_frac_low_mapq', 0):.2f}",
                f"{D.get_path(r, 'evidence.bam_frac_supplementary', 0):.2f}",
            ])
        _table(ax, ["isoform", "class", "n_span", "low_mapq", "suppl"], rows_r,
               f"Mapping-down-weighted calls ({len(mapping_recs)})",
               col_widths=[0.3, 0.26, 0.14, 0.15, 0.15])
    else:
        _placeholder(ax, "no mapping_or_repeat verdicts\n(or BAM axis not evaluable)")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# --- PAGE 4: rescue & firewall ------------------------------------------------
def page4(pdf, ctx, prefix):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Reference-bias rescue & circularity firewall — defensibility",
                 fontsize=13, weight="bold", x=0.07, ha="left")
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.28,
                          left=0.09, right=0.96, top=0.9, bottom=0.12)

    # rescue flow grouped bar
    ax = fig.add_subplot(gs[0, 0])
    fw = D.rescue_firewall(ctx.records)
    import numpy as np
    groups, fired, held = [], [], []
    if ctx.has_variant:
        groups.append("variant"); fired.append(fw["variant"]["rescue_fired"]); held.append(fw["variant"]["held_circular"])
    if ctx.has_pangenome:
        groups.append("pangenome"); fired.append(fw["pangenome"]["rescue_fired"]); held.append(fw["pangenome"]["held_circular"])
    if groups:
        x = np.arange(len(groups))
        ax.bar(x - 0.2, fired, 0.4, label="rescue fired", color="#1A7F37", edgecolor="white")
        ax.bar(x + 0.2, held, 0.4, label="held (circular-risk)", color="#2C7FB8", edgecolor="white")
        ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=8)
        ax.set_title("Rescue axis: promoted vs held", fontsize=10, loc="left", weight="bold")
        ax.legend(fontsize=7, frameon=False)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    else:
        _placeholder(ax, "no variant / pangenome axis supplied\n(reference-bias rescue not exercised)")

    # rescue outcome donut: promoted vs held vs none
    ax = fig.add_subplot(gs[0, 1])
    promoted = fw["promoted_total"]; held_total = fw["held_by_firewall_total"]
    other = len(ctx.records) - promoted - held_total
    if promoted or held_total:
        vals = [promoted, held_total, other]
        labs = [f"promoted\n{promoted}", f"held by firewall\n{held_total}", f"no rescue\n{other}"]
        ax.pie(vals, labels=labs, colors=["#1A7F37", "#2C7FB8", "#ECECEC"],
               textprops={"fontsize": 7.5}, wedgeprops={"width": 0.42, "edgecolor": "white"})
        ax.set_title("Rescue outcomes (firewall is the thesis)", fontsize=10, loc="left", weight="bold")
    else:
        _placeholder(ax, "0 reference-bias rescues this run")

    # held-rescue ledger (conditional)
    ax = fig.add_subplot(gs[1, :])
    held_recs = [r for r in ctx.records if D.get_path(r, "circularity_flag") is True]
    if held_recs:
        rows_r = []
        for r in held_recs:
            trace = D.get_path(r, "rule_trace", []) or []
            why = next((t for t in trace if "circular" in t.lower()), trace[-1] if trace else "")
            rows_r.append([
                str(D.get_path(r, "isoform"))[:20],
                D.get_path(r, "structural_category", "")[:20],
                D.get_path(r, "confidence_class", "")[:12],
                why[:70],
            ])
        _table(ax, ["isoform", "structural_category", "class", "why held (rule_trace)"],
               rows_r, f"Held-rescue ledger — circularity firewall ({len(held_recs)})",
               col_widths=[0.16, 0.22, 0.14, 0.48], max_rows=12)
    else:
        _placeholder(ax, "0 circular-risk rescues — firewall inactive this run (nothing held back)")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# --- PAGE 5: consensus (conditional) ------------------------------------------
def page5(pdf, ctx, prefix):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Multi-caller consensus axis", fontsize=13, weight="bold", x=0.07, ha="left")
    fig.text(0.07, 0.935, "Consensus is methodological corroboration — capped at "
             "MEDIUM_CONF_NOVEL, never HIGH.", fontsize=8, style="italic", color="#666666")
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.28,
                          left=0.09, right=0.96, top=0.88, bottom=0.12)
    cev = [r for r in ctx.records if D.get_path(r, "evidence.consensus_evaluable") is True]

    # n_callers distribution
    ax = fig.add_subplot(gs[0, 0])
    import numpy as np
    ncs = [D.get_path(r, "evidence.n_callers", 0) or 0 for r in cev]
    if ncs:
        mx = max(ncs)
        bins = np.arange(0.5, mx + 1.5, 1)
        ax.hist(ncs, bins=bins, color="#2A9D8F", edgecolor="white")
        ax.set_xticks(range(1, mx + 1))
        ax.set_xlabel("# callers supporting the chain", fontsize=8)
        ax.set_title("Caller support distribution", fontsize=10, loc="left", weight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    else:
        _placeholder(ax, "no consensus-evaluable records")

    # n_callers x confidence class
    ax = fig.add_subplot(gs[0, 1])
    kset = sorted({D.get_path(r, "evidence.n_callers", 0) or 0 for r in cev})
    cols = [c for c in CLASS_ORDER if any(D.get_path(r, "confidence_class") == c for r in cev)]
    mat = [[sum(1 for r in cev if (D.get_path(r, "evidence.n_callers", 0) or 0) == k
                and D.get_path(r, "confidence_class") == c) for c in cols] for k in kset]
    _heatmap(ax, mat, [f"{k} callers" for k in kset], cols, "Caller support → confidence class")

    # consensus-driven promotions of UNKNOWN-support novels
    ax = fig.add_subplot(gs[1, :])
    unk = [r for r in cev if D.get_path(r, "novelty_support") == "UNKNOWN"]
    classes = [c for c in CLASS_ORDER if any(D.get_path(r, "confidence_class") == c for r in unk)]
    vals = [sum(1 for r in unk if D.get_path(r, "confidence_class") == c) for c in classes]
    if unk:
        _hbar(ax, classes, vals, [CLASS_COLORS[c] for c in classes],
              "Verdicts of UNKNOWN-support novels under the consensus axis "
              "(promotions to MEDIUM/LOW come from caller agreement)")
    else:
        _placeholder(ax, "no UNKNOWN-support records in the consensus-evaluable set")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# --- PAGE 6: bio_flags pass-through (conditional) -----------------------------
def page6(pdf, ctx, prefix):
    def nonnull(field):
        return D.numeric_values(ctx.records, field)
    cage = nonnull("bio_flags.dist_to_CAGE_peak")
    polya = nonnull("bio_flags.dist_to_polyA_site")
    nmd_vals = [D.get_path(r, "bio_flags.predicted_NMD") for r in ctx.records]
    motif_vals = [D.get_path(r, "bio_flags.polyA_motif_found") for r in ctx.records]
    has_nmd = any(v not in (None, "NA") for v in nmd_vals)
    has_motif = any(v not in (None, "NA") for v in motif_vals)
    if not (cage or polya or has_nmd or has_motif):
        return  # nothing supplied; skip the page entirely

    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("SQANTI3 biological priors (pass-through)", fontsize=13, weight="bold", x=0.07, ha="left")
    fig.text(0.07, 0.935, "These are SQANTI3 annotations passed through verbatim and are "
             "NOT used in the PanIsoGuard verdict.", fontsize=8, style="italic", color="#C1272D")
    gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.28,
                          left=0.09, right=0.96, top=0.88, bottom=0.12)

    ax = fig.add_subplot(gs[0, 0])
    if cage:
        ax.hist([min(max(v, -500), 500) for v in cage], bins=30, color="#3FAE5A", edgecolor="white")
        ax.set_title("dist_to_CAGE_peak (clipped ±500)", fontsize=10, loc="left", weight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    else:
        _placeholder(ax, "no CAGE prior supplied")

    ax = fig.add_subplot(gs[0, 1])
    if polya:
        ax.hist([min(max(v, -500), 500) for v in polya], bins=30, color="#8CC63F", edgecolor="white")
        ax.set_title("dist_to_polyA_site (clipped ±500)", fontsize=10, loc="left", weight="bold")
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    else:
        _placeholder(ax, "no polyA-site prior supplied")

    ax = fig.add_subplot(gs[1, 0])
    if has_nmd:
        from collections import Counter
        c = Counter(v for v in nmd_vals if v not in (None, "NA"))
        labs = list(c.keys())
        _hbar(ax, labs, [c[k] for k in labs], ["#5D4037"] * len(labs), "predicted_NMD")
    else:
        _placeholder(ax, "no predicted_NMD prior")

    ax = fig.add_subplot(gs[1, 1])
    if has_motif:
        from collections import Counter
        c = Counter(v for v in motif_vals if v not in (None, "NA"))
        labs = list(c.keys())
        _hbar(ax, labs, [c[k] for k in labs], ["#D08C3F"] * len(labs), "polyA_motif_found")
    else:
        _placeholder(ax, "no polyA_motif prior")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# --- PAGE 7: provenance & defensibility appendix ------------------------------
def page7(pdf, ctx, prefix):
    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle("Methods, provenance & defensibility appendix",
                 fontsize=13, weight="bold", x=0.07, ha="left")
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.25,
                          left=0.07, right=0.96, top=0.9, bottom=0.08)
    m = ctx.prov["meta"]; ax_states = ctx.prov["axes"]

    # provenance table
    ax = fig.add_subplot(gs[:, 0])
    keys = ["tool_version", "ruleset_version", "sqanti3_version_target", "config",
            "classification", "isoforms", "ref_gtf", "sj_tab", "bam"]
    rows_r = [[k, str(m.get(k, "—"))[:46]] for k in keys if k in m]
    for a, st in ax_states.items():
        rows_r.append([f"axis.{a}", st[:46]])
    _table(ax, ["key", "value"], rows_r, "Run provenance", col_widths=[0.32, 0.68], max_rows=30)

    # rule-trace frequency (top reasons)
    ax = fig.add_subplot(gs[0, 1])
    from collections import Counter
    rc = Counter()
    for r in ctx.records:
        for t in (D.get_path(r, "rule_trace", []) or []):
            key = t.split(" reason=")[0].split(" -> ")[-1] if "->" in t else t
            rc[t[:42]] += 1
    top = rc.most_common(8)
    if top:
        _hbar(ax, [t for t, _ in top][::-1], [c for _, c in top][::-1],
              ["#7BA6C9"] * len(top), "Most frequent rule-trace lines")

    # defensibility checklist
    ax = fig.add_subplot(gs[1, 1])
    fw = D.rescue_firewall(ctx.records)
    checks = [
        ("✓" if not ctx.mismatch else "✗",
         "class tally reconciles (provenance ↔ JSONL)",
         "#1A7F37" if not ctx.mismatch else "#C1272D"),
        ("✓", f"firewall held {fw['held_by_firewall_total']} circular-risk rescue(s)", "#1A7F37"),
        ("ℹ", f"consensus axis {'used' if ctx.has_consensus else 'not used'} this run", "#666666"),
        ("ℹ", f"axes in play: " + ", ".join(
            n for n, on in [("SR", ctx.has_sr), ("BAM", ctx.has_bam),
                            ("variant", ctx.has_variant), ("pangenome", ctx.has_pangenome),
                            ("consensus", ctx.has_consensus)] if on) or "catalog only",
         "#666666"),
        ("ℹ", "SQANTI3 priors consumed, never recomputed", "#666666"),
    ]
    _callout(ax, "Defensibility checklist", [(f"{mark}  {txt}", col) for mark, txt, col in checks])

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


def build_report(prefix, out_pdf=None):
    jsonl = prefix + ".attribution.jsonl"
    prov_log = prefix + ".provenance.log"
    if not os.path.exists(jsonl):
        raise FileNotFoundError(f"missing {jsonl} (run `panisoguard adjudicate --out-prefix {prefix}` first)")
    records = D.load_records(jsonl)
    prov = D.load_provenance(prov_log)
    ctx = Ctx(records, prov)
    out_pdf = out_pdf or (prefix + ".report.pdf")
    with PdfPages(out_pdf) as pdf:
        page1(pdf, ctx, prefix)
        page2(pdf, ctx, prefix)
        page3(pdf, ctx, prefix)
        page4(pdf, ctx, prefix)
        if ctx.has_consensus:
            page5(pdf, ctx, prefix)
        page6(pdf, ctx, prefix)  # self-skips if no bio priors
        page7(pdf, ctx, prefix)
    return out_pdf, len(records)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Render a PDF report from PanIsoGuard adjudication output.")
    ap.add_argument("--prefix", required=True,
                    help="adjudicate --out-prefix value (reads <prefix>.attribution.jsonl + .provenance.log)")
    ap.add_argument("--out", default=None, help="output PDF (default <prefix>.report.pdf)")
    args = ap.parse_args(argv)
    out, n = build_report(args.prefix, args.out)
    sys.stderr.write(f"wrote {out} ({n} isoforms)\n")


if __name__ == "__main__":
    main()
