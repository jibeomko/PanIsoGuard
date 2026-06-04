#!/usr/bin/env python3
"""PanIsoGuard adjudication report -- a publication-grade multi-page PDF summary.

Companion tool (keeps the C++ core dependency-free): an optional Python post-process, in
the spirit of -- but independent from, and copying no design or code of -- SQANTI3's
separate report step. Consumes the files `panisoguard adjudicate` already writes -- it
recomputes nothing about the isoforms, it only visualizes the verdicts and the evidence
behind them:

    <prefix>.attribution.jsonl   per-isoform verdict + evidence + rule_trace (richest)
    <prefix>.provenance.log      run metadata, axis on/off, authoritative class tally

Output: <prefix>.report.pdf (landscape, one figure per page).

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
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D
    from matplotlib.patches import FancyBboxPatch, Rectangle
except ImportError:  # pragma: no cover
    sys.stderr.write(
        "panisoguard.report needs matplotlib (optional dependency).\n"
        "  pip install matplotlib   (or: conda install matplotlib)\n")
    raise

# ---------------------------------------------------------------------------
# Design system: one semantic palette + human-readable labels, learned once.
# ---------------------------------------------------------------------------
INK = "#1F2933"          # near-black for titles
INK2 = "#52606D"         # secondary text
MUTE = "#9AA5B1"         # tertiary / captions
HAIR = "#E4E7EB"         # hairlines / gridlines
PANEL_BG = "#FFFFFF"
ACCENT = "#2D7DD2"       # brand accent (rule under titles)

CLASS_COLORS = {
    "HIGH_CONF_KNOWN": "#0E7C66",            # teal — trusted catalog pass-through
    "HIGH_CONF_NOVEL": "#2E9E5B",            # green — confident novel
    "MEDIUM_CONF_NOVEL": "#86BC4C",          # light green — confident novel
    "PAN_REF_RESCUED_FALSE_NOVEL": "#2D7DD2",  # blue — reclassification, not a grade
    "LOW_CONF_PARTIAL": "#E8A33D",           # amber — tentative
    "AMBIGUOUS": "#9AA5B1",                  # grey — honest abstention
    "ARTIFACT": "#D1495B",                   # red — reject
}
CLASS_ORDER = [
    "HIGH_CONF_KNOWN", "HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL",
    "PAN_REF_RESCUED_FALSE_NOVEL", "LOW_CONF_PARTIAL", "AMBIGUOUS", "ARTIFACT",
]
CLASS_LABEL = {
    "HIGH_CONF_KNOWN": "High-confidence known",
    "HIGH_CONF_NOVEL": "High-confidence novel",
    "MEDIUM_CONF_NOVEL": "Medium-confidence novel",
    "PAN_REF_RESCUED_FALSE_NOVEL": "Reference-bias rescue",
    "LOW_CONF_PARTIAL": "Low-confidence / partial",
    "AMBIGUOUS": "Ambiguous",
    "ARTIFACT": "Likely artifact",
}
MECH_COLORS = {
    "none": "#CBD2D9", "mapping_or_repeat": "#6D4C41", "noncanonical": "#C77F33",
    "rt_switch": "#E0A93B", "degradation": "#EBD08A",
    "variant_created": "#2D7DD2", "population_known": "#5BA3D9",
}
MECH_LABEL = {
    "none": "No artifact mechanism", "mapping_or_repeat": "Mapping / repeat",
    "noncanonical": "Non-canonical motif", "rt_switch": "RT switching",
    "degradation": "Intra-priming", "variant_created": "Variant-created",
    "population_known": "Population-known",
}
SUPPORT_LABEL = {
    "SUPPORTED": "Supported", "PARTIAL": "Partial",
    "UNSUPPORTED": "Unsupported", "UNKNOWN": "Not evaluable",
}
SUPPORT_COLORS = {"SUPPORTED": "#2E9E5B", "PARTIAL": "#86BC4C",
                  "UNSUPPORTED": "#D1495B", "UNKNOWN": "#CBD2D9"}
AXIS_ON = "#17A2A2"
AXIS_OFF = "#CBD2D9"
NOVEL_CATS = {"novel_in_catalog", "novel_not_in_catalog"}
# The (support × mechanism) projection grid covers only ARTIFACT mechanisms. The two
# rescue mechanisms (variant_created / population_known) are set in the engine's early
# reference-bias-rescue branch and bypass the grid entirely, so they are shown on the
# rescue/firewall page, not here (including them would mis-colour the grid borders).
ARTIFACT_MECHS = ["none", "mapping_or_repeat", "noncanonical", "rt_switch", "degradation"]
SEQ_CMAP = LinearSegmentedColormap.from_list("pig_seq", ["#F4F8FB", "#2D7DD2", "#16456E"])


def _setup_style():
    plt.rcParams.update({
        "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"],
        "font.size": 9.0,
        "axes.titlesize": 11.5, "axes.titleweight": "bold",
        "axes.labelsize": 8.5, "axes.labelcolor": INK2,
        "axes.edgecolor": HAIR, "axes.linewidth": 1.0,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "text.color": INK, "figure.facecolor": PANEL_BG, "axes.facecolor": PANEL_BG,
        "axes.grid": False, "savefig.facecolor": PANEL_BG, "figure.dpi": 150,
    })


def clabel(c):
    return CLASS_LABEL.get(c, c)


# ---------------------------------------------------------------------------
# Drawing primitives (clean, consistent, despined)
# ---------------------------------------------------------------------------
def _despine(ax, keep=("left", "bottom")):
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(s in keep)
    ax.tick_params(length=0)


def _axtitle(ax, text, sub=None):
    # point-based offsets so title/subtitle never overlap regardless of panel height
    ax.set_title(text, loc="left", color=INK, pad=21 if sub else 7)
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords="axes fraction", xytext=(0, 3),
                    textcoords="offset points", ha="left", va="bottom",
                    fontsize=7.6, color=MUTE)


def _page_header(fig, title, n, subtitle=None):
    fig.text(0.045, 0.945, title, fontsize=16, weight="bold", color=INK, ha="left", va="center")
    if subtitle:
        fig.text(0.045, 0.905, subtitle, fontsize=8.6, color=MUTE, ha="left", va="center")
    fig.add_artist(Line2D([0.045, 0.955], [0.882, 0.882], color=ACCENT, lw=2.2,
                          solid_capstyle="round"))
    fig.text(0.955, 0.945, f"{n}", fontsize=11, color=MUTE, ha="right", va="center", weight="bold")
    fig.text(0.955, 0.92, "/ 7", fontsize=7.5, color="#C3CAD2", ha="right", va="center")


def _footer(fig, prov, prefix):
    m = prov["meta"]
    fig.add_artist(Line2D([0.045, 0.955], [0.045, 0.045], color=HAIR, lw=0.8))
    fig.text(0.045, 0.028, "PanIsoGuard report", fontsize=7, color=MUTE, ha="left", weight="bold")
    txt = (f"{m.get('tool_version', '?')} · ruleset {m.get('ruleset_version', '?')} · "
           f"SQANTI3 target {m.get('sqanti3_version_target', '?')} · {os.path.basename(prefix)}")
    fig.text(0.955, 0.028, txt, fontsize=7, color=MUTE, ha="right")


def _placeholder(ax, msg):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    box = FancyBboxPatch((0.06, 0.22), 0.88, 0.56, transform=ax.transAxes,
                         boxstyle="round,pad=0.01,rounding_size=0.03",
                         facecolor="#F7F9FB", edgecolor=HAIR, lw=1.0)
    ax.add_patch(box)
    ax.text(0.5, 0.56, "—", transform=ax.transAxes, ha="center", va="center",
            fontsize=15, color="#C3CAD2", weight="bold")
    ax.text(0.5, 0.40, msg, transform=ax.transAxes, ha="center", va="center",
            fontsize=8.2, color=MUTE, wrap=True)


def _grid_x(ax):
    ax.set_axisbelow(True)
    ax.xaxis.grid(True, color=HAIR, lw=0.8)


def _hbar(ax, labels, values, colors, title, sub=None, fmt_total=None, value_fmt="{:,}"):
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, height=0.66, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.invert_yaxis()
    vmax = max(values) if values and max(values) > 0 else 1
    for yi, v in zip(y, values):
        lab = value_fmt.format(v)
        if fmt_total:
            lab += f"  ({100 * v / fmt_total:.0f}%)"
        ax.text(v + vmax * 0.015, yi, lab, va="center", ha="left", fontsize=8, color=INK2)
    ax.set_xlim(0, vmax * 1.16)
    _despine(ax, keep=("left",)); _grid_x(ax)
    ax.set_xticks([])
    _axtitle(ax, title, sub)


def _lollipop(ax, rows, title, sub=None):
    """rows: list of (label, frac, k, n, on). A 0..1 track with a marker = cleaner gauge."""
    y = np.arange(len(rows))
    for i, (lab, frac, k, n, on) in enumerate(rows):
        ax.plot([0, 1], [i, i], color=HAIR, lw=3, solid_capstyle="round", zorder=1)
        col = AXIS_ON if on else AXIS_OFF
        ax.plot([0, frac], [i, i], color=col, lw=3, solid_capstyle="round", zorder=2)
        ax.scatter([frac], [i], s=70, color=col, zorder=3, edgecolor="white", lw=1.2)
        ax.text(1.06, i, f"{k:,}/{n:,}", va="center", ha="left", fontsize=8, color=INK2)
        if not on:
            ax.text(frac, i + 0.32, "not evaluable", va="bottom", ha="left",
                    fontsize=6.6, color=MUTE, style="italic")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows])
    ax.invert_yaxis(); ax.set_xlim(-0.02, 1.32); ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xticks([0, 0.5, 1.0]); ax.set_xticklabels(["0%", "50%", "100%"], fontsize=7)
    _despine(ax, keep=("bottom",))
    _axtitle(ax, title, sub)


def _ribbon(ax, segments, title, sub=None):
    """A single 100%-stacked bar with each segment labelled directly beneath it —
    at-a-glance composition. segments: [(label, count, color)]."""
    total = sum(c for _, c, _ in segments) or 1
    left = 0.0
    for lab, count, color in segments:
        w = count / total
        if w <= 0:
            continue
        ax.barh(0, w, left=left, color=color, height=0.55, zorder=3)
        ax.text(left + w / 2, 0.06, f"{count}", va="center", ha="center",
                fontsize=12, color="white", weight="bold", zorder=4)
        if w > 0.045:
            ax.text(left + w / 2, -0.52, f"{lab}\n{100 * count / total:.0f}%", va="top",
                    ha="center", fontsize=8, color=INK2, linespacing=1.3)
        left += w
    ax.set_xlim(0, 1); ax.set_ylim(-1.25, 0.45)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    _axtitle(ax, title, sub)


def _stat_cards(ax, cards, title=None):
    """cards: list of (value, label, color). Modern KPI cards in a row."""
    ax.axis("off")
    if title:
        ax.text(0, 1.06, title, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=11.5, weight="bold", color=INK)
    n = len(cards)
    gap = 0.035
    w = (1 - gap * (n - 1)) / n
    for i, (val, label, color) in enumerate(cards):
        x0 = i * (w + gap)
        card = FancyBboxPatch((x0, 0.12), w, 0.74, transform=ax.transAxes,
                              boxstyle="round,pad=0,rounding_size=0.04",
                              facecolor=color, alpha=0.10, edgecolor=color, lw=1.3)
        ax.add_patch(card)
        ax.add_patch(Rectangle((x0, 0.12), 0.012, 0.74, transform=ax.transAxes, color=color))
        ax.text(x0 + w / 2, 0.60, str(val), transform=ax.transAxes, ha="center", va="center",
                fontsize=23, weight="bold", color=color)
        for j, line in enumerate(label.split("\n")):
            ax.text(x0 + w / 2, 0.34 - j * 0.11, line, transform=ax.transAxes, ha="center",
                    va="center", fontsize=7.8, color=INK2)


def _callout(ax, title, lines):
    ax.axis("off")
    ax.text(0, 1.0, title, transform=ax.transAxes, ha="left", va="top",
            fontsize=11.5, weight="bold", color=INK)
    y = 0.80
    for txt, color, big in lines:
        ax.text(0.01, y, txt, transform=ax.transAxes, ha="left", va="top",
                fontsize=10.5 if big else 9, color=color, weight="bold" if big else "normal")
        y -= 0.15 if big else 0.125


def _heatmap(ax, matrix, row_labels, col_labels, title, border_colors=None, sub=None):
    arr = np.array(matrix, dtype=float)
    vmax = arr.max() if arr.size and arr.max() > 0 else 1
    ax.imshow(arr, cmap=SEQ_CMAP, aspect="auto", vmin=0, vmax=vmax)
    # white cell separators for a crisp grid
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", lw=2.2)
    ax.tick_params(which="minor", length=0)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, fontsize=7.4, rotation=28, ha="right")
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            v = int(arr[i, j])
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=8,
                        color="white" if arr[i, j] > vmax * 0.5 else INK,
                        weight="bold" if arr[i, j] > vmax * 0.5 else "normal")
            if border_colors is not None and border_colors[i][j] is not None:
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=border_colors[i][j], lw=2.6, zorder=5))
    ax.tick_params(length=0)
    _axtitle(ax, title, sub)


def _table(ax, col_labels, rows, title, col_widths=None, max_rows=12, sub=None):
    ax.axis("off")
    _axtitle(ax, title, sub)
    shown = rows[:max_rows]
    tbl = ax.table(cellText=shown, colLabels=col_labels, loc="upper center",
                   cellLoc="left", colWidths=col_widths, bbox=[0, 0, 1, 0.9])
    tbl.auto_set_font_size(False); tbl.set_fontsize(7.4)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("white"); cell.set_linewidth(1.0)
        cell.set_height(0.10)
        if r == 0:
            cell.set_facecolor(ACCENT); cell.set_text_props(weight="bold", color="white")
        else:
            cell.set_facecolor("#F4F8FB" if r % 2 else "#FFFFFF")
            cell.set_text_props(color=INK2)
    if len(rows) > max_rows:
        ax.text(0.5, -0.04, f"… {len(rows) - max_rows} more (see attribution.jsonl)",
                transform=ax.transAxes, ha="center", fontsize=6.8, color=MUTE)


# ---------------------------------------------------------------------------
# Projection rulebook (Page 2 heatmap borders)
# ---------------------------------------------------------------------------
def projected_class(support, mech):
    has_mech = mech not in (None, "none")
    if support == "SUPPORTED":
        return "HIGH_CONF_NOVEL" if not has_mech else "MEDIUM_CONF_NOVEL"
    if support == "PARTIAL":
        return "MEDIUM_CONF_NOVEL" if not has_mech else "LOW_CONF_PARTIAL"
    if support == "UNSUPPORTED":
        return "LOW_CONF_PARTIAL" if not has_mech else "ARTIFACT"
    if support == "UNKNOWN":
        return "ARTIFACT" if mech == "mapping_or_repeat" else "AMBIGUOUS"
    return "AMBIGUOUS"


# ---------------------------------------------------------------------------
class Ctx:
    def __init__(self, records, prov):
        self.records = records
        self.prov = prov
        self.novel = [r for r in records
                      if D.get_path(r, "structural_category") in NOVEL_CATS]
        self.has_consensus = any(D.get_path(r, "evidence.consensus_evaluable") is True for r in records)
        self.has_bam = any(D.get_path(r, "evidence.bam_evaluable") is True for r in records)
        self.has_variant = any(D.get_path(r, "evidence.variant_evaluable") is True for r in records)
        self.has_pangenome = any(D.get_path(r, "evidence.pangenome_evaluable") is True for r in records)
        self.has_sr = any(D.get_path(r, "evidence.sj_evaluable") is True for r in records)
        cc = D.count_field(records, "confidence_class")
        self.jsonl_counts = {k: cc.get(k, 0) for k in CLASS_ORDER}
        pc = prov["class_counts"]
        self.mismatch = bool(pc) and any(
            pc.get(k, 0) != self.jsonl_counts.get(k, 0) for k in set(pc) | set(self.jsonl_counts))


def _newpage():
    return plt.figure(figsize=(11.69, 8.27))  # A4 landscape


# ---------------------------------------------------------------------------
# PAGE 1 — executive summary
# ---------------------------------------------------------------------------
def page1(pdf, ctx, prefix):
    fig = _newpage()
    _page_header(fig, "Adjudication summary", 1,
                 "Orthogonal confidence layer — SQANTI3 features are consumed as priors, never recomputed.")
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.30, left=0.16, right=0.95, top=0.80, bottom=0.10)

    nov = ctx.novel
    def cc_in(rs, names):
        return sum(1 for r in rs if D.get_path(r, "confidence_class") in names)
    trust = cc_in(nov, {"HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL"})
    look = cc_in(nov, {"LOW_CONF_PARTIAL", "AMBIGUOUS"})
    art = cc_in(nov, {"ARTIFACT"})
    resc = cc_in(nov, {"PAN_REF_RESCUED_FALSE_NOVEL"})

    # headline ribbon (top, full width feel)
    ax = fig.add_subplot(gs[0, :])
    _ribbon(ax, [("Trustworthy", trust, "#2E9E5B"), ("Needs a look", look, "#E8A33D"),
                 ("Artifact", art, "#D1495B"), ("Rescue", resc, "#2D7DD2")],
            f"Verdict on {len(nov):,} novel isoforms",
            "Trustworthy = high/medium-confidence novel · the share you can act on at a glance")

    # axis availability lollipop
    ax = fig.add_subplot(gs[1, 0])
    rows = [(lab, (k / n if n else 0), k, n, k > 0) for lab, k, n in D.axis_evaluability(ctx.records)]
    _lollipop(ax, rows, "Evidence axes in play",
              "fraction of isoforms for which each axis could be evaluated")

    # confidence-class distribution (clean, full labels)
    ax = fig.add_subplot(gs[1, 1])
    counts = ctx.prov["class_counts"] or ctx.jsonl_counts
    present = [c for c in CLASS_ORDER if counts.get(c, 0) > 0]
    total = ctx.prov["total"] or sum(counts.values()) or 1
    _hbar(ax, [clabel(c) for c in present], [counts.get(c, 0) for c in present],
          [CLASS_COLORS[c] for c in present], "Confidence-class distribution",
          "every isoform · count (share of all)", fmt_total=total)
    if ctx.mismatch:
        ax.text(1.0, 1.10, "⚠ provenance tally ≠ JSONL", transform=ax.transAxes,
                ha="right", color="#D1495B", fontsize=8, weight="bold")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# PAGE 2 — two-axis decision space
# ---------------------------------------------------------------------------
def page2(pdf, ctx, prefix):
    fig = _newpage()
    _page_header(fig, "How each verdict was reached", 2,
                 "Every novel isoform is scored on two axes; the pair projects to one confidence class.")
    gs = fig.add_gridspec(2, 2, hspace=0.72, wspace=0.34, left=0.18, right=0.95, top=0.80, bottom=0.13)
    nov = ctx.novel

    ax = fig.add_subplot(gs[0, 0])
    # Rescue records (variant_created / population_known) bypass the projection grid in
    # the engine, so they are excluded here (shown on the rescue page instead).
    grid_recs = [r for r in nov if D.get_path(r, "primary_mechanism") in ARTIFACT_MECHS]
    rows, cols = D.NOVELTY_ORDER, ARTIFACT_MECHS
    ct = D.crosstab(grid_recs, "novelty_support", "primary_mechanism")
    mat = [[ct.get((rl, cl), 0) for cl in cols] for rl in rows]
    borders = [[CLASS_COLORS.get(projected_class(rl, cl)) for cl in cols] for rl in rows]
    _heatmap(ax, mat, [SUPPORT_LABEL[r] for r in rows], [MECH_LABEL[c] for c in cols],
             "Novelty support  ×  artifact mechanism",
             border_colors=borders, sub="cell = isoform count · coloured border = the class this pair projects to")

    ax = fig.add_subplot(gs[0, 1])
    cats = sorted({D.get_path(r, "structural_category") for r in nov})
    bottoms = np.zeros(len(cats))
    for sup in D.NOVELTY_ORDER:
        vals = [sum(1 for r in nov if D.get_path(r, "structural_category") == c
                    and D.get_path(r, "novelty_support") == sup) for c in cats]
        ax.bar(range(len(cats)), vals, bottom=bottoms, label=SUPPORT_LABEL[sup],
               color=SUPPORT_COLORS[sup], width=0.62, zorder=3)
        bottoms += np.array(vals)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([c.replace("novel_in_catalog", "NIC").replace("novel_not_in_catalog", "NNC")
                        for c in cats], fontsize=8)
    _despine(ax); ax.set_axisbelow(True); ax.yaxis.grid(True, color=HAIR, lw=0.8)
    ax.legend(fontsize=7, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.10), frameon=False)
    _axtitle(ax, "Novelty support by category", "is novel-junction support better for NIC or NNC?")

    ax = fig.add_subplot(gs[1, 0])
    mc = D.count_field(nov, "primary_mechanism")
    present = [m for m in ARTIFACT_MECHS if mc.get(m, 0) > 0]
    _hbar(ax, [MECH_LABEL[m] for m in present], [mc.get(m, 0) for m in present],
          [MECH_COLORS[m] for m in present], "Dominant artifact mechanism",
          "engine priority: mapping > non-canonical > RT-switch > intra-priming")

    ax = fig.add_subplot(gs[1, 1])
    ratios = []
    for r in nov:
        if D.get_path(r, "evidence.sj_evaluable") is True:
            n = D.get_path(r, "evidence.n_novel_junctions", 0) or 0
            k = D.get_path(r, "evidence.n_novel_jx_sr_supported", 0) or 0
            if n > 0:
                ratios.append(k / n)
    if ctx.has_sr and ratios:
        ax.hist(ratios, bins=11, range=(0, 1), color="#2E9E5B", zorder=3, rwidth=0.92)
        ax.axvline(0.5, ls=(0, (4, 3)), color="#D1495B", lw=1.3)
        ax.set_xlabel("fraction of an isoform's novel junctions backed by short reads")
        _despine(ax); ax.set_axisbelow(True); ax.yaxis.grid(True, color=HAIR, lw=0.8)
        _axtitle(ax, "Short-read corroboration", "to the right of the dashed line = majority backed")
    else:
        _placeholder(ax, "Short-read axis not evaluable\n(no --sj-tab / --ref-gtf supplied)")
        _axtitle(ax, "Short-read corroboration")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# PAGE 3 — mechanism diagnostics
# ---------------------------------------------------------------------------
def page3(pdf, ctx, prefix):
    fig = _newpage()
    _page_header(fig, "Why calls were down-weighted", 3,
                 "Which artifact mechanism drove each verdict, and whether the evidence supports it.")
    gs = fig.add_gridspec(2, 2, hspace=0.85, wspace=0.34, left=0.17, right=0.95, top=0.80, bottom=0.16)
    nov = ctx.novel

    ax = fig.add_subplot(gs[0, 0])
    rows = list(D.MECHANISM_ORDER); cols = list(CLASS_ORDER)
    ct = D.crosstab(nov, "primary_mechanism", "confidence_class")
    mat = [[ct.get((rl, cl), 0) for cl in cols] for rl in rows]
    keep_r = [i for i in range(len(rows)) if sum(mat[i]) > 0]
    keep_c = [j for j in range(len(cols)) if sum(mat[i][j] for i in range(len(rows))) > 0]
    mat2 = [[mat[i][j] for j in keep_c] for i in keep_r]
    _heatmap(ax, mat2, [MECH_LABEL[rows[i]] for i in keep_r], [clabel(cols[j]) for j in keep_c],
             "Mechanism → resulting class", sub="where each artifact mechanism ends up")

    ax = fig.add_subplot(gs[0, 1])
    if ctx.has_bam:
        gating, reported = [], []
        for r in ctx.records:
            if D.get_path(r, "evidence.bam_evaluable") is True:
                gating.append(D.get_path(r, "evidence.bam_frac_low_mapq", 0) or 0)
                reported.append(D.get_path(r, "evidence.bam_frac_softclip", 0) or 0)
        ax.hist(gating, bins=12, range=(0, 1), color="#6D4C41", alpha=0.9,
                label="low-MAPQ fraction — gates the verdict", zorder=3, rwidth=0.95)
        ax.hist(reported, bins=12, range=(0, 1), color="#E0A93B", alpha=0.55, hatch="///",
                label="soft-clip fraction — reported only, not used", zorder=2, rwidth=0.95)
        ax.set_xlabel("fraction of junction-spanning reads")
        _despine(ax); ax.set_axisbelow(True); ax.yaxis.grid(True, color=HAIR, lw=0.8)
        ax.legend(fontsize=6.6, frameon=False, loc="upper right")
        _axtitle(ax, "BAM mapping-evidence quality", "are mapping artifacts backed by poor alignment?")
    else:
        _placeholder(ax, "BAM axis not evaluable\nmapping-artifact mechanism could not be assessed")
        _axtitle(ax, "BAM mapping-evidence quality")

    ax = fig.add_subplot(gs[1, 0])
    def tri(field):
        t = sum(1 for r in ctx.records if D.get_path(r, field) is True)
        f = sum(1 for r in ctx.records if D.get_path(r, field) is False)
        return [str(t), str(f), str(len(ctx.records) - t - f)]
    rows_t = [["Non-canonical motif"] + tri("evidence.noncanonical"),
              ["RT-switch stage"] + tri("evidence.rts_stage"),
              ["Caller chain available"] + tri("evidence.chain_available")]
    _table(ax, ["Structural prior", "true", "false", "n/a"], rows_t,
           "Structural-signal availability", col_widths=[0.46, 0.18, 0.18, 0.18],
           sub="null is shown as n/a — never silently counted as false")

    ax = fig.add_subplot(gs[1, 1])
    mapping_recs = [r for r in ctx.records if D.get_path(r, "primary_mechanism") == "mapping_or_repeat"]
    if ctx.has_bam and mapping_recs:
        rows_r = [[str(D.get_path(r, "isoform"))[:16], clabel(D.get_path(r, "confidence_class"))[:16],
                   str(D.get_path(r, "evidence.bam_n_spanning", "")),
                   f"{D.get_path(r, 'evidence.bam_frac_low_mapq', 0):.2f}",
                   f"{D.get_path(r, 'evidence.bam_frac_supplementary', 0):.2f}"] for r in mapping_recs]
        _table(ax, ["isoform", "class", "reads", "low-MAPQ", "suppl."], rows_r,
               f"Mapping-down-weighted calls ({len(mapping_recs)})",
               col_widths=[0.26, 0.30, 0.13, 0.16, 0.15])
    else:
        _placeholder(ax, "No mapping / repeat verdicts\n(or BAM axis not evaluable)")
        _axtitle(ax, "Mapping-down-weighted calls")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# PAGE 4 — rescue & firewall
# ---------------------------------------------------------------------------
def page4(pdf, ctx, prefix):
    fig = _newpage()
    _page_header(fig, "Reference-bias rescue & circularity firewall", 4,
                 "A novel junction explained by reference bias is reclassified — but only on independent evidence.")
    gs = fig.add_gridspec(2, 2, hspace=0.75, wspace=0.30, left=0.10, right=0.95, top=0.80, bottom=0.14,
                          height_ratios=[1.0, 0.9])
    fw = D.rescue_firewall(ctx.records)

    # scorecards (top, full width)
    ax = fig.add_subplot(gs[0, :])
    held, promoted = fw["held_by_firewall_total"], fw["promoted_total"]
    rate = held / (held + promoted) if (held + promoted) else float("nan")
    _stat_cards(ax, [
        (promoted, "rescues PROMOTED\nreclassified as reference-bias", "#2E9E5B"),
        (held, "rescues HELD by firewall\ncircular-risk provenance — not promoted", "#2D7DD2"),
        ("n/a" if rate != rate else f"{100 * rate:.0f}%", "firewall hold-rate\nheld / (held + promoted)", "#6D4C41"),
    ], title="The thesis: nothing is promoted on circular evidence")

    ax = fig.add_subplot(gs[1, 0])
    groups, promotedb, heldb = [], [], []
    if ctx.has_variant:
        groups.append("Variant"); promotedb.append(fw["variant"]["promoted"]); heldb.append(fw["variant"]["held_circular"])
    if ctx.has_pangenome:
        groups.append("Pangenome"); promotedb.append(fw["pangenome"]["promoted"]); heldb.append(fw["pangenome"]["held_circular"])
    if groups:
        x = np.arange(len(groups))
        ax.bar(x - 0.19, promotedb, 0.36, label="promoted", color="#2E9E5B", zorder=3)
        ax.bar(x + 0.19, heldb, 0.36, label="held (circular-risk)", color="#2D7DD2", zorder=3)
        ax.set_xticks(x); ax.set_xticklabels(groups)
        _despine(ax); ax.set_axisbelow(True); ax.yaxis.grid(True, color=HAIR, lw=0.8)
        ax.legend(fontsize=7.5, frameon=False)
        _axtitle(ax, "Rescue axis: promoted vs held")
    else:
        _placeholder(ax, "No variant / pangenome axis supplied\n(reference-bias rescue not exercised)")
        _axtitle(ax, "Rescue axis: promoted vs held")

    ax = fig.add_subplot(gs[1, 1])
    held_recs = [r for r in ctx.records if D.get_path(r, "circularity_flag") is True]
    if held_recs:
        rows_r = []
        for r in held_recs:
            trace = D.get_path(r, "rule_trace", []) or []
            why = next((t for t in trace if "circular" in t.lower()), trace[-1] if trace else "")
            rows_r.append([str(D.get_path(r, "isoform"))[:14],
                           D.get_path(r, "structural_category", "").replace("novel_", "")[:14],
                           why[:54]])
        _table(ax, ["isoform", "category", "why held (rule trace)"], rows_r,
               f"Held-rescue ledger ({len(held_recs)})", col_widths=[0.18, 0.20, 0.62], max_rows=9)
    else:
        _placeholder(ax, "0 circular-risk rescues\nfirewall inactive this run — nothing held back")
        _axtitle(ax, "Held-rescue ledger")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# PAGE 5 — consensus (conditional)
# ---------------------------------------------------------------------------
def page5(pdf, ctx, prefix):
    fig = _newpage()
    _page_header(fig, "Multi-caller consensus", 5,
                 "Caller agreement is methodological corroboration — it promotes to MEDIUM at most, never HIGH.")
    gs = fig.add_gridspec(2, 2, hspace=0.78, wspace=0.34, left=0.16, right=0.95, top=0.80, bottom=0.14)
    cev = [r for r in ctx.records if D.get_path(r, "evidence.consensus_evaluable") is True]

    ax = fig.add_subplot(gs[0, 0])
    ncs = [D.get_path(r, "evidence.n_callers", 0) or 0 for r in cev]
    if ncs:
        mx = max(ncs)
        ax.hist(ncs, bins=np.arange(0.5, mx + 1.5, 1), color="#17A2A2", zorder=3, rwidth=0.86)
        ax.set_xticks(range(1, mx + 1)); ax.set_xlabel("number of callers recovering the chain")
        _despine(ax); ax.set_axisbelow(True); ax.yaxis.grid(True, color=HAIR, lw=0.8)
        _axtitle(ax, "Caller-support distribution", "how much do callers agree?")
    else:
        _placeholder(ax, "no consensus-evaluable records"); _axtitle(ax, "Caller-support distribution")

    ax = fig.add_subplot(gs[0, 1])
    kset = sorted({D.get_path(r, "evidence.n_callers", 0) or 0 for r in cev})
    cols = [c for c in CLASS_ORDER if any(D.get_path(r, "confidence_class") == c for r in cev)]
    mat = [[sum(1 for r in cev if (D.get_path(r, "evidence.n_callers", 0) or 0) == k
                and D.get_path(r, "confidence_class") == c) for c in cols] for k in kset]
    _heatmap(ax, mat, [f"{k} caller{'s' if k != 1 else ''}" for k in kset],
             [clabel(c) for c in cols], "Caller support → confidence class")

    ax = fig.add_subplot(gs[1, :])
    unk = [r for r in cev if D.get_path(r, "novelty_support") == "UNKNOWN"]
    classes = [c for c in CLASS_ORDER if any(D.get_path(r, "confidence_class") == c for r in unk)]
    vals = [sum(1 for r in unk if D.get_path(r, "confidence_class") == c) for c in classes]
    if unk:
        _hbar(ax, [clabel(c) for c in classes], vals, [CLASS_COLORS[c] for c in classes],
              "Verdicts of short-read-unknown novels under the consensus axis",
              "promotions to MEDIUM / LOW come purely from cross-caller agreement")
    else:
        _placeholder(ax, "no short-read-unknown records in the consensus set")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# PAGE 6 — bio priors (conditional)
# ---------------------------------------------------------------------------
def page6(pdf, ctx, prefix):
    cage = D.numeric_values(ctx.records, "bio_flags.dist_to_CAGE_peak")
    polya = D.numeric_values(ctx.records, "bio_flags.dist_to_polyA_site")
    nmd_vals = [D.get_path(r, "bio_flags.predicted_NMD") for r in ctx.records]
    motif_vals = [D.get_path(r, "bio_flags.polyA_motif_found") for r in ctx.records]
    has_nmd = any(v not in (None, "NA") for v in nmd_vals)
    has_motif = any(v not in (None, "NA") for v in motif_vals)
    if not (cage or polya or has_nmd or has_motif):
        return

    fig = _newpage()
    _page_header(fig, "SQANTI3 biological priors", 6,
                 "Passed through verbatim — these are NOT used in the PanIsoGuard verdict.")
    gs = fig.add_gridspec(2, 2, hspace=0.6, wspace=0.32, left=0.14, right=0.95, top=0.80, bottom=0.13)

    ax = fig.add_subplot(gs[0, 0])
    if cage:
        ax.hist([min(max(v, -500), 500) for v in cage], bins=30, color="#2E9E5B", zorder=3)
        _despine(ax); ax.yaxis.grid(True, color=HAIR, lw=0.8); ax.set_axisbelow(True)
        _axtitle(ax, "Distance to CAGE peak", "clipped to ±500 bp")
    else:
        _placeholder(ax, "no CAGE prior supplied"); _axtitle(ax, "Distance to CAGE peak")

    ax = fig.add_subplot(gs[0, 1])
    if polya:
        ax.hist([min(max(v, -500), 500) for v in polya], bins=30, color="#86BC4C", zorder=3)
        _despine(ax); ax.yaxis.grid(True, color=HAIR, lw=0.8); ax.set_axisbelow(True)
        _axtitle(ax, "Distance to polyA site", "clipped to ±500 bp")
    else:
        _placeholder(ax, "no polyA-site prior supplied"); _axtitle(ax, "Distance to polyA site")

    from collections import Counter
    ax = fig.add_subplot(gs[1, 0])
    if has_nmd:
        c = Counter(v for v in nmd_vals if v not in (None, "NA"))
        _hbar(ax, list(c.keys()), [c[k] for k in c], ["#6D4C41"] * len(c), "Predicted NMD")
    else:
        _placeholder(ax, "no predicted-NMD prior"); _axtitle(ax, "Predicted NMD")

    ax = fig.add_subplot(gs[1, 1])
    if has_motif:
        c = Counter(v for v in motif_vals if v not in (None, "NA"))
        _hbar(ax, list(c.keys()), [c[k] for k in c], ["#C77F33"] * len(c), "PolyA motif found")
    else:
        _placeholder(ax, "no polyA-motif prior"); _axtitle(ax, "PolyA motif found")

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------
# PAGE 7 — provenance appendix
# ---------------------------------------------------------------------------
def page7(pdf, ctx, prefix):
    fig = _newpage()
    _page_header(fig, "Provenance & defensibility appendix", 7,
                 "Exactly how this report was produced, and the checks that make it auditable.")
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.95, left=0.06, right=0.96, top=0.80, bottom=0.10,
                          width_ratios=[0.78, 1.0])
    m = ctx.prov["meta"]; ax_states = ctx.prov["axes"]

    ax = fig.add_subplot(gs[:, 0])
    keys = ["tool_version", "ruleset_version", "sqanti3_version_target", "config",
            "classification", "isoforms", "ref_gtf", "sj_tab", "bam", "caller_support"]
    rows_r = [[k, str(m.get(k, "—"))[:36]] for k in keys if k in m]
    for a, st in ax_states.items():
        rows_r.append([f"axis.{a}", st[:36]])
    _table(ax, ["key", "value"], rows_r, "Run provenance", col_widths=[0.36, 0.64], max_rows=26)

    ax = fig.add_subplot(gs[0, 1])
    from collections import Counter
    rc = Counter()
    for r in ctx.records:
        for t in (D.get_path(r, "rule_trace", []) or []):
            # compact label: drop the verbose "project(...)" wrapper for the legend
            key = t.replace("project(", "").replace(") -> ", " → ")
            rc[key[:30]] += 1
    top = rc.most_common(7)
    if top:
        _hbar(ax, [t for t, _ in top][::-1], [c for _, c in top][::-1],
              ["#7BA6C9"] * len(top), "Most frequent rule-trace lines")

    ax = fig.add_subplot(gs[1, 1])
    fw = D.rescue_firewall(ctx.records)
    in_play = ", ".join(n for n, on in [("short-read", ctx.has_sr), ("BAM", ctx.has_bam),
                        ("variant", ctx.has_variant), ("pangenome", ctx.has_pangenome),
                        ("consensus", ctx.has_consensus)] if on) or "catalog only"
    checks = [
        (("✓ " if not ctx.mismatch else "✗ ") + "class tally reconciles (provenance ↔ JSONL)",
         "#2E9E5B" if not ctx.mismatch else "#D1495B", False),
        (f"✓ firewall held {fw['held_by_firewall_total']} circular-risk rescue(s)", "#2E9E5B", False),
        (f"• consensus axis {'used' if ctx.has_consensus else 'not used'} this run", INK2, False),
        (f"• axes in play: {in_play}", INK2, False),
        ("• SQANTI3 priors consumed, never recomputed", INK2, False),
    ]
    _callout(ax, "Defensibility checklist", checks)

    _footer(fig, ctx.prov, prefix)
    pdf.savefig(fig); plt.close(fig)


def build_report(prefix, out_pdf=None):
    jsonl = prefix + ".attribution.jsonl"
    prov_log = prefix + ".provenance.log"
    if not os.path.exists(jsonl):
        raise FileNotFoundError(f"missing {jsonl} (run `panisoguard adjudicate --out-prefix {prefix}` first)")
    _setup_style()
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
        page6(pdf, ctx, prefix)
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
