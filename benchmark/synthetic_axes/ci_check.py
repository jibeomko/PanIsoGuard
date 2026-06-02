#!/usr/bin/env python3
"""CTest integration driver for the synthetic-axes protocol.

Validates the evidence axes end-to-end through `panisoguard adjudicate` on a fully
synthetic contig (gen.py), against known per-category truth. Tiered by external-tool
availability so it stays green in a minimal CI (htslib only, no samtools):

  * short-read axis  (cat A genuine -> HIGH_CONF_NOVEL, cat B noncanon -> ARTIFACT)
        -> always run (pure text).
  * variant axis     (cat D reference-bias -> PAN_REF_RESCUED_FALSE_NOVEL)
        -> always run; htslib faidx auto-builds the .fai, so no samtools needed.
  * BAM mapping axis (cat C mapping-artifact -> ARTIFACT/mapping)
        -> run ONLY if `samtools` is on PATH (needs a sorted+indexed BAM); otherwise
           skipped with a notice (not a failure).

Exit code is non-zero on any assertion failure. Usage:  ci_check.py /path/to/panisoguard
"""
import json
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEN = HERE / "gen.py"
NPC = 10  # isoforms per category (40 total) -- fast, still unambiguous


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.stderr.write(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n")
        sys.stderr.write(r.stdout + "\n" + r.stderr + "\n")
        sys.exit(1)
    return r


def load_classes(adj_tsv: Path):
    cls, header = {}, None
    for ln in open(adj_tsv):
        f = ln.rstrip("\n").split("\t")
        if header is None:
            header = {n: k for k, n in enumerate(f)}
            continue
        cls[f[header["isoform_id"]]] = f[header["confidence_class"]]
    return cls


def by_category(cat_of, cls):
    table = defaultdict(Counter)
    for iso, c in cat_of.items():
        table[c][cls.get(iso, "MISSING")] += 1
    return table


def expect(table, category, want, failures, label):
    """Assert every isoform in `category` got class `want`."""
    got = table[category]
    if set(got) != {want}:
        failures.append(f"{label}: category {category} expected all={want}, got {dict(got)}")


def main() -> int:
    args = [a for a in sys.argv[1:]]
    emit = None
    if "--emit-metrics" in args:
        k = args.index("--emit-metrics")
        emit = args[k + 1]
        del args[k:k + 2]
    if not args:
        sys.stderr.write("usage: ci_check.py <panisoguard-binary> [--emit-metrics PATH]\n")
        return 2
    pig = args[0]
    samtools = shutil.which("samtools")

    with tempfile.TemporaryDirectory(prefix="pig_sa_") as tmp:
        work = Path(tmp)
        run([sys.executable, str(GEN), str(work), str(NPC)])

        cat_of = {}
        for ln in open(work / "truth.tsv"):
            if ln.startswith("isoform"):
                continue
            iso, category, _truth = ln.rstrip("\n").split("\t")
            cat_of[iso] = category

        base = ["--classification", str(work / "classification.tsv"),
                "--isoforms-gtf", str(work / "caller.gtf"),
                "--ref-gtf", str(work / "catalog.gtf"),
                "--sj-tab", str(work / "real.SJ.tab")]
        failures = []

        # --- short-read axis (text only) -------------------------------------
        run([pig, "adjudicate", *base, "--out-prefix", str(work / "base")])
        t_base = by_category(cat_of, load_classes(work / "base.adjudicated.tsv"))
        expect(t_base, "A", "HIGH_CONF_NOVEL", failures, "base")
        expect(t_base, "B", "ARTIFACT", failures, "base")

        # --- variant axis (faidx auto-built; no samtools) --------------------
        run([pig, "adjudicate", *base,
             "--reference", str(work / "ref.fa"),
             "--reference-haplotype", str(work / "hap.fa"),
             "--haplotype-provenance", "external",
             "--out-prefix", str(work / "var")])
        t_var = by_category(cat_of, load_classes(work / "var.adjudicated.tsv"))
        expect(t_var, "D", "PAN_REF_RESCUED_FALSE_NOVEL", failures, "variant")
        expect(t_var, "A", "HIGH_CONF_NOVEL", failures, "variant")  # rescue must not perturb A
        expect(t_var, "B", "ARTIFACT", failures, "variant")        # ...or B

        # --- BAM mapping axis (needs samtools to build sorted+indexed BAM) ----
        if samtools:
            run([samtools, "sort", "-O", "bam", "-o", str(work / "syn.bam"), str(work / "syn.sam")])
            run([samtools, "index", str(work / "syn.bam")])
            run([pig, "adjudicate", *base,
                 "--bam", str(work / "syn.bam"),
                 "--reference", str(work / "ref.fa"),
                 "--out-prefix", str(work / "bam")])
            t_bam = by_category(cat_of, load_classes(work / "bam.adjudicated.tsv"))
            expect(t_bam, "C", "ARTIFACT", failures, "bam")  # low-MAPQ spanning -> mapping artifact
            bam_note = "C->ARTIFACT(mapping) checked"
        else:
            bam_note = "SKIPPED (samtools not found) -- BAM/mapping axis (category C) not exercised"

        print(f"[synthetic_axes] {NPC}/category  short-read A/B + variant D checked; {bam_note}")
        if failures:
            for f in failures:
                sys.stderr.write(f"[synthetic_axes] FAIL: {f}\n")
            return 1

        if emit:
            def summ(table):  # category -> dominant class (categories are uniform here)
                return {c: (table[c].most_common(1)[0][0] if table[c] else "MISSING")
                        for c in ("A", "B", "C", "D")}
            metrics = {
                "n_per_category": NPC,
                "samtools_used": bool(samtools),
                "axes_checked": ["short_read", "variant"] + (["mapping"] if samtools else []),
                "base": summ(t_base),
                "variant": summ(t_var),
                "bam": summ(t_bam) if samtools else None,
            }
            Path(emit).parent.mkdir(parents=True, exist_ok=True)
            with open(emit, "w") as fh:
                json.dump(metrics, fh, indent=2, sort_keys=True)
                fh.write("\n")
            print(f"[synthetic_axes] wrote metrics -> {emit}")
        print("[synthetic_axes] PASS")
        return 0


if __name__ == "__main__":
    sys.exit(main())
