#!/usr/bin/env python3
"""CTest integration driver for the controlled-truth (incomplete-reference) protocol.

Self-contained: needs only python3 + the panisoguard binary (no samtools, no real
genome, no network). It synthesizes a small disjoint-window reference GTF, runs the
protocol's own `gen.py` to derive a labelled truth set, adjudicates it through the
freshly-built binary, and ASSERTS that the short-read corroboration axis is correct
on decisive calls (precision(genuine)=1.0, specificity(false)=1.0). Exit code is
non-zero on any failure, so `ctest` treats a regression as a test failure.

Usage (CTest passes the binary path):  ci_check.py /path/to/panisoguard
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
GEN = HERE / "gen.py"


def make_reference_gtf(path: Path, n_tx: int = 60) -> None:
    """Write a synthetic reference GTF: n_tx multi-exon transcripts in disjoint
    2 kb windows (no shared junctions), 4 exons each so gen.py can both hide some
    (genuine) and fabricate a shifted-exon junction in others (false)."""
    chrom = "chrSYN"
    with open(path, "w") as g:
        for t in range(n_tx):
            b = t * 2000 + 1000  # 1-based window start; windows never overlap
            exons = [(b, b + 79), (b + 200, b + 279), (b + 400, b + 479), (b + 600, b + 679)]
            tid = f"TX{t:02d}"  # zero-padded -> string sort == numeric sort
            for (s, e) in exons:
                g.write(f'{chrom}\tsyn\texon\t{s}\t{e}\t.\t+\t.\t'
                        f'gene_id "G{t:02d}"; transcript_id "{tid}";\n')


def run(cmd, **kw):
    r = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.stderr.write(f"command failed ({r.returncode}): {' '.join(map(str, cmd))}\n")
        sys.stderr.write(r.stdout + "\n" + r.stderr + "\n")
        sys.exit(1)
    return r


def pred(cls: str) -> str:
    if cls in ("HIGH_CONF_NOVEL", "MEDIUM_CONF_NOVEL"):
        return "genuine"
    if cls in ("LOW_CONF_PARTIAL", "ARTIFACT"):
        return "false"
    return "abstain"  # AMBIGUOUS / HIGH_CONF_KNOWN


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

    with tempfile.TemporaryDirectory(prefix="pig_ct_") as tmp:
        work = Path(tmp)
        ref_gtf = work / "synthetic_ref.gtf"
        make_reference_gtf(ref_gtf)

        run([sys.executable, str(GEN), str(ref_gtf), str(work)])
        run([pig, "adjudicate",
             "--classification", str(work / "classification.tsv"),
             "--isoforms-gtf", str(work / "caller.gtf"),
             "--ref-gtf", str(work / "reduced_catalog.gtf"),
             "--sj-tab", str(work / "real.SJ.tab"),
             "--out-prefix", str(work / "e1")])

        truth = {}
        for ln in open(work / "truth.tsv"):
            if ln.startswith("isoform"):
                continue
            i, t = ln.rstrip("\n").split("\t")
            truth[i] = t

        cls, header = {}, None
        for ln in open(work / "e1.adjudicated.tsv"):
            f = ln.rstrip("\n").split("\t")
            if header is None:
                header = {n: k for k, n in enumerate(f)}
                continue
            cls[f[header["isoform_id"]]] = f[header["confidence_class"]]

        TP = FP = FN = TN = ab_g = 0
        # A truth isoform absent from the adjudicated output is a COVERAGE failure
        # (the tool emitted no verdict), not a benign abstention -- flag it explicitly.
        missing = [i for i in truth if i not in cls]

        for i, t in truth.items():
            p = pred(cls.get(i, "MISSING"))
            if t == "genuine":
                TP += p == "genuine"; FN += p == "false"; ab_g += p == "abstain"
            else:
                FP += p == "genuine"; TN += p == "false"

        prec = TP / (TP + FP) if (TP + FP) else 0.0
        spec = TN / (TN + FP) if (TN + FP) else 0.0
        print(f"[controlled_truth] genuine={TP+FN+ab_g} false={FP+TN}  "
              f"TP={TP} FP={FP} TN={TN} FN={FN} abstain(genuine)={ab_g}  "
              f"precision={prec:.3f} specificity={spec:.3f}")

        ok = True
        for name, cond in [
            ("every truth isoform adjudicated (no missing output)", not missing),
            ("at least one decisive genuine call", TP > 0),
            ("at least one decisive false call", TN > 0),
            ("precision(genuine) == 1.0", prec == 1.0),
            ("specificity(false) == 1.0", spec == 1.0),
            ("no genuine misclassified as false", FN == 0),
        ]:
            if not cond:
                sys.stderr.write(f"[controlled_truth] FAIL: {name}\n")
                ok = False
        if not ok:
            return 1
        if emit:
            metrics = {
                "TP": TP, "FP": FP, "TN": TN, "FN": FN,
                "n_genuine": TP + FN + ab_g, "n_false": FP + TN,
                "abstain_genuine": ab_g,
                "precision_genuine": round(prec, 4),
                "specificity_false": round(spec, 4),
            }
            Path(emit).parent.mkdir(parents=True, exist_ok=True)
            with open(emit, "w") as fh:
                json.dump(metrics, fh, indent=2, sort_keys=True)
                fh.write("\n")
            print(f"[controlled_truth] wrote metrics -> {emit}")
        print("[controlled_truth] PASS")
        return 0


if __name__ == "__main__":
    sys.exit(main())
