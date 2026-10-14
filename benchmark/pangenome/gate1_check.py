#!/usr/bin/env python3
"""GATE-1 sensitivity + circularity-firewall check for the pangenome rescue axis.

Specificity (that genuine novel junctions are NOT falsely rescued) is shown on real
caller output in benchmark/results/pangenome (run a real chr22 FLAIR set through the
HPRC-derived junctions -> 0 rescues). This script is the complementary constructed
check on REAL graph-derived junction coordinates:

  POS       novel intron == a real HPRC graph junction  -> REFERENCE_BIAS
            (only when --pangenome-provenance is population/external)
  CIRCULAR  same isoform, provenance=unknown (default)   -> held UNRESOLVED (firewall),
            circularity_flag=true  (NOT promoted)
  NEARMISS  novel intron 5 bp off a graph junction       -> held (exact match required)
  CONTROL   novel intron absent from the graph           -> held (no support)
  PATHOK    two graph junctions that share a haplotype path (6th column) -> REFERENCE_BIAS
  PATHBAD   two graph junctions on disjoint paths        -> held UNRESOLVED, reason path_inconsistent
            (the two path cases need the haplotype_ids column; skipped when it is absent)

Usage: gate1_check.py <panisoguard-binary> <graph_junctions.tsv>
"""
import subprocess
import sys
import tempfile
from pathlib import Path

CHROM = "chr22"


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"FAILED: {' '.join(map(str, cmd))}\n{r.stdout}\n{r.stderr}\n")
        sys.exit(1)
    return r


def load_junctions(path, n=6):
    """Pick n graph junctions (strand +, intron-sized) to build POS isoforms from."""
    return [(s, e) for s, e, _ in load_junctions_with_paths(path, n)]


def load_junctions_with_paths(path, n=None):
    """(start, end, set(path ids)) for the usable + strand junctions; ids empty when the column is absent."""
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or f[0] != CHROM or f[3] != "+":
            continue
        s, e = int(f[1]), int(f[2])
        if 60 <= (e - s + 1) <= 2000 and s > 200:
            ids = set(f[5].split(",")) if len(f) >= 6 and f[5] not in ("", ".", "NA") else set()
            out.append((s, e, ids))
            if n is not None and len(out) >= n:
                break
    return out


def path_pairs(junctions):
    """One (shared-path) and one (disjoint-path) pair of non-overlapping junctions, if the file has path ids."""
    ok = bad = None
    withids = [j for j in junctions if j[2]]
    for i, a in enumerate(withids):
        for b in withids[i + 1:]:
            if b[0] <= a[1] + 200:  # need room for an exon between them
                continue
            if ok is None and a[2] & b[2]:
                ok = (a, b)
            if bad is None and not (a[2] & b[2]):
                bad = (a, b)
            if ok and bad:
                return ok, bad
    return ok, bad


def verdicts(adj_tsv):
    h, v = None, {}
    for line in open(adj_tsv):
        f = line.rstrip("\n").split("\t")
        if h is None:
            h = {n: i for i, n in enumerate(f)}
            continue
        v[f[h["isoform_id"]]] = (f[h["confidence_class"]], f[h["circularity_flag"]], f[h["verdict_reason"]])
    return v


def main():
    if len(sys.argv) < 3:
        sys.stderr.write("usage: gate1_check.py <panisoguard-binary> <graph_junctions.tsv>\n")
        return 2
    pig, pj = sys.argv[1], sys.argv[2]
    jx = load_junctions(pj)
    if len(jx) < 3:
        sys.stderr.write(f"need >=3 usable junctions in {pj}, got {len(jx)}\n")
        return 2

    with tempfile.TemporaryDirectory(prefix="pig_gate1_") as t:
        d = Path(t)
        HDR = ("isoform\tchrom\tstrand\tstructural_category\tassociated_gene\tassociated_transcript\t"
               "subcategory\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\tn_indels_junc\t"
               "dist_to_CAGE_peak\tdist_to_polyA_site\tfilter_result\n")
        gtf = open(d / "caller.gtf", "w")
        cls = open(d / "cls.tsv", "w"); cls.write(HDR)
        cat = open(d / "catalog.gtf", "w")
        # decoy catalog transcript far from the test introns (keeps test introns novel)
        for (a, b) in [(1000, 1100), (1300, 1400)]:
            cat.write(f'{CHROM}\tdecoy\texon\t{a}\t{b}\t.\t+\t.\tgene_id "D"; transcript_id "DECOY";\n')
        cat.close()

        def emit(iso, *introns):
            # transcript with the given introns (1-based inclusive), 100-bp exons around them
            bounds = [introns[0][0] - 100]
            for (s_, e_) in introns:
                bounds += [s_ - 1, e_ + 1]
            bounds.append(introns[-1][1] + 100)
            for a, b in zip(bounds[0::2], bounds[1::2]):
                gtf.write(f'{CHROM}\tsyn\texon\t{a}\t{b}\t.\t+\t.\tgene_id "{iso}"; transcript_id "{iso}";\n')
            cls.write(f"{iso}\t{CHROM}\t+\tnovel_not_in_catalog\t{iso}\tnovel\tmulti-exon\tFALSE\t"
                      f"non_canonical\t10.0\t0\tNA\tNA\tIsoform\n")

        truth = {}
        for i, (s, e) in enumerate(jx):
            emit(f"POS_{i}", (s, e));          truth[f"POS_{i}"] = "POS"           # exact graph junction
            emit(f"NEARMISS_{i}", (s, e + 5)); truth[f"NEARMISS_{i}"] = "NEARMISS"  # 5 bp off
            emit(f"CONTROL_{i}", (s, s + 207)); truth[f"CONTROL_{i}"] = "CONTROL"   # not in graph
        pair_ok, pair_bad = path_pairs(load_junctions_with_paths(pj))
        if pair_ok:
            emit("PATHOK", pair_ok[0][:2], pair_ok[1][:2]); truth["PATHOK"] = "PATHOK"
        if pair_bad:
            emit("PATHBAD", pair_bad[0][:2], pair_bad[1][:2]); truth["PATHBAD"] = "PATHBAD"
        gtf.close(); cls.close()

        common = ["--classification", str(d / "cls.tsv"), "--isoforms-gtf", str(d / "caller.gtf"),
                  "--ref-gtf", str(d / "catalog.gtf"), "--pangenome-junctions", pj]

        # Run 1: population provenance (promotes)
        run([pig, "adjudicate", *common, "--pangenome-provenance", "population",
             "--out-prefix", str(d / "pop")])
        vp = verdicts(d / "pop.adjudicated.tsv")
        # Run 2: unknown provenance (circularity firewall holds)
        run([pig, "adjudicate", *common, "--pangenome-provenance", "unknown",
             "--out-prefix", str(d / "unk")])
        vu = verdicts(d / "unk.adjudicated.tsv")

        fails = []
        npos = sum(1 for k in truth if truth[k] == "POS")
        rescued_pop = held_unk = circ_unk = 0
        for iso, cat_ in truth.items():
            cls_pop, _, reason_pop = vp.get(iso, ("MISSING", "", ""))
            cls_unk, circ_unk_flag, _ = vu.get(iso, ("MISSING", "", ""))
            if cat_ == "PATHOK":
                if cls_pop != "REFERENCE_BIAS":
                    fails.append(f"{iso}: two junctions on a shared path -> {cls_pop} (expected REFERENCE_BIAS)")
            elif cat_ == "PATHBAD":
                if cls_pop == "REFERENCE_BIAS" or reason_pop != "path_inconsistent":
                    fails.append(f"{iso}: two junctions on disjoint paths -> {cls_pop}/{reason_pop} "
                                 "(expected UNRESOLVED/path_inconsistent)")
            elif cat_ == "POS":
                if cls_pop != "REFERENCE_BIAS":
                    fails.append(f"{iso}: population -> {cls_pop} (expected REFERENCE_BIAS)")
                else:
                    rescued_pop += 1
                if cls_unk == "REFERENCE_BIAS":
                    fails.append(f"{iso}: unknown provenance was PROMOTED (firewall breach!)")
                else:
                    held_unk += 1
                    if circ_unk_flag == "true":
                        circ_unk += 1
            else:  # NEARMISS / CONTROL must never be rescued
                if cls_pop == "REFERENCE_BIAS":
                    fails.append(f"{iso} ({cat_}): falsely rescued (expected held)")

        print(f"[gate1] POS={npos}  rescued(population)={rescued_pop}/{npos}  "
              f"held+flagged(unknown firewall)={circ_unk}/{npos}  "
              f"NEARMISS/CONTROL false rescues={sum(1 for f in fails if 'falsely rescued' in f)}  "
              f"path-consistency: shared-path pair {'rescued' if pair_ok and 'PATHOK' not in ' '.join(fails) else ('n/a' if not pair_ok else 'FAIL')}, "
              f"disjoint-path pair {'held' if pair_bad and 'PATHBAD' not in ' '.join(fails) else ('n/a' if not pair_bad else 'FAIL')}")
        if fails:
            for f in fails:
                sys.stderr.write(f"[gate1] FAIL: {f}\n")
            return 1
        print("[gate1] PASS (rescue fires on real graph junctions; firewall holds circular-risk; "
              "near-miss/control not rescued)")
        return 0


if __name__ == "__main__":
    sys.exit(main())
