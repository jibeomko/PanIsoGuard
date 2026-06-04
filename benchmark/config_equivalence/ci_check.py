#!/usr/bin/env python3
"""Regression guard: the shipped config/rules.default.toml must encode EXACTLY the
built-in engine defaults.

A `consensus_min_callers` sentinel in the default TOML once silently DISABLED the
consensus axis when a user passed `--config config/rules.default.toml` (the TOML said
999 while the built-in was 2). This test would have caught it: it adjudicates the same
fixture WITH and WITHOUT the default config and asserts the verdicts are byte-identical,
AND that the consensus axis actually fires (a >=2-caller novel reaches MEDIUM_CONF_NOVEL,
not AMBIGUOUS). Any future default that drifts between the TOML and the built-ins fails.

Usage: ci_check.py <path-to-panisoguard>   (run by the `integration_config_equivalence` CTest)
"""
import os
import pathlib
import subprocess
import sys
import tempfile

PIG = sys.argv[1]
REPO = pathlib.Path(__file__).resolve().parent.parent.parent
DEFAULT_CONFIG = REPO / "config" / "rules.default.toml"

CLS = (
    "isoform\tchrom\tstrand\tstructural_category\tassociated_gene\tassociated_transcript\t"
    "subcategory\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\tn_indels_junc\t"
    "dist_to_CAGE_peak\tdist_to_polyA_site\tfilter_result\n"
    "isoX\tchr1\t+\tnovel_not_in_catalog\tG1\tnovel\tat_least_one_novel_splicesite\t"
    "FALSE\tcanonical\t10.0\t0\tNA\tNA\tIsoform\n")
# isoX: 2 exons -> one novel intron [201,399]. (No --sj-tab => novelty-support UNKNOWN,
# which is exactly when the consensus axis is allowed to corroborate.)
ISO_GTF = (
    'chr1\tt\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "isoX";\n'
    'chr1\tt\texon\t400\t600\t.\t+\t.\tgene_id "G1"; transcript_id "isoX";\n')
# reference catalog WITHOUT isoX's intron (so isoX is genuinely novel)
REF_GTF = (
    'chr1\tt\texon\t100\t200\t.\t+\t.\tgene_id "REF"; transcript_id "REF";\n'
    'chr1\tt\texon\t250\t350\t.\t+\t.\tgene_id "REF"; transcript_id "REF";\n')
# combine-style caller-support matrix: isoX recovered by 3 callers
SUPPORT = (
    "pig_id\tchrom\tstrand\tn_introns\tn_callers\tn_isoforms\tcallers\tnovelty\tnative_ids\n"
    "PIG.1\tchr1\t+\t1\t3\t3\tflair,isoquant,bambu\tnovel\tflair=isoX;isoquant=x2;bambu=x3\n")


def adjudicate(workdir, prefix, config=None):
    cmd = [PIG, "adjudicate",
           "--classification", str(workdir / "cls.tsv"),
           "--isoforms-gtf", str(workdir / "iso.gtf"),
           "--ref-gtf", str(workdir / "ref.gtf"),
           "--caller-support", str(workdir / "support.tsv"),
           "--out-prefix", str(workdir / prefix)]
    if config:
        cmd += ["--config", str(config)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stderr.write(f"adjudicate failed ({prefix}):\n{r.stderr}\n")
        sys.exit(1)
    return (workdir / f"{prefix}.adjudicated.tsv").read_text()


def verdict_of(tsv, iso):
    for line in tsv.splitlines()[1:]:
        f = line.split("\t")
        if f and f[0] == iso:
            return f[6]  # confidence_class column
    return None


def main():
    if not DEFAULT_CONFIG.exists():
        sys.stderr.write(f"missing {DEFAULT_CONFIG}\n"); return 1
    with tempfile.TemporaryDirectory() as td:
        w = pathlib.Path(td)
        (w / "cls.tsv").write_text(CLS)
        (w / "iso.gtf").write_text(ISO_GTF)
        (w / "ref.gtf").write_text(REF_GTF)
        (w / "support.tsv").write_text(SUPPORT)

        builtin = adjudicate(w, "builtin")            # built-in defaults
        shipped = adjudicate(w, "shipped", DEFAULT_CONFIG)  # config/rules.default.toml

        # (1) shipped default config must reproduce the built-in defaults exactly
        if builtin != shipped:
            sys.stderr.write("FAIL: config/rules.default.toml does NOT match the built-in "
                             "defaults (adjudicated output differs).\n--- builtin ---\n"
                             f"{builtin}\n--- shipped ---\n{shipped}\n")
            return 1
        # (2) the consensus axis must actually fire (guards the silent-disable bug class)
        v = verdict_of(builtin, "isoX")
        if v != "MEDIUM_CONF_NOVEL":
            sys.stderr.write(f"FAIL: consensus axis did not fire — isoX = {v}, expected "
                             "MEDIUM_CONF_NOVEL (>=2 callers, no short-read support).\n")
            return 1
        print(f"OK: shipped default config == built-in defaults; consensus fires (isoX -> {v})")
        return 0


if __name__ == "__main__":
    sys.exit(main())
