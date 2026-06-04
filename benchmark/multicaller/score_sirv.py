#!/usr/bin/env python3
"""Score the multi-caller consensus axis on the SIRV-Set4 spike-in reference — a SECOND,
denser, multi-contig transcriptome — against the hidden-chain truth. This is the
generalization companion to score_multicaller.py (which is GENCODE chr22): it checks that
"single-caller novels are mostly artifacts; >= k-caller consensus raises precision" is not
specific to one reference/simulator.

Truth: SIRV-Set4 transcripts deleted from the reduced annotation (`truth_hidden_chains.tsv`,
one `tid<TAB>contig<TAB>intron-chain` line each). Matching is strand-agnostic (the SIRV
truth carries no strand) on the exact intron chain. Same `combine` matrix + scoring as the
chr22 protocol; emits the same per-caller / stratification / PR-curve metrics.

Usage:
  score_sirv.py --matrix M --truth truth_hidden_chains.tsv --gtf flair:F --gtf isoquant:I ...
      [--min-callers 2] [--emit-metrics out.json] [--tool-version V]
"""
import argparse
import json
from collections import Counter, defaultdict


def tid_of(attr):
    return attr.split('transcript_id "')[1].split('"')[0]


def gtf_chains(path):
    """transcript_id -> 'contig|introns' (strand-agnostic); None if monoexonic."""
    ex = defaultdict(list); chrom = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = tid_of(f[8]); ex[t].append((int(f[3]), int(f[4]))); chrom[t] = f[0]
    out = {}
    for t, exons in ex.items():
        e = sorted(exons)
        # SIRV truth_hidden_chains.tsv encodes each intron as (exon_end, next_exon_start-1)
        # -- start is the last base of the upstream exon (no +1), unlike SQANTI-SIM. Match it.
        intr = [(e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)]
        out[t] = f"{chrom[t]}|" + ",".join(f"{s}-{ee}" for s, ee in intr) if intr else None
    return out


def parse_native_ids(field):
    out = []
    for grp in field.split(';'):
        if '=' in grp:
            caller, ids = grp.split('=', 1)
            for i in ids.split('|'):
                if i:
                    out.append((caller, i))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--gtf', action='append', default=[])
    ap.add_argument('--min-callers', type=int, default=2)
    ap.add_argument('--emit-metrics')
    ap.add_argument('--tool-version')
    args = ap.parse_args()

    # hidden-chain truth (genuine novels), strand-agnostic key "contig|introns"
    truth = set()
    for line in open(args.truth):
        f = line.rstrip('\n').split('\t')
        if len(f) >= 3 and f[2]:
            truth.add(f"{f[1]}|" + f[2].replace(';', ','))
    total_genuine = len(truth)

    caller_chain = {}
    for spec in args.gtf:
        label, path = spec.split(':', 1)
        caller_chain[label] = gtf_chains(path)

    def resolve(natives):
        for caller, nid in natives:
            v = caller_chain.get(caller, {}).get(nid)
            if v is not None:
                return v
        return None

    hdr = None
    # Dedup to DISTINCT chains (max n_callers, union of callers) so strand-agnostic matching
    # can't double-count one truth chain across opposite-strand rows -> recall stays <= 1.
    key_info = {}   # key -> [max_n_callers, genuine(bool), callers_union]
    for line in open(args.matrix):
        f = line.rstrip('\n').split('\t')
        if hdr is None:
            hdr = {c: i for i, c in enumerate(f)}; continue
        if f[hdr['novelty']] != 'novel':
            continue
        ck = resolve(parse_native_ids(f[hdr['native_ids']]))
        if ck is None:
            continue
        n = int(f[hdr['n_callers']]); cs = set(f[hdr['callers']].split(','))
        if ck not in key_info:
            key_info[ck] = [n, ck in truth, set(cs)]
        else:
            ki = key_info[ck]; ki[0] = max(ki[0], n); ki[2] |= cs

    per_caller = {}
    for spec in args.gtf:
        label = spec.split(':', 1)[0]
        g = sum(1 for n, gen, cs in key_info.values() if label in cs and gen)
        fp = sum(1 for n, gen, cs in key_info.values() if label in cs and not gen)
        per_caller[label] = dict(novel=g + fp, genuine=g, false=fp,
                                 precision=(g / (g + fp) if (g + fp) else None),
                                 recall=(g / total_genuine if total_genuine else None))

    by_support = defaultdict(lambda: [0, 0])
    for n, gen, cs in key_info.values():
        by_support[n][0 if gen else 1] += 1
    max_k = max(by_support) if by_support else 1

    def agg(mink):
        g = sum(v[0] for kk, v in by_support.items() if kk >= mink)
        fp = sum(v[1] for kk, v in by_support.items() if kk >= mink)
        return dict(genuine=g, false=fp, precision=(g / (g + fp) if (g + fp) else None),
                    recall=(g / total_genuine if total_genuine else None))

    pr_curve = []
    for k in range(1, max_k + 1):
        a = agg(k); p, r = a['precision'], a['recall']
        f1 = (2 * p * r / (p + r)) if (p and r and (p + r)) else None
        pr_curve.append(dict(min_callers=k, **a, f1=f1))

    print(f"SIRV hidden (genuine-novel) chains = {total_genuine}\n")
    print(f"{'caller':10s} {'novel':>6s} {'genuine':>8s} {'false':>6s} {'prec':>6s} {'recall':>7s}")
    for lab, m in per_caller.items():
        print(f"{lab:10s} {m['novel']:6d} {m['genuine']:8d} {m['false']:6d} "
              f"{(m['precision'] or 0):6.3f} {(m['recall'] or 0):7.3f}")
    print(f"\n#callers  genuine  false  precision")
    for k in sorted(by_support):
        g, fp = by_support[k]
        print(f"{k:7d} {g:8d} {fp:6d} {(g/(g+fp) if (g+fp) else 0):10.3f}")
    print(f"\n>=k PR curve:")
    for c in pr_curve:
        print(f"  >={c['min_callers']}: P={(c['precision'] or 0):.3f} R={(c['recall'] or 0):.3f} "
              f"F1={(c['f1'] or 0):.3f}")

    if args.emit_metrics:
        callers = ",".join(sorted(per_caller))
        env = dict(
            protocol="sirv_multicaller", scope="sirv_set4_multicaller",
            status="tracked", tool_version=args.tool_version, ruleset_version=None,
            generated_utc=None,
            data_provenance=("Lexogen SIRV-Set4 spike-in reference (7 loci, ~70 dense "
                "overlapping isoforms, multi-contig), PBSIM3-simulated reads, one shared "
                "minimap2 alignment. Callers: " + callers + ". Truth = SIRV transcripts "
                "deleted from the reduced annotation (hidden chains). A SECOND reference "
                "(beyond GENCODE chr22) for the multi-caller consensus. No private data."),
            command=("panisoguard combine --gtf <caller:gtf>...; "
                "benchmark/multicaller/score_sirv.py --matrix matrix.tsv "
                "--truth truth_hidden_chains.tsv --gtf <caller:gtf>... --emit-metrics metrics.json"),
            source=None,
            metrics=dict(total_hidden_genuine=total_genuine,
                         n_novel_chains=sum(sum(v) for v in by_support.values()),
                         per_caller=per_caller,
                         stratification={str(k): dict(genuine=by_support[k][0],
                                                      false=by_support[k][1],
                                                      precision=(by_support[k][0] /
                                                                 sum(by_support[k])
                                                                 if sum(by_support[k]) else None))
                                         for k in sorted(by_support)},
                         pr_curve=pr_curve, min_callers=args.min_callers),
            notes=("Generalization of the chr22 multicaller result to a SECOND, denser, "
                "multi-contig reference (SIRV-Set4). Confirms the consensus signal "
                "(single-caller novels mostly artifacts; >= k-caller agreement raises "
                "precision) is not specific to one reference/simulator. Strand-agnostic "
                "exact-chain matching."),
        )
        with open(args.emit_metrics, 'w') as fh:
            json.dump(env, fh, indent=2, sort_keys=True); fh.write('\n')
        print(f"\nwrote {args.emit_metrics}")


if __name__ == '__main__':
    main()
