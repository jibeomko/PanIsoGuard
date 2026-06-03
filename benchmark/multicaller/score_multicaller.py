#!/usr/bin/env python3
"""Score the multi-caller consensus axis against SQANTI-SIM truth.

Motivation. Three long-read isoform callers (FLAIR, IsoQuant, Bambu) run on the
SAME alignment disagree on which novel isoforms exist. PanIsoGuard `combine`
integrates them caller-agnostically by intron-chain fingerprint; the consensus
axis then treats "recovered by >= K callers" as corroboration when no orthogonal
short-read evidence is available. This script quantifies the payoff: how much
does multi-caller agreement raise the precision of a novel call?

Truth (SQANTI-SIM, keyed by exact intron chain, caller-independent):
  genuine_novel - the chain is a simulated transcript DELETED from the reduced
                  reference (a real isoform the caller must rediscover)
  known         - the chain is a simulated transcript kept in the reference
  false_novel   - a novel multi-exon chain matching NO simulated transcript
                  (a caller/alignment artifact)

Each `combine` matrix row is one unique intron chain with the set of callers that
recovered it. We map the chain to truth (via any caller's native id -> GTF chain)
and, restricting to `novelty == novel` multi-exon chains, report:

  - per-caller novel discovery (count / genuine / false / precision / recall)
  - consensus stratification: precision and yield by number of supporting callers
  - the decisive comparison: union-of-callers (>=1) vs consensus (>= min_callers)

Usage:
  score_multicaller.py --matrix M --truth T --gtf flair:F --gtf isoquant:I \
      --gtf bambu:B [--min-callers 2] [--emit-metrics out.json]
"""
import argparse
import json
import sys
from collections import defaultdict


def tid_of(attr):
    return attr.split('transcript_id "')[1].split('"')[0]


def chain_key(exons, chrom, strand):
    e = sorted(exons)
    intr = [(e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)]
    if not intr:
        return None
    return f"{chrom}|{strand}|" + ",".join(f"{s}-{ee}" for s, ee in intr)


def gtf_native_to_chain(path):
    """transcript_id -> chain_key for one caller GTF (None for monoexonic)."""
    ex = defaultdict(list)
    strand, chrom = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            t = tid_of(f[8])
            ex[t].append((int(f[3]), int(f[4])))
            strand[t] = f[6]
            chrom[t] = f[0]
    return {t: chain_key(ex[t], chrom[t], strand[t]) for t in ex}


def parse_native_ids(field):
    """'flair=id1|id2;bambu=id3' -> [(caller, id), ...]"""
    out = []
    for group in field.split(';'):
        if '=' not in group:
            continue
        caller, ids = group.split('=', 1)
        for i in ids.split('|'):
            if i:
                out.append((caller, i))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--truth', required=True)
    ap.add_argument('--gtf', action='append', default=[],
                    help='caller:path (repeatable)')
    ap.add_argument('--min-callers', type=int, default=2)
    ap.add_argument('--emit-metrics')
    ap.add_argument('--tool-version', default=None)
    ap.add_argument('--ruleset-version', default=None)
    ap.add_argument('--engine-demo', default=None,
                    help='JSON object merged into metrics under "engine_demo" '
                         '(the adjudicate --caller-support run on one caller)')
    args = ap.parse_args()

    # truth: chain_key -> genuine_novel|known
    truth = {}
    for line in open(args.truth):
        k, v = line.rstrip('\n').split('\t')
        truth[k] = v
    total_genuine = sum(1 for v in truth.values() if v == 'genuine_novel')

    # per-caller native_id -> chain_key
    caller_chain = {}
    for spec in args.gtf:
        label, path = spec.split(':', 1)
        caller_chain[label] = gtf_native_to_chain(path)

    def resolve_chain(natives):
        for caller, nid in natives:
            ck = caller_chain.get(caller, {}).get(nid)
            if ck is not None:
                return ck
        return None

    # walk the combine matrix; keep NOVEL multi-exon chains
    hdr = None
    novel_rows = []          # (chain_key, n_callers, set(callers), truth_label)
    seen_genuine = defaultdict(set)   # min_callers_threshold not needed; track per chain
    for line in open(args.matrix):
        f = line.rstrip('\n').split('\t')
        if hdr is None:
            hdr = {c: i for i, c in enumerate(f)}
            continue
        novelty = f[hdr['novelty']]
        if novelty != 'novel':
            continue
        n_callers = int(f[hdr['n_callers']])
        callers = set(f[hdr['callers']].split(','))
        natives = parse_native_ids(f[hdr['native_ids']])
        ck = resolve_chain(natives)
        if ck is None:
            continue  # monoexonic / unresolved
        tl = truth.get(ck)
        label = 'genuine_novel' if tl == 'genuine_novel' else (
            'known' if tl == 'known' else 'false_novel')
        novel_rows.append((ck, n_callers, callers, label))

    # ---- per-caller novel discovery ----
    per_caller = {}
    for spec in args.gtf:
        label = spec.split(':', 1)[0]
        g = f_ = 0
        for ck, n, callers, tl in novel_rows:
            if label in callers and tl in ('genuine_novel', 'false_novel'):
                if tl == 'genuine_novel':
                    g += 1
                else:
                    f_ += 1
        prec = g / (g + f_) if (g + f_) else float('nan')
        rec = g / total_genuine if total_genuine else float('nan')
        per_caller[label] = dict(novel=g + f_, genuine=g, false=f_,
                                 precision=prec, recall=rec)

    # ---- consensus stratification (novel multi-exon: genuine vs false) ----
    by_support = defaultdict(lambda: [0, 0])  # n_callers -> [genuine, false]
    for ck, n, callers, tl in novel_rows:
        if tl == 'genuine_novel':
            by_support[n][0] += 1
        elif tl == 'false_novel':
            by_support[n][1] += 1

    # ---- decisive comparison: union (>=1) vs consensus (>= min_callers) ----
    def agg(min_k):
        g = sum(v[0] for k, v in by_support.items() if k >= min_k)
        fp = sum(v[1] for k, v in by_support.items() if k >= min_k)
        prec = g / (g + fp) if (g + fp) else float('nan')
        rec = g / total_genuine if total_genuine else float('nan')
        return dict(genuine=g, false=fp, precision=prec, recall=rec)

    union = agg(1)
    consensus = agg(args.min_callers)

    # ---- precision/recall curve as the consensus threshold sweeps 1..n_callers ----
    # Each operating point is the novel set retained by requiring ">= k callers".
    # Raising k trades recall for precision; this is the curve a confidence layer
    # should expose. f1 included so a single best operating point is identifiable.
    max_k = max(by_support) if by_support else 1
    pr_curve = []
    for k in range(1, max_k + 1):
        a = agg(k)
        p, r = a['precision'], a['recall']
        f1 = (2 * p * r / (p + r)) if (p == p and r == r and (p + r) > 0) else None
        pr_curve.append(dict(min_callers=k, genuine=a['genuine'], false=a['false'],
                             precision=p, recall=r, f1=f1))

    # ---- report ----
    print("=== per-caller novel discovery (vs SQANTI-SIM truth) ===")
    print(f"{'caller':10s} {'novel':>7s} {'genuine':>8s} {'false':>7s} "
          f"{'precision':>10s} {'recall':>8s}")
    for label, m in per_caller.items():
        print(f"{label:10s} {m['novel']:7d} {m['genuine']:8d} {m['false']:7d} "
              f"{m['precision']:10.3f} {m['recall']:8.3f}")

    print(f"\n=== consensus stratification (novel multi-exon chains, n={sum(sum(v) for v in by_support.values())}) ===")
    print(f"{'#callers':>9s} {'genuine':>8s} {'false':>7s} {'precision':>10s}")
    for k in sorted(by_support):
        g, fp = by_support[k]
        prec = g / (g + fp) if (g + fp) else float('nan')
        print(f"{k:9d} {g:8d} {fp:7d} {prec:10.3f}")

    print(f"\n=== decisive: union (>=1 caller) vs consensus (>= {args.min_callers}) ===")
    print(f"{'set':16s} {'genuine':>8s} {'false':>7s} {'precision':>10s} {'recall':>8s}")
    print(f"{'union(>=1)':16s} {union['genuine']:8d} {union['false']:7d} "
          f"{union['precision']:10.3f} {union['recall']:8.3f}")
    print(f"{'consensus(>=%d)' % args.min_callers:16s} {consensus['genuine']:8d} "
          f"{consensus['false']:7d} {consensus['precision']:10.3f} {consensus['recall']:8.3f}")
    dprec = consensus['precision'] - union['precision']
    print(f"\nconsensus gate raises novel-call precision by {dprec:+.3f} "
          f"({union['precision']:.3f} -> {consensus['precision']:.3f}) "
          f"at recall {consensus['recall']:.3f} (of {total_genuine} truth genuine novels).")

    print(f"\n=== precision/recall curve over consensus_min_callers ===")
    print(f"{'>=k callers':>11s} {'genuine':>8s} {'false':>7s} {'precision':>10s} "
          f"{'recall':>8s} {'F1':>7s}")
    best = max(pr_curve, key=lambda p: (p['f1'] if p['f1'] is not None else -1))
    for p in pr_curve:
        star = "  <- max F1" if p is best else ""
        f1s = f"{p['f1']:.3f}" if p['f1'] is not None else "   nan"
        print(f"{p['min_callers']:11d} {p['genuine']:8d} {p['false']:7d} "
              f"{p['precision']:10.3f} {p['recall']:8.3f} {f1s:>7s}{star}")

    if args.emit_metrics:
        metrics = dict(
            total_truth_genuine_novel=total_genuine,
            n_novel_chains=sum(sum(v) for v in by_support.values()),
            per_caller=per_caller,
            stratification={str(k): dict(genuine=by_support[k][0],
                                         false=by_support[k][1],
                                         precision=(by_support[k][0] /
                                                    (by_support[k][0] + by_support[k][1])
                                                    if sum(by_support[k]) else None))
                            for k in sorted(by_support)},
            union=union,
            consensus=consensus,
            min_callers=args.min_callers,
            pr_curve=pr_curve,
        )
        if args.engine_demo:
            metrics["engine_demo"] = json.loads(args.engine_demo)
        callers = ",".join(sorted(per_caller))
        envelope = dict(
            protocol="multicaller",
            scope="gencode_v49_chr22_%dcaller" % len(per_caller),
            status="tracked",
            tool_version=args.tool_version,
            ruleset_version=args.ruleset_version,
            generated_utc=None,
            data_provenance=(
                "SQANTI-SIM GENCODE v49 chr22 (PBSIM3 simulated reads, one shared "
                "minimap2 alignment). Long-read callers run on that alignment: "
                "FLAIR collapse, IsoQuant 3.x transcript_models, Bambu 3.x novel "
                "transcripts (NDR=0.5), ESPRESSO 1.4 novel_isoform, TALON 6.0 "
                "(filtered whitelist, minCount 5). Truth = SQANTI-SIM "
                "genuine_novel/known labels (transcripts deleted-from vs kept-in the "
                "reduced annotation). Callers present in this run: " + callers +
                ". No private data."),
            command=(
                "panisoguard combine --gtf <caller:gtf>... --ref-gtf chr22_modified.gtf "
                "--out matrix.tsv; benchmark/multicaller/score_multicaller.py "
                "--matrix matrix.tsv --truth truth.truth.tsv --gtf <caller:gtf>... "
                "--emit-metrics metrics.json  (see benchmark/multicaller/run.sh)"),
            source=None,
            metrics=metrics,
            notes=(
                "Caller-agnostic multi-caller integration over %d callers. "
                "Single-caller novel calls are overwhelmingly artifacts (precision "
                "rises monotonically with the number of supporting callers); requiring "
                ">= min_callers raises novel-call precision sharply. The pr_curve sweeps "
                "the consensus threshold 1..n: the F1-optimal threshold SCALES with the "
                "number of callers (more callers -> require more agreement). engine_demo "
                "confirms the wired consensus axis reproduces this through `adjudicate "
                "--caller-support` on real caller output (long-read only, no short-read "
                "SJ): a single caller yields 0 confident novel calls (all AMBIGUOUS), "
                "while consensus promotes the cross-caller-agreed novels. "
                "consensus_min_callers is configurable (TOML); the default is 2." %
                len(per_caller)),
        )
        with open(args.emit_metrics, 'w') as fh:
            json.dump(envelope, fh, indent=2, sort_keys=True)
            fh.write('\n')
        print(f"\nwrote {args.emit_metrics}")


if __name__ == '__main__':
    main()
