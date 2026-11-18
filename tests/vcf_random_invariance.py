#!/usr/bin/env python3
"""Randomized check of the VCF haplotype route (--haplotype-vcf) on a locus of one or two novel
junctions. Each random file holds 1-5 records at or next to the splice dinucleotides: single
bases, MNPs, multi-allelic sites, insertions, deletions, a deletion of the whole intron, symbolic
and `*` alleles; phased in sets, unphased, partly or fully missing, haploid, without GT, with an
allele index the record does not define, with a REF that is not the reference; rows repeated.

The call must not depend on how the same information is written:
  order      rows that share a POS, in any order
  relabel    every genotype with its alleles swapped (haplotype 1 <-> 2)
  group      the alleles swapped inside one phase set, or in one unphased record
  duplicate  one row written twice
  distractor a record added where no motif is read
  rewrite    an indel moved through a repeat, or a record written with flanking reference bases
and it must not claim more than the file supports:
  fasta      substitutions, fully phased: two haplotype FASTAs give the same call
  motif      every haplotype motif the tool reports equals the one read off the haplotype sequence
  sound      REFERENCE_BIAS only when, under every orientation of the groups the file does not
             link and with every unknown base free, one haplotype explains every junction

The haplotype sequence is built by applying the alleles as written. Its motif is found without
coordinates: the ten exon bases on either side of the intron are located in the sequence, and the
two bases after the left ten and before the right ten are the motif (none between them: the
haplotype lacks the intron).

Usage: vcf_random_invariance.py <panisoguard> [--n 300] [--seed 1] [--jobs 4] [--bgzip PATH]
                                [--tabix PATH] [--work DIR]
Exit 0 when every check holds, 1 otherwise; prints the failing files.
"""
import argparse, collections, csv, itertools, json, os, random, shutil, subprocess, sys, tempfile
import multiprocessing as mp

BASES = "ACGT"
REF_MOTIFS = [("GG", "AA"), ("GG", "AG"), ("GA", "AG"), ("CT", "AG"), ("GT", "AA"), ("GT", "CG"), ("AT", "AG"), ("GC", "AA"), ("TT", "AG"), ("AT", "AA")]
CANON = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
LEN = 420
JUNC = [(60, 100), (160, 200)]       # introns, 0-based half-open; the last 180 bases hold no motif
FAR = 400
HEADER = (f'##fileformat=VCFv4.2\n##contig=<ID=chrT2,length={LEN}>\n##ALT=<ID=DEL,Description="Deletion">\n'
          '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
          '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n')
PIG = BGZIP = TABIX = None


def ref_seq(mi, nj, rep):
    """A sequence without two equal neighbours, the motif pair planted at each intron; with `rep`
    the motifs repeat into the exons, so an indel can slide across the boundary."""
    rng = random.Random(7919 + mi); seq = []
    while len(seq) < LEN:
        b = rng.choice(BASES)
        if not seq or seq[-1] != b: seq.append(b)
    don, acc = REF_MOTIFS[mi]
    for s, e in JUNC[:nj]:
        seq[s:s + 2] = don; seq[e - 2:e] = acc
        if rep: seq[s - 4:s] = don * 2; seq[e:e + 4] = acc * 2
    return "".join(seq)


def make_locus(root, mi, nj, rep):
    d = os.path.join(root, f"locus_m{mi}_j{nj}_r{rep}"); os.makedirs(d, exist_ok=True)
    with open(f"{d}/ref.fa", "w") as f: f.write(">chrT2\n" + ref_seq(mi, nj, rep) + "\n")
    with open(f"{d}/ref.fa.fai", "w") as f: f.write(f"chrT2\t{LEN}\t7\t{LEN}\t{LEN + 1}\n")
    ex = [(21, 60), (101, 160)] + ([(201, 240)] if nj == 2 else [])
    with open(f"{d}/caller.gtf", "w") as f:
        f.write("".join(f'chrT2\tx\texon\t{a}\t{b}\t.\t+\t.\tgene_id "G1"; transcript_id "iso1";\n' for a, b in ex))
    with open(f"{d}/ref.gtf", "w") as f:
        f.write('chrT2\tx\texon\t1\t10\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\nchrT2\tx\texon\t410\t420\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    with open(f"{d}/sq.tsv", "w") as f:
        f.write("isoform\tchrom\tstrand\tstructural_category\tall_canonical\tRTS_stage\tperc_A_downstream_TTS\niso1\tchrT2\t+\tnovel_not_in_catalog\tnon_canonical\tFALSE\t0\n")
    return d


# ---- random records (pos is 0-based) ------------------------------------------------------------
def rand_alt(rng, ref, keep=(), pos=0, target=None):
    """A same-length ALT; at the motif positions it takes the base of the file's canonical target 7 times in 10."""
    while True:
        alt = "".join(ref[i] if i in keep else (target[pos + i] if target and pos + i in target and rng.random() < 0.7 else rng.choice(BASES)) for i in range(len(ref)))
        if alt != ref: return alt


def gen_genotype(rng, n_alts):
    pool = ["."] + ["0"] * 3 + ["1"] * 5 + (["2"] * 2 if n_alts == 2 else ["2"] if rng.random() < 0.03 else [])
    u = rng.random()
    if u < 0.04: gt, sep = None, ""
    elif u < 0.08: gt, sep = (rng.choice(pool),), ""
    else: gt, sep = (rng.choice(pool), rng.choice(pool)), ("|" if rng.random() < 0.65 else "/")
    ps = None
    if gt is not None and len(gt) == 2:
        ps = rng.choice([None, None, 1, 1, 1, 2, 2]) if sep == "|" else rng.choice([None] * 8 + [1, 2])
    return gt, sep, ps


def gen_record(rng, seq, nj, target=None):
    s, e = rng.choice(JUNC[:nj]); d, a = s, e - 2; r = rng.random(); info = "."
    if r < 0.66:                                                     # same-length alleles
        t = rng.random()
        if t < 0.66: pos = rng.choice([d, d + 1, a, a + 1]); n = 1
        elif t < 0.84: pos = rng.choice([d, a]); n = 2
        elif t < 0.91: pos = rng.choice([d - 1, a - 1]); n = 3
        elif t < 0.96: pos = rng.choice([d + 1, a + 1]); n = 2
        else: pos = d; n = a + 2 - d                                 # one MNP over both windows of the junction
        ref = seq[pos:pos + n]
        keep = tuple(sorted(set(range(n)) - {d - pos, d + 1 - pos, a - pos, a + 1 - pos})) if n > 3 else ()
        if n == 3: keep = (0,) if rng.random() < 0.5 else ()
        alts = [rand_alt(rng, ref, keep, pos, target)]
        if rng.random() < 0.18:
            a2 = rand_alt(rng, ref, keep, pos, target if rng.random() < 0.3 else None)
            if a2 != alts[0]: alts.append(a2)
    elif r < 0.78:                                                   # insertion after pos, from the boundary inwards
        left = rng.random() < 0.5
        pos = rng.choice([s - 1, s - 1, s, s + 1, s + 2] if left else [e - 1, e - 1, e - 2, e - 3, e - 4])
        want = (target[d] + target[d + 1]) if target and left else (target[a] + target[a + 1]) if target else "GT"
        ins = rng.choice([want, want, want[0], want[1], want + rng.choice(BASES), rng.choice(BASES) + want,
                          "".join(rng.choice(BASES) for _ in range(rng.choice([1, 2, 3])))])
        ref = seq[pos]; alts = [ref + ins]
        if rng.random() < 0.15: alts.append(rand_alt(rng, ref))
    elif r < 0.89:                                                   # deletion of 1-3 intron bases next to a boundary
        k = rng.choice([1, 1, 2, 3])
        q0 = rng.choice([s, s, s + 1, s + 2]) if rng.random() < 0.5 else rng.choice([e, e, e - 1, e - 2]) - k
        pos = q0 - 1; ref = seq[pos:pos + 1 + k]; alts = [ref[0]]
        if rng.random() < 0.15: alts.append(rand_alt(rng, ref))
    elif r < 0.93:                                                   # the whole intron
        pos = s - 1; ref = seq[pos:e]; alts = [ref[0]]
    else:                                                            # no sequence given
        pos = rng.choice([d, d + 1, a, a + 1]); ref = seq[pos]; alts = [rng.choice(["<DEL>", "<DEL>", "*"])]
        if alts[0] == "<DEL>" and rng.random() < 0.5: info = f"END={pos + 1 + rng.choice([1, 2, e - pos - 1])}"
        if rng.random() < 0.3: alts.append(rand_alt(rng, ref))
    gt, sep, ps = gen_genotype(rng, len(alts))
    return {"pos": pos, "ref": ref, "alts": alts, "gt": gt, "sep": sep, "ps": ps, "info": info}


def gen_file(rng, seq, nj):
    tgt = rng.choice(sorted(CANON)); target = {}
    for s, e in JUNC[:nj]: target.update({s: tgt[0][0], s + 1: tgt[0][1], e - 2: tgt[1][0], e - 1: tgt[1][1]})
    recs = [gen_record(rng, seq, nj, target) for _ in range(rng.choice([1, 2, 2, 3, 3, 4, 5]))]
    if len(recs) < 5 and rng.random() < 0.15:                         # the same site written again, with the same or another genotype
        r = dict(rng.choice(recs))
        if rng.random() < 0.6:
            g = gen_record(rng, seq, nj, target)
            if g["gt"] is not None and len(g["alts"]) >= len(r["alts"]): r.update(gt=g["gt"], sep=g["sep"], ps=g["ps"])
        recs.append(r)
    rng.shuffle(recs)
    return sorted(recs, key=lambda r: r["pos"])


def gen_file_near(rng, seq, nj):
    """A file that puts a canonical motif on one haplotype (by substitutions, or by an insertion or
    deletion where one fits), then 0-3 changes to it: an added record (also at an occupied
    position), another phase set, a slash, swapped or missing alleles, a repeat with another
    genotype. These sit where the way of writing could tip the call."""
    recs = []
    for s, e in JUNC[:nj]:
        d, a = s, e - 2
        cands = sorted(sorted(CANON), key=lambda t: (sum(x != y for x, y in zip(t[0] + t[1], seq[d:d + 2] + seq[a:a + 2])), rng.random()))
        tgt = cands[0] if rng.random() < 0.8 else rng.choice(cands)
        for p0, want, left in ((d, tgt[0], True), (a, tgt[1], False)):
            cur = seq[p0:p0 + 2]; diff = [i for i in (0, 1) if cur[i] != want[i]]
            if not diff: continue
            u = rng.random()
            if u < 0.12 and left and seq[s + 1:s + 3] == want:        # dropping the first intron base leaves the motif
                recs.append({"pos": s - 1, "ref": seq[s - 1:s + 1], "alts": [seq[s - 1]]})
            elif u < 0.12 and not left and seq[e - 3:e - 1] == want:
                recs.append({"pos": e - 2, "ref": seq[e - 2:e], "alts": [seq[e - 2]]})
            elif u < 0.24:                                            # the motif inserted one base inside the intron
                p = s if left else e - 2
                recs.append({"pos": p, "ref": seq[p], "alts": [seq[p] + (want[1] + rng.choice(BASES) if left else rng.choice(BASES) + want[0])]}
                            if (seq[p] == want[0] if left else seq[e - 1] == want[1]) else {"pos": p0, "ref": cur, "alts": [want]})
            elif len(diff) == 2 and u < 0.62 or len(diff) == 1 and u < 0.35: recs.append({"pos": p0, "ref": cur, "alts": [want]})
            else: recs.extend({"pos": p0 + i, "ref": cur[i], "alts": [want[i]]} for i in diff)
    scheme = rng.choice(["het1", "het1", "het0", "hom", "flip2", "mixed", "mixed"])
    for r in recs:
        g = {"het1": (("1", "0"), "|", 1), "het0": (("1", "0"), "|", None), "hom": (("1", "1"), rng.choice("|/"), None), "flip2": (("0", "1"), "|", 2),
             "mixed": rng.choice([(("1", "1"), "|", None), (("1", "0"), "|", 1), (("1", "0"), "|", 1)])}[scheme]
        r.update(gt=g[0], sep=g[1], ps=g[2], info=".")
    target = {}
    for s, e in JUNC[:nj]:
        t = rng.choice(sorted(CANON)); target.update({s: t[0][0], s + 1: t[0][1], e - 2: t[1][0], e - 1: t[1][1]})
    for _ in range(rng.choice([0, 1, 1, 2, 2, 3])):
        u = rng.random(); j = rng.randrange(len(recs)) if recs else 0
        two = bool(recs) and recs[j]["gt"] is not None and len(recs[j]["gt"]) == 2
        if u < 0.25 or not recs: recs.append(gen_record(rng, seq, nj, target))
        elif u < 0.45:                                              # another record on an occupied position
            b = recs[j]; alt = rand_alt(rng, b["ref"])
            if alt not in b["alts"]:
                g = gen_record(rng, seq, nj, target)
                src = g if g["gt"] is not None and len(g["alts"]) == 1 else b
                recs.append({"pos": b["pos"], "ref": b["ref"], "alts": [alt], "gt": src["gt"], "sep": src["sep"], "ps": src["ps"], "info": "."})
        elif u < 0.60 and two: recs[j] = dict(recs[j], ps=rng.choice([None, 1, 2]))
        elif u < 0.70 and two: recs[j] = dict(recs[j], sep="/", ps=None)
        elif u < 0.80 and two: recs[j] = swapped(recs[j])
        elif u < 0.88 and two:
            g = list(recs[j]["gt"]); g[rng.randrange(2)] = "."; recs[j] = dict(recs[j], gt=tuple(g))
        elif u < 0.96 and two:                                      # the same site again with another genotype
            g = rng.choice([("0", "1"), ("1", "0"), ("1", "1"), ("0", "0")]); recs.append(dict(recs[j], gt=g, ps=rng.choice([recs[j]["ps"], 1, 2, None])))
        elif recs: recs[j] = dict(recs[j], gt=rng.choice([None, ("1",)]), sep="", ps=None)
    rng.shuffle(recs)
    return sorted(recs, key=lambda r: r["pos"])


def row(r):
    head = f"chrT2\t{r['pos'] + 1}\t.\t{r['ref']}\t{','.join(r['alts'])}\t.\tPASS\t{r.get('info', '.')}\t"
    if r["gt"] is None: return head + "DP\t10\n"
    g = r["sep"].join(r["gt"]) if len(r["gt"]) == 2 else r["gt"][0]
    return head + (f"GT:PS\t{g}:{r['ps']}\n" if r["ps"] is not None else f"GT\t{g}\n")


def key(r): return (r["pos"], r["ref"], tuple(r["alts"]))
def gkey(r):
    """What the genotype states: a phased one its alleles in order and its phase set, an unphased one its alleles in no order."""
    if r["gt"] is None: return None
    if r["sep"] == "/": return ("/", tuple(sorted(r["gt"])))
    return (r["sep"], r["gt"], r["ps"])
def span(r): return range(r["pos"], r["pos"] + len(r["ref"]))
def is_sub(a, ref): return len(a) == len(ref) and a[0] not in "<*"
def defined(r, a): return a == "." or int(a) <= len(r["alts"])
def swapped(r): return dict(r, gt=(r["gt"][1], r["gt"][0]))
def group_of(r):
    if r["gt"] is None or len(r["gt"]) != 2: return None
    if r["sep"] == "/": return ("U", key(r), gkey(r))
    return ("PS", r["ps"])


def features(recs, seq):
    f = {"haploid": any(r["gt"] is not None and len(r["gt"]) == 1 for r in recs),
         "gtless": any(r["gt"] is None for r in recs),
         "missing": any(r["gt"] is not None and any(a == "." or not defined(r, a) for a in r["gt"]) for r in recs),
         "nonsub": any(not is_sub(a, r["ref"]) for r in recs for a in r["alts"]),
         "symbolic": any(a[0] in "<*" for r in recs for a in r["alts"]),
         "refbad": any(seq[r["pos"]:r["pos"] + len(r["ref"])] != r["ref"] for r in recs)}
    uniq = {}
    for r in recs: uniq.setdefault(key(r), set()).add(gkey(r))
    f["conflict_dup"] = any(len(v) > 1 for v in uniq.values())
    ks = list(uniq)
    f["overlap"] = any(set(range(a[0], a[0] + len(a[1]))) & set(range(b[0], b[0] + len(b[1]))) for a, b in itertools.combinations(ks, 2))
    return f


# ---- running the tool ---------------------------------------------------------------------------
def run_vcf(locus, work, name, recs):
    vcf = f"{work}/{name}.vcf"
    with open(vcf, "w") as f: f.write(HEADER + "".join(row(r) for r in recs))
    with open(vcf + ".gz", "wb") as f: subprocess.run([BGZIP, "-c", vcf], stdout=f, check=True)
    subprocess.run([TABIX, "-f", "-p", "vcf", vcf + ".gz"], check=True)
    return run_pig(locus, work, name, ["--haplotype-vcf", vcf + ".gz"])


def run_pig(locus, work, name, hap_args):
    p = subprocess.run([PIG, "adjudicate", "--classification", f"{locus}/sq.tsv", "--isoforms-gtf", f"{locus}/caller.gtf", "--ref-gtf", f"{locus}/ref.gtf",
                        "--reference", f"{locus}/ref.fa", *hap_args, "--haplotype-provenance", "wgs", "--out-prefix", f"{work}/{name}"], capture_output=True, text=True)
    if p.returncode: return ["ERROR", p.stderr.strip().splitlines()[-1][:120] if p.stderr.strip() else str(p.returncode), ""]
    with open(f"{work}/{name}.attribution.jsonl") as f: r = json.loads(f.readline())
    with open(f"{work}/{name}.junctions.tsv") as f: js = list(csv.DictReader(f, delimiter="\t"))
    return [r["confidence_class"], r["verdict_reason"], ";".join(j["hap_motifs"] + "/" + j["variant_phase"] for j in js)]


# ---- what the file supports ---------------------------------------------------------------------
def orientations(items):
    groups = sorted({str(group_of(r)) for r, mode in items if mode == "gt" and r["gt"][0] != r["gt"][1]})
    for bits in itertools.product((0, 1), repeat=len(groups)): yield dict(zip(groups, bits))


def stated(recs, seq):
    """One (record, mode) per site: "gt" when the genotype can be applied, "unknown" when the site
    is written again with another genotype, has no GT, or its REF is not the reference."""
    uniq = {}
    for r in recs: uniq.setdefault(key(r), []).append(r)
    return [(rs[0], "unknown" if len({gkey(r) for r in rs}) > 1 or rs[0]["gt"] is None or seq[k[0]:k[0] + len(k[1])] != k[1] else "gt")
            for k, rs in uniq.items()]


def guaranteed_cells(recs, seq, nj):
    """Substitutions only. True when under every orientation one haplotype is known and canonical at
    every junction; a base two records disagree on, or one the file leaves open, is unknown."""
    items = stated(recs, seq)
    for o in orientations(items):
        cells = [collections.defaultdict(set), collections.defaultdict(set)]
        for r, mode in items:
            if mode == "unknown":
                for h in (0, 1):
                    for p in span(r): cells[h][p].add("?")
                continue
            flip = o.get(str(group_of(r)), 0)
            for idx, a in enumerate(r["gt"]):
                h = idx ^ flip
                if a == "." or not defined(r, a):
                    for p in span(r): cells[h][p].add("?")
                elif a != "0":
                    for i, p in enumerate(span(r)): cells[h][p].add(r["alts"][int(a) - 1][i])
        ok = False
        for h in (0, 1):
            def g(p):
                c = cells[h].get(p)
                return seq[p] if not c else next(iter(c)) if len(c) == 1 else "?"
            if all((g(s) + g(s + 1), g(e - 2) + g(e - 1)) in CANON for s, e in JUNC[:nj]): ok = True
        if not ok: return False
    return True


def haplotypes(items, seq, o):
    """The two haplotype sequences with every allele applied as written, each with the sequence
    that carries its same-length alleles only (the exon flanks as that haplotype has them)."""
    out = []
    for h in (0, 1):
        H = list(seq); S = list(seq)
        for r, _ in sorted(items, key=lambda x: -x[0]["pos"]):
            a = int(r["gt"][h ^ o.get(str(group_of(r)), 0)])
            if a == 0: continue
            alt = r["alts"][a - 1]
            H[r["pos"]:r["pos"] + len(r["ref"])] = alt
            if len(alt) == len(r["ref"]): S[r["pos"]:r["pos"] + len(r["ref"])] = alt
        out.append(("".join(H), "".join(S)))
    return out


def motif_of(H, S, s, e):
    """The motif of intron [s, e) on haplotype sequence H, DEL when H lacks the intron, None when
    the exon flanks are not found once each."""
    fl, fr = S[s - 10:s], S[e:e + 10]
    if H.count(fl) != 1 or H.count(fr) != 1: return None
    x, y = H.find(fl) + 10, H.find(fr)
    if y == x: return "DEL"
    if y - x < 4: return None
    return H[x:x + 2] + "-" + H[y - 2:y]


def explains(m): return m == "DEL" or tuple(m.split("-")) in CANON


def guaranteed_sequences(items, seq, nj):
    """Files in which every genotype is complete and no two records touch. True / False as
    guaranteed_cells; None when a flank cannot be placed and the answer hangs on it."""
    for o in orientations(items):
        ok = open_ = False
        for H, S in haplotypes(items, seq, o):
            ms = [motif_of(H, S, s, e) for s, e in JUNC[:nj]]
            if any(m is None for m in ms): open_ = True
            elif all(explains(m) for m in ms): ok = True
        if not ok: return None if open_ else False
    return True


# ---- the same haplotypes written another way ----------------------------------------------------
def rewrites(rng, recs, seq):
    """Records rewritten without changing a haplotype: an insertion or deletion moved as far left
    and as far right as the reference allows, and a record padded with the reference bases around
    it. Only for a record with a complete diploid genotype that no other record comes near: an
    allele the file leaves open makes the REF span of its record unknown, so there the span counts."""
    out = []
    for j, r in enumerate(recs):
        if r["gt"] is None or len(r["gt"]) != 2 or len(r["alts"]) != 1 or r["alts"][0][0] in "<*" or any(a not in "01" for a in r["gt"]): continue
        if seq[r["pos"]:r["pos"] + len(r["ref"])] != r["ref"]: continue
        ref, alt = r["ref"], r["alts"][0]
        others = [x for i, x in enumerate(recs) if i != j]
        def free(lo, hi): return lo >= 12 and hi <= LEN - 12 and all(x["pos"] + len(x["ref"]) + 2 <= lo or x["pos"] - 2 >= hi for x in others)
        if len(ref) != len(alt) and (len(alt) == 1 or len(ref) == 1) and ref[0] == alt[0]:
            q0 = r["pos"] + 1; q1 = q0 + len(ref) - 1; ins = alt[1:]
            for step in (-1, 1):
                a, b, s = q0, q1, ins
                while 12 < a and b < LEN - 12:
                    if step < 0:
                        if b > a:
                            if seq[a - 1] != seq[b - 1]: break
                        else:
                            if seq[a - 1] != s[-1]: break
                            s = s[-1] + s[:-1]
                    else:
                        if b > a:
                            if seq[a] != seq[b]: break
                        else:
                            if seq[a] != s[0]: break
                            s = s[1:] + s[0]
                    a += step; b += step
                if a != q0 and free(min(a, q0) - 1, max(b, q1) + 1):
                    out.append(("shift", j, dict(r, pos=a - 1, ref=seq[a - 1:b], alts=[seq[a - 1] + s])))
        lp, rp = rng.choice([0, 1, 2]), rng.choice([0, 1, 2])
        if lp + rp and free(r["pos"] - lp, r["pos"] + len(ref) + rp):
            L, R = seq[r["pos"] - lp:r["pos"]], seq[r["pos"] + len(ref):r["pos"] + len(ref) + rp]
            out.append(("pad", j, dict(r, pos=r["pos"] - lp, ref=L + ref + R, alts=[L + alt + R])))
    return out


def one_file(args):
    i, seed, root = args
    rng = random.Random(seed * 1000003 + i)
    mi = rng.randrange(len(REF_MOTIFS)); nj = 1 if rng.random() < 0.6 else 2; rep = int(rng.random() < 0.35)
    locus = os.path.join(root, f"locus_m{mi}_j{nj}_r{rep}"); seq = ref_seq(mi, nj, rep)
    work = os.path.join(root, "w", f"{i % 97:02d}"); os.makedirs(work, exist_ok=True)
    recs = gen_file(rng, seq, nj) if rng.random() < 0.3 else gen_file_near(rng, seq, nj)
    if not recs: recs = gen_file(rng, seq, nj)
    if rng.random() < 0.04:                                           # a REF that is not the reference
        j = rng.randrange(len(recs)); r = recs[j]
        recs[j] = dict(r, ref={"A": "C", "C": "G", "G": "T", "T": "A"}[r["ref"][0]] + r["ref"][1:])
    feat = features(recs, seq); n = f"f{i}"
    base = run_vcf(locus, work, n, recs)
    out = {"i": i, "locus": os.path.basename(locus), "rows": [row(r).strip() for r in recs], "feat": feat, "base": base, "fail": {}, "ran": [], "stat": {}}
    def unordered(res):   # the haplotype motifs of each junction in no order: which index carries one is a label
        return res[:2] + [";".join(",".join(sorted(x.split("/")[0].split(","))) + "/" + x.split("/")[1] for x in res[2].split(";"))] if res[2] else res
    def check(name, tag, alt_recs, upto=3, view=lambda r: r):
        res = run_vcf(locus, work, f"{n}_{tag}", alt_recs if name == "order" else sorted(alt_recs, key=lambda r: r["pos"])); out["ran"].append(name)
        if view(res)[:upto] != view(base)[:upto]: out["fail"].setdefault(name, []).append({"rows": [row(r).strip() for r in alt_recs], "got": res})
    # order: rows that share a POS
    idx = {}
    for j, r in enumerate(recs): idx.setdefault(r["pos"], []).append(j)
    if any(len(v) > 1 for v in idx.values()):
        seen = 0
        for combo in itertools.product(*[list(itertools.permutations(v)) for v in idx.values()]):
            order = [j for grp in combo for j in grp]
            if order == list(range(len(recs))): continue
            check("order", f"o{seen}", [recs[j] for j in order], view=unordered); seen += 1
            if seen >= 4: break
    if not feat["haploid"]:
        dip = [r for r in recs if r["gt"] is not None]
        if dip:
            check("relabel", "rl", [swapped(r) if r["gt"] is not None else r for r in recs], view=unordered)
            gs = sorted({str(group_of(r)) for r in dip})
            if len(gs) > 1:
                for gi, g in enumerate(gs[:3]):
                    check("group", f"g{gi}", [swapped(r) if r["gt"] is not None and str(group_of(r)) == g else r for r in recs], view=unordered)
    j = rng.randrange(len(recs)); check("duplicate", "dup", recs[:j + 1] + [recs[j]] + recs[j + 1:])
    # a record where no motif is read: beyond every window, and in a window in the middle of an
    # intron (clear of every record; an allele without a sequence may reach anywhere). Its genotype
    # has the file's ploidy, which a wider genotype would change.
    s0, e0 = JUNC[0]
    shaped = [r for r in recs if r["gt"] is not None]
    like = rng.choice(shaped) if shaped else {"gt": ("1", "1"), "sep": "|", "ps": None}
    gt = tuple(rng.choice("01") for _ in like["gt"])
    for tag, pos in (("far", FAR), ("mid", (s0 + e0) // 2)):
        if tag == "mid" and (feat["symbolic"] or any(r["pos"] - 14 <= pos < r["pos"] + len(r["ref"]) + 14 for r in recs)): continue
        alt = rng.choice([b for b in BASES if b != seq[pos]])
        check("distractor", tag, recs + [{"pos": pos, "ref": seq[pos], "alts": [alt], "gt": gt, "sep": like["sep"], "ps": like["ps"], "info": "."}])
    for k, (kind, j, new) in enumerate(rewrites(rng, recs, seq)[:4]):
        check("rewrite", f"{kind}{k}", recs[:j] + [new] + recs[j + 1:])
    subs = not (feat["haploid"] or feat["nonsub"])
    plain = not (feat["gtless"] or feat["missing"] or feat["conflict_dup"] or feat["refbad"] or feat["overlap"] or feat["haploid"])
    if subs and plain and all(r["sep"] == "|" for r in recs) and len({r["ps"] for r in recs}) == 1:
        haps = [list(seq), list(seq)]
        for r in recs:
            for h, a in enumerate(r["gt"]):
                if a != "0": haps[h][r["pos"]:r["pos"] + len(r["ref"])] = r["alts"][int(a) - 1]
        for nm, s in zip("AB", haps):
            with open(f"{work}/{n}_hap{nm}.fa", "w") as f: f.write(">chrT2\n" + "".join(s) + "\n")
            with open(f"{work}/{n}_hap{nm}.fa.fai", "w") as f: f.write(f"chrT2\t{LEN}\t7\t{LEN}\t{LEN + 1}\n")
        res = run_pig(locus, work, f"{n}_fa", ["--reference-haplotype", f"{work}/{n}_hapA.fa", "--reference-haplotype", f"{work}/{n}_hapB.fa"]); out["ran"].append("fasta")
        if res[:2] != base[:2]: out["fail"].setdefault("fasta", []).append({"got": res})
    g = None
    if subs: g = guaranteed_cells(recs, seq, nj)
    elif plain and not feat["symbolic"]:
        items = stated(recs, seq); g = guaranteed_sequences(items, seq, nj)
        # the motifs the tool reports, haplotype by haplotype, against the sequences as written
        out["ran"].append("motif")
        told = [x.split("/")[0] for x in base[2].split(";")] if base[2] else []
        for (s, e), t in zip(JUNC[:nj], told):
            for (H, S), m in zip(haplotypes(items, seq, {}), t.split(",") if t != "NA" else ["NA", "NA"]):
                want = motif_of(H, S, s, e)
                state = "open" if want is None else "held" if m == "NA" else "same" if m == want else "differs"
                out["stat"][state] = out["stat"].get(state, 0) + 1
                if state == "differs": out["fail"].setdefault("motif", []).append({"junction": [s, e], "tool": m, "sequence": want})
    if g is not None:
        out["ran"].append("sound"); out["guaranteed"] = g
        if base[0] == "REFERENCE_BIAS" and not g: out["fail"].setdefault("sound", []).append({"got": base})
    return out


def _init(pig, bgzip, tabix):
    global PIG, BGZIP, TABIX
    PIG, BGZIP, TABIX = pig, bgzip, tabix


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("panisoguard"); ap.add_argument("--n", type=int, default=300); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--jobs", type=int, default=4); ap.add_argument("--bgzip", default="bgzip"); ap.add_argument("--tabix", default="tabix")
    ap.add_argument("--work"); ap.add_argument("--keep", action="store_true", help="keep the work directory")
    a = ap.parse_args()
    work = a.work or tempfile.mkdtemp(prefix="pig_random_")
    os.makedirs(work, exist_ok=True)
    for mi in range(len(REF_MOTIFS)):
        for nj in (1, 2):
            for rep in (0, 1): make_locus(work, mi, nj, rep)
    with mp.Pool(a.jobs, initializer=_init, initargs=(os.path.abspath(a.panisoguard), a.bgzip, a.tabix)) as pool:
        res = pool.map(one_file, [(i, a.seed, work) for i in range(a.n)], chunksize=4)
    with open(os.path.join(work, "results.jsonl"), "w") as f:
        for r in res: f.write(json.dumps(r) + "\n")
    checks = ("order", "relabel", "group", "duplicate", "distractor", "rewrite", "fasta", "motif", "sound")
    print(f"files {len(res)}  tool runs {sum(1 + len(r['ran']) for r in res)}  classes {dict(sorted(collections.Counter(r['base'][0] for r in res).items()))}")
    failed = 0
    for c in checks:
        ran = [r for r in res if c in r["ran"]]; bad = [r for r in ran if c in r["fail"]]; failed += len(bad)
        print(f"{c:10s} files checked {len(ran):5d}  failing {len(bad):4d}")
    stat = collections.Counter()
    for r in res: stat.update(r["stat"])
    print(f"motif: haplotype motifs equal to the sequence's {stat['same']}, not reported by the tool {stat['held']}, flank not placed {stat['open']}, different {stat['differs']}")
    g = [r for r in res if "guaranteed" in r]
    print(f"sound: guaranteed and REFERENCE_BIAS {sum(1 for r in g if r['guaranteed'] and r['base'][0] == 'REFERENCE_BIAS')}, "
          f"guaranteed but held {sum(1 for r in g if r['guaranteed'] and r['base'][0] != 'REFERENCE_BIAS')}, not guaranteed {sum(1 for r in g if not r['guaranteed'])}")
    errors = [r for r in res if r["base"][0] == "ERROR"]
    print(f"errors {len(errors)}")
    for r in [r for r in res if r["fail"]][:10] + errors[:5]:
        print(f"\nfile {r['i']} ({r['locus']}): {r['base']}")
        for x in r["rows"]: print("  " + x)
        for c, fs in r["fail"].items():
            for f in fs[:2]: print(f"  {c}: {json.dumps(f)}")
    if not a.keep and not a.work: shutil.rmtree(work, ignore_errors=True)
    sys.exit(1 if failed or errors else 0)


if __name__ == "__main__":
    main()
