# 05. 참조 게놈이 이 사람과 다르면 무엇이 달라질까?

[02](02_short_read_support.md)에서 `iso_A`의 novel junction은 short read 9개가 걸쳤는데도 지지를 못 받았다. 참조 게놈(GRCh38 같은 표준 게놈, 이 예제에서는 `toy.fa`)에서 이 junction이 non-canonical이기 때문이었다. 그런데 이 사람의 게놈에서는 염기 하나가 달라서 canonical이다. "주석에 없고 motif도 이상한 junction"이 사실은 이 사람에게는 평범한 junction일 수 있다는 것이다. 참조 게놈이 모든 사람을 대표하지 못해서 생기는 이런 오판을 reference bias라고 부른다. 이 노트에서는 PanIsoGuard가 개인 게놈 서열(haplotype FASTA)이나 pangenome junction 목록으로 이런 경우를 찾아내는 방법, 그리고 그 판정을 언제 믿고 언제 보류하는지를 본다. 마지막으로 이 기능의 검증이 왜 아직 충분하지 않은지도 적어 둔다.

> 관련 문서: [docs/decision_engine.md](../docs/decision_engine.md) "Projection to confidence classes" 1번, "Provenance and circularity" · 코드: `src/evidence/variant_motif.cpp`, `src/io/pangenome_reader.cpp`, `build_evidence()`

## 1. 이 사람의 게놈에서 `iso_A`의 junction은 어떻게 생겼을까?

예제에는 이 사람의 haplotype 두 개가 들어 있다. 사람은 염색체를 두 벌 가지니, 부모에게서 하나씩 받은 두 서열이다. 참조 게놈과 비교해서 어디가 다른지, 그리고 `iso_A`의 novel intron(0-based `[500, 730)`) 양 끝이 세 서열에서 어떻게 보이는지 확인한다.

```python
import pysam
ref, h1, h2 = (pysam.FastaFile(f"data/{n}.fa") for n in ("toy", "hap1", "hap2"))
diff = [i + 1 for i, (a, b) in enumerate(zip(ref.fetch("chrT"), h1.fetch("chrT"))) if a != b]
print("positions where hap1 differs from the reference (1-based):", diff,
      "| hap2 differs at:", [i + 1 for i, (a, b) in enumerate(zip(ref.fetch("chrT"), h2.fetch("chrT"))) if a != b])
s, e = 500, 730                                    # iso_A's novel intron, 0-based half-open
for name, fa in (("reference", ref), ("hap1", h1), ("hap2", h2)):
    print(f"{name:9s} donor {fa.fetch('chrT', s, s + 2)}  acceptor {fa.fetch('chrT', e - 2, e)}  (intron end context ...{fa.fetch('chrT', e - 8, e)}|{fa.fetch('chrT', e, e + 6)}...)")
```

```text
positions where hap1 differs from the reference (1-based): [730] | hap2 differs at: []
reference donor GT  acceptor AC  (intron end context ...GATAGGAC|TAGATC...)
hap1      donor GT  acceptor AG  (intron end context ...GATAGGAG|TAGATC...)
hap2      donor GT  acceptor AC  (intron end context ...GATAGGAC|TAGATC...)
```

haplotype 1은 730번 염기 하나만 참조와 다르다(C → G). 이 사람은 이 자리가 이형접합(heterozygous)인 셈이다. 그 한 염기 때문에 haplotype 1에서 `iso_A`의 intron은 GT…AG, 곧 canonical이다. 괄호 안의 `|`는 intron과 다음 exon의 경계인데, 바뀐 염기가 `|` 왼쪽, 곧 intron 안쪽에 있다는 점을 기억해 두자(4절).

## 2. PanIsoGuard는 haplotype으로 무엇을 판정할까?

`--reference`(참조 FASTA)와 `--reference-haplotype`(haplotype FASTA, 여러 번 줄 수 있음)를 주면, PanIsoGuard는 novel junction마다 참조와 haplotype에서 motif를 따로 읽어 세 가지 중 하나로 분류한다(`HaplotypeProvider::classify()`).

| 분류 | 조건 | 쓰임 |
|---|---|---|
| `CREATED` | 참조에서는 non-canonical, haplotype **하나 이상**에서 canonical | rescue 후보 |
| `DISRUPTED` | 참조에서는 canonical, **모든** haplotype에서 non-canonical | 계산만 하고 아직 쓰지 않음 |
| `NONE` | 그 밖의 경우 | 없음 |

그리고 isoform의 **모든** novel junction이 `CREATED`일 때만 isoform을 reference bias로 본다(`variant_rescue`). 규칙을 Python으로 옮겨서 isoform마다 계산하고 PanIsoGuard의 결과(00의 step5)와 비교한다. novel junction 목록은 [01](01_intron_chain_and_novelty.md) 3절의 결과를 그대로 옮겨 적었다.

```python
import json, pysam

COMP = str.maketrans("ACGT", "TGCA")
CANONICAL = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}

def motif(fa, chrom, s, e, strand):
    """Transcript-oriented donor/acceptor of intron [s, e) and whether it is canonical."""
    left, right = fa.fetch(chrom, s, s + 2).upper(), fa.fetch(chrom, e - 2, e).upper()
    if strand == "-":                                    # the donor is at the genomic right end
        left, right = right[::-1].translate(COMP), left[::-1].translate(COMP)
    return f"{left}-{right}", (left, right) in CANONICAL

def classify(s, e, strand="+"):
    ref_m, ref_ok = motif(REF, "chrT", s, e, strand)
    haps = [motif(h, "chrT", s, e, strand) for h in HAPS]
    if not ref_ok and any(ok for _, ok in haps):
        return "CREATED", ref_m, [m for m, _ in haps]
    if ref_ok and not any(ok for _, ok in haps):
        return "DISRUPTED", ref_m, [m for m, _ in haps]
    return "NONE", ref_m, [m for m, _ in haps]

REF = pysam.FastaFile("data/toy.fa")
HAPS = [pysam.FastaFile("data/hap1.fa"), pysam.FastaFile("data/hap2.fa")]
novel = {"iso_A": [(500, 730)], "iso_B": [(800, 1020)], "iso_C": [(800, 860), (940, 1000)],
         "iso_E": [(500, 1000)], "iso_F": [(200, 450)], "iso_G": [(800, 1004)]}   # from note 01
pig = {json.loads(l)["isoform"]: json.loads(l)["evidence"] for l in open("out/step5.attribution.jsonl")}
for iso, js in novel.items():
    verdicts = [classify(s, e) for s, e in js]
    rescue = all(v[0] == "CREATED" for v in verdicts)
    print(f"{iso}: {[(v[0], v[1], v[2]) for v in verdicts]} -> rescue {rescue} | PanIsoGuard variant_rescue {pig[iso]['variant_rescue']}")
print("same intron read as '-' strand:", motif(REF, "chrT", 800, 1020, "+"), "vs", motif(REF, "chrT", 800, 1020, "-"))
```

```text
iso_A: [('CREATED', 'GT-AC', ['GT-AG', 'GT-AC'])] -> rescue True | PanIsoGuard variant_rescue True
iso_B: [('NONE', 'GT-AG', ['GT-AG', 'GT-AG'])] -> rescue False | PanIsoGuard variant_rescue False
iso_C: [('NONE', 'GT-AG', ['GT-AG', 'GT-AG']), ('NONE', 'GT-AG', ['GT-AG', 'GT-AG'])] -> rescue False | PanIsoGuard variant_rescue False
iso_E: [('NONE', 'GT-AG', ['GT-AG', 'GT-AG'])] -> rescue False | PanIsoGuard variant_rescue False
iso_F: [('NONE', 'GT-TT', ['GT-TT', 'GT-TT'])] -> rescue False | PanIsoGuard variant_rescue False
iso_G: [('NONE', 'GT-AG', ['GT-AG', 'GT-AG'])] -> rescue False | PanIsoGuard variant_rescue False
same intron read as '-' strand: ('GT-AG', True) vs ('CT-AC', False)
```

`iso_A`만 `CREATED`이고, haplotype 두 개 중 하나(hap1)에서만 canonical이어도 된다. `iso_F`는 두 haplotype에서도 GT-TT라서 구제되지 않는다. 이 사람에게서도 non-canonical이니, 여전히 가짜일 가능성이 크다.

마지막 줄은 strand 처리를 보여 준다. 같은 좌표라도 `-` strand 전사체라면 전사 방향으로 읽어야 하니, 오른쪽 끝을 역상보로 바꾼 것이 donor, 왼쪽 끝을 역상보로 바꾼 것이 acceptor다. `+`로 읽어 GT-AG인 intron을 `-`로 읽으면 CT-AC가 되어 non-canonical이다. strand를 틀리게 주면 motif 판정도 틀린다.

## 3. 출처가 불확실한 haplotype은 왜 믿지 않을까?

haplotype FASTA가 어디서 왔는지가 중요하다. 이 사람의 DNA(WGS나 HiFi 조립)에서 만든 서열이라면 RNA와 독립된 증거다. 그런데 **판정하려는 바로 그 RNA read**에서 변이를 불러 만든 서열이라면 이야기가 다르다. RNA read가 non-canonical junction을 억지로 설명하려고 정렬되면서 생긴 mismatch가 "변이"로 불리고, 그 변이가 다시 그 junction을 "canonical"로 만들어 줄 수 있다. 자기 자신을 근거로 자기를 구제하는 순환 논리다.

PanIsoGuard는 haplotype이 어디서 왔는지 확인할 방법이 없으니, 사용자가 `--haplotype-provenance`로 알려 줘야 한다. 이 값이 정확히 `wgs`나 `external`일 때만 구제하고, 나머지는 모두 `AMBIGUOUS`로 보류하면서 `circularity_flag`를 켠다. 이것을 circularity firewall이라고 부른다. 값을 바꿔 가며 돌려 봤다.

```bash
PIG=../build/panisoguard; mkdir -p out
A="--classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab --reference data/toy.fa --reference-haplotype data/hap1.fa --reference-haplotype data/hap2.fa"
for p in "" rna_derived unknown WGS wgs external; do
  $PIG adjudicate $A ${p:+--haplotype-provenance $p} --out-prefix out/prov 2>/dev/null
  printf '%-14s %-50s %s\n' "${p:-(not given)}" "$(grep '^axis.variant' out/prov.provenance.log | cut -f2)" \
    "$(grep -P '^iso_A\t' out/prov.adjudicated.tsv | cut -f6,7,10 | tr '\t' ' ')"
done
```

```text
(not given)    on (circular-risk: held, not promoted)             variant_created AMBIGUOUS true
rna_derived    on (circular-risk: held, not promoted)             variant_created AMBIGUOUS true
unknown        on (circular-risk: held, not promoted)             variant_created AMBIGUOUS true
WGS            on (circular-risk: held, not promoted)             variant_created AMBIGUOUS true
wgs            on                                                 variant_created PAN_REF_RESCUED_FALSE_NOVEL false
external       on                                                 variant_created PAN_REF_RESCUED_FALSE_NOVEL false
```

값을 주지 않으면 `unknown`으로 보고 보류한다. 대문자 `WGS`도 보류된다. 문자열이 정확히 같아야 하고 오류도 내지 않는다. 오타가 나면 구제 쪽이 아니라 보류 쪽으로 빠지게 만든 것이다. 결과를 읽을 때는 `provenance.log`의 `axis.variant` 줄에서 `circular-risk`가 붙었는지 보면 된다.

rescue는 isoform의 **모든** novel junction이 `CREATED`여야 한다. 하나만 설명되고 나머지가 진짜 새로운 junction일 수도 있기 때문이다. `iso_A`의 novel junction과 `iso_B`의 novel junction을 둘 다 가진 isoform `iso_AB`를 만들어서 확인했다. SQANTI3 값은 `iso_A`의 줄을 그대로 복사했다.

```python
# iso_AB = iso_A's novel acceptor (730) followed by iso_B's novel acceptor (1020): two novel junctions.
exons = [(101, 200), (401, 500), (731, 800), (1021, 1100), (1301, 1400)]
attr = 'gene_id "GENE_T"; transcript_id "iso_AB";'
with open("out/caller_ab.gtf", "w") as f:
    f.write(open("data/caller1.gtf").read())
    for a, b in exons:
        f.write(f"chrT\tcaller1\texon\t{a}\t{b}\t.\t+\t.\t{attr}\n")
rows = open("data/caller1_classification.txt").read().splitlines()
row_a = next(r for r in rows if r.startswith("iso_A\t"))
with open("out/classification_ab.txt", "w") as f:
    f.write("\n".join(rows + ["iso_AB\t" + row_a.split("\t", 1)[1]]) + "\n")   # same SQANTI3 values as iso_A
print("added iso_AB:", exons)
```

```text
added iso_AB: [(101, 200), (401, 500), (731, 800), (1021, 1100), (1301, 1400)]
```

```bash
$PIG adjudicate --classification out/classification_ab.txt --isoforms-gtf out/caller_ab.gtf \
  --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab --reference data/toy.fa \
  --reference-haplotype data/hap1.fa --reference-haplotype data/hap2.fa --haplotype-provenance wgs \
  --out-prefix out/ab 2>/dev/null
grep '"iso_AB"' out/ab.attribution.jsonl | python3 -c 'import json,sys; r=json.loads(sys.stdin.read()); e=r["evidence"]; print(r["confidence_class"], "| n_novel", e["n_novel_junctions"], "sr", e["n_novel_jx_sr_supported"], "variant_rescue", e["variant_rescue"]); print(r["rule_trace"])'
```

```text
LOW_CONF_PARTIAL | n_novel 2 sr 1 variant_rescue False
['sj_support=1/2 -> PARTIAL', 'all_canonical=non_canonical -> mechanism=noncanonical', 'project(PARTIAL,noncanonical) -> LOW_CONF_PARTIAL']
```

두 junction 중 하나(730)만 `CREATED`라서 rescue가 되지 않고, 조합표로 넘어갔다. 1020 쪽 junction은 short read 12개가 확인했으니 지지 수준은 `PARTIAL`이고, 결과는 `LOW_CONF_PARTIAL`이다. `wgs`를 줬는데도 이렇다.

## 4. short read로는 왜 이 변이를 볼 수 없을까?

1절에서 바뀐 염기(730)가 intron 안쪽에 있다는 점을 짚어 두었다. canonical motif의 두 염기는 정의상 intron의 첫 두 염기와 마지막 두 염기라서, 그 motif를 만든 변이는 늘 intron 안에 있다. intron은 splicing 때 잘려 나가니, 그 junction을 쓰는 mRNA에서 나온 read에는 이 염기가 없다. [02](02_short_read_support.md) 5절의 STAR 정렬 결과에서 730번을 덮는 read를 찾아봤다.

```python
import collections, pysam
# Which short reads actually contain genome position 730 (0-based 729), and which base do they carry?
seen, iso_a = collections.Counter(), 0
for r in pysam.AlignmentFile("out/star/default/Aligned.out.sam"):
    origin = r.query_name.split(".")[0]
    pairs = dict((ref, q) for q, ref in r.get_aligned_pairs(matches_only=True))
    if 729 in pairs:
        seen[(origin, r.query_sequence[pairs[729]])] += 1
    if origin == "iso_A":
        iso_a += 1
print("reads covering chrT:730, by source isoform and base:", dict(seen))
print("iso_A reads in total:", iso_a, "| iso_A reads covering 730:", sum(v for (o, _), v in seen.items() if o == "iso_A"))
```

```text
reads covering chrT:730, by source isoform and base: {('iso_known', 'G'): 25, ('iso_D', 'G'): 10}
iso_A reads in total: 9 | iso_A reads covering 730: 0
```

730번을 덮는 read는 모두 `iso_known`과 `iso_D`에서 왔다. 이 두 isoform에서는 730번이 exon 3 한가운데라 read에 남는다. 정작 `iso_A`의 read 9개는 하나도 이 자리를 덮지 않는다. 예제의 short read는 모두 haplotype 1에서 뽑았기 때문에 G만 보이는데, 실제 이형접합 시료라면 G와 C가 섞여 보였을 것이다.

그러니 "이 junction을 쓰는 read가 변이를 갖고 있는가"는 확인할 수 없는 질문이다. 변이를 확인하려면 DNA를 보거나(haplotype FASTA), 이 junction을 쓰는 read가 **같은 haplotype의 다른 표지**를 갖고 있는지 봐야 한다. 이 점은 6절의 검증 문제와 이어진다.

## 5. pangenome 목록으로도 구제할 수 있을까?

개인 게놈이 없어도, 여러 사람의 조립 게놈을 모은 pangenome graph에서 "이 junction은 누군가의 haplotype 경로 위에 실제로 있다"는 목록을 뽑아 줄 수 있다. PanIsoGuard는 그 목록을 파일(`--pangenome-junctions`)로 받는다. 형식은 탭으로 나눈 `chrom`, `start`, `end`(1-based, `SJ.tab`과 같은 규칙), `strand`, 그리고 선택 열 `n_haplotypes`다. `iso_A`의 junction이 12개 haplotype에 있다는 목록을 만들어 넣어 봤다.

```bash
printf '# chrom\tstart(1-based)\tend\tstrand\tn_haplotypes\nchrT\t501\t730\t+\t12\n' > out/pangenome.tsv
cat out/pangenome.tsv
for p in unknown population; do
  $PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
    --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab \
    --pangenome-junctions out/pangenome.tsv --pangenome-provenance $p --out-prefix out/pan 2>/dev/null
  printf '%-11s %s\n' $p "$(grep -P '^iso_A\t' out/pan.adjudicated.tsv | cut -f6,7,10 | tr '\t' ' ')"
done
grep '"iso_A"' out/pan.attribution.jsonl | python3 -c 'import json,sys; print(json.loads(sys.stdin.read())["rule_trace"])'
```

```text
# chrom	start(1-based)	end	strand	n_haplotypes
chrT	501	730	+	12
unknown     population_known AMBIGUOUS true
population  population_known PAN_REF_RESCUED_FALSE_NOVEL false
['all 1/1 novel junctions realizable on a pangenome haplotype path (independent population data) -> PAN_REF_RESCUED_FALSE_NOVEL (reference bias, population_known)']
```

규칙은 haplotype 축과 같다. isoform의 모든 novel junction이 목록에 있어야 하고, 목록의 출처가 `population`이나 `external`일 때만 구제한다. 목록이 이 시료 자신의 read로 만든 graph에서 나왔다면 다시 순환 논리가 되기 때문이다. 기전 이름만 `population_known`으로 다르다. pangenome 축은 haplotype 축보다 먼저 확인한다. 그래서 개인 haplotype의 출처가 불확실해도, 독립적인 pangenome 목록이 있으면 구제된다.

## 6. 이 rescue는 얼마나 검증되었을까?

이 저장소에는 reference-bias rescue를 실제 데이터에 돌린 벤치마크가 여럿 있다. GIAB 네 사람(137개 junction), HPRC의 서아프리카계 두 사람(86개 junction) 등에서 "N/N 구제, 거짓 0개"라는 결과가 나온다([docs/validation.md](../docs/validation.md)). 그런데 이 숫자가 무엇을 보여 주는지는 따져 봐야 한다. 벤치마크가 정답 라벨 `CREATED`를 붙이는 코드는 이렇다.

```python
# benchmark/hg002/scan_variant_axis.py 발췌 (실행하지 않음)
    cr, ch = canonical(ref, a, b, st), canonical(hap, a, b, st)
    if not cr and ch:
        created.append((a, b, st))
```

참조에서 canonical이 아니고(`not cr`) haplotype에서 canonical이면(`ch`) `CREATED`다. 2절의 `classify()`, 곧 PanIsoGuard의 rescue 규칙과 **같은 조건**이다. 그러니 "`CREATED` 86개를 모두 구제했다"는 것은 코드가 명세대로 동작한다는 확인이다. 이 junction들이 정말 reference bias로 생긴 가짜 novelty인지, 이 사람에게서 실제로 쓰이는 splice site인지는 보여 주지 못한다. "firewall이 86개를 모두 보류했다"도 같다. provenance 값 하나로 갈리는 분기가 그대로 동작했다는 뜻이다.

독립적으로 확인하려면 정답 라벨이 규칙과 다른 곳에서 와야 한다. 4절에서 봤듯이 motif를 만든 변이 자체는 read로 볼 수 없으니, 방법은 두 가지 정도다.

- **phasing**: long read는 exon 여러 개를 한 번에 읽는다. 같은 유전자의 다른 exon에 이형접합 SNP가 있으면, novel junction을 쓰는 long read가 모두 motif를 만든 haplotype(여기서는 hap1) 쪽 allele을 갖고 있는지, 참조 junction을 쓰는 read는 반대쪽(hap2) allele을 갖고 있는지 볼 수 있다. 그렇다면 이 junction은 정말 그 haplotype에서만 쓰인다는 독립적인 증거가 된다.
- **여러 사람 비교**: 같은 변이를 가진 사람에게서만 이 junction이 나타나고, 없는 사람에게서는 나타나지 않는지 본다.

이 노트를 쓰고 나서 첫 번째 방법을 실제 데이터에 해 봤다([benchmark/refbias_phasing](../benchmark/refbias_phasing)). HPRC의 서아프리카계 두 사람(HG03516, HG02717)에서 rescue 대상 junction 가운데 이형접합인 것만 골라, 그 junction을 쓰는 long read가 어느 haplotype에서 왔는지 유전자의 다른 이형접합 SNP로 가렸다. phase는 RNA가 아니라 부모 정보로 나눈 HiFi 조립에서 가져왔다. 결과는 저장소에 들어 있어서 데이터 없이 읽을 수 있다.

```python
import json
m = json.load(open("../benchmark/results/refbias_phasing/metrics.json"))["metrics"]
for who, r in m.items():
    print(f"{who}: {r['junctions']} junctions, {r['heterozygous']} heterozygous | {r['verdicts']}")
    print(f"   junction reads from the motif haplotype / the other: "
          f"{r['tested_junction_reads_motif_haplotype']} / {r['tested_junction_reads_other_haplotype']}"
          f" | other reads at the same loci: {r['tested_locus_reads_motif_haplotype']} / {r['tested_locus_reads_other_haplotype']}")
    for row in r["rows"]:
        if row["verdict"].startswith("CONTRADICTED"):
            print("   contradicted:", {k: row[k] for k in ("junction", "haplotypes", "J_H", "J_O", "J_indel_near")})
```

```text
HG02717: 47 junctions, 31 heterozygous | {'consistent': 16, 'uninformative': 15, 'untestable': 16}
   junction reads from the motif haplotype / the other: 400 / 0 | other reads at the same loci: 1839 / 1078
HG03516: 41 junctions, 23 heterozygous | {'CONTRADICTED': 1, 'consistent': 11, 'uninformative': 11, 'untestable': 18}
   junction reads from the motif haplotype / the other: 172 / 3 | other reads at the same loci: 669 / 654
   contradicted: {'junction': 'chr3:184709997-184735943:+', 'haplotypes': 'mat', 'J_H': 0, 'J_O': 3, 'J_indel_near': 1.0}
```

판정할 수 있었던 이형접합 junction 28개 가운데 27개에서, junction을 쓰는 read는 모두 motif가 canonical인 haplotype에서 왔다. 같은 자리의 다른 read를 보면 두 haplotype이 모두 발현되고 있으니, 우연히 한쪽만 보인 것이 아니다. rescue의 설명("이 사람의 이 haplotype에서는 평범한 splice site다")이 규칙과 독립된 증거로 처음 확인된 셈이다.

예외 하나(HG03516 chr3:184709997)는 junction read 3개가 모두 반대쪽 haplotype, 곧 motif가 non-canonical인 쪽에서 왔다. 세 read 모두 donor 바로 앞에 2 bp 삽입이 있어서(`J_indel_near` 1.0) 정렬이 만든 가짜 junction으로 보인다. [03](03_artifact_mechanisms.md)의 BAM 축이라면 잡았을 흔적이지만, rescue가 조합표보다 먼저 판정을 끝내기 때문에 보지 않는다([08](08_one_isoform_end_to_end.md) 연습문제). 그렇다고 mapping 흔적이 rescue를 막게 바꾸면 되는 것도 아니다. HG02717의 한 junction은 read 81개 모두 같은 흔적이 있는데 phasing 검사를 통과했다. 동형접합 junction 34개는 두 haplotype이 같으니 이 방법으로 검사할 수 없다.

얼마나 자주 일어나는지도 적어 둘 만하다. 벤치마크에서 이런 junction은 서아프리카계 개인 한 명당 40개 남짓이었다. caller가 부른 novel junction은 한 사람당 9만 개쯤이었으니 0.05% 정도다. 드물지만 이 사람에게 고유한 junction이라는 점에서 중요할 수 있는 경우를 가려 주는 안전장치로 보는 것이 맞다.

## 정리

- haplotype 축은 novel junction의 motif를 참조와 haplotype에서 따로 읽어 `CREATED`(참조 non-canonical, haplotype canonical), `DISRUPTED`, `NONE`으로 나눔. isoform의 모든 novel junction이 `CREATED`일 때만 rescue 후보가 됨.
- rescue는 haplotype의 출처가 정확히 `wgs`나 `external`일 때만 `PAN_REF_RESCUED_FALSE_NOVEL`이 되고, 나머지는 `AMBIGUOUS` + `circularity_flag`로 보류됨(대문자 `WGS`도 보류). 판정하려는 RNA에서 만든 haplotype으로 그 RNA를 구제하는 순환을 막으려는 것임.
- pangenome 목록도 같은 규칙(모든 junction, 출처 `population`/`external`)으로 구제하고, haplotype 축보다 먼저 확인함.
- motif를 만드는 변이는 늘 intron 안에 있어서 그 junction을 쓰는 read에는 나타나지 않음. short read가 이 junction을 확인하지 못하는 이유이기도 함.
- 벤치마크의 "N/N 구제, 거짓 0개"는 정답 라벨이 rescue 규칙과 같은 조건이라 구현 확인일 뿐임. 규칙과 독립된 phasing 검사에서는 판정 가능한 이형접합 junction 28개 중 27개가 설명과 맞았고, 1개는 정렬 artifact로 보임. 동형접합 34개는 검사하지 못함.

다음 노트 [06](06_multi_caller_consensus.md)에서는 caller 여러 개의 결과를 합쳐서 short read 대신 쓰는 방법을 본다.

## 연습문제

### 문제 1

> haplotype FASTA를 하나만(`hap2.fa`, 참조와 같은 서열) 줬다면 `iso_A`는 어떻게 될까? 반대로 `hap1.fa`만 줬다면?

<details>
<summary>풀이</summary>

2절의 `classify()`를 haplotype 목록만 바꿔서 다시 부르면 된다. 앞 블록의 함수를 이어서 쓴다.

```python
for haps in (["hap2"], ["hap1"]):
    HAPS = [pysam.FastaFile(f"data/{h}.fa") for h in haps]
    print(haps, classify(500, 730))
```

```text
['hap2'] ('NONE', 'GT-AC', ['GT-AC'])
['hap1'] ('CREATED', 'GT-AC', ['GT-AG'])
```

`hap2`만 주면 참조와 같으니 `NONE`이고 rescue는 없다. 이형접합 변이는 두 haplotype을 모두 줘야 빠짐없이 잡힌다. 반대로 `hap1`만 줘도 `CREATED`가 된다. 조건이 "haplotype 하나 이상에서 canonical"이기 때문이다.

한 가지 조심할 점이 있다. SNV만 반영해서 참조와 길이가 같은 haplotype이어야 좌표가 맞는다. 삽입이나 결실이 들어간 haplotype FASTA는 그 뒤로 좌표가 밀려서, 같은 좌표에서 엉뚱한 두 염기를 읽게 된다. 저장소의 실제 데이터 벤치마크도 SNV만 반영한 haplotype을 만들어 썼다([benchmark/refbias_cohort](../benchmark/refbias_cohort)).

</details>

### 문제 2

> `iso_A`는 step2(`SJ.out.tab`까지)에서 `ARTIFACT`였고 step5(haplotype + `wgs`)에서 `PAN_REF_RESCUED_FALSE_NOVEL`이 되었다. 이 등급 이름의 "FALSE_NOVEL"은 `iso_A`가 가짜라는 뜻일까?

<details>
<summary>풀이</summary>

"novelty가 가짜"라는 뜻이지 "isoform이 가짜"라는 뜻이 아니다. `iso_A`가 non-canonical이라 수상해 보였던 이유가 참조 게놈 쪽에 있었다는 뜻이다. 이 사람의 haplotype 1에서 `iso_A`는 평범한 canonical splicing으로 만들어질 수 있는 isoform이다. 그러니 이 등급은 "버려라"가 아니라 "참조 게놈 기준의 artifact 판정을 거둬라"로 읽어야 한다. 이 사람에게서 실제로 쓰이는 isoform일 수 있고, 그렇다면 오히려 이 사람의 유전형과 연결된 흥미로운 isoform이다. 다만 그것이 정말 쓰이는지는 6절의 독립 검증이 있어야 말할 수 있다.

</details>

## 더 깊이 보기

<details>
<summary><code>DISRUPTED</code>는 왜 쓰지 않나</summary>

`DISRUPTED`는 참조에서는 canonical인데 이 사람의 모든 haplotype에서 non-canonical이 된 junction이다. 이론상 이 사람에게서는 쓰일 수 없는 junction이니, 그런 novel junction을 short read가 확인했다면 오히려 이상한 일이다. `HaplotypeProvider`는 이 경우를 계산하지만, 판정에서 어떻게 다룰지(`ARTIFACT` 쪽으로 볼지, 보류할지) 정하지 않아서 `build_evidence()`가 쓰지 않는다([docs/decision_engine.md](../docs/decision_engine.md) "Future work"). 이 예제에는 `DISRUPTED`인 novel junction이 없다.

</details>

<details>
<summary>pangenome 목록의 <code>n_haplotypes</code>는 무엇에 쓰나</summary>

`[axis_pangenome] min_haplotypes`(기본 1) 이상의 haplotype에서 나온 junction만 인정한다. 한두 사람에게만 있는 드문 junction을 빼고 싶을 때 올린다. 같은 junction이 목록에 여러 줄 있으면 `SJ.tab`처럼 더하지 않고 가장 큰 값을 쓴다(`PangenomeJunctions::add()`). 이미 여러 사람을 센 값이라 더하면 중복으로 세게 되기 때문이다. strand가 `+`, `-`, `1`, `2`가 아니면 오류를 내고 멈춘다. `SJ.tab`의 strand 0을 조용히 "모름"으로 읽는 것과 다르게 동작한다.

</details>
