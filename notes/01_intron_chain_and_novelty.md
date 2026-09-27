# 01. novel isoform은 정확히 무엇이 새로울까?

[00](00_overview.md)의 표에서 `iso_A`의 novel junction 수는 참조 GTF를 넣기 전에는 0, 넣은 뒤에는 1이었다. PanIsoGuard는 무엇과 무엇을 비교해서 "새롭다"고 셀까? 이 노트에서는 isoform을 intron chain으로 바꾸는 방법, 파일 형식마다 다른 좌표를 하나로 맞추는 방법, 그리고 참조 주석과 비교해서 novel junction을 세는 방법을 차례로 따라가 본다. 마지막으로 SQANTI3의 구조 범주(NIC, NNC)와 PanIsoGuard의 novel junction 수가 어떻게 다른지 실제 시뮬레이션 데이터로 확인한다.

> 관련 문서: [docs/input_formats.md](../docs/input_formats.md), [types.hpp](../include/panisoguard/types.hpp)의 좌표 규칙 주석 · 코드: `src/io/gtf_reader.cpp`, `src/io/bed12_reader.cpp`, `src/io/sj_tab_reader.cpp`

## 1. isoform을 intron chain으로 바꾸면 무엇이 남을까?

`iso_A`는 caller GTF에 exon 다섯 줄로 적혀 있다.

```bash
PIG=../build/panisoguard; mkdir -p out
grep -P '\texon\t' data/caller1.gtf | grep '"iso_A"'
```

```text
chrT	caller1	exon	101	200	.	+	.	gene_id "GENE_T"; transcript_id "iso_A";
chrT	caller1	exon	401	500	.	+	.	gene_id "GENE_T"; transcript_id "iso_A";
chrT	caller1	exon	731	800	.	+	.	gene_id "GENE_T"; transcript_id "iso_A";
chrT	caller1	exon	1001	1100	.	+	.	gene_id "GENE_T"; transcript_id "iso_A";
chrT	caller1	exon	1301	1400	.	+	.	gene_id "GENE_T"; transcript_id "iso_A";
```

PanIsoGuard는 이 exon 목록에서 exon 사이의 간격, 곧 intron만 뽑아 쓴다. 정렬된 intron 목록을 **intron chain**이라고 부른다. exon 다섯 개면 intron은 네 개다.

```python
def exons_of(path):
    tx = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            tid = f[8].split('transcript_id "')[1].split('"')[0]
            tx.setdefault(tid, (f[0], f[6], []))[2].append((int(f[3]), int(f[4])))
    return tx

def intron_chain(exons):
    """GTF exons (1-based inclusive) -> introns, 0-based half-open [start, end)."""
    e = sorted(exons)
    return [(e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)]

iso = exons_of("data/caller1.gtf")
for t in ("iso_known", "iso_A"):
    chrom, strand, ex = iso[t]
    print(f"{t:9s} exons {ex}")
    print(f"{'':9s} introns (0-based half-open) {intron_chain(ex)}")
```

```text
iso_known exons [(101, 200), (401, 500), (701, 800), (1001, 1100), (1301, 1400)]
          introns (0-based half-open) [(200, 400), (500, 700), (800, 1000), (1100, 1300)]
iso_A     exons [(101, 200), (401, 500), (731, 800), (1001, 1100), (1301, 1400)]
          introns (0-based half-open) [(200, 400), (500, 730), (800, 1000), (1100, 1300)]
```

두 isoform의 intron chain은 두 번째 intron 하나만 다르다. `iso_known`은 (500, 700), `iso_A`는 (500, 730)이다. 좌표를 읽는 법은 2절에서 설명한다.

intron chain으로 바꾸면 사라지는 정보가 두 가지 있다. 첫 exon이 어디서 시작하는지(TSS)와 마지막 exon이 어디서 끝나는지(TES)다. 이것은 일부러 버리는 것이다. long read는 5' 끝이 잘려 나가기 쉽고 poly-A 위치도 조금씩 흔들려서, 끝 좌표가 몇 bp 다르다고 다른 isoform으로 보면 같은 isoform이 여러 개로 쪼개진다. SQANTI3도 FSM을 판정할 때 양 끝을 보지 않는다. caller 여러 개의 결과를 합칠 때도 같은 이유로 intron chain만 비교한다([06](06_multi_caller_consensus.md)).

## 2. 세 파일은 같은 intron을 같은 좌표로 적을까?

같은 intron이라도 파일 형식마다 좌표를 적는 방식이 다르다.

| 형식 | 무엇을 적나 | 좌표 규칙 | `iso_B`의 novel intron |
|---|---|---|---|
| GTF | exon | 1-based, 양 끝 포함 | exon `701–800`과 `1021–1100` 사이 |
| STAR `SJ.out.tab` | intron | 1-based, 양 끝 포함 | `801 1020` |
| BED12 | exon(block) | 0-based, 끝 미포함(half-open), 시작점 기준 상대 좌표 | block 끝 800, 다음 block 시작 1020 |

PanIsoGuard는 읽는 순간 모든 intron을 **0-based half-open `[start, end)`**로 바꾼다. `start`는 intron의 첫 염기(0부터 센 위치), `end`는 intron 마지막 염기 바로 다음 위치다. 이렇게 하면 intron 길이가 `end - start`로 바로 나온다. [types.hpp](../include/panisoguard/types.hpp) 맨 위 주석에 형식별 변환 규칙이 정리되어 있다. 정말 세 형식이 같은 값이 되는지 `iso_B`의 novel intron으로 확인해 본다.

```python
# The same intron of iso_B, as written by three file formats, converted to 0-based half-open.
gtf_exons = [(701, 800), (1021, 1100)]                      # GTF: 1-based inclusive exons
from_gtf = (gtf_exons[0][1], gtf_exons[1][0] - 1)           # [e_i, s_{i+1} - 1)

sj = next(l.split("\t") for l in open("data/short_reads.SJ.out.tab") if l.split("\t")[1:3] == ["801", "1020"])
from_sj = (int(sj[1]) - 1, int(sj[2]))                      # SJ.tab: 1-based inclusive intron -> [s-1, e)

bed_start = 700                                             # BED: 0-based chromStart of the first block
sizes, starts = [100, 80], [0, 320]                         # the two blocks around the intron
from_bed = (bed_start + starts[0] + sizes[0], bed_start + starts[1])   # [block_i end, block_{i+1} start)

print("SJ.tab row:", " ".join(x.strip() for x in sj))
print("GTF ->", from_gtf, "| SJ.tab ->", from_sj, "| BED12 ->", from_bed)
print("1-based inclusive:", (from_gtf[0] + 1, from_gtf[1]), "| length", from_gtf[1] - from_gtf[0], "bp")
```

```text
SJ.tab row: chrT 801 1020 1 1 0 12 0 48
GTF -> (800, 1020) | SJ.tab -> (800, 1020) | BED12 -> (800, 1020)
1-based inclusive: (801, 1020) | length 220 bp
```

세 형식 모두 `(800, 1020)`이 된다. GTF 쪽 식이 조금 헷갈리는데, 풀어 보면 이렇다. exon이 800에서 끝나면(1-based) intron은 801에서 시작하고, 0-based로는 800이다. 그래서 exon 끝 좌표가 그대로 intron의 0-based 시작점이 된다. 다음 exon이 1021에서 시작하면 intron의 마지막 염기는 1020(1-based)이고, half-open 끝은 "마지막 염기의 0-based 위치 + 1"이라서 다시 1020이다.

좌표를 이렇게 하나로 맞춰 두는 이유는, PanIsoGuard가 junction을 비교할 때 **좌표가 정확히 같아야** 같은 것으로 보기 때문이다. 1 bp라도 어긋나면 다른 junction이다(연습문제 1). 실제로 caller 결과를 BED12로 줘도 판정이 똑같이 나오는지 확인해 본다.

```python
# Write caller 1 as BED12 (0-based chromStart, blocks relative to it).
tx = {}
for line in open("data/caller1.gtf"):
    f = line.rstrip("\n").split("\t")
    if f[2] == "exon":
        t = f[8].split('transcript_id "')[1].split('"')[0]
        tx.setdefault(t, (f[0], f[6], []))[2].append((int(f[3]) - 1, int(f[4])))   # -> 0-based half-open
with open("out/caller1.bed", "w") as out:
    for t, (chrom, strand, ex) in tx.items():
        ex.sort()
        s, e = ex[0][0], ex[-1][1]
        sizes = ",".join(str(b - a) for a, b in ex)
        starts = ",".join(str(a - s) for a, b in ex)
        out.write(f"{chrom}\t{s}\t{e}\t{t}\t0\t{strand}\t{s}\t{e}\t0\t{len(ex)}\t{sizes}\t{starts}\n")
print(open("out/caller1.bed").read().splitlines()[2])
```

```text
chrT	100	1400	iso_A	0	+	100	1400	0	5	100,100,70,100,100	0,300,630,900,1200
```

```bash
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-bed out/caller1.bed \
  --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab --bam data/long_reads.bam \
  --out-prefix out/toy_bed 2>/dev/null
diff out/toy.adjudicated.tsv out/toy_bed.adjudicated.tsv && echo "GTF and BED12 input give identical verdicts"
```

```text
GTF and BED12 input give identical verdicts
```

`out/toy.adjudicated.tsv`는 [00](00_overview.md) 3절에서 같은 입력을 GTF로 준 결과다. BED12 쪽 `iso_A` 줄을 보면 세 번째 block의 크기가 70이다. exon 3이 731–800으로 짧아진 것이 여기 보인다.

## 3. 무엇과 비교해서 "새롭다"고 할까?

PanIsoGuard는 참조 GTF(`--ref-gtf`)를 읽어서 **catalog**을 만든다. 참조 주석의 모든 전사체에서 intron을 뽑아 `(염색체, strand, start, end)` 네 값의 집합으로 모아 둔 것이다. 그다음 isoform의 intron을 하나씩 catalog에서 찾아보고, 없으면 novel junction으로 센다(`build_evidence()`의 `catalog->has_intron()`). 전사체 단위가 아니라 intron 하나하나를 본다는 점이 중요하다. 손으로 똑같이 세어서 PanIsoGuard 결과(00의 step1)와 비교해 본다.

```python
import csv

def chains(path):
    tx = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            t = f[8].split('transcript_id "')[1].split('"')[0]
            tx.setdefault(t, [f[0], f[6], []])[2].append((int(f[3]), int(f[4])))
    out = {}
    for t, (chrom, strand, ex) in tx.items():
        ex.sort()
        out[t] = (chrom, strand, [(ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1)])
    return out

catalog = {(c, s, j) for c, s, js in chains("data/reference.gtf").values() for j in js}
print("catalog:", len(catalog), "known introns ->", sorted(j for _, _, j in catalog))

pig = {r["isoform_id"]: r["n_novel_junctions"] for r in csv.DictReader(open("out/step1.adjudicated.tsv"), delimiter="\t")}
for t, (c, s, js) in chains("data/caller1.gtf").items():
    novel = [j for j in js if (c, s, j) not in catalog]
    print(f"{t:9s} novel introns {str(novel):24s} by hand {len(novel)} | PanIsoGuard {pig[t]}")
```

```text
catalog: 6 known introns -> [(200, 400), (200, 700), (500, 700), (800, 1000), (800, 1300), (1100, 1300)]
iso_known novel introns []                       by hand 0 | PanIsoGuard 0
iso_ism   novel introns []                       by hand 0 | PanIsoGuard 0
iso_A     novel introns [(500, 730)]             by hand 1 | PanIsoGuard 1
iso_B     novel introns [(800, 1020)]            by hand 1 | PanIsoGuard 1
iso_C     novel introns [(800, 860), (940, 1000)] by hand 2 | PanIsoGuard 2
iso_D     novel introns []                       by hand 0 | PanIsoGuard 0
iso_E     novel introns [(500, 1000)]            by hand 1 | PanIsoGuard 1
iso_F     novel introns [(200, 450)]             by hand 1 | PanIsoGuard 1
iso_G     novel introns [(800, 1004)]            by hand 1 | PanIsoGuard 1
```

아홉 개 모두 손계산과 같다. catalog에는 참조 전사체 세 개의 intron 여섯 개가 들어 있다. TX1의 네 개에 TX2의 (200, 700), TX3의 (800, 1300)이 더해졌다.

catalog의 열쇠에는 strand가 들어 있다. 같은 좌표라도 반대 strand의 intron은 다른 intron이다. 한 유전자 자리에 반대 방향 유전자가 겹쳐 있을 때(antisense) 서로의 intron을 "이미 알려진 것"으로 잘못 세지 않으려는 것이다. 참조 GTF의 strand를 전부 `-`로 바꿔서 돌려 보면 효과가 바로 보인다.

```bash
awk 'BEGIN{FS=OFS="\t"} {$7="-"; print}' data/reference.gtf > out/reference_minus.gtf
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
  --ref-gtf out/reference_minus.gtf --out-prefix out/minus 2>/dev/null
paste <(cut -f1,4,8 out/step1.adjudicated.tsv) <(cut -f8,7 out/minus.adjudicated.tsv) | column -t
```

```text
isoform_id  structural_category      n_novel_junctions  confidence_class  n_novel_junctions
iso_A       novel_not_in_catalog     1                  AMBIGUOUS         4
iso_B       novel_not_in_catalog     1                  AMBIGUOUS         4
iso_C       novel_not_in_catalog     2                  AMBIGUOUS         5
iso_D       novel_in_catalog         0                  AMBIGUOUS         2
iso_E       novel_in_catalog         1                  AMBIGUOUS         3
iso_F       novel_not_in_catalog     1                  AMBIGUOUS         4
iso_G       novel_not_in_catalog     1                  AMBIGUOUS         4
iso_ism     incomplete-splice_match  0                  LOW_CONF_PARTIAL  3
iso_known   full-splice_match        0                  HIGH_CONF_KNOWN   4
```

셋째 열이 원래 참조로 센 값, 마지막 열이 strand를 뒤집은 참조로 센 값이다. 모든 intron이 novel이 되었다. 그런데 `iso_known`은 intron 네 개가 전부 novel로 세어졌는데도 여전히 `HIGH_CONF_KNOWN`이다. `structural_category`는 SQANTI3가 원래 참조로 미리 계산해 둔 값이고, PanIsoGuard는 FSM이면 novel junction 수를 보지 않고 바로 `HIGH_CONF_KNOWN`을 주기 때문이다.

이 실험은 실제로 조심해야 할 점을 하나 보여 준다. PanIsoGuard는 **범주는 SQANTI3에게서 받고, novel junction은 직접 센다.** 그러니 `--ref-gtf`에는 SQANTI3를 돌릴 때 쓴 것과 같은 주석을 줘야 한다. 두 주석이 다르면 "SQANTI3는 NNC라는데 PanIsoGuard는 novel junction이 0개"처럼 앞뒤가 맞지 않는 증거가 조용히 만들어진다. 경고는 나오지 않는다.

## 4. SQANTI3의 범주와 novel junction 수는 같은 말일까?

같은 말이 아니다. 이 예제에서 둘을 나란히 놓으면 이렇다.

| isoform | SQANTI3 범주 / subcategory | novel junction 수 | 무엇이 새로운가 |
|---|---|---:|---|
| `iso_A`, `iso_B`, `iso_F`, `iso_G` | NNC / at_least_one_novel_splicesite | 1 | 새 splice site를 쓰는 intron 하나 |
| `iso_C` | NNC / at_least_one_novel_splicesite | 2 | 새 exon 하나가 intron 두 개를 새로 만듦 |
| `iso_E` | NIC / combination_of_known_splicesites | 1 | 아는 donor(501)와 아는 acceptor(1000)를 새로 이은 intron |
| `iso_D` | NIC / combination_of_known_junctions | 0 | intron은 모두 알려져 있고 조합만 새로움 |
| `iso_ism` | ISM / 3prime_fragment | 0 | TX1의 앞부분이 잘린 조각 |

NIC라고 해서 novel junction이 없는 것은 아니다. `iso_E`처럼 알려진 splice site 두 개를 새로 이으면 intron 자체는 catalog에 없으니 novel junction이 된다. 이런 isoform은 short read로 확인할 대상이 있다. 반면 `iso_D`처럼 알려진 intron만 새로 조합한 isoform은 확인할 novel junction이 없어서 `AMBIGUOUS`로 보류된다([00](00_overview.md) 연습문제 2).

실제 데이터에서는 이런 isoform이 얼마나 될까? SQANTI-SIM으로 만든 chr22 시뮬레이션(정답을 아는 데이터, [07](07_evaluation.md))을 PanIsoGuard 0.0.4로 다시 돌려서 NIC만 subcategory별로 세어 봤다. 이 블록은 [benchmark/sqanti_sim/run.sh](../benchmark/sqanti_sim/run.sh)가 만든 작업 폴더가 있어야 돌아간다. 폴더가 약 260 MB라 저장소에는 넣지 않았다.

```bash
W=/mnt/Data/pig_sqantisim/work        # benchmark/sqanti_sim/run.sh의 작업 폴더
$PIG adjudicate --classification $W/sqanti_out/flair_classification.txt \
  --isoforms-gtf $W/flair.isoforms.gtf --ref-gtf $W/chr22_modified.gtf --sj-tab $W/truth.SJ.tab \
  --out-prefix out/sqsim 2>/dev/null
```

```python
import csv, collections
W = "/mnt/Data/pig_sqantisim/work"
cls = {r["isoform"]: r for r in csv.DictReader(open(f"{W}/sqanti_out/flair_classification.txt"), delimiter="\t")}
tab = collections.Counter()
for r in csv.DictReader(open("out/sqsim.adjudicated.tsv"), delimiter="\t"):
    c = cls[r["isoform_id"]]
    if c["structural_category"] == "novel_in_catalog":
        has = "novel jx >= 1" if int(r["n_novel_junctions"]) else "novel jx = 0"
        tab[(c["subcategory"], has, r["confidence_class"])] += 1
for k, v in sorted(tab.items()):
    print(f"{k[0]:34s} {k[1]:14s} {k[2]:18s} {v:4d}")
print("NIC total:", sum(tab.values()))
```

```text
combination_of_known_junctions     novel jx = 0   AMBIGUOUS           124
combination_of_known_splicesites   novel jx = 0   AMBIGUOUS             1
combination_of_known_splicesites   novel jx >= 1  HIGH_CONF_NOVEL     119
combination_of_known_splicesites   novel jx >= 1  LOW_CONF_PARTIAL      3
combination_of_known_splicesites   novel jx >= 1  MEDIUM_CONF_NOVEL    12
intron_retention                   novel jx = 0   AMBIGUOUS            24
intron_retention                   novel jx >= 1  HIGH_CONF_NOVEL       3
intron_retention                   novel jx >= 1  MEDIUM_CONF_NOVEL     1
NIC total: 287
```

NIC 287개 가운데 149개(124 + 1 + 24)가 novel junction이 없어서 `AMBIGUOUS`로 보류된다. 124개는 `iso_D`처럼 알려진 junction의 조합이고, 24개는 intron retention(intron 하나가 exon 안에 그대로 남은 것, 곧 intron이 하나 사라진 것)이다. intron retention은 새 intron을 만들지 않으니 novel junction이 없는 것이 당연하다. 반대로 novel junction이 있는 NIC 138개는 대부분(119개) `HIGH_CONF_NOVEL`이 된다.

여기서 [07](07_evaluation.md)에 나올 숫자 하나가 설명된다. SQANTI-SIM에서 진짜 NIC의 recall이 0.477(279개 중 133개)밖에 안 되는 이유가 대부분 이 149개다. PanIsoGuard의 한계라기보다 설계상 선택이다. short read 하나로는 멀리 떨어진 두 junction이 같은 분자에 함께 있는지 확인할 수 없으니, 확인할 수 없는 것은 판정하지 않고 보류한다. 다만 이 선택 때문에 recall을 잃는다는 사실은 그대로 적어 둬야 한다.

그런데 표에 이상한 줄이 하나 있다. `combination_of_known_splicesites`, 곧 SQANTI3가 "알려진 splice site로 새 intron을 만들었다"고 분류했는데 PanIsoGuard는 novel junction을 하나도 찾지 못한 isoform이다. 두 intron이 참조 주석의 어느 유전자에 들어 있는지 찾아봤다.

```python
W = "/mnt/Data/pig_sqantisim/work"
iso = "ENST00000756367.1_PBSIM_simulated_read_3"   # the one NIC "known splice sites" with no novel intron
ex = sorted((int(f[3]), int(f[4])) for f in (l.split("\t") for l in open(f"{W}/flair.isoforms.gtf"))
            if len(f) > 8 and f[2] == "exon" and f'"{iso}"' in f[8])
introns = [(ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1)]     # 1-based inclusive
seen = {}
for l in open(f"{W}/chr22_modified.gtf"):
    f = l.split("\t")
    if len(f) > 8 and f[2] == "exon" and f[6] == "-":
        t = f[8].split('transcript_id "')[1].split('"')[0]
        g = f[8].split('gene_id "')[1].split('"')[0]
        seen.setdefault((t, g), []).append((int(f[3]), int(f[4])))
for a, b in introns:
    hits = sorted({g for (t, g), e in seen.items() for i, (x, y) in enumerate(sorted(e)[:-1])
                   if (y + 1, sorted(e)[i + 1][0] - 1) == (a, b)})
    print(f"intron {a}-{b} (- strand) is annotated in gene(s): {hits}")
```

```text
intron 37340970-37342600 (- strand) is annotated in gene(s): ['ENSG00000237862.4', 'ENSG00000243902.7']
intron 37342704-37342976 (- strand) is annotated in gene(s): ['ENSG00000237862.4']
```

```bash
grep -P "^ENST00000756367.1_PBSIM_simulated_read_3\t" /mnt/Data/pig_sqantisim/work/sqanti_out/flair_classification.txt | cut -f8,9,11
```

```text
novel_in_catalog	combination_of_known_splicesites	ENSG00000243902.7
```

SQANTI3는 이 isoform을 유전자 ENSG00000243902.7(ELFN2)에 연결했다. 그런데 두 번째 intron은 옆 유전자 ENSG00000237862.4의 전사체에만 있다. SQANTI3는 연결한 유전자 안에서만 junction을 찾기 때문에 이 intron을 새것으로 본 것으로 보인다. PanIsoGuard의 catalog은 유전자를 구분하지 않고 좌표와 strand만 보기 때문에 알려진 intron으로 본다. 코드 동작으로 보면 PanIsoGuard 쪽은 설명한 그대로다. SQANTI3 쪽은 결과에서 거꾸로 추정한 것이고 SQANTI3 소스로 확인하지는 않았다.

## 정리

- PanIsoGuard는 isoform을 intron chain(정렬된 intron 목록)으로 바꿔서 봄. 첫 exon의 시작과 마지막 exon의 끝은 일부러 버림.
- GTF(1-based 양 끝 포함 exon), `SJ.out.tab`(1-based 양 끝 포함 intron), BED12(0-based half-open block)의 intron은 읽는 순간 모두 0-based half-open `[start, end)`로 바뀜. `iso_B`의 novel intron은 세 형식 모두 `(800, 1020)`임.
- novel junction은 참조 GTF로 만든 catalog(염색체, strand, start, end의 집합)에 없는 intron임. 전사체가 아니라 intron 하나하나를 비교하고, strand가 다르면 다른 intron임.
- 구조 범주(FSM, ISM, NIC, NNC)는 SQANTI3에게서 받고, novel junction 수는 PanIsoGuard가 직접 셈. `--ref-gtf`는 SQANTI3에 준 것과 같은 주석이어야 함.
- NIC라도 novel junction이 있을 수 있음(`iso_E`). 알려진 junction만 조합한 NIC(`iso_D`)와 intron retention은 확인할 대상이 없어서 보류됨. SQANTI-SIM에서는 NIC 287개 중 149개가 이렇게 보류됨.

다음 노트 [02](02_short_read_support.md)에서는 찾아낸 novel junction을 short read로 어떻게 확인하는지 본다. `iso_A`의 novel junction을 short read 9개가 걸쳤는데도 왜 지지를 못 받는지가 거기서 나온다.

## 연습문제

### 문제 1

> 누군가 `SJ.out.tab`의 intron 시작 좌표를 0-based로 착각해서 1씩 빼서 적었다. PanIsoGuard는 오류를 낼까? 판정은 어떻게 될까?

<details>
<summary>풀이</summary>

오류는 나지 않고 판정만 조용히 바뀐다. 둘째 열에서 1을 뺀 파일로 돌려서 원래 결과(00의 step2)와 비교했다.

```bash
awk 'BEGIN{FS=OFS="\t"} {$2=$2-1; print}' data/short_reads.SJ.out.tab > out/shifted.SJ.tab
head -2 out/shifted.SJ.tab
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
  --ref-gtf data/reference.gtf --sj-tab out/shifted.SJ.tab --out-prefix out/shifted 2>/dev/null
paste <(cut -f1,5,7 out/step2.adjudicated.tsv) <(cut -f5,7 out/shifted.adjudicated.tsv) | grep -E 'iso_(B|C|E)|isoform' | column -t
```

```text
chrT	200	400	1	1	1	30	0	50
chrT	200	700	1	1	1	10	0	48
isoform_id  novelty_support  confidence_class  novelty_support  confidence_class
iso_B       SUPPORTED        HIGH_CONF_NOVEL   UNSUPPORTED      LOW_CONF_PARTIAL
iso_C       PARTIAL          LOW_CONF_PARTIAL  UNSUPPORTED      LOW_CONF_PARTIAL
iso_E       SUPPORTED        HIGH_CONF_NOVEL   UNSUPPORTED      LOW_CONF_PARTIAL
```

1 bp 어긋난 junction은 catalog에서도, short read 표에서도 다른 junction이다. 그래서 short read 지지가 전부 사라지고 `HIGH_CONF_NOVEL` 두 개가 `LOW_CONF_PARTIAL`로 떨어진다. `SJ.out.tab`은 여러 도구가 직접 만들거나 고쳐 쓰는 경우가 많으니, 결과에서 `SUPPORTED`가 이상하게 적으면 좌표 규칙부터 의심해 볼 만하다. 이런 실수를 알아채는 가장 쉬운 방법은 참조 주석에 있는 intron이 `SJ.out.tab`에서 얼마나 발견되는지 세어 보는 것이다. 정상이라면 발현되는 유전자의 알려진 intron은 대부분 있어야 한다.

</details>

### 문제 2

> caller 결과를 BED12로 줬는데, caller가 isoform 이름 앞에 도구 이름을 붙여서 `iso_A`가 `FLAIR_iso_A`로 적혀 있다. SQANTI3 분류표에는 `iso_A`로 되어 있다. 어떻게 될까?

<details>
<summary>풀이</summary>

PanIsoGuard는 SQANTI3 분류표의 `isoform` 열 값과 caller 파일의 이름이 **정확히 같은 것**끼리 짝을 짓는다. 짝이 없으면 intron chain이 없으니 novel junction을 셀 수 없다.

```bash
sed 's/\tiso_A\t/\tFLAIR_iso_A\t/' out/caller1.bed > out/renamed.bed
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-bed out/renamed.bed \
  --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab --out-prefix out/renamed 2>/dev/null
grep '"iso_A"' out/renamed.attribution.jsonl | python3 -c 'import json,sys; r=json.loads(sys.stdin.read()); print(r["confidence_class"], r["evidence"]["chain_available"], r["rule_trace"])'
```

```text
AMBIGUOUS False ['novelty-support UNKNOWN reason=caller_chain_absent', 'all_canonical=non_canonical -> mechanism=noncanonical', 'project(UNKNOWN,noncanonical) -> AMBIGUOUS']
```

`iso_A`는 `ARTIFACT`(step2) 대신 `AMBIGUOUS`가 되고, 이유는 `caller_chain_absent`로 남는다. 이번에도 경고는 없다. stderr에는 "read 9 caller isoform chains"만 찍혀서, 이름이 하나 어긋났다는 것을 거기서는 알 수 없다. 결과에 `caller_chain_absent`가 보이면 이름 짝짓기부터 확인해야 한다.

</details>

## 더 깊이 보기

<details>
<summary>catalog에는 intron 말고 무엇이 더 들어 있나</summary>

`Catalog` 클래스([gtf.hpp](../include/panisoguard/gtf.hpp))에는 intron 집합 말고도 splice site 집합(`has_site`)과 intron chain fingerprint 집합(`has_chain`)이 있다. 2026년 9월 기준 코드에서 `adjudicate`가 쓰는 것은 intron 집합뿐이다. fingerprint 집합은 `combine`이 "이 chain이 참조 전사체와 똑같은가"를 표시할 때만 쓴다([06](06_multi_caller_consensus.md)). splice site 집합은 헤더 주석에 NIC와 NNC를 가르는 데 쓴다고 적혀 있지만, 그 판단은 SQANTI3 범주를 그대로 받아 쓰기 때문에 실제로 호출하는 곳은 없다(`grep -rn has_site src`로 확인).

</details>

<details>
<summary>GTF의 strand는 어느 줄에서 읽나</summary>

GTF reader는 한 전사체의 exon 줄 가운데 **처음 나온 줄**의 strand를 그 전사체의 strand로 쓴다(`read_gtf_transcripts()`). 뒤쪽 exon 줄의 strand가 다르게 적혀 있어도 확인하지 않고 오류도 내지 않는다. strand가 `.`로 적힌 전사체는 `+`나 `-`인 catalog intron과 절대 같아지지 않으니, 모든 intron이 novel로 세어진다.

</details>
