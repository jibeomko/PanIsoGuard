# 03. artifact는 어떤 흔적을 남길까?

[02](02_short_read_support.md)의 질문은 "누가 이 junction을 확인해 주나"였다. 이 노트의 질문은 반대쪽이다. 이 isoform이 가짜라면 어떤 흔적이 남아 있을까? PanIsoGuard는 두 곳에서 흔적을 찾는다. 하나는 SQANTI3가 이미 계산해서 분류표에 적어 둔 QC 값이고, 다른 하나는 long-read 정렬 파일(BAM)이다. 흔적이 여러 개면 하나를 골라 `primary_mechanism`으로 보고하는데, 고르는 순서도 이 노트에서 확인한다.

> 관련 문서: [docs/decision_engine.md](../docs/decision_engine.md) "Mechanism priority and its rationale", [benchmark/bam_axis](../benchmark/bam_axis), [benchmark/gm12878_realdata](../benchmark/gm12878_realdata) · 코드: `RuleEngine::evaluate()`(`src/core/rules.cpp`), `src/evidence/bam_features.cpp`

## 1. SQANTI3는 어떤 흔적을 적어 줄까?

SQANTI3 분류표의 열 가운데 PanIsoGuard가 기전을 정하는 데 쓰는 것은 세 개다. 예제의 값부터 본다.

```bash
PIG=../build/panisoguard; mkdir -p out
cut -f1,8,21,22,45 data/caller1_classification.txt | column -t
```

```text
isoform    structural_category      RTS_stage  all_canonical  perc_A_downstream_TTS
iso_A      novel_not_in_catalog     FALSE      non_canonical  45.0
iso_B      novel_not_in_catalog     FALSE      canonical      45.0
iso_C      novel_not_in_catalog     FALSE      canonical      45.0
iso_D      novel_in_catalog         FALSE      canonical      45.0
iso_E      novel_in_catalog         FALSE      canonical      45.0
iso_F      novel_not_in_catalog     FALSE      non_canonical  45.0
iso_G      novel_not_in_catalog     FALSE      canonical      45.0
iso_ism    incomplete-splice_match  FALSE      canonical      45.0
iso_known  full-splice_match        FALSE      canonical      45.0
```

| 열 | SQANTI3가 계산하는 것 | PanIsoGuard의 기전 |
|---|---|---|
| `all_canonical` | isoform의 모든 junction이 canonical motif(GT-AG, GC-AG, AT-AC)인지 | `non_canonical`이면 `noncanonical` |
| `RTS_stage` | 역전사 효소가 반복 서열 사이를 건너뛰어 가짜 junction을 만들었을 가능성(RT switching) | `TRUE`면 `rt_switch` |
| `perc_A_downstream_TTS` | isoform 끝 바로 뒤 게놈 20 bp 중 A의 비율(%) | 60 이상이면 `degradation` |

세 번째 값은 설명이 조금 필요하다. poly-A 꼬리에 붙는 oligo-dT primer는 mRNA 중간의 A가 많은 구간에도 붙을 수 있다(intra-priming). 그러면 3' 끝이 실제보다 앞에서 잘린 가짜 isoform이 생긴다. isoform이 끝난 바로 뒤 게놈에 A가 많으면 이런 일을 의심한다. 기준 60%는 SQANTI3 기본 필터와 같다. 예제는 모든 isoform이 1400에서 끝나서 이 값이 45.0으로 같다.

PanIsoGuard는 이 세 값을 다시 계산하지 않고 그대로 받아 쓴다. `iso_A`와 `iso_F`의 `non_canonical`이 정말 게놈 서열과 맞는지만 확인해 봤다.

```python
import pysam
fa = pysam.FastaFile("data/toy.fa")
for name, (s, e) in {"iso_A": (500, 730), "iso_B": (800, 1020), "iso_F": (200, 450)}.items():
    donor, acceptor = fa.fetch("chrT", s, s + 2), fa.fetch("chrT", e - 2, e)   # first / last two intron bases
    print(f"{name} novel intron [{s}, {e}): {donor}...{acceptor}")
```

```text
iso_A novel intron [500, 730): GT...AC
iso_B novel intron [800, 1020): GT...AG
iso_F novel intron [200, 450): GT...TT
```

intron의 첫 두 염기(donor)와 마지막 두 염기(acceptor)다. `iso_B`는 GT…AG로 canonical이고, `iso_A`의 GT…AC와 `iso_F`의 GT…TT는 non-canonical이다. SQANTI3의 값과 맞는다.

그러면 이 세 흔적은 실제로 가짜를 얼마나 잘 가려낼까? SQANTI-SIM 시뮬레이션(정답을 아는 데이터, [07](07_evaluation.md))의 NIC와 NNC에서 정답별로 세어 봤다. 이 블록은 [benchmark/sqanti_sim/run.sh](../benchmark/sqanti_sim/run.sh)의 작업 폴더가 있어야 돌아간다.

```python
import csv, collections
W = "/mnt/Data/pig_sqantisim/work"
ex = collections.defaultdict(list); st = {}; ch = {}
for l in open(f"{W}/flair.isoforms.gtf"):
    f = l.rstrip("\n").split("\t")
    if len(f) > 8 and f[2] == "exon":
        t = f[8].split('transcript_id "')[1].split('"')[0]
        ex[t].append((int(f[3]), int(f[4]))); st[t] = f[6]; ch[t] = f[0]
def key(t):
    e = sorted(ex[t])
    intr = [(e[i][1] + 1, e[i + 1][0] - 1) for i in range(len(e) - 1)]
    return f"{ch[t]}|{st[t]}|" + ",".join(f"{a}-{b}" for a, b in intr) if intr else None
truth = dict(l.rstrip("\n").split("\t") for l in open(f"{W}/truth.truth.tsv"))
tab = collections.Counter(); n = collections.Counter()
for r in csv.DictReader(open(f"{W}/sqanti_out/flair_classification.txt"), delimiter="\t"):
    if r["structural_category"] not in ("novel_in_catalog", "novel_not_in_catalog"):
        continue
    k = key(r["isoform"]); lab = truth.get(k, "false_novel") if k else "monoexonic"
    if lab not in ("genuine_novel", "false_novel"):
        continue
    n[lab] += 1
    tab[(lab, "RTS_stage=TRUE")] += r["RTS_stage"] == "TRUE"
    tab[(lab, "non_canonical")] += r["all_canonical"] == "non_canonical"
    pa = r["perc_A_downstream_TTS"]
    tab[(lab, "perc_A>=60")] += pa not in ("NA", "") and float(pa) >= 60
print("NIC+NNC multi-exon:", dict(n))
for flag in ("non_canonical", "RTS_stage=TRUE", "perc_A>=60"):
    g, f = tab[("genuine_novel", flag)], tab[("false_novel", flag)]
    print(f"{flag:15s} genuine {g:3d}/{n['genuine_novel']} ({g / n['genuine_novel']:.1%}) | false {f:3d}/{n['false_novel']} ({f / n['false_novel']:.1%})")
```

```text
NIC+NNC multi-exon: {'false_novel': 124, 'genuine_novel': 576}
non_canonical   genuine   5/576 (0.9%) | false  65/124 (52.4%)
RTS_stage=TRUE  genuine  57/576 (9.9%) | false  25/124 (20.2%)
perc_A>=60      genuine  15/576 (2.6%) | false   4/124 (3.2%)
```

세 흔적의 쓸모가 크게 다르다. `non_canonical`은 가짜의 절반(52.4%)에 붙고 진짜에는 1%도 안 붙는다. 이 데이터에서 가장 강한 신호다. `RTS_stage`는 가짜에 두 배쯤 자주 붙지만 진짜에도 10%나 붙는다. `perc_A`는 진짜와 가짜에서 거의 차이가 없다. 이 시뮬레이션(PBSIM3)은 전사체 서열에서 read를 뽑을 뿐 intra-priming을 흉내 내지 않으니, 여기서 `perc_A`가 쓸모없게 나오는 것은 당연하다. 이 데이터로는 degradation 기전을 평가할 수 없다는 뜻이다.

## 2. long-read 정렬에서는 무엇을 읽을까?

`--bam`을 주면 PanIsoGuard는 novel junction마다 그 junction을 걸친 long read를 모아서 네 가지를 센다. "걸쳤다"는 read의 CIGAR에 **이 junction과 좌표가 정확히 같은 `N`(intron 건너뛰기) 연산**이 있다는 뜻이다. 근처에서 조금 다른 intron으로 정렬된 read는 세지 않는다.

| 특징 | read 하나가 해당되는 조건 | 무엇을 의심하나 |
|---|---|---|
| low MAPQ | MAPQ < 20 | 반복 서열이나 paralog라 정렬 위치가 불확실함 |
| supplementary | FLAG에 secondary(256)나 supplementary(2048)가 있음 | chimeric read, 한 read가 여러 조각으로 정렬됨 |
| indel near | junction 양 끝 ±10 bp 안에 삽입(`I`)이나 결실(`D`) | 작은 indel을 aligner가 intron으로 잘못 그림 |
| soft-clip | read 양 끝 중 한 곳에 20 bp 이상 soft-clip(`S`) | read 끝이 정렬되지 않음 |

저장소의 테스트용 BAM(`tests/data/tiny/mini.bam`)에는 이 네 경우를 하나씩 담은 read가 들어 있다.

```bash
samtools view ../tests/data/tiny/mini.bam | cut -f1-6
```

```text
r1_clean	0	chr1	1901	60	100M1000N100M
r2_lowmapq	0	chr1	1901	5	100M1000N100M
r3_supp	2048	chr1	1901	60	100M1000N100M
r4_softclip	0	chr1	1901	60	30S100M1000N100M
r5_indel	0	chr1	1901	60	100M1000N5I95M
r6_nonspan	0	chr1	5001	60	200M
```

앞의 다섯 read는 모두 chr1:2001–3000(0-based `[2000, 3000)`)의 intron을 걸친다. 규칙을 Python으로 옮겨서 read마다 어디에 해당하는지 따져 본다(`BamReader::features_at_junction()`과 같은 순서로 CIGAR를 읽는다).

```python
import pysam

JX = (2000, 3000)            # the intron to test, 0-based half-open (chr1:2001-3000)
W, MIN_MAPQ, MIN_CLIP = 10, 20, 20

def near(pos):
    return abs(pos - JX[0]) <= W or abs(pos - JX[1]) <= W

for r in pysam.AlignmentFile("../tests/data/tiny/mini.bam"):
    pos, spanning, clip, indel = r.reference_start, False, False, False
    for i, (op, n) in enumerate(r.cigartuples):
        if op == 4 and n >= MIN_CLIP and i in (0, len(r.cigartuples) - 1):  # S at either end
            clip = True
        elif op == 3:                                                        # N
            spanning |= (pos, pos + n) == JX
            pos += n
        elif op == 2:                                                        # D
            indel |= near(pos)
            pos += n
        elif op == 1:                                                        # I
            indel |= near(pos)
        elif op in (0, 7, 8):                                                # M = X
            pos += n
    flags = [name for name, hit in [("low_mapq", r.mapping_quality < MIN_MAPQ),
                                    ("supplementary", r.is_supplementary or r.is_secondary),
                                    ("softclip", clip), ("indel_near", indel)] if hit]
    print(f"{r.query_name:12s} {r.cigarstring:18s} MAPQ {r.mapping_quality:2d} spanning={spanning!s:5s} {flags}")
```

```text
r1_clean     100M1000N100M      MAPQ 60 spanning=True  []
r2_lowmapq   100M1000N100M      MAPQ  5 spanning=True  ['low_mapq']
r3_supp      100M1000N100M      MAPQ 60 spanning=True  ['supplementary']
r4_softclip  30S100M1000N100M   MAPQ 60 spanning=True  ['softclip']
r5_indel     100M1000N5I95M     MAPQ 60 spanning=True  ['indel_near']
r6_nonspan   200M               MAPQ 60 spanning=False []
```

`r5_indel`의 5 bp 삽입은 `N` 바로 뒤, 곧 intron 끝(3000)에서 일어났으니 ±10 bp 안이다. `r6_nonspan`은 intron을 걸치지 않아서 어떤 특징이 있어도 세지 않는다. 같은 BAM을 PanIsoGuard에 넣어서 확인한다. 이 intron 하나만 가진 isoform을 위한 GTF와 분류표를 즉석에서 만들었다. 참조 주석에는 다른 intron을 두어 이 intron이 novel로 세어지게 했다.

```bash
printf 'chr1\tt\texon\t1901\t2000\t.\t+\t.\ttranscript_id "nov1";\nchr1\tt\texon\t3001\t3100\t.\t+\t.\ttranscript_id "nov1";\n' > out/mini_iso.gtf
printf 'chr1\tt\texon\t1901\t2000\t.\t+\t.\ttranscript_id "ref1";\nchr1\tt\texon\t2501\t3100\t.\t+\t.\ttranscript_id "ref1";\n' > out/mini_ref.gtf
printf 'isoform\tchrom\tstrand\tstructural_category\tRTS_stage\tall_canonical\tperc_A_downstream_TTS\nnov1\tchr1\t+\tnovel_not_in_catalog\tFALSE\tcanonical\t10\n' > out/mini_class.txt
$PIG adjudicate --classification out/mini_class.txt --isoforms-gtf out/mini_iso.gtf --ref-gtf out/mini_ref.gtf \
  --bam ../tests/data/tiny/mini.bam --out-prefix out/mini 2>/dev/null
python3 -c '
import json
r = json.loads(open("out/mini.attribution.jsonl").readline())
print({k: v for k, v in r["evidence"].items() if k.startswith("bam")})
print(r["confidence_class"], r["rule_trace"])'
```

```text
{'bam_evaluable': True, 'bam_n_spanning': 5, 'bam_frac_low_mapq': 0.2, 'bam_frac_supplementary': 0.2, 'bam_frac_softclip': 0.2, 'bam_frac_indel_near': 0.2}
AMBIGUOUS ['novelty-support UNKNOWN reason=short_read/catalog_axis_absent', 'no artifact mechanism flagged -> mechanism=none', 'project(UNKNOWN,none) -> AMBIGUOUS']
```

걸친 read는 5개이고, 네 특징이 각각 5개 중 1개(0.2)씩이다. 분모가 걸친 read 수라는 점이 중요하다. 네 비율 모두 기준을 넘지 않아서 기전은 `none`이다. (`SJ.out.tab`을 주지 않았으니 지지 수준은 `UNKNOWN`이고, 그래서 `AMBIGUOUS`다.)

## 3. 흔적이 얼마나 많아야 기전으로 인정할까?

비율 하나가 기준을 **넘으면**(같으면 안 됨) `mapping_or_repeat` 기전이 된다. 기준은 [config/rules.default.toml](../config/rules.default.toml)의 `[axis_artifact.mapping]`에 있다.

| 특징 | 기준 | 기본 상태 |
|---|---|---|
| low MAPQ 비율 | > 0.5 | 켜짐 |
| supplementary 비율 | > 0.5 | 켜짐 |
| indel near 비율 | > 0.5 | 켜짐 (0.0.4부터) |
| soft-clip 비율 | > 1.01 | **꺼짐** (비율은 1을 넘을 수 없으니 절대 걸리지 않음) |

예제의 `iso_G`가 이 기전에 걸린 isoform이다. `iso_G`의 novel junction은 801–1004로, 알려진 intron 801–1000보다 acceptor가 4 bp 뒤에 있다. 이 junction을 걸친 long read 4개 중 3개에 junction 5 bp 뒤 3 bp 결실이 있다.

```bash
samtools view data/long_reads.bam | grep '^iso_G' | cut -f1-6
python3 -c '
import json
for line in open("out/toy.attribution.jsonl"):
    r = json.loads(line)
    if r["isoform"] == "iso_G":
        print({k: v for k, v in r["evidence"].items() if k.startswith("bam")})
        print(r["rule_trace"])'
```

```text
iso_G.r1	0	chrT	101	60	100M200N100M200N100M204N5M3D88M200N100M
iso_G.r2	0	chrT	101	60	100M200N100M200N100M204N5M3D88M200N100M
iso_G.r3	0	chrT	101	60	100M200N100M200N100M204N5M3D88M200N100M
iso_G.r4	0	chrT	101	60	100M200N100M200N100M204N96M200N100M
{'bam_evaluable': True, 'bam_n_spanning': 4, 'bam_frac_low_mapq': 0, 'bam_frac_supplementary': 0, 'bam_frac_softclip': 0, 'bam_frac_indel_near': 0.75}
['sj_support=0/1 -> UNSUPPORTED', 'bam mapping artifact (low_mapq_frac=0.000000, supplementary_frac=0.000000, indel_near_frac=0.750000, softclip_frac=0.000000) -> mechanism=mapping_or_repeat', 'project(UNSUPPORTED,mapping_or_repeat) -> ARTIFACT']
```

`out/toy`는 [00](00_overview.md) 3절의 실행 결과다. indel near 비율 0.75가 기준 0.5를 넘어서 `mapping_or_repeat`이 되었다.

왜 junction 근처의 indel이 가짜 junction의 흔적일까? read에 junction 근처의 작은 indel(시퀀싱 오류이거나 실제 변이)이 있으면, aligner는 indel을 그대로 두는 것보다 intron 경계를 몇 bp 옮기는 쪽이 점수가 더 좋을 때가 있다. 그러면 알려진 junction 옆에 몇 bp 어긋난 "새" junction이 생기고, 그 옆에 남은 indel이 흔적으로 보인다. `iso_G`가 딱 그런 모양이다. SQANTI-SIM에서 이 신호는 진짜 novel junction에서 가장 높아야 0.19였고 가짜에서는 평균 0.35였다. 기준 0.5로 가짜 124개 중 41개를 잡았고 진짜는 하나도 잡지 않았다([benchmark/bam_axis](../benchmark/bam_axis)).

soft-clip만 기본으로 꺼 둔 이유는, 이 특징이 junction 근처가 아니라 **read 양 끝**의 soft-clip을 보기 때문이다. 실제 long read의 끝에는 adapter나 poly-A가 정렬되지 않고 남는 일이 흔하다. 그래서 junction 자체는 멀쩡한 read에서도 이 특징이 켜진다. 실제 GM12878 ONT direct RNA 데이터에 이 특징들을 돌려 본 결과가 저장소에 들어 있다.

```python
import json
m = json.load(open("../benchmark/results/gm12878_realdata/metrics.json"))["metrics"]
for k, (hit, n) in sorted(m["bam_signal_fire_rate_on_real_ONT"].items()):
    print(f"{k:22s} {hit:3d}/{n} = {hit / n:.0%}")
print("mean indel_near fraction:", m["indel_near_mean_by_consensus"])
```

```text
indel_near_gt_0.5      175/191 = 92%
low_mapq_gt_0.5         11/191 = 6%
softclip_gt_0.5        187/191 = 98%
supplementary_gt_0.5     4/191 = 2%
mean indel_near fraction: {'ge2_caller': 0.732, 'single_caller': 0.925}
```

soft-clip을 켰다면 novel junction 191개 중 187개(98%)가 `mapping_or_repeat`이 되었을 것이다. 문제는 soft-clip만이 아니다. indel near도 92%에서 켜진다. ONT read는 시퀀싱 오류 자체가 대부분 indel이라, 진짜 junction 근처에도 indel이 흔하다. caller 두 개 이상이 찾은(더 믿을 만한) junction도 평균 비율이 0.73이다. 기준 0.5는 오류가 적은 PacBio HiFi read(SQANTI-SIM은 HiFi를 흉내 냄)에 맞춘 값이라, ONT에는 그대로 쓸 수 없다. `--bam`은 선택 입력이고, 쓴다면 시퀀싱 화학에 맞게 기준을 다시 정해야 한다.

## 4. 흔적이 여러 개면 무엇을 보고할까?

한 isoform에 흔적이 여러 개일 수 있다. PanIsoGuard는 **mapping > noncanonical > rt_switch > degradation** 순서로 보고, 처음 걸리는 하나를 `primary_mechanism`으로 정한다. 정렬이 틀렸다면 그 정렬로 계산한 motif나 QC 값도 믿을 수 없으니, 가장 앞 단계의 문제를 먼저 보고한다는 논리다. 분류표를 일부러 고쳐서 확인해 본다. `iso_B`에는 intra-priming 흔적을, `iso_E`에는 RT switching 흔적을, `iso_G`에는 세 흔적을 모두 더했다.

```python
import csv
rows = list(csv.reader(open("data/caller1_classification.txt"), delimiter="\t"))
col = {name: i for i, name in enumerate(rows[0])}
edits = {"iso_B": {"perc_A_downstream_TTS": "70.0"},
         "iso_E": {"RTS_stage": "TRUE"},
         "iso_G": {"RTS_stage": "TRUE", "all_canonical": "non_canonical", "perc_A_downstream_TTS": "70.0"}}
for r in rows[1:]:
    for k, v in edits.get(r[0], {}).items():
        r[col[k]] = v
csv.writer(open("out/edited_classification.txt", "w"), delimiter="\t", lineterminator="\n").writerows(rows)
print("edited:", edits)
```

```text
edited: {'iso_B': {'perc_A_downstream_TTS': '70.0'}, 'iso_E': {'RTS_stage': 'TRUE'}, 'iso_G': {'RTS_stage': 'TRUE', 'all_canonical': 'non_canonical', 'perc_A_downstream_TTS': '70.0'}}
```

```bash
$PIG adjudicate --classification out/edited_classification.txt --isoforms-gtf data/caller1.gtf \
  --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab --bam data/long_reads.bam \
  --out-prefix out/edited 2>/dev/null
```

```python
import json
for line in open("out/edited.attribution.jsonl"):
    r = json.loads(line)
    if r["isoform"] in ("iso_B", "iso_E", "iso_G"):
        print(f'{r["isoform"]}: {r["novelty_support"]} x {r["primary_mechanism"]} -> {r["confidence_class"]}')
        for t in r["rule_trace"][1:]:
            print("    ", t)
```

```text
iso_B: SUPPORTED x degradation -> MEDIUM_CONF_NOVEL
     perc_A_downstream_TTS=70.000000 >= 60.000000 -> mechanism=degradation
     project(SUPPORTED,degradation) -> MEDIUM_CONF_NOVEL
iso_E: SUPPORTED x rt_switch -> MEDIUM_CONF_NOVEL
     RTS_stage=TRUE -> mechanism=rt_switch
     project(SUPPORTED,rt_switch) -> MEDIUM_CONF_NOVEL
iso_G: UNSUPPORTED x mapping_or_repeat -> ARTIFACT
     bam mapping artifact (low_mapq_frac=0.000000, supplementary_frac=0.000000, indel_near_frac=0.750000, softclip_frac=0.000000) -> mechanism=mapping_or_repeat
     project(UNSUPPORTED,mapping_or_repeat) -> ARTIFACT
```

`iso_B`와 `iso_E`는 short read로 확인된 isoform이라 흔적이 하나 생겨도 `MEDIUM_CONF_NOVEL`에 머문다. `iso_G`는 흔적이 넷인데 `mapping_or_repeat`만 보고된다.

여기서 문서와 다른 점을 하나 발견했다. [docs/decision_engine.md](../docs/decision_engine.md)는 "기전이 여러 개면 첫 번째만 primary이고, 걸린 조건은 모두 `rule_trace`에 나열된다"고 적고 있었다. 그런데 위 출력의 `iso_G` `rule_trace`에는 mapping 한 줄뿐이다. 코드(`RuleEngine::evaluate()`)가 기전을 `if … else if …`로 고르기 때문에, 처음 걸린 조건만 기록하고 나머지는 보지도 않는다. 나머지 흔적은 `rule_trace`가 아니라 JSON의 `evidence` 쪽(`noncanonical: true`, `rts_stage: true`)에 남는다. 문서는 이 노트를 쓰면서 코드에 맞게 고쳤다.

## 정리

- 기전은 두 곳에서 옴. SQANTI3 분류표의 `all_canonical`, `RTS_stage`, `perc_A_downstream_TTS`(다시 계산하지 않음)와, `--bam`으로 준 long read에서 센 네 가지 비율임.
- BAM 비율의 분모는 novel junction과 좌표가 정확히 같은 `N`을 가진 read 수임. 비율이 0.5를 넘으면 `mapping_or_repeat`임. soft-clip은 기본으로 꺼져 있음.
- SQANTI-SIM에서 가장 강한 흔적은 `non_canonical`(가짜 52%, 진짜 1%)이었고, `RTS_stage`는 약했으며, `perc_A`는 intra-priming을 흉내 내지 않는 시뮬레이션이라 평가할 수 없었음.
- indel near 기준 0.5는 HiFi용임. ONT 실제 데이터에서는 novel junction의 92%에서 켜져서 그대로 쓸 수 없음.
- 흔적이 여러 개면 mapping > noncanonical > rt_switch > degradation 순서로 하나만 보고하고, `rule_trace`에도 그 하나만 남음. 나머지는 JSON의 `evidence`에서 봐야 함.

다음 노트 [04](04_projection.md)에서는 [02](02_short_read_support.md)의 지지 수준과 이 노트의 기전을 조합해서 등급 하나를 만드는 규칙을 본다.

## 연습문제

### 문제 1

> `iso_C`의 `bam_n_spanning`은 10이다. 그런데 BAM에 `iso_C` read는 5개뿐이다. 10은 어디서 나왔을까?

<details>
<summary>풀이</summary>

```bash
python3 -c '
import json
for line in open("out/toy.attribution.jsonl"):
    r = json.loads(line)
    if r["isoform"] == "iso_C":
        print("n_novel_junctions", r["evidence"]["n_novel_junctions"], "| bam_n_spanning", r["evidence"]["bam_n_spanning"])'
samtools view -c data/long_reads.bam
samtools view data/long_reads.bam | grep -c '^iso_C'
```

```text
n_novel_junctions 2 | bam_n_spanning 10
34
5
```

`bam_n_spanning`은 read 수가 아니라 **novel junction마다 센 spanning read 수의 합**이다(`build_evidence()`의 `bam_n_spanning_total += f.n_spanning`). `iso_C`는 novel junction이 두 개이고, read 5개가 두 junction을 모두 걸치니 5 + 5 = 10이다. 반면 네 비율은 junction별로 계산한 뒤 **가장 큰 값**을 쓴다. 그래서 JSON의 `bam_frac_*`는 "가장 나쁜 junction의 비율"로 읽어야 한다.

</details>

### 문제 2

> `iso_G`의 read 4개 중 결실이 있는 read가 3개가 아니라 2개였다면?

<details>
<summary>풀이</summary>

`iso_G.r3`의 결실을 없앤 BAM을 만들어 돌렸다(samtools가 필요하다).

```bash
awk 'BEGIN{FS=OFS="\t"} $1=="iso_G.r3" {sub("5M3D88M", "96M", $6)} 1' data/long_reads.sam | samtools view -b -o out/half.bam -
samtools index out/half.bam
samtools view out/half.bam | grep '^iso_G' | cut -f1,6
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
  --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab --bam out/half.bam --out-prefix out/half 2>/dev/null
python3 -c '
import json
for line in open("out/half.attribution.jsonl"):
    r = json.loads(line)
    if r["isoform"] == "iso_G":
        print("indel_near", r["evidence"]["bam_frac_indel_near"], "->", r["primary_mechanism"], r["confidence_class"])'
```

```text
iso_G.r1	100M200N100M200N100M204N5M3D88M200N100M
iso_G.r2	100M200N100M200N100M204N5M3D88M200N100M
iso_G.r3	100M200N100M200N100M204N96M200N100M
iso_G.r4	100M200N100M200N100M204N96M200N100M
indel_near 0.5 -> none LOW_CONF_PARTIAL
```

비율이 정확히 0.5라서 기준(0.5 **초과**)을 넘지 못한다. 기전은 `none`이 되고, short read 지지도 없으니 `UNSUPPORTED × none → LOW_CONF_PARTIAL`이다. read가 4개뿐일 때는 read 하나 차이로 `ARTIFACT`와 `LOW_CONF_PARTIAL`이 갈린다. 비율 기준은 spanning read가 적은 junction에서 이렇게 불안정하다. 지금 코드에는 spanning read 수의 최솟값 조건이 없어서, read 1개짜리 junction도 그 read 하나에 indel이 있으면 비율 1.0으로 바로 걸린다.

</details>

## 더 깊이 보기

<details>
<summary>무작위 서열에서도 RT switching 흔적이 나온다</summary>

예제를 처음 만들 때 seed 7로 만든 서열로 SQANTI3를 돌렸더니 isoform 아홉 개 중 여덟 개에 `RTS_stage = TRUE`가 붙었다. SQANTI3 6.0.1의 `rt_switching.py`는 junction 주변에서 exon 쪽과 intron 쪽 서열에 같은 8-mer(`min_match = 8`, 기본값)가 있는지 찾는데, 위치를 조금씩 옮겨 가며(`wiggle_count = 1`) 찾고, QC가 `-a` 옵션으로 부르기 때문에 mismatch 1개까지 허용한다. 무작위 서열에서도 이 조건이 우연히 맞을 수 있다. seed 1–20을 모두 돌려 보니 3과 7에서는 여덟 개, 14와 16에서는 한 개에 `TRUE`가 붙었다. junction 하나가 걸리면 그 junction을 쓰는 isoform 전부에 붙는다. seed 7에서 `FALSE`로 남은 하나가 exon 4–5 사이 intron(1101–1300)을 쓰지 않는 `iso_D`뿐이었으니, 그 intron 하나가 걸린 것으로 보인다(SQANTI3의 junction별 출력으로 확인하지는 않았다). 1절에서 진짜 novel isoform의 9.9%에 `RTS_stage = TRUE`가 붙은 데에도 이런 우연이 섞여 있을 수 있다. 이 예제는 RT switching을 다루지 않으려고 seed 1을 썼다.

</details>

<details>
<summary>기전이 없는 축(`not_evaluable`)은 어떻게 다루나</summary>

`--bam`이 없으면 BAM 비율은 모두 0으로 남고 `bam_evaluable = false`가 된다. 이때 mapping 조건은 아예 검사하지 않는다(`ev.bam_evaluable`이 조건에 들어 있다). SQANTI3 분류표에 `RTS_stage`나 `all_canonical` 열이 없거나 값이 비어 있으면 그 기전도 검사하지 않는다(`rts_evaluable`, `canon_evaluable`). `perc_A_downstream_TTS`가 `NA`면 degradation을 검사하지 않는다. 어느 경우든 "흔적이 없다(`none`)"와 구별되지 않고 같은 결과가 되니, 입력이 빠졌는지는 `provenance.log`와 JSON의 `*_evaluable`로 확인해야 한다.

</details>
