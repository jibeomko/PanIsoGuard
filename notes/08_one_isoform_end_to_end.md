# 08. isoform 하나는 어떻게 등급 하나가 될까?

시리즈의 마지막 노트다. 주인공 `iso_A`를 입력 파일에서 최종 등급까지 한 번에 따라간다. 이번에는 PanIsoGuard를 쓰지 않고, Python 표준 라이브러리만으로 파일을 직접 읽어서 계산한다. 앞 노트들에서 한 조각씩 확인한 규칙을 이어 붙이면 [00](00_overview.md) 1절의 표가 그대로 나오는지 보는 것이 목표다. 마지막에는 예제의 isoform 아홉 개 전부를 입력 조합 아홉 가지로 같은 방식으로 다시 계산하는 스크립트([check_notes.py](check_notes.py))를 돌린다.

> 관련 문서: [docs/algorithm.md](../docs/algorithm.md) · 코드: [check_notes.py](check_notes.py)

## 1. `iso_A`의 증거를 하나씩 모으면?

[check_notes.py](check_notes.py)에는 파일을 읽고 증거를 계산하는 짧은 함수들이 들어 있다. 모두 앞 노트에서 한 번씩 본 계산이다. pysam이나 numpy 없이 표준 라이브러리만 쓰기 때문에, BAM 대신 같은 내용의 SAM 텍스트(`data/long_reads.sam`)를 읽는다. 이 함수들로 `iso_A`의 증거를 순서대로 모은다. 코드는 `notes/` 폴더에서 실행한다.

```python
import check_notes as cn        # the stdlib-only helpers in notes/check_notes.py

chains = cn.gtf_chains("data/caller1.gtf")
chrom, strand, introns = chains["iso_A"]
print("1. intron chain:", chrom, strand, introns)

catalog = {(c, s, j) for c, s, js in cn.gtf_chains("data/reference.gtf").values() for j in js}
novel = [j for j in introns if (chrom, strand, j) not in catalog]
print("2. novel junctions:", novel)
jx = novel[0]

sj = cn.read_sj("data/short_reads.SJ.out.tab")
print("3. SJ.out.tab row for it:", sj.get((chrom, strand, jx)), "-> short-read supported 0/1 -> UNSUPPORTED")

row = cn.read_classification("data/caller1_classification.txt")["iso_A"]
print("4. SQANTI3:", {k: row[k] for k in ("structural_category", "all_canonical", "RTS_stage", "perc_A_downstream_TTS")})

n, hits = cn.bam_features(cn.read_sam("data/long_reads.sam"), jx)
print(f"5. long reads spanning it: {n}, flagged: {hits}")

ref = cn.read_fasta("data/toy.fa")
haps = [cn.read_fasta("data/hap1.fa"), cn.read_fasta("data/hap2.fa")]
print("6. canonical on reference:", cn.canonical(ref, *jx, strand),
      "| on hap1, hap2:", [cn.canonical(h, *jx, strand) for h in haps], "-> CREATED")
```

```text
1. intron chain: chrT + [(200, 400), (500, 730), (800, 1000), (1100, 1300)]
2. novel junctions: [(500, 730)]
3. SJ.out.tab row for it: None -> short-read supported 0/1 -> UNSUPPORTED
4. SQANTI3: {'structural_category': 'novel_not_in_catalog', 'all_canonical': 'non_canonical', 'RTS_stage': 'FALSE', 'perc_A_downstream_TTS': '45.0'}
5. long reads spanning it: 5, flagged: {'low_mapq': 0, 'supplementary': 0, 'softclip': 0, 'indel_near': 0}
6. canonical on reference: False | on hap1, hap2: [True, False] -> CREATED
```

한 줄씩 앞 노트와 이어 보면 이렇다.

1. `iso_A`의 intron 네 개. 0-based half-open 좌표다([01](01_intron_chain_and_novelty.md) 1–2절).
2. 참조 주석의 intron 여섯 개와 비교하면 `(500, 730)` 하나가 새롭다. 1-based로 501–730, 길이 230 bp다([01](01_intron_chain_and_novelty.md) 3절).
3. `SJ.out.tab`에 이 junction이 없다. short read 9개가 정렬됐지만 STAR가 non-canonical junction 필터로 뺐다([02](02_short_read_support.md) 5절).
4. SQANTI3는 `non_canonical`을 적었다. 그래서 기전은 `noncanonical`이다([03](03_artifact_mechanisms.md) 1절).
5. long read 5개가 이 junction을 걸쳤고 네 흔적 모두 0개다. mapping 기전은 없다([03](03_artifact_mechanisms.md) 2–3절).
6. 참조에서는 non-canonical(GT-AC), haplotype 1에서는 canonical(GT-AG)이다. 유일한 novel junction이 `CREATED`이니 rescue 후보다([05](05_reference_bias.md) 2절).

## 2. 입력을 하나씩 더하면 등급은 어떻게 바뀔까?

같은 함수의 `verdict()`는 [04](04_projection.md)의 규칙 전체(FSM/ISM 통과, rescue, 지지 수준, 기전 순서, 조합표)를 담고 있다. 어떤 입력을 넣었는지를 사전 하나로 넘기면, 그 입력만 보고 판정한다. [00](00_overview.md) 1절의 여섯 단계를 그대로 재현해서, 그때 PanIsoGuard가 쓴 `out/step*.adjudicated.tsv`와 비교한다.

```python
import csv
import check_notes as cn

chains = cn.gtf_chains("data/caller1.gtf")
row = cn.read_classification("data/caller1_classification.txt")["iso_A"]
catalog = {(c, s, j) for c, s, js in cn.gtf_chains("data/reference.gtf").values() for j in js}
sj, reads = cn.read_sj("data/short_reads.SJ.out.tab"), cn.read_sam("data/long_reads.sam")
hap = dict(ref_seq=cn.read_fasta("data/toy.fa"), haps=[cn.read_fasta("data/hap1.fa"), cn.read_fasta("data/hap2.fa")])
steps = [{}, dict(catalog=catalog), dict(catalog=catalog, sj=sj), dict(catalog=catalog, sj=sj, reads=reads),
         dict(catalog=catalog, sj=sj, reads=reads, ref_seq_circular=True, **hap),
         dict(catalog=catalog, sj=sj, reads=reads, ref_seq_circular=False, **hap)]
for i, inp in enumerate(steps):
    mine = cn.verdict("iso_A", row, chains["iso_A"], inp)
    pig = next(r for r in csv.DictReader(open(f"out/step{i}.adjudicated.tsv"), delimiter="\t") if r["isoform_id"] == "iso_A")
    same = mine[:3] == (pig["novelty_support"], pig["primary_mechanism"], pig["confidence_class"])
    print(f"step{i}: by hand {mine[0]:11s} {mine[1]:15s} {mine[2]:28s} | PanIsoGuard {'same' if same else 'DIFFERENT'}")
```

```text
step0: by hand UNKNOWN     noncanonical    AMBIGUOUS                    | PanIsoGuard same
step1: by hand UNKNOWN     noncanonical    AMBIGUOUS                    | PanIsoGuard same
step2: by hand UNSUPPORTED noncanonical    ARTIFACT                     | PanIsoGuard same
step3: by hand UNSUPPORTED noncanonical    ARTIFACT                     | PanIsoGuard same
step4: by hand UNKNOWN     variant_created AMBIGUOUS                    | PanIsoGuard same
step5: by hand UNKNOWN     variant_created PAN_REF_RESCUED_FALSE_NOVEL  | PanIsoGuard same
```

여섯 단계 모두 같다. 단계마다 무엇이 판정을 바꿨는지 정리하면 이렇다.

| 단계 | 더한 입력 | 바뀐 증거 | 적용된 규칙 | 등급 |
|---|---|---|---|---|
| step0 | 분류표 + caller GTF | 기전 `noncanonical`(SQANTI3) | `UNKNOWN` × 흔적 → 보류 | `AMBIGUOUS` |
| step1 | + 참조 GTF | novel junction 1개 | 확인할 입력 없음, 여전히 `UNKNOWN` | `AMBIGUOUS` |
| step2 | + `SJ.out.tab` | short read 지지 0/1 | `UNSUPPORTED` × `noncanonical` | `ARTIFACT` |
| step3 | + BAM | spanning 5개, 흔적 0 | 기전 그대로 | `ARTIFACT` |
| step4 | + haplotype (출처 unknown) | `CREATED` 1/1 | rescue 후보지만 firewall이 보류 | `AMBIGUOUS` + `circularity_flag` |
| step5 | + 출처 `wgs` | 같음 | rescue | `PAN_REF_RESCUED_FALSE_NOVEL` |

step0에서 step1로 갈 때는 novel junction을 찾았지만 등급은 그대로다. step3에서 BAM을 더해도 흔적이 없어 등급은 그대로다. 판정을 실제로 움직인 입력은 short read 표와 haplotype 두 가지였다. step4와 step5의 지지 수준이 `UNKNOWN`으로 돌아가는 것은 rescue가 조합표보다 먼저 판정을 끝내기 때문이다([04](04_projection.md) 2절).

## 3. 예제 전체를 다시 계산하면?

[check_notes.py](check_notes.py)는 같은 계산을 isoform 아홉 개 전부, 입력 조합 아홉 가지(위의 여섯 단계, pangenome, caller 합의, caller 합의 + BAM)에 대해 한다. PanIsoGuard도 같은 조합으로 돌려서 isoform마다 지지 수준, 기전, 등급, novel junction 수, 확인된 junction 수, circularity flag가 모두 같은지 비교하고, 노트에 나온 주요 숫자도 확인한다. 하나라도 다르면 오류를 내며 멈춘다.

```bash
python3 check_notes.py
```

```text
step0          9 isoforms, iso_A -> AMBIGUOUS
step1          9 isoforms, iso_A -> AMBIGUOUS
step2          9 isoforms, iso_A -> ARTIFACT
step3          9 isoforms, iso_A -> ARTIFACT
step4          9 isoforms, iso_A -> AMBIGUOUS
step5          9 isoforms, iso_A -> PAN_REF_RESCUED_FALSE_NOVEL
pangenome      9 isoforms, iso_A -> PAN_REF_RESCUED_FALSE_NOVEL
consensus      9 isoforms, iso_A -> AMBIGUOUS
consensus_bam  9 isoforms, iso_A -> AMBIGUOUS
ok   iso_A across the six steps (00)
ok   iso_C is PARTIAL -> LOW_CONF_PARTIAL (02)
ok   iso_G is a mapping artifact (03)
ok   iso_A novel junction absent from SJ.out.tab (02)
ok   fingerprint of iso_known (06)
ok   consensus lifts iso_G, BAM drops it (06)
ok   AUPRC 0.9712 from the SQANTI-SIM confusion table (07)
All checks passed.
```

81개(9 × 9) 판정이 모두 같다. 이 스크립트는 CTest(`integration_study_notes`)에도 들어 있어서, 판정 규칙이 바뀌면 이 노트들과 어긋났다는 것을 CI가 알려 준다. 스크립트가 PanIsoGuard의 규칙을 그대로 옮겨 쓴 것이라, 규칙이 둘 다에서 똑같이 틀려도 이 검사는 통과한다는 점은 기억해 둬야 한다. 이 스크립트가 보장하는 것은 "노트에 적힌 설명이 실제 동작과 같다"이지 "판정이 맞다"가 아니다. 판정이 맞는지는 [07](07_evaluation.md)의 정답 데이터로만 잴 수 있다.

## 정리

- `iso_A`의 판정은 여섯 가지 증거로 정해짐. intron chain, novel junction 하나(501–730), short read 지지 없음, SQANTI3의 `non_canonical`, long read 흔적 없음, haplotype 1에서 canonical이 되는 SNV임.
- 판정을 실제로 움직인 입력은 short read 표(`AMBIGUOUS` → `ARTIFACT`)와 haplotype과 그 출처(`ARTIFACT` → 보류 → rescue)였음.
- 이 과정을 표준 라이브러리만으로 다시 계산해서 예제의 81개 판정을 모두 재현했고, CTest로 묶어 두었음. 이것은 설명과 동작이 같다는 확인이지, 판정이 맞다는 증거는 아님.

시리즈를 다 따라왔다면 [README](README.md)의 "문서와 다르게 나온 것들"에서 이 공부 중에 발견해서 고친 문서 오류들을 한 번에 볼 수 있다.

## 연습문제

### 문제 1

> step3에서 `iso_A`의 long read 5개 중 3개가 low MAPQ였다면 등급은 어떻게 될까? step5(`wgs`)에서는?

<details>
<summary>풀이</summary>

`verdict()`에 넘기는 read 목록에서 `iso_A` read 세 개의 MAPQ만 5로 바꿔 보면 된다. 앞 블록의 변수를 이어서 쓴다. read는 SAM의 이름으로 고른다.

```python
names = [l.split("\t")[0] for l in open("data/long_reads.sam") if not l.startswith("@")]
iso_a = [i for i, n in enumerate(names) if n.startswith("iso_A.")][:3]
low = [(f, p, 5 if i in iso_a else q, c) for i, (f, p, q, c) in enumerate(reads)]
print("MAPQ set to 5 for:", [names[i] for i in iso_a])
for i in (3, 5):
    print(f"step{i}:", cn.verdict("iso_A", row, chains["iso_A"], dict(steps[i], reads=low))[:3])
```

```text
MAPQ set to 5 for: ['iso_A.r1', 'iso_A.r2', 'iso_A.r3']
step3: ('UNSUPPORTED', 'mapping_or_repeat', 'ARTIFACT')
step5: ('UNKNOWN', 'variant_created', 'PAN_REF_RESCUED_FALSE_NOVEL')
```

`iso_A`의 spanning read 5개 중 3개(0.6)가 low MAPQ가 되면 기준 0.5를 넘으니, step3의 기전은 `noncanonical`에서 `mapping_or_repeat`로 바뀐다. 등급은 그대로 `ARTIFACT`다.

step5는 여전히 rescue다. reference-bias rescue는 조합표보다 먼저 판정을 끝내기 때문에 mapping 흔적을 보지 않는다. 정렬 위치부터 의심스러운 junction도 "reference bias"로 구제될 수 있다는 뜻이다. rescue된 isoform은 JSON의 `evidence`에서 `bam_frac_*` 값을 따로 확인해 보는 것이 좋다.

</details>
