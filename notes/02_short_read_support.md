# 02. short read가 novel junction을 봤다는 것을 어떻게 확인할까?

[01](01_intron_chain_and_novelty.md)에서 novel junction을 찾았다. 그다음 질문은 이 junction이 진짜인지다. long read는 길어서 isoform 전체를 한 번에 읽지만 오류도 많다. 그래서 같은 시료의 short-read RNA-seq에서도 이 junction을 걸친 read가 나왔다면, long read의 오류로 생긴 junction일 가능성이 크게 줄어든다. PanIsoGuard는 이 확인을 STAR가 만든 `SJ.out.tab` 하나로 한다. 이 노트에서는 junction 하나가 언제 "지지받았다"고 인정되는지, 그리고 junction 여러 개의 결과를 isoform 하나의 지지 수준으로 어떻게 묶는지 따라가 본다. 주인공 `iso_A`의 novel junction은 short read 9개가 걸쳤는데도 지지를 못 받는데, 그 이유도 여기서 찾는다.

> 관련 문서: [docs/decision_engine.md](../docs/decision_engine.md) "Two-axis evidence space", [docs/input_formats.md](../docs/input_formats.md) · 코드: `build_evidence()`(`src/core/adjudicator.cpp`), `src/io/sj_tab_reader.cpp`

## 1. `SJ.out.tab`에는 무엇이 적혀 있을까?

STAR는 short read를 정렬하면서 발견한 splice junction을 한 줄에 하나씩 `SJ.out.tab`에 모아 둔다. 예제의 short read 160개를 STAR로 정렬해서 얻은 표다([00](00_overview.md) "더 깊이 보기").

```bash
PIG=../build/panisoguard; mkdir -p out
cat data/short_reads.SJ.out.tab
```

```text
chrT	201	400	1	1	1	30	0	50
chrT	201	700	1	1	1	10	0	48
chrT	501	700	1	1	1	25	0	49
chrT	501	1000	1	1	0	6	0	49
chrT	801	860	1	1	0	7	0	49
chrT	801	1000	1	1	1	20	0	50
chrT	801	1020	1	1	0	12	0	48
chrT	801	1300	1	1	1	10	0	47
chrT	941	1000	1	1	0	1	0	50
chrT	1101	1300	1	1	1	30	0	50
```

열은 아홉 개다.

| 열 | 이름 | 뜻 | PanIsoGuard가 쓰나 |
|---|---|---|---|
| 1 | chrom | 염색체 | 씀 |
| 2, 3 | intron start, end | intron의 첫 염기와 마지막 염기, 1-based | 씀 (0-based half-open으로 바꿔서) |
| 4 | strand | 0 = 모름, 1 = +, 2 = - | 씀 |
| 5 | motif | 0 = non-canonical, 1 = GT/AG, 2 = CT/AC, 3 = GC/AG, 4 = CT/GC, 5 = AT/AC, 6 = GT/AT | 씀 (1 이상인지만) |
| 6 | annotated | STAR에 준 주석(GTF)에 있던 junction이면 1 | 안 씀 |
| 7 | n_uniq | 이 junction을 걸친 uniquely mapped read 수 | 씀 |
| 8 | n_multi | 여러 곳에 정렬된(multi-mapped) read 수 | 안 씀 |
| 9 | max overhang | junction 양쪽으로 read가 뻗은 길이 중 가장 긴 값 | 안 씀 |

6열의 `annotated`를 쓰지 않는 이유는, 무엇이 새로운지는 PanIsoGuard가 자기 catalog으로 따로 정하기 때문이다([01](01_intron_chain_and_novelty.md)). STAR에 준 주석과 `--ref-gtf`가 다를 수도 있다.

표에서 `iso_B`의 novel junction(801–1020)은 12개, `iso_C`의 두 novel junction(801–860, 941–1000)은 7개와 1개, `iso_E`의 novel junction(501–1000)은 6개 read가 걸쳤다. `iso_A`의 novel junction(501–730)은 표에 아예 없다. 5절에서 이유를 찾는다.

## 2. junction 하나는 언제 "지지받았다"고 볼까?

`build_evidence()`에서 novel junction 하나를 확인하는 조건은 세 가지다.

1. `SJ.out.tab`에 **염색체, strand, 좌표가 정확히 같은** 줄이 있다. 1 bp라도 다르면 없는 것이다.
2. 그 줄의 `n_uniq`가 **3 이상**이다(`sj_min_uniq_reads`).
3. 그 줄의 `motif`가 **1 이상**, 곧 STAR가 canonical motif로 분류했다(`sj_require_canonical_motif`).

`n_multi`는 세지 않는다. 여러 곳에 정렬된 read는 어느 자리에서 왔는지 모르니, 이 junction을 확인해 주는 증거로 쓰기 어렵다(연습문제 1). 이 규칙을 Python으로 그대로 옮겨서 isoform마다 계산하고, PanIsoGuard의 결과(00의 step2)와 비교해 본다. `chains()`는 [01](01_intron_chain_and_novelty.md)과 같은 함수다.

```python
import csv

def chains(path):
    """GTF -> {transcript: (chrom, strand, introns as 0-based half-open)}"""
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
STRAND = {"1": "+", "2": "-"}                          # STAR strand code; 0 = undefined
sj = {}
for f in csv.reader(open("data/short_reads.SJ.out.tab"), delimiter="\t"):
    key = (f[0], STRAND.get(f[3], "."), (int(f[1]) - 1, int(f[2])))
    sj[key] = {"motif": int(f[4]), "n_uniq": int(f[6])}

def supported(key, min_uniq=3, need_canonical=True):
    r = sj.get(key)
    return r is not None and r["n_uniq"] >= min_uniq and (not need_canonical or r["motif"] >= 1)

pig = {r["isoform_id"]: r["novelty_support"] for r in csv.DictReader(open("out/step2.adjudicated.tsv"), delimiter="\t")}
for t, (c, s, js) in chains("data/caller1.gtf").items():
    if t in ("iso_known", "iso_ism"):
        continue                                        # FSM / ISM never reach this axis
    novel = [j for j in js if (c, s, j) not in catalog]
    reads = [sj.get((c, s, j), {}).get("n_uniq", 0) for j in novel]
    k = sum(supported((c, s, j)) for j in novel)
    level = ("UNKNOWN" if not novel else "SUPPORTED" if k == len(novel) else
             "UNSUPPORTED" if k == 0 else "PARTIAL")
    print(f"{t}: novel {novel} unique reads {reads} -> {k}/{len(novel)} {level:11s} | PanIsoGuard {pig[t]}")
```

```text
iso_A: novel [(500, 730)] unique reads [0] -> 0/1 UNSUPPORTED | PanIsoGuard UNSUPPORTED
iso_B: novel [(800, 1020)] unique reads [12] -> 1/1 SUPPORTED   | PanIsoGuard SUPPORTED
iso_C: novel [(800, 860), (940, 1000)] unique reads [7, 1] -> 1/2 PARTIAL     | PanIsoGuard PARTIAL
iso_D: novel [] unique reads [] -> 0/0 UNKNOWN     | PanIsoGuard UNKNOWN
iso_E: novel [(500, 1000)] unique reads [6] -> 1/1 SUPPORTED   | PanIsoGuard SUPPORTED
iso_F: novel [(200, 450)] unique reads [0] -> 0/1 UNSUPPORTED | PanIsoGuard UNSUPPORTED
iso_G: novel [(800, 1004)] unique reads [0] -> 0/1 UNSUPPORTED | PanIsoGuard UNSUPPORTED
```

일곱 개 모두 같다. `out/step2.adjudicated.tsv`는 [00](00_overview.md)의 step2 실행 결과다(참조 GTF + `SJ.out.tab`).

## 3. junction 여러 개의 결과를 isoform 하나로 어떻게 묶을까?

위 코드의 `level` 줄이 그 규칙이다. isoform의 novel junction이 `n`개, 그중 확인된 것이 `k`개라면 이렇다.

| 조건 | 지지 수준 | 예제 |
|---|---|---|
| `k = n` (모두 확인) | `SUPPORTED` | `iso_B`, `iso_E` |
| `0 < k < n` (일부만) | `PARTIAL` | `iso_C` (1/2) |
| `k = 0` | `UNSUPPORTED` | `iso_A`, `iso_F`, `iso_G` |
| 판정할 수 없음 | `UNKNOWN` | `iso_D` |

`UNKNOWN`에는 이유가 세 가지 있고, `rule_trace` 첫 줄에 어느 것인지 적힌다. 모두 앞 노트들에서 한 번씩 만났다.

| `rule_trace`의 이유 | 뜻 | 만난 곳 |
|---|---|---|
| `caller_chain_absent` | 분류표의 isoform 이름과 같은 intron chain이 caller 파일에 없음 | [01](01_intron_chain_and_novelty.md) 연습문제 2 |
| `short_read/catalog_axis_absent` | `--sj-tab`이나 `--ref-gtf`가 없음 | [00](00_overview.md) 연습문제 1 |
| `no_novel_junctions` | 확인할 novel junction이 없음 (알려진 junction의 조합 등) | `iso_D` |

`UNSUPPORTED`와 `UNKNOWN`은 다르다는 점을 기억해 두면 좋다. `UNSUPPORTED`는 "short read를 봤는데 없었다", `UNKNOWN`은 "볼 수 없었다"이다. [04](04_projection.md)에서 보겠지만 같은 artifact 흔적이라도 `UNSUPPORTED`와 만나면 `ARTIFACT`, `UNKNOWN`과 만나면 대개 `AMBIGUOUS`가 된다.

## 4. `PARTIAL`은 왜 더 이상 `MEDIUM_CONF_NOVEL`이 아닐까?

`iso_C`는 새 exon 하나가 끼어들어 novel junction이 두 개 생긴 isoform이다. 앞쪽(801–860)은 read 7개가 확인했지만, 뒤쪽(941–1000)은 read가 1개뿐이라 기준(3개)에 못 미친다. 0.0.3까지는 이런 `PARTIAL`을 `MEDIUM_CONF_NOVEL`로 올렸고, 0.0.4부터는 `LOW_CONF_PARTIAL`에 둔다. 같은 입력을 두 버전에 넣어 봤다. `$OLD`는 `v0.0.3` 태그를 따로 빌드한 바이너리다.

```bash
A="--classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf --ref-gtf data/reference.gtf --sj-tab data/short_reads.SJ.out.tab"
$OLD adjudicate $A --out-prefix out/v003 2>/dev/null
$PIG adjudicate $A --out-prefix out/v004 2>/dev/null
for v in v003 v004; do
  printf '%s  %s\n' $v "$(grep '"iso_C"' out/$v.attribution.jsonl | python3 -c 'import json,sys; r=json.loads(sys.stdin.read()); print(r["confidence_class"], "|", " ; ".join(r["rule_trace"]))')"
done
```

```text
v003  MEDIUM_CONF_NOVEL | sj_support=1/2 -> PARTIAL ; no artifact mechanism flagged -> mechanism=none ; project(PARTIAL,none) -> MEDIUM_CONF_NOVEL
v004  LOW_CONF_PARTIAL | sj_support=1/2 -> PARTIAL ; no artifact mechanism flagged -> mechanism=none ; project(PARTIAL,none) -> LOW_CONF_PARTIAL
```

증거는 똑같고 마지막 조합 규칙만 다르다. 왜 바꿨을까? novel junction이 두 개인 isoform에서 하나만 확인됐다면, 나머지 하나는 여전히 long read만의 주장이다. isoform은 junction 전부가 맞아야 맞는 것이니, 확인 안 된 junction이 하나라도 남으면 isoform 전체가 확인됐다고 할 수 없다.

규칙만 따져서가 아니라 실제로 틀리는 사례가 나와서 바꿨다. SQANTI-SIM 시뮬레이션 정답으로 PanIsoGuard를 SQANTI3 필터와 비교했을 때([07](07_evaluation.md)), NNC에서 PanIsoGuard만 진짜로 부른 가짜 isoform이 셋 있었다. 셋 다 `rule_trace`가 `sj_support=1/2 -> PARTIAL`이었다. 0.0.4에서 규칙을 바꾸자 이 셋이 빠지면서 precision이 0.982에서 0.989로, 가짜를 걸러 내는 비율(specificity)이 0.946에서 0.966으로 올랐고, 진짜를 찾는 비율(recall)은 그대로였다([benchmark/sqanti3_filter_h2h](../benchmark/sqanti3_filter_h2h)).

## 5. `iso_A`의 junction은 read 9개가 걸쳤는데 왜 지지를 못 받을까?

예제를 만들 때 `iso_A`의 novel junction(501–730)을 걸치는 short read를 9개 넣었다. 기준(3개)의 세 배다. 그런데 1절의 `SJ.out.tab`에는 이 junction이 없다. read가 정렬되지 않은 걸까? STAR를 다시 돌려서 정렬 결과(SAM)를 직접 봤다. 이 블록은 STAR가 있어야 돌아간다(SQANTI3 conda 환경에 STAR 2.7.11b가 들어 있다).

```bash
rm -rf out/star && mkdir -p out/star/idx && cd out/star
STAR --runMode genomeGenerate --genomeDir idx --genomeFastaFiles ../../data/toy.fa \
     --sjdbGTFfile ../../data/reference.gtf --sjdbOverhang 99 --genomeSAindexNbases 4 --outFileNamePrefix idx/ > /dev/null
STAR --genomeDir idx --readFilesIn ../../data/short_reads.fq --outSAMtype SAM --outFileNamePrefix default/ > /dev/null
grep '^iso_A' default/Aligned.out.sam | cut -f1-6 | head -3
grep -c $'\t[0-9]*M230N[0-9]*M\t' default/Aligned.out.sam
grep -P '\t501\t730\t' default/SJ.out.tab || echo "501-730 is not in SJ.out.tab"
STAR --genomeDir idx --readFilesIn ../../data/short_reads.fq --outSAMtype None \
     --outSJfilterDistToOtherSJmin 0 0 0 0 --outFileNamePrefix nofilter/ > /dev/null
grep -P '\t501\t730\t' nofilter/SJ.out.tab
cd ../..
```

```text
iso_A.j2.106	0	chrT	444	255	57M230N43M
iso_A.j2.107	0	chrT	454	255	47M230N53M
iso_A.j2.108	0	chrT	453	255	48M230N52M
9
501-730 is not in SJ.out.tab
chrT	501	730	0	0	0	9	0	48
```

read 9개는 모두 제자리에 정렬되었다. CIGAR의 `230N`이 바로 501–730, 230 bp짜리 intron이다. 그런데 기본 설정의 `SJ.out.tab`에는 이 junction이 빠졌다. STAR는 정렬과 별개로 `SJ.out.tab`에 적을 junction을 한 번 더 거르는데, 그중 `--outSJfilterDistToOtherSJmin`의 기본값이 `10 0 5 10`이다. 순서대로 non-canonical, GT/AG, GC/AG, AT/AC junction에 대한 값이고, non-canonical junction은 donor나 acceptor가 다른 junction의 donor나 acceptor에서 10 bp 이상 떨어져 있어야 남긴다는 뜻이다. `iso_A`의 junction은 donor(501)가 알려진 intron 501–700의 donor와 같은 자리라서 거리가 0이고, 그래서 빠졌다. 이 필터를 끄자(`0 0 0 0`) 마지막 줄처럼 나타난다.

그런데 필터를 꺼서 나타난 줄을 보면 strand가 0(모름), motif가 0(non-canonical)이다. 참조 게놈에서 이 intron의 양 끝은 `GT…AC`이고, STAR는 canonical motif가 아니면 strand를 정하지 못한다(unstranded library에서는 motif로 strand를 추론하기 때문이다). 그러니 이 줄을 PanIsoGuard에 줘도 두 관문이 더 남는다. strand 0은 `+`인 `iso_A`와 열쇠가 맞지 않고, motif 0은 3번 조건에 걸린다. 관문을 하나씩 열어 가며 돌려 봤다.

```bash
sed 's/^sj_require_canonical_motif *= *true/sj_require_canonical_motif = false/' ../config/rules.default.toml > out/no_motif.toml
grep '^sj_require_canonical_motif' out/no_motif.toml
awk 'BEGIN{FS=OFS="\t"} $2==501 && $3==730 {$4=1} 1' out/star/nofilter/SJ.out.tab > out/nofilter_plus.SJ.tab
grep -P '\t501\t730\t' out/nofilter_plus.SJ.tab
run() {
  $PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
    --ref-gtf data/reference.gtf --sj-tab "$1" ${2:+--config $2} --out-prefix out/barrier 2>/dev/null
  printf '%-34s %-22s %s\n' "$1" "${2:-default config}" "$(grep -P '^iso_A\t' out/barrier.adjudicated.tsv | cut -f5-7 | tr '\t' ' ')"
}
run data/short_reads.SJ.out.tab
run out/star/nofilter/SJ.out.tab
run out/star/nofilter/SJ.out.tab out/no_motif.toml
run out/nofilter_plus.SJ.tab
run out/nofilter_plus.SJ.tab out/no_motif.toml
```

```text
sj_require_canonical_motif = false   # require SJ.tab motif code >= 1
chrT	501	730	1	0	0	9	0	48
data/short_reads.SJ.out.tab        default config         UNSUPPORTED noncanonical ARTIFACT
out/star/nofilter/SJ.out.tab       default config         UNSUPPORTED noncanonical ARTIFACT
out/star/nofilter/SJ.out.tab       out/no_motif.toml      UNSUPPORTED noncanonical ARTIFACT
out/nofilter_plus.SJ.tab           default config         UNSUPPORTED noncanonical ARTIFACT
out/nofilter_plus.SJ.tab           out/no_motif.toml      SUPPORTED noncanonical MEDIUM_CONF_NOVEL
```

관문 세 개(STAR의 필터, strand 0, motif 조건)를 모두 열어야 비로소 `SUPPORTED`가 된다. 이때도 SQANTI3가 `all_canonical = non_canonical`을 적어 두었으니 기전은 `noncanonical`이고, 등급은 `MEDIUM_CONF_NOVEL`에 그친다([04](04_projection.md)).

이 세 관문은 모두 같은 사실 하나에서 나온다. **short read는 참조 게놈에 정렬되고, 참조 게놈에서 이 junction은 non-canonical이다.** non-canonical junction은 정렬 오류로 생기는 경우가 훨씬 많아서, STAR도 PanIsoGuard도 이런 junction을 기본적으로 믿지 않는다. 대개는 옳은 판단이다(`iso_F`가 그런 경우다).

`iso_A`의 사정은 다르다. 이 사람의 haplotype 1에서는 730번 염기 하나가 달라서 이 junction의 acceptor가 `AG`, 곧 canonical이다. 그런데 730번은 intron의 마지막 염기라 splicing 때 잘려 나간다. mRNA에서 온 read에는 이 염기가 아예 없다. 그러니 short read를 아무리 많이 모아도 "이 사람에게서는 canonical"이라는 사실을 read 안에서 찾을 수 없다. 이것을 알려면 이 사람의 게놈 서열을 따로 봐야 한다. [05](05_reference_bias.md)의 주제다.

## 6. 기준을 바꾸면 어떻게 될까?

기준값은 코드에 고정되어 있지 않고 TOML 설정 파일로 바꿀 수 있다. 기본값은 [config/rules.default.toml](../config/rules.default.toml)에 있다. read 1개만 있어도 확인된 것으로 치게 바꿔 본다.

```bash
sed 's/^sj_min_uniq_reads *= *3/sj_min_uniq_reads = 1/' ../config/rules.default.toml > out/min1.toml
grep '^sj_min_uniq_reads' out/min1.toml
$PIG adjudicate $A --config out/min1.toml --out-prefix out/min1 2>/dev/null
paste <(cut -f1,5,7 out/v004.adjudicated.tsv) <(cut -f5,7 out/min1.adjudicated.tsv) | grep -E 'isoform|iso_C' | column -t
grep -E '^(ruleset|config|thresholds)' out/min1.provenance.log
```

```text
sj_min_uniq_reads = 1      # STAR SJ.tab n_uniq to corroborate a novel junction
isoform_id  novelty_support  confidence_class  novelty_support  confidence_class
iso_C       PARTIAL          LOW_CONF_PARTIAL  SUPPORTED        HIGH_CONF_NOVEL
ruleset_version	default-0.0.2
config	out/min1.toml
thresholds	sj_min_uniq_reads=1 sj_require_canonical_motif=true consensus_min_callers=2 max_perc_A_downstream_TTS=60 min_mapq=20 softclip_min_bp=20 junction_window_bp=10 max_low_mapq_frac=0.5 max_supplementary_frac=0.5 max_indel_near_frac=0.5 max_softclip_frac=1.01 pangenome_min_haplotypes=1
```

`iso_C`의 두 번째 junction(read 1개)이 기준을 넘으면서 `iso_C`가 `HIGH_CONF_NOVEL`로 올라간다. 기준 하나로 등급이 두 단계 뛴다.

`provenance.log`의 마지막 세 줄을 같이 보면 조심할 점이 하나 있다. `ruleset_version`은 TOML의 `[meta]`에 적힌 이름을 그대로 옮긴 것이라, 파일을 복사해서 기준값만 바꾸면 이름은 `default-0.0.2` 그대로다. 그래서 실제로 쓴 기준값은 `thresholds` 줄에 따로 적힌다. 여기서는 `sj_min_uniq_reads=1`이 보인다. 설정 파일이 나중에 바뀌거나 사라져도 이 줄로 판정을 재현할 수 있다.

이 노트를 처음 쓸 때는 `thresholds` 줄이 없었고, 0.0.4의 `PARTIAL` 규칙 변경(4절)도 `ruleset_version`에 드러나지 않았다(0.0.3과 0.0.4가 모두 `builtin-0.0.1`). 조합 규칙은 TOML이 아니라 코드에 있기 때문이다. 이 노트를 쓰고 나서 기본 규칙 이름을 `0.0.2`로 올리고, 앞으로 기본 기준값이나 조합 규칙이 바뀌면 함께 올리도록 코드에 적어 두었다([rules.hpp](../include/panisoguard/rules.hpp)).

기준을 어디에 둘지는 short read의 깊이에 달렸다. SQANTI-SIM에서 기준을 1, 3, 5, 10, 60으로 바꿔 가며 돌린 결과는 [benchmark/results/sqanti_sim/sweep.tsv](../benchmark/results/sqanti_sim/sweep.tsv)에 있다. 그 데이터는 정답에서 만든 `SJ.tab`에 모든 junction의 read 수를 50으로 적어 두었기 때문에, 기준이 50 이하일 때는 결과가 거의 변하지 않고 50을 넘는 순간 한꺼번에 무너진다. 실제 short read의 깊이로 다시 따져 본 결과는 [07](07_evaluation.md)에 있다.

## 정리

- novel junction 하나는 `SJ.out.tab`에 염색체, strand, 좌표가 정확히 같은 줄이 있고, `n_uniq ≥ 3`, `motif ≥ 1`일 때 확인된 것으로 침. `n_multi`와 `annotated`는 쓰지 않음.
- isoform의 지지 수준은 novel junction `n`개 중 확인된 `k`개로 정함. 모두면 `SUPPORTED`, 일부면 `PARTIAL`, 없으면 `UNSUPPORTED`, 볼 수 없으면 `UNKNOWN`(이유 세 가지)임.
- 0.0.4부터 `PARTIAL`은 `LOW_CONF_PARTIAL`에 머묾. 확인 안 된 junction이 남은 isoform을 confident novel로 부르던 것이 SQANTI-SIM에서 가짜 3개를 통과시켰기 때문임.
- `iso_A`의 junction은 short read 9개가 정렬됐지만 STAR의 non-canonical 필터, strand 0, motif 조건이라는 세 관문에 모두 걸림. 셋 다 "참조 게놈에서 non-canonical"이라는 한 사실에서 나옴. canonical로 만드는 염기는 intron 안에 있어서 read로는 볼 수 없음.
- 기준은 TOML로 바꿀 수 있음. 실제로 쓴 기준값은 `provenance.log`의 `thresholds` 줄에 남음. `ruleset_version`은 설정 파일에 적힌 이름일 뿐이라 기준값을 바꿔도 따라 바뀌지 않음.

다음 노트 [03](03_artifact_mechanisms.md)에서는 두 번째 축인 artifact 흔적을 본다. SQANTI3가 적어 준 QC 값과 long-read BAM에서 무엇을 읽는지다.

## 연습문제

### 문제 1

> `iso_C`의 두 번째 junction(941–1000)에 uniquely mapped read는 1개지만 multi-mapped read가 40개 있다고 해 보자. 확인된 것으로 칠까?

<details>
<summary>풀이</summary>

치지 않는다. 8열(`n_multi`)을 40으로 바꿔서 돌렸다.

```bash
A="--classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf --ref-gtf data/reference.gtf"
awk 'BEGIN{FS=OFS="\t"} $2==941 && $3==1000 {$7=1; $8=40} 1' data/short_reads.SJ.out.tab > out/multi.SJ.tab
grep -P '\t941\t1000\t' out/multi.SJ.tab
$PIG adjudicate $A --sj-tab out/multi.SJ.tab --out-prefix out/multi 2>/dev/null
grep -P '^iso_C\t' out/multi.adjudicated.tsv | cut -f1,5,7,9
```

```text
chrT	941	1000	1	1	0	1	40	50
iso_C	PARTIAL	LOW_CONF_PARTIAL	1
```

그대로 `PARTIAL`이다. multi-mapped read는 반복 서열이나 paralog처럼 비슷한 서열이 여러 곳에 있을 때 생긴다. 그런 read가 이 자리의 junction을 걸쳤다는 보장이 없고, 오히려 이 junction 자체가 반복 서열 때문에 생긴 정렬 artifact일 가능성을 시사한다. 그래서 확인 증거로 쓰지 않는다.

</details>

### 문제 2

> 시료 두 개의 `SJ.out.tab`을 그냥 이어 붙여서 하나로 줬더니, 같은 junction이 두 줄(read 1개, 2개)로 들어갔다. PanIsoGuard는 어느 줄을 쓸까?

<details>
<summary>풀이</summary>

어느 한 줄을 고르지 않고 합친다. 같은 염색체, strand, 좌표의 줄이 다시 나오면 `n_uniq`와 `n_multi`를 더한다(`SjTable::add()`).

```bash
(cat data/short_reads.SJ.out.tab; printf 'chrT\t941\t1000\t1\t1\t0\t2\t0\t44\n') > out/merged.SJ.tab
grep -P '\t941\t1000\t' out/merged.SJ.tab
$PIG adjudicate $A --sj-tab out/merged.SJ.tab --out-prefix out/merged 2>&1 | grep 'short-read'
grep -P '^iso_C\t' out/merged.adjudicated.tsv | cut -f1,5,7,9
```

```text
chrT	941	1000	1	1	0	1	0	50
chrT	941	1000	1	1	0	2	0	44
read 10 short-read junctions
iso_C	SUPPORTED	HIGH_CONF_NOVEL	2
```

11줄을 읽었지만 junction은 10개로 세고, 941–1000의 read는 1 + 2 = 3개가 되어 기준을 넘는다. 여러 시료를 합칠 때 편리한 동작이지만, 거꾸로 보면 "시료 하나에서 3개"와 "시료 셋에서 1개씩"을 구별하지 않는다는 뜻이기도 하다.

</details>

## 더 깊이 보기

<details>
<summary>SQANTI-SIM의 short-read 표는 왜 현실적이지 않나</summary>

[benchmark/sqanti_sim](../benchmark/sqanti_sim)의 `truth.SJ.tab`은 short read를 정렬해서 만든 것이 아니라 정답에서 바로 만든 표다. 시뮬레이션한 전사체의 intron을 모두 적고, read 수는 전부 50, motif는 전부 1(canonical), strand는 전사체의 strand로 채웠다(`prep_truth_sj.py`). 그래서 5절 같은 일이 생기지 않는다. 진짜 junction은 빠짐없이 확인되고, 가짜 junction은 하나도 확인되지 않는다. 이런 표로 잰 정확도는 short read 쪽이 완벽하다는 가정 위의 숫자다. [07](07_evaluation.md)에서는 short read를 실제로 시뮬레이션하고 STAR로 정렬해서 다시 쟀는데, 결론은 바뀌지 않았다.

</details>

<details>
<summary>STAR의 strand 0은 언제 나오나</summary>

STAR 매뉴얼은 4열을 "strand (0: undefined, 1: +, 2: -)"로 설명한다. stranded library가 아니면 STAR는 intron 양 끝의 motif로 strand를 추론한다. `GT/AG`면 +, 그 역상보인 `CT/AC`면 -다. non-canonical motif는 어느 쪽인지 정할 근거가 없어서 0이 된다. 이 노트의 예제에서도 `501 730` 줄이 strand 0으로 나왔다. PanIsoGuard의 `SJ.tab` reader는 0을 `.`(모름)으로 읽고, `.`인 줄은 `+`나 `-`인 isoform과 절대 짝지어지지 않는다.

</details>
