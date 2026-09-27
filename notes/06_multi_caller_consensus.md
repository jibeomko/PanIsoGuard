# 06. caller 여러 개가 같은 isoform을 찾았다면 믿어도 될까?

short read가 없으면 PanIsoGuard의 novel isoform은 거의 모두 `AMBIGUOUS`다([04](04_projection.md)). 그런데 같은 long-read 데이터에 isoform caller를 여러 개(FLAIR, IsoQuant, Bambu 등) 돌려 보는 일은 흔하다. 서로 다른 알고리즘이 같은 isoform을 찾았다면 그만큼 믿을 만하지 않을까? 이 노트에서는 caller들의 결과를 intron chain으로 맞춰 합치는 `panisoguard combine`, 그 결과를 판정에 넣는 `--caller-support`를 따라가 본다. 그리고 caller의 합의가 어떤 가짜는 잡아내지 못하는지도 예제로 확인한다.

> 관련 문서: [docs/relationship_to_merge_tools.md](../docs/relationship_to_merge_tools.md), [benchmark/multicaller](../benchmark/multicaller), [benchmark/merge_comparison](../benchmark/merge_comparison) · 코드: `src/core/fingerprint.cpp`, `src/core/consensus.cpp`, `src/cli/cmd_combine.cpp`

## 1. 서로 다른 caller의 isoform을 어떻게 같은 것으로 알아볼까?

예제에는 두 번째 caller의 결과(`data/caller2.gtf`)가 있다. isoform 네 개인데, 셋은 caller 1의 isoform과 intron chain이 같고 첫 exon의 시작이나 마지막 exon의 끝만 다르다. 나머지 하나(`c2_4`)는 caller 2만 찾았다.

`combine`은 isoform마다 **fingerprint**라는 64비트 숫자를 계산해서, 이 숫자가 같은 isoform끼리 묶는다. fingerprint는 염색체 이름, strand, 정렬한 intron 좌표를 이어 붙인 바이트열에 FNV-1a라는 간단한 해시 함수를 적용한 값이다. intron만 넣으니 양 끝(TSS, TES)이 달라도 intron chain이 같으면 fingerprint가 같다([01](01_intron_chain_and_novelty.md) 1절). 같은 계산을 Python으로 옮겨서 두 caller의 isoform을 묶어 봤다.

```python
import struct

def chains(path):
    tx = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            t = f[8].split('transcript_id "')[1].split('"')[0]
            tx.setdefault(t, [f[0], f[6], []])[2].append((int(f[3]), int(f[4])))
    return {t: (c, s, sorted(e)) for t, (c, s, e) in tx.items()}

def fingerprint(chrom, strand, introns):
    """FNV-1a 64-bit over chrom, strand and the sorted introns (as in src/core/fingerprint.cpp)."""
    h = 1469598103934665603
    data = chrom.encode() + strand.encode() + b"".join(struct.pack("<qq", s, e) for s, e in sorted(introns))
    for byte in data:
        h = ((h ^ byte) * 1099511628211) % 2**64
    return h

groups = {}
for caller in ("caller1", "caller2"):
    for t, (c, s, ex) in chains(f"data/{caller}.gtf").items():
        introns = [(ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1)]
        groups.setdefault(fingerprint(c, s, introns), []).append(f"{caller}:{t} {ex[0][0]}-{ex[-1][1]}")
for fp, members in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    print(f"{fp:020d}  {members}")
```

```text
04486241254994410266  ['caller1:iso_known 101-1400', 'caller2:c2_1 121-1380']
14203394687731003782  ['caller1:iso_B 101-1400', 'caller2:c2_2 131-1400']
13332283987921940502  ['caller1:iso_G 101-1400', 'caller2:c2_3 101-1400']
12649984371697548413  ['caller1:iso_ism 401-1400']
10280539423540787832  ['caller1:iso_A 101-1400']
10954041290435748490  ['caller1:iso_C 101-1400']
01762008240599116803  ['caller1:iso_D 101-1400']
10304000854924771737  ['caller1:iso_E 101-1400']
16185630266840774612  ['caller1:iso_F 101-1400']
16905696038402022381  ['caller2:c2_4 101-1400']
```

`iso_known`(101–1400)과 `c2_1`(121–1380)은 양 끝이 다르지만 같은 묶음이다. 묶음은 모두 10개다. intron은 0-based 64비트 정수 두 개(`start`, `end`)를 little-endian으로 이어 붙였는데, C++ 코드가 `int64_t` 값을 메모리에 있는 그대로 해시하기 때문이다. 그래서 C++ 쪽 값은 이 코드를 x86 같은 little-endian 기계에서 돌렸을 때만 이 숫자와 같다("더 깊이 보기"에서 C++로 직접 확인했다).

64비트 해시도 아주 드물게 서로 다른 chain에 같은 값을 줄 수 있다(충돌). `combine`은 fingerprint가 같을 때 intron 목록까지 비교해서, 다르면 따로 묶는다(`ConsensusBuilder::add()`). intron이 없는 단일 exon isoform은 intron chain으로 구별할 수 없으니 하나씩 따로 둔다.

## 2. `combine`은 무엇을 적어 줄까?

```bash
PIG=../build/panisoguard; mkdir -p out
$PIG combine --gtf caller1:data/caller1.gtf --gtf caller2:data/caller2.gtf --ref-gtf data/reference.gtf \
  --out out/matrix.tsv 2>&1 | tail -2
column -t out/matrix.tsv
```

```text
integrated 10 isoforms  (multi-caller: 3, known: 1, novel: 9)
wrote out/matrix.tsv
pig_id      chrom  strand  n_introns  n_callers  n_isoforms  callers          novelty  native_ids
PIG.000001  chrT   +       5          1          1           caller1          novel    caller1=iso_C
PIG.000002  chrT   +       4          2          2           caller1,caller2  known    caller1=iso_known;caller2=c2_1
PIG.000003  chrT   +       4          2          2           caller1,caller2  novel    caller1=iso_G;caller2=c2_3
PIG.000004  chrT   +       4          2          2           caller1,caller2  novel    caller1=iso_B;caller2=c2_2
PIG.000005  chrT   +       4          1          1           caller2          novel    caller2=c2_4
PIG.000006  chrT   +       4          1          1           caller1          novel    caller1=iso_A
PIG.000007  chrT   +       3          1          1           caller1          novel    caller1=iso_E
PIG.000008  chrT   +       4          1          1           caller1          novel    caller1=iso_F
PIG.000009  chrT   +       2          1          1           caller1          novel    caller1=iso_D
PIG.000010  chrT   +       3          1          1           caller1          novel    caller1=iso_ism
```

1절의 Python 묶음과 똑같이 10개로 묶였다. 줄마다 묶음 하나이고, `n_callers`는 그 intron chain을 찾은 caller 수, `native_ids`는 각 caller에서의 원래 이름이다.

`novelty` 열은 조심해서 읽어야 한다. 이 열은 "참조 전사체와 intron chain 전체가 **똑같은가**"만 본다. 그래서 TX1과 같은 `iso_known`만 `known`이고, TX1의 앞부분이 잘린 `iso_ism`(SQANTI3의 ISM)도, 알려진 intron만 조합한 `iso_D`(NIC)도 모두 `novel`이다. SQANTI3의 범주나 [01](01_intron_chain_and_novelty.md)의 novel junction 수와 다른 기준이니, 이 열로 NIC와 NNC를 세면 안 된다.

## 3. caller 합의는 판정을 어떻게 바꿀까?

`combine`의 표를 `--caller-support`로 주면, PanIsoGuard는 SQANTI3 분류표의 isoform 이름으로 `native_ids`에서 자기 줄을 찾아 `n_callers`를 가져온다. 그리고 short read 지지 수준이 `UNKNOWN`일 때만 이 값을 쓴다. 2 이상이면 흔적이 없을 때 `MEDIUM_CONF_NOVEL`, 흔적이 있을 때 `LOW_CONF_PARTIAL`이다([04](04_projection.md) 3절). short read 없이, BAM도 없이 돌린 결과와 BAM을 더한 결과를 나란히 놓았다.

```bash
B="--classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf --ref-gtf data/reference.gtf --caller-support out/matrix.tsv"
$PIG adjudicate $B --out-prefix out/cs 2>/dev/null
$PIG adjudicate $B --bam data/long_reads.bam --out-prefix out/csb 2>/dev/null
paste <(cut -f1,6,7 out/cs.adjudicated.tsv) <(cut -f6,7 out/csb.adjudicated.tsv) | column -t
grep '"iso_B"' out/cs.attribution.jsonl | python3 -c 'import json,sys; r=json.loads(sys.stdin.read()); print(r["evidence"]["n_callers"], r["rule_trace"])'
```

```text
isoform_id  primary_mechanism  confidence_class   primary_mechanism  confidence_class
iso_A       noncanonical       AMBIGUOUS          noncanonical       AMBIGUOUS
iso_B       none               MEDIUM_CONF_NOVEL  none               MEDIUM_CONF_NOVEL
iso_C       none               AMBIGUOUS          none               AMBIGUOUS
iso_D       none               AMBIGUOUS          none               AMBIGUOUS
iso_E       none               AMBIGUOUS          none               AMBIGUOUS
iso_F       noncanonical       AMBIGUOUS          noncanonical       AMBIGUOUS
iso_G       none               MEDIUM_CONF_NOVEL  mapping_or_repeat  ARTIFACT
iso_ism     none               LOW_CONF_PARTIAL   none               LOW_CONF_PARTIAL
iso_known   none               HIGH_CONF_KNOWN    none               HIGH_CONF_KNOWN
2 ['novelty-support UNKNOWN reason=short_read/catalog_axis_absent', 'no artifact mechanism flagged -> mechanism=none', 'caller-consensus n_callers=2 >= 2 corroborates the novel chain (no short-read axis) -> consensus-supported', 'project(UNKNOWN,none) -> MEDIUM_CONF_NOVEL']
```

caller 두 개가 찾은 `iso_B`와 `iso_G`가 `MEDIUM_CONF_NOVEL`로 올라갔다. caller 하나만 찾은 나머지는 그대로 `AMBIGUOUS`다.

문제는 `iso_G`다. [03](03_artifact_mechanisms.md)에서 본 대로 `iso_G`는 junction 근처의 indel을 aligner가 intron으로 잘못 그린 가짜다. 그런데 caller 두 개가 모두 이 가짜를 찾았다. 두 caller가 **같은 정렬 파일**에서 isoform을 만들었다면 당연한 일이다. 정렬에서 생긴 오류는 그 정렬을 읽는 caller 모두에게 똑같이 보인다. caller의 합의는 알고리즘의 독립성을 보여 줄 뿐, 데이터의 독립성을 보여 주지 않는다. BAM을 더하자(오른쪽 두 열) mapping 흔적이 잡히면서 `iso_G`만 `ARTIFACT`로 떨어졌다. mapping 흔적이 합의보다 앞서는 이유가 이것이다.

같은 이유로 합의는 short read 확인(`HIGH`까지 가능)보다 약하게, `MEDIUM`까지만 올린다. short read는 다른 시퀀싱 방법, 다른 분자에서 온 증거지만, caller 합의는 같은 read를 다르게 읽은 결과다.

## 4. 실제로는 caller 몇 개가 동의해야 할까?

저장소의 [benchmark/multicaller](../benchmark/multicaller)에는 SQANTI-SIM chr22 데이터(정답을 아는 시뮬레이션)에 caller 다섯 개(FLAIR, IsoQuant, Bambu, ESPRESSO, TALON)를 **한 정렬 파일로** 돌리고 `combine`으로 합친 결과가 있다. 저장소에 들어 있는 결과 파일을 읽었다.

```python
import json
m = json.load(open("../benchmark/results/multicaller/metrics.json"))["metrics"]
print("isoforms found by exactly k of 5 callers (SQANTI-SIM chr22):")
for k, s in m["stratification"].items():
    print(f"  k={k}: genuine {s['genuine']:4d}  false {s['false']:4d}  precision {s['precision']:.3f}")
print("keep isoforms found by >= k callers:")
for p in m["pr_curve"]:
    print(f"  >={p['min_callers']}: precision {p['precision']:.3f}  recall {p['recall']:.3f}  F1 {p['f1']:.3f}")
```

```text
isoforms found by exactly k of 5 callers (SQANTI-SIM chr22):
  k=1: genuine   13  false 1092  precision 0.012
  k=2: genuine   47  false  116  precision 0.288
  k=3: genuine  245  false    9  precision 0.965
  k=4: genuine  126  false    2  precision 0.984
  k=5: genuine  329  false    3  precision 0.991
keep isoforms found by >= k callers:
  >=1: precision 0.383  recall 0.988  F1 0.553
  >=2: precision 0.852  recall 0.971  F1 0.908
  >=3: precision 0.980  recall 0.910  F1 0.944
  >=4: precision 0.989  recall 0.592  F1 0.740
  >=5: precision 0.991  recall 0.428  F1 0.598
```

caller 하나만 찾은 novel isoform 1,105개 중 진짜는 13개(1.2%)뿐이다. 둘이 찾으면 29%, 셋 이상이면 96% 이상이 진짜다. 기준을 "k개 이상"으로 두고 F1을 보면 3개 이상에서 가장 좋다. 이 기준은 caller 수에 따라 달라진다. PanIsoGuard의 기본값 2(`consensus_min_callers`)는 caller가 다섯 개라면 precision 0.852로 느슨한 편이고, caller가 두 개뿐이면 3으로 올릴 수 없다(연습문제 2).

정렬 파일을 같이 썼는데도 합의가 이렇게 잘 가려내는 이유는, 이 데이터에서 가짜의 대부분이 정렬 오류가 아니라 caller마다 다른 판단(어떤 read 묶음을 isoform으로 볼지)에서 나왔기 때문으로 보인다. 가짜가 어디서 생겼는지 하나하나 확인하지는 않았다. `iso_G` 같은 정렬 오류형 가짜는 이 표의 k=3–5 줄에 섞여 남을 수 있다.

이 방법은 PanIsoGuard가 새로 만든 것이 아니다. [benchmark/merge_comparison](../benchmark/merge_comparison)에서 `combine`의 묶음이 `gffcompare -i`와 정확히 같다는 것을 확인했고, 여러 caller의 합의로 거르는 것도 이미 쓰이는 방법이다. PanIsoGuard는 그 합의를 판정 규칙의 한 축으로 넣었을 뿐이다.

## 정리

- `combine`은 isoform을 fingerprint(염색체, strand, 정렬한 intron에 대한 64비트 FNV-1a 해시)로 묶음. 양 끝은 보지 않고, 해시가 같아도 intron 목록이 다르면 따로 묶음.
- `combine` 표의 `novelty`는 참조 전사체와 intron chain 전체가 똑같은지만 봄. ISM과 NIC도 `novel`로 나옴.
- `--caller-support`는 short read 지지 수준이 `UNKNOWN`일 때만 쓰이고, caller 2개 이상이면 `MEDIUM_CONF_NOVEL`(흔적이 있으면 `LOW_CONF_PARTIAL`)까지 올림.
- 같은 정렬 파일을 쓴 caller들은 정렬 오류를 함께 부름. `iso_G`는 caller 두 개가 동의했지만 BAM의 mapping 흔적으로만 잡혔음.
- SQANTI-SIM에서 caller 5개 중 하나만 찾은 novel isoform의 precision은 0.012, 셋 이상은 0.98이었음.

다음 노트 [07](07_evaluation.md)에서는 지금까지 본 판정이 정답과 얼마나 맞는지 재는 방법, 그리고 그 숫자를 읽을 때 조심할 점을 본다.

## 연습문제

### 문제 1

> `combine`에 준 caller 1 GTF의 전사체 이름이 SQANTI3 분류표의 이름과 다르면(예: `iso_B` 대신 `c1_B`) 어떻게 될까?

<details>
<summary>풀이</summary>

전사체 이름을 바꾼 GTF로 `combine`을 돌리고, 원래 분류표로 판정했다.

```bash
sed 's/transcript_id "iso_/transcript_id "c1_/' data/caller1.gtf > out/caller1_renamed.gtf
$PIG combine --gtf caller1:out/caller1_renamed.gtf --gtf caller2:data/caller2.gtf --ref-gtf data/reference.gtf \
  --out out/matrix_renamed.tsv 2>/dev/null
grep -P '\tcaller1,caller2\t' out/matrix_renamed.tsv | cut -f1,5,9
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
  --ref-gtf data/reference.gtf --caller-support out/matrix_renamed.tsv --out-prefix out/cs_renamed 2>&1 | grep 'caller-support'
grep '"iso_B"' out/cs_renamed.attribution.jsonl | python3 -c 'import json,sys; r=json.loads(sys.stdin.read()); e=r["evidence"]; print(r["confidence_class"], "consensus_evaluable", e["consensus_evaluable"], "n_callers", e["n_callers"])'
```

```text
PIG.000002	2	caller1=c1_known;caller2=c2_1
PIG.000003	2	caller1=c1_G;caller2=c2_3
PIG.000004	2	caller1=c1_B;caller2=c2_2
read caller-support for 13 native isoform id(s) (consensus axis)
WARNING: 7 of 7 novel SQANTI3 isoform id(s) (e.g. "iso_A") are not in the --caller-support matrix; check that the ids match exactly.
AMBIGUOUS consensus_evaluable False n_callers None
```

`combine`의 묶음은 그대로다. 하지만 PanIsoGuard는 분류표의 이름 `iso_B`로 `native_ids`를 찾는데, 표에는 `c1_B`만 있다. 그래서 `iso_B`는 합의 정보가 없는 것(`consensus_evaluable = false`)으로 처리되어 `AMBIGUOUS`로 남는다. "13 native isoform id(s)"는 표에서 읽은 이름 수일 뿐이라, 이 노트를 처음 쓸 때는 짝이 맞지 않는다는 것을 알 수 없었다. 지금은 novel isoform 7개 모두 표에 없다고 경고한다. `combine`에는 SQANTI3에 넣은 것과 같은 caller GTF를 줘야 한다.

</details>

### 문제 2

> 합의 기준(`consensus_min_callers`)을 3으로 올리면 caller가 두 개뿐인 이 예제는 어떻게 될까?

<details>
<summary>풀이</summary>

```bash
sed 's/^consensus_min_callers *= *2/consensus_min_callers      = 3/' ../config/rules.default.toml > out/min3.toml
grep '^consensus_min_callers' out/min3.toml
$PIG adjudicate --classification data/caller1_classification.txt --isoforms-gtf data/caller1.gtf \
  --ref-gtf data/reference.gtf --caller-support out/matrix.tsv --config out/min3.toml --out-prefix out/cs3 2>/dev/null
grep -P '^iso_(B|G)\t' out/cs3.adjudicated.tsv | cut -f1,7
```

```text
consensus_min_callers      = 3      # a novel chain recovered by >= this many callers is
iso_B	AMBIGUOUS
iso_G	AMBIGUOUS
```

caller가 두 개뿐이니 어떤 isoform도 기준을 넘지 못하고, 합의 축이 사실상 꺼진다. 4절에서 caller 다섯 개일 때 3이 가장 좋았다고 해서 기준을 그대로 옮겨 쓰면 이렇게 된다. 기준은 쓰는 caller 수에 맞춰 정해야 한다.

</details>

## 더 깊이 보기

<details>
<summary>fingerprint 값을 C++로 직접 확인하기</summary>

`combine`은 fingerprint 값을 출력하지 않아서, 1절의 Python 값이 C++ 값과 정말 같은지는 묶음 결과로만 알 수 있다. 확실히 하려고 `fingerprint.cpp`만 따로 컴파일해서 `iso_known`과 `iso_A`의 intron chain을 넣어 봤다(PanIsoGuard를 빌드한 conda 환경의 `g++`).

```bash
cat > out/fp_check.cpp <<'EOF'
#include <cstdio>
#include "panisoguard/fingerprint.hpp"
using namespace panisoguard;
int main() {
  IntronChain a{"chrT", Strand::kPlus, {{200, 400}, {500, 700}, {800, 1000}, {1100, 1300}}};  // iso_known
  IntronChain b{"chrT", Strand::kPlus, {{200, 400}, {500, 730}, {800, 1000}, {1100, 1300}}};  // iso_A
  std::printf("iso_known %020llu\niso_A     %020llu\n",
              (unsigned long long)fingerprint_intron_chain(a), (unsigned long long)fingerprint_intron_chain(b));
}
EOF
g++ -std=c++17 -I../include out/fp_check.cpp ../src/core/fingerprint.cpp -o out/fp_check && ./out/fp_check
```

```text
iso_known 04486241254994410266
iso_A     10280539423540787832
```

1절의 Python 값과 같다.

</details>

<details>
<summary><code>combine</code>의 known 표시와 충돌</summary>

`combine`의 `novelty = known`은 catalog에 같은 fingerprint가 있는지로만 정한다(`Catalog::has_chain()`). 여기서는 intron 목록을 다시 비교하지 않으니, 아주 드물게 해시 충돌이 나면 새 isoform이 `known`으로 표시될 수 있다. 코드 주석은 실제 catalog 크기에서 이 확률을 약 1e-9로 적고 있다. 이 표시는 `combine` 표에만 쓰이고, `adjudicate`의 novel junction 판정은 intron 하나하나를 좌표로 비교하니 영향을 받지 않는다.

</details>
