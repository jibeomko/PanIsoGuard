# PanIsoGuard 공부 노트

PanIsoGuard를 돌리면 novel isoform마다 `confidence_class`, `primary_mechanism`, `rule_trace`가 나온다. 그런데 이 값들이 입력 파일에서 어떤 규칙을 거쳐 나오는지 막상 설명하려니 막막했다. 그래서 작은 예제 유전자 하나를 만들어서, 입력을 하나씩 더해 가며 판정이 어떻게 바뀌는지 직접 계산해 봤다. 이 폴더는 그 기록이다. 저장소의 문서([docs/](../docs/))와 코드([src/](../src/))를 따라가며 적었고, 문서에 적힌 설명은 실제로 돌려서 확인했다.

long-read RNA-seq과 isoform이라는 말은 알지만, SQANTI3의 범주나 splice junction의 좌표 규칙은 처음 보는 사람을 떠올리며 썼다.

## 읽는 순서

처음이라면 00부터 순서대로 읽는 것을 추천한다. 노트마다 질문 하나로 시작해서 예제로 먼저 보여 주고, 규칙은 그다음에 나온다. 코드 세부나 부수적인 확인은 노트 끝 "더 깊이 보기"에 접어 두었으니 처음 읽을 때는 건너뛰어도 괜찮다.

모든 노트가 같은 예제를 쓴다. 1,700 bp짜리 가상 염색체 `chrT`의 유전자 하나(GENE_T)이고, 참조 전사체 세 개와 caller가 부른 isoform 아홉 개가 있다([00](00_overview.md) 1절, 입력은 [data/](data/)). 주인공은 novel junction 하나를 가진 `iso_A`인데, 입력을 하나씩 더할 때마다 판정이 이렇게 바뀐다.

| 더한 입력 | `iso_A`의 판정 | 노트 |
|---|---|---|
| SQANTI3 분류표 + caller GTF | `AMBIGUOUS` | [01](01_intron_chain_and_novelty.md) |
| + 참조 GTF | `AMBIGUOUS` (novel junction 1개) | [01](01_intron_chain_and_novelty.md) |
| + short-read `SJ.out.tab` | `ARTIFACT` (지지 없음 + non-canonical) | [02](02_short_read_support.md), [04](04_projection.md) |
| + long-read BAM | `ARTIFACT` (mapping 흔적 없음) | [03](03_artifact_mechanisms.md) |
| + 개인 haplotype FASTA (출처 unknown) | `AMBIGUOUS` (rescue 보류) | [05](05_reference_bias.md) |
| + 출처 `wgs` | `PAN_REF_RESCUED_FALSE_NOVEL` | [05](05_reference_bias.md) |

[08](08_one_isoform_end_to_end.md)에서는 이 표를 PanIsoGuard 없이 Python 표준 라이브러리만으로 처음부터 다시 계산한다.

## 노트 목록

| 노트 | 궁금했던 것 | 관련 문서 |
|---|---|---|
| [00. PanIsoGuard는 novel isoform 하나에 무엇을 물을까?](00_overview.md) | 무엇을 판정하고, 입력에서 판정까지 어떤 단계를 거칠까 | README, architecture |
| [01. novel isoform은 정확히 무엇이 새로울까?](01_intron_chain_and_novelty.md) | intron chain, 파일마다 다른 좌표 규칙, novel junction을 세는 방법 | input_formats |
| [02. short read가 novel junction을 봤다는 것을 어떻게 확인할까?](02_short_read_support.md) | `SJ.out.tab`, 지지 수준 네 가지, `PARTIAL` 규칙 변경, `iso_A`가 지지를 못 받는 이유 | decision_engine |
| [03. artifact는 어떤 흔적을 남길까?](03_artifact_mechanisms.md) | SQANTI3의 QC 값, long-read BAM의 네 흔적, 기전 우선순위 | decision_engine, bam_axis |
| [04. 지지 수준과 artifact 흔적을 어떻게 등급 하나로 합칠까?](04_projection.md) | 2축 조합표, 조합표 앞의 예외들, `ablate`, 규칙 전체의 Python 재현 | decision_engine |
| [05. 참조 게놈이 이 사람과 다르면 무엇이 달라질까?](05_reference_bias.md) | haplotype과 pangenome으로 reference bias를 찾는 방법, circularity firewall, 검증의 한계 | decision_engine, validation |
| [06. caller 여러 개가 같은 isoform을 찾았다면 믿어도 될까?](06_multi_caller_consensus.md) | intron chain fingerprint, `combine`, caller 합의의 강점과 약점 | relationship_to_merge_tools |
| [07. 이 판정이 맞는지는 어떻게 잴까?](07_evaluation.md) | 시뮬레이션 정답, precision/recall/AUPRC, "baseline"의 뜻, 한 줄 규칙과의 비교 | validation, sqanti3_filter_h2h |
| [08. isoform 하나는 어떻게 등급 하나가 될까?](08_one_isoform_end_to_end.md) | `iso_A`를 입력에서 등급까지 손으로 계산하고 PanIsoGuard와 맞춰 보기 | algorithm |

[check_notes.py](check_notes.py)는 예제의 isoform 아홉 개를 입력 조합 아홉 가지로 PanIsoGuard 없이 다시 판정해서, PanIsoGuard의 결과와 하나라도 다르면 오류를 내며 멈춘다. 노트에 나온 주요 숫자도 함께 확인한다. CTest(`integration_study_notes`)로 CI에서도 돌아가서, 판정 규칙이 바뀌면 노트가 낡았다는 것을 알려 준다.

## 확인 방법

PanIsoGuard 0.0.4(이 저장소에서 빌드)로 돌렸다. Python 쪽은 SQANTI3 6.0.1이 설치하는 conda 환경(Python 3.11.15, numpy 1.26.4, scikit-learn 1.5.2, matplotlib 3.10.9, pysam 0.22.1, samtools 1.23.1, STAR 2.7.11b)을 그대로 썼다. `check_notes.py`는 표준 라이브러리만 쓴다.

노트의 코드 블록은 모두 실제로 실행했고, 바로 아래 출력은 그때 나온 그대로다. bash 블록은 `notes/` 폴더에서 `PIG=../build/panisoguard`로 두고 돌렸고, 출력 파일은 `notes/out/`에 쓴다(저장소에는 넣지 않음). 코드는 노트마다 위에서부터 한 세션에서 이어서 돌렸다. 다 쓴 뒤에는 `notes/out/`을 지우고 모든 노트의 코드 블록을 처음부터 다시 돌려서, 출력이 노트에 적힌 것과 글자 하나까지 같은지 비교했다.

두 종류의 블록은 따로 준비가 필요하다.

- **STAR를 쓰는 블록**([02](02_short_read_support.md) 5절): STAR가 있어야 한다.
- **SQANTI-SIM 블록**([01](01_intron_chain_and_novelty.md) 4절, [03](03_artifact_mechanisms.md) 1절, [04](04_projection.md) 5절, [07](07_evaluation.md) 1절): [benchmark/sqanti_sim/run.sh](../benchmark/sqanti_sim/run.sh)가 만드는 작업 폴더(약 260 MB)가 있어야 한다. 저장소에는 넣지 않았다. 이 폴더에서 계산한 결과 중 저장소에 들어 있는 것([benchmark/results/](../benchmark/results/))은 데이터 없이 읽을 수 있다.

[02](02_short_read_support.md) 4절의 `$OLD`는 `v0.0.3` 태그를 따로 빌드한 바이너리다.

## 문서와 다르게 나온 것들

저장소 문서의 설명은 대부분 실행 결과와 맞았다. 다르게 나온 것 중 고친 것과 아직 남은 것을 모아 둔다.

**이 노트를 쓰면서 고친 것**

- [docs/decision_engine.md](../docs/decision_engine.md)는 BAM의 indel-near 비율이 판정에 쓰이지 않는다고 적고 있었음. 0.0.4부터 쓰임([03](03_artifact_mechanisms.md) 3절). `src/core/adjudicator.cpp`의 주석도 같은 내용이 낡아 있었음.
- 같은 문서가 "걸린 기전 조건은 모두 `rule_trace`에 나열된다"고 적고 있었음. 실제로는 첫 번째 기전만 남고, 나머지는 JSON의 `evidence`에만 있음([03](03_artifact_mechanisms.md) 4절).
- 같은 문서의 `ablate --without consensus`는 없는 옵션임. `ablate --axes consensus`가 맞음. `rule_trace`의 `UNKNOWN` 이유 이름(`axis_absent` → `short_read/catalog_axis_absent`)과 circularity 표시 이름(`circularity_risk` → `circularity_flag`)도 코드와 달랐음.
- [docs/validation.md](../docs/validation.md)에는 "junction을 쓰는 read가 alt allele을 갖고 있는지"로 reference-bias rescue를 독립 검증하자는 제안이 있었음. motif를 만드는 염기는 늘 intron 안에 있어서 그 read에는 나타나지 않으니 불가능함. phasing으로 바꿔서 실제로 해 봤고, 판정 가능한 28개 중 27개가 맞았음([05](05_reference_bias.md) 4, 6절).

**노트를 쓴 뒤 코드에서 고친 것**

- 0.0.4의 `PARTIAL` 규칙 변경에도 `ruleset_version`이 `builtin-0.0.1` 그대로였음. `0.0.2`로 올리고, 실제로 쓴 기준값 전부를 `provenance.log`의 `thresholds` 줄에 남기게 함([02](02_short_read_support.md) 6절).
- isoform 이름이 분류표와 caller 파일에서 어긋나거나([01](01_intron_chain_and_novelty.md) 연습문제 2), `--caller-support` 표의 이름이 어긋나도([06](06_multi_caller_consensus.md) 연습문제 1) 경고가 없었음. 이제 짝을 찾지 못한 수를 경고함.
- `examples/tiny`의 기대 출력이 0.0.4 버전 표기 때문에 어긋나 있었는데 아무 테스트도 돌리지 않아 몰랐음. 버전 줄은 비교에서 빼고 CTest(`example_tiny`)에 넣음.

**0.0.4에서 이미 고친 것**(이 노트의 [07](07_evaluation.md)과 같은 분석에서 나옴)

- README의 "AUPRC 0.970 vs 0.831 baseline"에서 0.831은 무작위 순위의 AUPRC(양성 비율)였음.
- [benchmark/sqanti_sim](../benchmark/sqanti_sim)의 README는 ISM 154개를 놓친 것을 "정답 라벨의 문제"라고 적고 있었음. 이 154개는 지운 전사체와 chain이 정확히 같은 진짜임.
- `PARTIAL` 지지를 `MEDIUM_CONF_NOVEL`로 올리던 규칙이 가짜 3개를 통과시켰음([02](02_short_read_support.md) 4절).

**아직 남아 있는 것**(고치지 않고 기록만 함)

- `SJ.tab` 좌표가 1 bp 밀려도 경고 없이 판정만 바뀜([01](01_intron_chain_and_novelty.md) 연습문제 1).
- BAM 비율 기준에 spanning read 수의 최솟값이 없어서, read가 적은 junction에서는 read 하나로 판정이 갈림([03](03_artifact_mechanisms.md) 연습문제 2).
- reference-bias rescue는 조합표보다 먼저 판정을 끝내서 mapping 흔적을 보지 않음([08](08_one_isoform_end_to_end.md) 연습문제 1). 실제 데이터에서도 이 때문에 정렬 artifact 하나가 구제됨([05](05_reference_bias.md) 6절). 다만 mapping 흔적으로 rescue를 막으면 맞는 rescue 하나가 막혀서, 규칙은 바꾸지 않음.
- SQANTI3는 junction을 연결한 유전자 안에서만 찾고 PanIsoGuard의 catalog은 유전자를 구분하지 않아서, SQANTI-SIM에서 한 건의 NIC 분류가 어긋났음([01](01_intron_chain_and_novelty.md) 4절).
