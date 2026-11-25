# QC pipeline reproducibility contract

## Confirmed transcriptomics outlier-removal branch (2026-11-23)

Two biological subjects are confirmed manual exclusions: `Ng_MOG1` and
`Ng_HD3`. Every CSF/PB row and every `__CELL__...` row derived from those
subjects is excluded. The matrix entering the rerun keeps the selected feature
set and stored expression values unchanged. PyWGCNA may then remove features
through its default `goodSamplesGenes` checks because the retained sample set
has changed; those removals belong to preprocess, not outlier removal.

The pre-removal variance-filtered results remain immutable. Each network gains
an additive `outlier_removal/` branch containing sample-selection records, the
ten before-preprocess QC outputs, `preprocess(cut=inf)` diagnostics, and an
explicitly pending selected-cut gate.

The dendrogram line is the exact midpoint of the global largest difference
between consecutive sample-linkage heights. This intentionally simple rule does
not auto-apply the cut. A warning is attached when retaining a unique largest
cluster is ambiguous or when that interpretation would remove more than 5% of
rows. The suggestion is review evidence only: no finite cut or downstream
analysis is run before approval.

Review decisions are recorded separately as `Use selected cut`, `Hold`, or
`Do not use`. The proposed value can be edited before confirmation and exported
to CSV; browser-local decisions do not alter expression matrices.

```powershell
.\.venv\Scripts\python.exe results_20261119\code\src\run_transcriptomics_outlier_removed_diagnostic.py --run-root results_20261119
```

## Cross-omics module-pair overlapping genes (2026-11-21)

모든 transcriptomics–proteomics module pair에 대해 weighted Jaccard 값뿐 아니라 실제
공통 유전자와 gene별 기여값을 저장한다.

- 결과 파일: `all_pair_gene_contributions.csv.gz`
- 식별 열: `transcript_module`, `proteomics_module`, `gene`
- membership 열: `transcript_abs_kME`, `proteomics_abs_kME`
- gene별 numerator 기여: `weighted_intersection = min(RNA |kME|, protein |kME|)`
- gene별 shared-union 값: `weighted_union = max(RNA |kME|, protein |kME|)`
- proteomics gene mapping: ambiguous multi-symbol probe 제외 후, proper module 전체에서
  signed kME가 가장 큰 representative probe 하나를 선택
- unassigned feature: 양쪽 모두 제외

각 pair의 저장된 gene 행 수는 `module_pair_similarity.csv`의 `shared_genes`와 같아야
하며, `weighted_intersection` 합은 pair의 전체 weighted intersection과 같아야 한다.
단, weighted Jaccard denominator에는 한쪽 module에만 존재하는 gene도 포함되므로
공통 gene table의 `weighted_union` 합은 전체 pair union과 일반적으로 같지 않다.

Dashboard의 `Cross-omics Gene Overlap`은 network, RNA module, protein module을 선택해
공통 gene과 양쪽 |kME|를 확인하고 검색·정렬·CSV export할 수 있다. 이는 descriptive
module membership overlap이며 통계적 유의성 또는 질병 연관성 검정이 아니다.

## Dashboard artifact stage navigation (2026-11-19)

Cell type comparison과 Individual run의 산출물 선택은 두 단계로 구성한다. 첫 번째
단계에서 분석 stage를 선택하고, 두 번째 단계에는 해당 stage의 산출물만 표시한다.

- `Before Preprocess`: feature selection, before heatmap/PCA/PCA 3D
- `Preprocess`: cut diagnostic, selected cut
- `After Preprocess`: after heatmap/PCA/PCA 3D
- `Find Modules`: soft threshold, gene dendrogram, ME dendrogram, module sizes
- `Downstream`: pairwise-HD diagnosis, one-hot diagnosis audit, cross-omics similarity

특정 실험에서 생성되지 않은 산출물은 해당 stage와 세부 탭 모두 표시하지 않는다.
Stage를 바꾸면 그 stage에서 사용 가능한 첫 번째 산출물을 선택한다.

## Module member explorer contract (2026-11-19)

Inspection dashboard는 각 완료된 network의 공식
`06_kme_module_diagnosis/.../feature_module_kme.csv.gz`를 source of truth로 사용한다.
Dashboard builder는 이를 `dashboard/module-members/` 아래의 module별 JavaScript shard로
변환한다. 전체 feature table을 한 번에 browser memory로 읽지 않고, 사용자가 선택한
module 하나만 지연 로딩하기 위한 구조이다.

- 표시 식별자: `feature_id`, `gene_id`, `SYMBOL`, `PROBEID`, `UNIPROT`
- 표시 membership 값: signed kME, absolute kME, `hub_rank_within_module`
- 기본 정렬: module 내부 hub rank 오름차순
- 검색 범위: feature/gene/protein identifier 전체
- export: 현재 module과 검색 조건에 해당하는 전체 행을 CSV로 저장
- raw-unassigned 규칙: feature는 목록에 보존하지만 kME와 hub rank는 `NA`로 유지
- 결측 규칙: 정의되지 않은 kME를 0으로 대체하지 않음
- 해석 규칙: kME는 module membership의 크기와 방향이며 p-value나 질병 연관성이 아님

진입점은 Cell type comparison의 각 condition panel, Individual run 상단,
`Final proper modules`의 각 module 행, comparison card이다. `file://` 사용을 보장하기
위해 `fetch()`나 외부 CDN을 사용하지 않고 local `<script>` shard를 로드한다.

이 문서는 `results_20261119`에서 현재 확정된 파이프라인 순서, 계산 규칙, 저장 구조의 기준 문서다. 코드나 과거 결과와 설명이 충돌하면 먼저 이 문서를 기준으로 차이를 확인하고, 실제 재실행 전에는 해당 단계의 `run_manifest.json`과 현재 코드를 함께 검증한다.

- 마지막 정리일: 2026-11-19
- 현재 구현 범위: proteomics raw QC부터 `findModules` parameter sensitivity study까지, transcriptomics raw/batch-corrected QC, seed-0 terminal replicate filtering, gene filtering, 세 가지 variance/CPM feature-selection 비교, selected-cut preprocess·after-preprocess QC, 기본 `findModules`, ME 병합 후보 진단과 `MEDissThres=0.10` 실제 재계산까지
- 현재 중단 지점: transcriptomics variance-filter의 기본값과 `MEDissThres=0.10` 재계산 결과 비교·선택 전
- 실행하지 않은 단계: parameter-tuned 최종 network, module-trait association, enrichment 및 기타 downstream 해석
- raw proteomics 정책: QC까지만 수행하며 PyWGCNA `preprocess()`는 실행하지 않는다.
- raw transcriptomics 정책: cell-separated/cell-integrated QC까지만 수행하며 PyWGCNA `preprocess()`는 실행하지 않는다.

## 1. 기록과 변경 관리

재현성은 파일 이름만으로 보장하지 않는다. 다음 네 요소를 함께 보존한다.

1. `code/src/`: 실제 생성 코드
2. `PIPELINE_SPEC.md`: 현재 유효한 규칙과 저장 구조
3. `CHANGELOG.md`: 규칙이 바뀐 이유와 영향
4. 각 단계의 `run_manifest.json`: 입력·출력·코드 SHA-256, parameter, Python 및 package version

생성 규칙을 바꿀 때는 코드, 이 문서, `CHANGELOG.md`를 같은 작업에서 갱신한다. 결과를 다시 생성한 경우에는 새 `run_manifest.json`도 함께 생성한다. 문서가 같아도 입력 또는 코드 hash가 다르면 동일 실행으로 간주하지 않는다.

## 2. 실행 환경과 현재 재현 명령

프로젝트 루트는 `C:\Research\ProjectCode\gene-analysis`, 실행 결과 루트는 `results_20261119`다. 현재 확정된 범위는 아래 순서로 실행한다.

```powershell
.\.venv\Scripts\python.exe results_20261119\code\src\run_proteomics_raw_qc.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_proteomics_batch_corrected_qc.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_proteomics_preprocess_diagnostic.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_proteomics_selected_cut.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_proteomics_find_modules.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_transcriptomics_stage1.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_transcriptomics_selected_cuts.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_transcriptomics_replicate_filtered_qc.py --run-root results_20261119 --seed 0
.\.venv\Scripts\python.exe results_20261119\code\src\run_transcriptomics_gene_filtered_qc.py --run-root results_20261119
.\.venv\Scripts\python.exe results_20261119\code\src\run_transcriptomics_variance_filtered.py --run-root results_20261119
```

세 번째 명령은 `01_sample_selection`, `02_before_preprocess`, `03_preprocess/<specimen>/01_cut_diagnostic`까지 만든다. 네 번째 명령은 확정된 cut을 새 PyWGCNA 객체에 적용하고 `02_selected_cut`과 `04_after_preprocess`를 만든다. 다섯 번째 명령은 CSF/PB selected-cut 객체를 검증한 후 기본 조건에서 만든 power·TOM·gene tree를 고정하고 네 가지 module-cut 조건과 signed kME를 생성한다. `analyseWGCNA()`는 실행하지 않으며 최종 condition도 자동 선택하지 않는다.

의존성 선언은 프로젝트 루트의 `pyproject.toml`을 사용하고, 실제 실행 환경은 각 `run_manifest.json`에 기록한다.

## 3. 공통 데이터 계약

### 3.1 Expression과 metadata

- expression matrix는 sample이 행, `PROBEID`가 열인 numeric matrix다.
- assay-level `PROBEID`를 gene symbol로 합치지 않는다.
- metadata의 문자형 열은 expression에 섞지 않고 sample annotation 또는 PyWGCNA `sampleInfo`로만 전달한다.
- expression 영역의 실제 비결측 문자열은 숫자로 강제 대체하지 않고 오류로 처리한다.
- `NaN`은 결측값으로 유지하고, 무한대는 오류로 처리한다.
- 입력 파일 자체는 수정하지 않는다.

### 3.2 값 변환

- raw 단계는 저장된 raw intensity를 그대로 사용한다.
- batch-corrected 단계는 전달받은 ComBat-corrected 값을 그대로 사용한다.
- 두 단계 모두 추가 `log2`, normalization, 결측 대체를 하지 않는다.
- figure의 색상 표시 범위를 제한하더라도 저장된 matrix나 통계값은 절단하지 않는다.

### 3.3 Outlier의 의미

- QC에서 표시한 outlier는 검토 후보이며 자동 제거 대상이 아니다.
- 실제 sample 제거는 별도 의사결정과 기록이 있어야 한다.
- CIS 제외는 QC outlier 제거가 아니라 분석 대상 정의다.

## 4. Formatted data

PyWGCNA 입력용 CSV는 sample이 행, feature가 열인 형태로 `formatted/`에 둔다. 기존 `results_20261102/outputs/1_formatted_files`는 삭제하지 않았고, 값 변환 없이 아래 구조로 복사했다.

```text
formatted/
├── 00_raw/
│   ├── proteomics/pywgcna_input/
│   └── transcriptomics/
│       ├── pywgcna_input/
│       └── pywgcna_input_cell_integrated/
├── 01_batch_corrected/
│   ├── proteomics/pywgcna_input/
│   └── transcriptomics/
│       ├── threshold_5_cells/
│       └── threshold_10_cells/
└── manifests/
    ├── formatted_files.csv
    └── formatting_summary.json
```

`formatted/manifests/formatted_files.csv`는 파일 경로·크기·SHA-256을, `formatting_summary.json`은 경로 매핑과 변환 코드 hash를 기록한다.

변환 코드의 역할은 다음과 같다.

- `code/src/format_files.py`: proteomics와 cell type별 transcriptomics 변환
- `code/src/format_corrected_transcriptomics_all_fast.py`: corrected transcriptomics ALL matrix 변환
- `code/src/format_celltype_integrated.py`: `SampleID + specimen + cell type` 통합 matrix와 metadata 생성
- `code/src/write_formatted_manifest.py`: formatted 파일 manifest 생성

Formatted transcriptomics는 cell type별 matrix와 `SampleID + specimen + cell type` 통합 matrix를 모두 보존한다. 통합 ALL matrix는 metadata의 specimen 값으로 CSF/PB를 나누어 같은 디렉토리 아래에 생성하며, 값·gene 순서·row 순서는 바꾸지 않는다. 이 split과 downstream QC/preprocess 규칙은 아래 10절에서 확정한다.

## 5. Raw proteomics

### 5.1 입력과 feature 순서

- expression: `data/proteomics/cohorts/1_raw/*_DL_Feature_Matrix.csv`
- metadata: `data/metadata/Proteomics_metadata.csv`
- cohort 순서: `EMORY_7K`, `KIST_11K`, `KIST_7K`, `SU_5K7K`
- 각 cohort 파일에서 `ProteinColumn == Original_Sample_ID`인 alias 행은 expression feature에서 제외한다.
- union feature 순서는 위 cohort 순서대로 처음 등장한 `PROBEID` 순서를 보존한다.
- common feature는 네 cohort에 모두 존재하는 `PROBEID`다.

### 5.2 Descriptive statistics와 outlier

- sample/protein 통계는 raw intensity로 계산한다.
- descriptive outlier는 `cohort + specimen` 그룹 안에서 계산한다.
- 지표는 sample median, sample IQR, common-panel missing rate의 robust z-score다.
- robust z threshold는 `5.0`이다.
- 세 지표 중 두 개 이상이 threshold를 넘으면 `descriptive_outlier=True`다.
- 최종 `outlier_candidate = descriptive_outlier OR pca_outlier`다.

### 5.3 Heatmap

- 열은 union `PROBEID` 저장 순서를 사용한다.
- 행은 `specimen → cohort → sample_id`로 정렬하며 specimen 순서는 `CSF`, `PB`, `Unknown`이다.
- 위 패널은 관측/결측 상태, 아래 패널은 raw intensity를 표시한다.
- intensity colormap은 `magma`, 결측 cell은 회색이다.
- 색상 범위만 전체 관측값의 1st–99th percentile로 제한한다.
- 두 패널에 `specimen / cohort` 라벨을 표시한다.
- cohort 경계는 흰색, specimen 경계는 청록색 선으로 표시한다.
- outlier는 heatmap 위에 덮지 않고 오른쪽 별도 annotation 영역에 빨간 사각형과 sample ID로 표시한다.
- ALL heatmap과 CSF/PB 전용 heatmap을 각각 저장한다.

### 5.4 PCA

- `ALL`, `CSF`, `PB`를 각각 독립적으로 적합한다.
- 각 범위에서 네 cohort 공통 panel을 출발점으로 사용한다.
- 해당 범위에서 결측값이 하나라도 있는 feature는 feature 전체를 제외한다.
- 결측값 대체는 하지 않는다.
- 그다음 variance가 0인 constant feature만 제외한다.
- 남은 complete-case nonconstant feature를 모두 사용한다.
- MAD/variance ranking이나 top-N feature 선택을 하지 않는다.
- raw intensity에 feature별 평균 중심화만 적용하고 variance scaling 또는 z-score scaling은 하지 않는다.
- 최대 10개 PC score를 저장한다.
- PCA outlier distance는 `cohort + specimen` 내 PC1–PC5 robust z-score의 Euclidean distance다.
- 전체 distance의 `0.999` quantile 이상을 `pca_outlier=True`로 표시한다.

### 5.5 저장 구조

```text
experiments/00_raw/proteomics/
├── 01_descriptive_statistics/
│   ├── data_summary.csv
│   ├── intensity_distribution.png
│   ├── sample_statistics.csv
│   ├── protein_statistics.csv
│   └── outlier_candidates.csv
├── 02_data_distribution/
│   ├── heatmap.png
│   ├── pca_metadata.png
│   ├── pca_3d.png
│   ├── pca_3d_interactive.html
│   ├── pca_scores.csv
│   └── by_specimen/
│       ├── CSF/  # heatmap + PCA 4개 산출물
│       └── PB/   # heatmap + PCA 4개 산출물
└── run_manifest.json
```

Raw proteomics는 여기서 종료한다. 이 실험에 PyWGCNA `preprocess()` 또는 `findModules()`를 추가하지 않는다.

## 6. Batch-corrected proteomics

### 6.1 입력

- expression: `data/proteomics/cohorts/2_batch_corrected/*_CombatCorrected_DL_Feature_Matrix.csv`
- formatted ALL matrix: `formatted/01_batch_corrected/proteomics/pywgcna_input/ALL/PROT_ALL_expression_matrix.csv`
- metadata: `data/metadata/Proteomics_metadata.csv`
- 값은 전달받은 corrected value이며 batch correction 자체를 이 파이프라인에서 다시 적합하지 않는다.

### 6.2 전체 저장 구조와 상태

```text
experiments/01_batch_corrected/proteomics/
├── 00_feature_provenance/                         # 완료
├── 01_descriptive_statistics/                     # 보존된 CIS 제외 전 전체 QC
├── 02_data_distribution/                          # 보존된 CIS 제외 전 전체 QC
├── 01_sample_selection/                           # 완료: CIS 제외
│   └── matrices/{CSF,PB}/
├── 02_before_preprocess/                          # 완료: CIS 제외 후 QC
│   ├── ALL/{01_descriptive_statistics,02_data_distribution}/
│   ├── CSF/{01_descriptive_statistics,02_data_distribution}/
│   └── PB/{01_descriptive_statistics,02_data_distribution}/
├── 03_preprocess/
│   ├── CSF/
│   │   ├── 01_cut_diagnostic/                     # 완료: cut=inf
│   │   └── 02_selected_cut/                       # 완료: cut=50
│   └── PB/
│       ├── 01_cut_diagnostic/                     # 완료: cut=inf
│       └── 02_selected_cut/                       # 완료: cut=80
├── 04_after_preprocess/                           # 완료: selected-cut 이후 QC
│   ├── CSF/{01_descriptive_statistics,02_data_distribution}/
│   └── PB/{01_descriptive_statistics,02_data_distribution}/
├── 05_find_modules/                               # 완료: CSF/PB 독립 network
│   ├── CSF/
│   └── PB/
└── run_manifest.json
```

`01_descriptive_statistics/`와 `02_data_distribution/`은 CIS 제외 전 385개 sample 전체를 본 기존 baseline 결과로 보존한다. Downstream preprocess 효과를 비교할 때 사용하는 canonical before 결과는 `02_before_preprocess/<specimen>/`이며, preprocess 입력의 기준 matrix는 `01_sample_selection/matrices/<specimen>/`다.

### 6.3 Feature provenance

`00_feature_provenance/excluded_features_due_to_missingness.csv`에는 전달받은 corrected matrix에서 제외된 feature를 기록한다.

- raw union feature: 11,083개
- corrected matrix에 남은 feature: 5,284개
- 결측 때문에 제외된 feature: 5,799개
- panel 부재로 결측이 생긴 feature: 3,487개
- panel 내부 결측이 있는 feature: 2,312개
- 두 원인의 중복: 0개

이 목록은 raw union matrix에서 결측이 하나라도 있는 feature 집합과 corrected output에서 빠진 feature 집합이 정확히 일치함을 확인한 provenance 기록이다. Batch correction을 재실행한 결과가 아니다.

### 6.4 `01_sample_selection`: CIS 제외

- 규칙: metadata의 `Diagnosis`가 정확히 `CIS`인 sample을 downstream 분석에서 제외한다.
- 입력: 385 samples × 5,284 features
- 제외: 13 samples(CSF 7, PB 6)
- 유지: 372 samples(CSF 257, PB 115)
- 전달받은 batch correction은 CIS 제외 후 다시 적합하지 않았다. 따라서 CIS가 correction 적합에 포함되었을 가능성을 제한 사항으로 유지한다.

산출물:

- `sample_selection.csv`
- `excluded_cis_samples.csv`
- `retained_samples.csv`
- `selection_summary.json`
- `matrices/CSF/PROT_CSF_after_cis_exclusion.csv`
- `matrices/PB/PROT_PB_after_cis_exclusion.csv`
- `run_manifest.json`

### 6.5 `02_before_preprocess`: CIS 제외 후 QC

ALL, CSF, PB 각각에 대해 raw 단계와 같은 형식의 10개 산출물을 저장한다.

- `ALL`: formatted ALL matrix에서 CIS를 제외한 372 samples × 5,284 features. CSF/PB로 분리하기 전 specimen 차이를 포함한 전체 분포를 독립적으로 계산한다.
- `CSF`: CIS 제외 후 CSF 257 samples × 5,284 features
- `PB`: CIS 제외 후 PB 115 samples × 5,284 features

- `01_descriptive_statistics/`: `data_summary.csv`, `intensity_distribution.png`, `sample_statistics.csv`, `protein_statistics.csv`, `outlier_candidates.csv`
- `02_data_distribution/`: `heatmap.png`, `pca_metadata.png`, `pca_3d.png`, `pca_3d_interactive.html`, `pca_scores.csv`

CSF/PB 입력은 `01_sample_selection/matrices/<specimen>/`이다. ALL은 중복 checkpoint matrix를 저장하지 않고 formatted ALL matrix에 `sample_selection.csv`의 유지 목록을 적용하며, 두 입력의 hash를 ALL `run_manifest.json`에 기록한다. 모든 범위에서 추가 변환이나 결측 대체를 하지 않는다. 이 결과가 selected-cut 이후 QC와 비교할 기준이다.

### 6.6 `03_preprocess/<specimen>/01_cut_diagnostic`

CSF와 PB 각각에서 PyWGCNA `preprocess()` 전체를 `cut=inf`로 실행한다.

- `TPMcutoff=1`: PyWGCNA 기본값
- `goodSamplesGenes()`: monkey patch하지 않은 PyWGCNA 기본 구현
- 외부 결측 대체 또는 complete-case 재필터 없음
- distance: Euclidean
- linkage: average
- `cut=inf`: dendrogram cut에 의한 sample 제거 없음
- `findModules()`: 실행하지 않음

산출물:

- `matrix_before_cut.csv`: TPMcutoff와 `goodSamplesGenes()` 적용 후, cut 적용 직전 matrix
- `PROT_<specimen>_cut_diagnostic.p`: 재사용 가능한 PyWGCNA 객체
- `dendrogram_cut_diagnostic.png`
- `dendrogram_leaf_order.csv`
- `sample_linkage.csv`
- `sample_distance_metrics.csv`
- `outlier_candidates_before_cut.csv`
- `removed_features.csv`
- `removed_samples.csv`
- `preprocess_summary.json`
- `run_manifest.json`

`matrix_before_cut.csv`와 `.p` 객체의 `datExpr`는 값, sample 순서, feature 순서가 정확히 같은지 코드에서 검증한다.

현재 결과:

- CSF: 257 samples × 5,284 features, TPMcutoff 제거 0, `goodSamplesGenes()` 제거 0, diagnostic 후보 6
- PB: 115 samples × 5,284 features, TPMcutoff 제거 0, `goodSamplesGenes()` 제거 0, diagnostic 후보 5
- diagnostic 후보 기준: cohort 내 median pairwise Euclidean distance robust z ≥ 5
- 후보 action: `flagged_only`; sample을 제거하지 않음

### 6.7 `03_preprocess/<specimen>/02_selected_cut`

동료의 기존 proteomics 분석을 재현하기 위해 CSF는 `cut=50`, PB는 `cut=80`으로 확정했다. 이 숫자는 모든 데이터에 공통인 절대 기준이 아니다. 현재 corrected matrix, Euclidean distance, average linkage 조합에서 diagnostic dendrogram의 높이 범위가 동료 결과와 같고, 유지 sample 수도 정확히 재현되는 것을 확인했기 때문에 사용한다.

- CSF linkage 최대 높이: 92.345; `cut=50`
- PB linkage 최대 높이: 95.917; `cut=80`
- PyWGCNA `cut_tree`의 cluster label 0을 유지한다.
- 이번 실행에서는 label 0이 각 specimen의 가장 큰 cluster인지 별도로 검증했다.
- diagnostic 객체를 수정하지 않고 CIS 제외 checkpoint에서 새 PyWGCNA 객체로 `preprocess()` 전체를 다시 실행한다.
- selected-cut 직전 matrix는 `01_cut_diagnostic/matrix_before_cut.csv`와 sample, feature, 값 및 hash가 같아야 한다.
- `TPMcutoff=1`, patch하지 않은 `goodSamplesGenes()` 및 추가 변환 없음 규칙은 diagnostic과 동일하다.

결과:

- CSF: 257 → 254 samples; 제거 `PRJ-BBB-0014_CSF`, `PRJ-BBB-0372_CSF`, `PRJ-BBB-0404_CSF`
- PB: 115 → 113 samples; 제거 `PRJ-BBB-0358_PB`, `PRJ-BBB-0383_PB`
- 두 범위 모두 feature는 5,284개 유지
- TPMcutoff 또는 `goodSamplesGenes()`에 의한 추가 sample/feature 제거는 0개

산출물:

- `matrix_after_cut.csv`
- `PROT_<specimen>_selected_cut.p`
- `dendrogram_selected_cut.png`
- `dendrogram_leaf_order.csv`
- `sample_linkage.csv`
- `sample_cluster_assignments.csv`
- `removed_samples.csv`
- `removed_features.csv`
- `preprocess_summary.json`
- `run_manifest.json`

### 6.8 `04_after_preprocess`: selected-cut 이후 QC

CSF와 PB 각각에 대해 `matrix_after_cut.csv`를 입력으로 before-preprocess와 같은 10개 QC 산출물을 다시 계산한다.

- `01_descriptive_statistics/`: `data_summary.csv`, `intensity_distribution.png`, `sample_statistics.csv`, `protein_statistics.csv`, `outlier_candidates.csv`
- `02_data_distribution/`: `heatmap.png`, `pca_metadata.png`, `pca_3d.png`, `pca_3d_interactive.html`, `pca_scores.csv`
- 추가 log 변환, 결측 대체, top-N feature 선택을 하지 않는다.
- PCA는 각 specimen의 cut 이후 matrix에서 독립적으로 다시 적합한다.
- CSF 254개 및 PB 113개 sample과 5,284개 feature를 사용한다.

### 6.9 `05_find_modules`: CSF/PB parameter sensitivity study

CSF와 PB를 합치지 않고 각 specimen에서 다음 네 조건을 순차 비교한다. 첫 기본 조건에서 선택된 power, adjacency, TOM, gene tree를 같은 specimen의 모든 조건에서 고정하며 `deepSplit`, `minModuleSize`, `MEDissThres`만 변경한다.

```text
05_find_modules/<specimen>/
├── 01_default_ds2_mms50_med0p20/
├── 02_ds3_mms50_med0p20/
├── 03_ds3_mms30_med0p20/
├── 04_ds3_mms30_med0p10/
├── condition_summary.csv
├── condition_comparison.png
├── colleague_reference_comparison.json
└── selected_condition.json
```

- 입력 feature: 각 specimen의 5,284개 전체 protein
- 추가 filtering, 결측 대체, ranking, top-N 선택, 값 변환: 없음
- 고정 power: CSF 4, PB 3
- 고정 network: `networkType="signed hybrid"`, `TOMType="signed"`
- 고정 PAM 설정: `pamRespectsDendro=False`
- 변경 순서: `2/50/0.2` → `3/50/0.2` → `3/30/0.2` → `3/30/0.1`
- `04`는 동료 결과 재현 조건이지만 자동 선택하지 않는다.
- PyWGCNA color 이름은 unassigned 판정에 사용하지 않는다. `cutreeHybrid()` 반환값의 color 변환 전 `Name == 0`을 raw unassigned로 저장한다.
- 각 condition에서 dynamic label, 최종 module, ME, signed kME와 `.p` 객체를 독립적으로 저장한다.

결과:

- CSF 최종 color group: 7 → 10 → 15 → 19
- PB 최종 color group: 4 → 7 → 11 → 13
- `04`는 동료 결과와 module color 및 protein 수가 정확히 일치했다.
- CSF `04`: 19개 color group, raw unassigned `orangered` 6개, minModuleSize 이상 해석 대상 17개
- PB `04`: 13개 color group, raw unassigned `maroon` 5개, minModuleSize 이상 해석 대상 11개
- `selected_condition.json`은 `pending_review`이며 최종 condition을 아직 선택하지 않았다.

각 condition 산출물:

- `PROT_<specimen>_modules.p`
- `soft_threshold_table.csv`
- `dynamic_tree_labels.csv`
- `feature_module_kme.csv`
- `module_summary.csv`
- `module_eigengenes.csv`
- `analysis_summary.json`
- `run_manifest.json`
- `figures/soft_threshold.png`
- `figures/module_eigengene_dendrogram.png`
- `figures/gene_dendrogram_modules.png`
- `figures/module_sizes.png`
- `figures/module_eigengene_heatmap.png`
- `figures/assigned_kme_summary.png`

## 7. 다음 의사결정 gate

현재 파이프라인은 네 조건의 sensitivity study까지 완료하고 여기서 멈춘다. 다음 단계는 자동으로 이어서 실행하지 않는다.

1. `condition_summary.csv`와 condition별 module size, kME, ME 중복성을 검토한다.
2. 필요하면 sample resampling을 사용해 새로 분리된 module의 안정성을 검토한다.
3. 선정 근거를 `selected_condition.json`에 기록한다.
4. condition이 확정된 뒤에만 module-trait association과 enrichment로 진행한다.

## 8. Matrix checkpoint 규칙

모든 중간 계산을 중복 저장하지 않고, sample 또는 feature 집합이 실제로 바뀔 수 있는 경계만 저장한다.

1. CIS 제외 직후: `01_sample_selection/matrices/<specimen>/`
2. TPMcutoff와 `goodSamplesGenes()` 이후, cut 직전: `01_cut_diagnostic/matrix_before_cut.csv`
3. selected finite cut 이후: `02_selected_cut/matrix_after_cut.csv`와 `PROT_<specimen>_selected_cut.p`

Checkpoint는 CSF/PB별로 저장한다. 동일 내용을 다시 합친 ALL matrix는 만들지 않는다. CSV는 독립 확인과 이식성을 위한 것이고, `.p`는 PyWGCNA 내부 상태 재사용을 위한 것이다.

## 9. 아직 확정하지 않은 규칙

다음 항목은 현재 문서가 미리 정하지 않는다.

- tuned module parameter의 선택 기준과 최종값
- raw unassigned branch의 저장·보고 규칙
- module-trait association model과 contrast 규칙
- enrichment 및 기타 downstream 분석 규칙

이 항목은 의사결정 후 코드와 함께 이 문서에 추가하고 실행한다.

## 10. Transcriptomics stage-1 QC와 cut diagnostic

### 10.1 분석 단위

- cell type은 `B_cell`, `CD4_T_cell`, `CD8_T_cell`, `cDCs`, `Classical_Monocytes`, `Macrophages`, `MAIT_gDT_cell`, `NK_cell`, `Non_classical_Monocytes`, `other`, `pDCs`, `Plasmablast`, `Proliferating_T_cell`의 13종이다.
- `cell_separate`는 specimen × cell type을 독립 network로 취급한다.
- `cell_integrated`는 `SampleID + specimen + cell type` 행을 유지한 채 specimen별로 하나의 통합 network를 만든다. 통합 network 내부의 cell type은 expression feature로 넣지 않고 metadata/PCA annotation으로만 사용한다.
- raw는 13개 cell-separated matrix와 1개 cell-integrated matrix의 descriptive/distribution만 생성한다. 각 raw 분석에는 ALL distribution과 `by_specimen/CSF`, `by_specimen/PB` distribution을 함께 저장한다.
- batch-corrected는 `threshold_5_cells`, `threshold_10_cells`를 독립 조건으로 보존하며 서로 합치지 않는다.

### 10.2 QC 계산 규칙

- 저장된 값을 그대로 사용하고 추가 log, scaling, normalization, 결측 대체, MAD/top-N feature 선택을 하지 않는다.
- heatmap은 gene stored order를 유지하고 sample을 specimen → cohort → cell type → sample ID 순으로 정렬한다. 색 범위만 관측값의 1st–99th percentile로 제한하며 값은 절단하지 않는다.
- outlier 표시는 heatmap 오른쪽의 독립 band에 두어 expression 색을 가리지 않는다.
- PCA는 분석 범위에서 하나라도 비유한 값이 있는 gene 전체와 zero-variance gene을 제외한 뒤 mean-centering만 적용한다. feature 일부를 임의 선택하지 않는다.
- sample이 3개 미만이면 PCA를 추정하지 않는다. `pca_scores.csv`, 정적 PNG, interactive HTML에는 `not_run` 사유를 명시하고 descriptive/heatmap은 계속 생성한다.
- descriptive outlier와 PCA-distance outlier는 robust z ≥ 5 후보로만 표시하며 자동 제거하지 않는다.

각 descriptive root에는 `data_summary.csv`, `intensity_distribution.png`, `sample_statistics.csv`, `gene_statistics.csv`, `outlier_candidates.csv`를 저장한다. 각 distribution root에는 `heatmap.png`, `pca_metadata.png`, `pca_3d.png`, `pca_3d_interactive.html`, `pca_scores.csv`를 저장한다.

### 10.3 Batch-corrected 순서와 저장 구조

두 threshold와 두 mode에 아래 순서를 동일하게 적용한다.

```text
<threshold>/<mode>/
├── 00_feature_provenance/
├── 01_sample_selection/
├── 02_before_preprocess/
│   ├── CSF/[<cell-type>/]{01_descriptive_statistics,02_data_distribution}/
│   └── PB/[<cell-type>/]{01_descriptive_statistics,02_data_distribution}/
└── 03_preprocess/
    ├── CSF/[<cell-type>/]01_cut_diagnostic/
    └── PB/[<cell-type>/]01_cut_diagnostic/
```

대괄호의 `<cell-type>/` 하위 경로는 `cell_separate`에만 있고, `cell_integrated`는 specimen 바로 아래에 산출물을 둔다.

- `00_feature_provenance`는 raw 대비 corrected gene 누락/추가 목록을 저장한다. 현재 두 threshold 모두 누락·추가 gene은 0개다.
- `01_sample_selection`은 diagnosis가 정확히 `CIS`인 행을 제외하고 specimen/network별 matrix checkpoint와 전체 선택 표를 저장한다. 두 threshold 모두 고유 source sample 4개가 제외된다.
- `02_before_preprocess`는 CIS 제외 matrix에 10개 QC 산출물을 생성한다.
- `03_preprocess/.../01_cut_diagnostic`는 PyWGCNA 기본 `TPMcutoff=1`, 기본 `goodSamplesGenes()` 기준, `cut=inf`로 실행한다. finite cut과 `findModules`는 실행하지 않는다.
- PyWGCNA 2.2.1이 pandas 2.2 이상에서 반복 filtering 중 bool Series에 정수 count를 대입해 실패하는 부분은 integer count Series를 만드는 호환 shim으로만 교체한다. `minFraction=0.5`, `minNSamples=4`, `minNGenes=4`, variance 판정식과 반복 조건은 바꾸지 않는다. 기존 threshold-5 실행은 TPM cutoff 이후 `goodSamplesGenes()` 제거 gene이 0개여서 원본 구현으로 완료됐고, threshold-10에서 실제 반복 제거 경로가 발생해 shim이 필요했다.
- sample이 4개 미만이면 기본값을 완화하지 않고 `not_run_insufficient_samples`를 기록한다. 현재 `threshold_10_cells / cell_separate / CSF / Non_classical_Monocytes`가 CIS 제외 후 2개 sample이라 유일하게 skip되었다.

완료된 diagnostic에는 PyWGCNA `.p`, `matrix_before_cut.csv`, `removed_features.csv`, `removed_samples.csv`, `sample_linkage.csv`, `dendrogram_leaf_order.csv`, `dendrogram_cut_diagnostic.png`, `sample_distance_metrics.csv`, `outlier_candidates_before_cut.csv`, `preprocess_summary.json`, `run_manifest.json`을 저장한다. 현재 56개 대상 중 55개가 완료됐고 1개가 명시적으로 skip되었다.

Heatmap 표시 규칙만 바뀐 경우 `code/src/refresh_transcriptomics_heatmaps.py`를 사용한다. 이 유틸리티는 matrix와 기존 `pca_scores.csv`를 읽어 heatmap만 다시 만들고, QC manifest에 `refresh_scope=heatmap only`를 기록해 다른 산출물이 재계산된 것으로 오해하지 않게 한다.

### 10.3 Reviewed selected-cut preprocess

cut 의사결정의 단일 입력은 `experiments/01_batch_corrected/transcriptomics/dendrogram/selected_cuts_review.csv`다. 모든 review 행은 `use` 또는 `reject`로 확정되어야 하며, `use` 행에는 유한한 양의 `selected_cut`이 있어야 한다.

- `use`: review된 cut 아래에서 유일한 최대 sample cluster를 유지한다.
- `reject`: sample cut을 적용하지 않고 `cut=inf`로 다시 preprocess한다.
- review CSV에 없는 network: diagnostic flagged count가 0이거나 표본 수 부족으로 diagnostic이 skip된 경우에만 no-cut으로 허용한다.
- 새 PyWGCNA 객체를 CIS 제외 checkpoint에서 만들고 `TPMcutoff=1` 및 기본 `goodSamplesGenes()`를 다시 실행한다.
- cut 직전 sample ID, gene ID, 값은 `01_cut_diagnostic/matrix_before_cut.csv`와 완전히 같아야 한다.
- 추가 log, scaling, normalization, 결측 대체, MAD/top-N 선택은 하지 않는다.
- finite cut에서 scipy cluster label 0이 최대 군집이 아닐 수 있으므로, PyWGCNA 2.2.1의 고정 label-0 보존이 정상 군집을 제거하지 않게 최대 군집 label과 0을 `preprocess()` 내부에서만 교환한다. 원래 label과 교환 여부를 결과에 기록한다.
- 산출물은 각 network의 `03_preprocess/.../02_selected_cut/`에 `matrix_after_cut.csv`, PyWGCNA `.p`, selected-cut dendrogram, linkage/leaf order, cluster assignment, 제거 sample/gene 목록, summary와 manifest로 저장한다.
- 전체 실행 감사표는 transcriptomics root의 `selected_cut_preflight.csv`, `selected_cut_preprocess_summary.csv`, `selected_cut_preprocess_manifest.json`에 저장한다.
- 현재 실행 결과는 56개 중 55개 완료, 표본 2개인 1개 network 명시적 skip이다. finite cut 20개에서 총 21개 row를 제거했고, 35개 완료 network는 no-cut으로 유지했다. PB/CD4 T cell의 두 threshold에서만 최대 군집 label이 1이어서 label 0과 교환했다.
- 이 단계에서는 after-preprocess QC와 `findModules`를 실행하지 않는다.

## 11. Transcriptomics terminal replicate filtering QC

`experiments/02_sample_filtered/transcriptomics`는 batch-corrected CIS 제외 checkpoint에서 시작한다. 기존 selected-cut matrix를 입력으로 사용하지 않는다. terminal replicate 선택은 dendrogram cut보다 먼저 수행되어야 하고, 표본 구성이 바뀌면 sample tree도 다시 계산해야 하기 때문이다.

### 11.1 Replicate 선택 규칙

- source sample ID 끝의 `CSF1`, `CSF2`, `CSF3` 또는 `PB1`, `PB2`처럼 specimen 뒤에 붙은 숫자만 제거해 replicate group을 만든다.
- `HD1`, `MS2`, `MOG1`, `NMO3` 등 diagnosis/sample 번호는 제거하지 않는다.
- group과 member ID를 정렬한 뒤 local `random.Random(seed)`를 사용하며 seed는 0으로 고정한다.
- 선택은 source sample 수준에서 한 번만 수행하고 두 threshold, cell-separated, cell-integrated에 공통 적용한다.
- cell-integrated에서는 선택되지 않은 source sample에 속한 모든 cell row를 제거한다.
- expression 값, gene 순서 및 유지된 row 값은 변경하지 않는다.

현재 결과는 CIS 제외 후 131개 source sample에서 반복 그룹 16개를 확인해 20개를 제외하고 111개를 유지했다. 제외 ID는 `results_20261022/.../2_random_sampled/removed_samples.csv`와 정확히 일치한다.

### 11.2 저장 구조와 중단 지점

```text
02_sample_filtered/transcriptomics/
├── 00_replicate_selection/
│   ├── selection_manifest.csv
│   ├── removed_samples.csv
│   ├── retained_samples.csv
│   ├── sampling_summary.json
│   └── run_manifest.json
├── threshold_5_cells/
│   ├── cell_separate/
│   └── cell_integrated/
├── threshold_10_cells/
│   ├── cell_separate/
│   └── cell_integrated/
├── qc_network_summary.csv
├── run_summary.json
└── run_manifest.json
```

각 threshold/mode 아래에는 `01_sample_selection/matrices`, `02_before_preprocess/<specimen>/[<cell-type>/]{01_descriptive_statistics,02_data_distribution}`를 둔다. 56개 네트워크 모두 filtered matrix와 10개 QC 산출물을 생성한다.

- stored expression value를 그대로 사용한다.
- 추가 log, scaling, normalization, 결측 대체, MAD/top-N 선택을 하지 않는다.
- PCA는 complete-case nonconstant gene 전체를 mean-center해 계산한다.
- `threshold_10_cells / cell_separate / CSF / Non_classical_Monocytes`는 표본이 2개라 PCA를 실행하지 않고 명시적 placeholder와 상태를 저장한다.
- 이 단계에서는 PyWGCNA `preprocess()`, dendrogram, sample cut, `.p` 객체, `findModules`를 생성하지 않는다.

## 12. Transcriptomics gene-filtered QC

`experiments/03_gene_filtered/transcriptomics`는 `02_sample_filtered`의 seed-0 matrix를 입력으로 사용한다. 따라서 적용 순서는 `CIS 제외 → terminal replicate sampling(seed=0) → gene filtering`이다. 이전 selected-cut preprocess 결과를 입력으로 사용하지 않는다.

### 12.1 필터 정의와 순서

2026-10-22의 `outputs/3_gene_filtering/transcriptomics`에 실제 저장된 다음 다섯 필터만 순서대로 적용한다.

1. `1_housekeeping`: `data/reference/gene_filtering/housekeeping/HK_genes.txt`의 첫 symbol token과 일치하는 gene
2. `2_mitocondrial`: reference workbook의 `10 Mitochondrial Genes` 열에 있는 10개 gene
3. `3_ribosomal`: `^RP[SL]\d+[A-Z]?$`
4. `4_non-coding`: 과거 실행과 동일한 non-coding family 정규식
5. `5_human_ig`: `^(IGHV|IGHJ|IGHD[0-9]|IGKV|IGKJ|IGLV|IGLJ)`

이 단계에는 TCR V/J/D 추가 필터를 넣지 않는다. 해당 필터는 2026-10-22의 저장된 5단계 결과에 포함되지 않았기 때문이다. 모든 네트워크가 동일한 gene ID와 순서를 가지는지 먼저 확인한 뒤, 각 필터를 직전 단계의 retained gene에 순차 적용한다.

확정된 gene 수는 다음과 같다.

- 입력: 36,617
- housekeeping 제거: 3,564 → 33,053
- mitochondrial 제거: 10 → 33,043
- ribosomal 제거: 58 → 32,985
- non-coding 제거: 14,554 → 18,431
- human immunoglobulin 제거: 386 → 최종 18,045

최종 retained gene 순서와 단계별 removed gene 행은 2026-10-22 결과와 exact match여야 한다. 일치하지 않으면 실행을 중단한다.

### 12.2 값 보존과 QC 규칙

- gene column만 제거하고 sample row, row 순서, retained gene 순서, stored expression value는 변경하지 않는다.
- 추가 log, normalization, scaling, 결측 대체, variance/MAD ranking, top-N 선택을 하지 않는다.
- 최종 matrix에 대해 기존 transcriptomics QC engine의 10개 산출물을 다시 생성한다.
- PCA는 필터 후 matrix에서 complete-case nonconstant gene 전체를 사용하고 mean-centering만 수행한다.
- outlier는 robust z 5 이상인 descriptive 또는 PCA-distance 후보로만 표시하며 자동 제거하지 않는다.
- 이 단계에서는 PyWGCNA `preprocess()`, dendrogram cut, `.p`, `findModules`를 실행하지 않는다.

### 12.3 저장 구조

```text
03_gene_filtered/transcriptomics/
├── 00_gene_filtering/
│   ├── 1_housekeeping/ ... 5_human_ig/
│   │   ├── removed_genes.csv
│   │   ├── retained_genes.csv
│   │   └── rule.json
│   ├── gene_filter_summary.csv
│   ├── all_removed_genes.csv
│   ├── final_retained_genes.csv
│   ├── gene_filter_assignments.csv
│   └── results_20261022_alignment.json
├── threshold_5_cells/
│   ├── cell_separate/
│   └── cell_integrated/
├── threshold_10_cells/
│   ├── cell_separate/
│   └── cell_integrated/
├── qc_network_summary.csv
├── run_summary.json
└── run_manifest.json
```

각 threshold/mode 아래에는 `01_gene_filtering/matrices`와 `02_after_gene_filtering/<specimen>/[<cell-type>/]{01_descriptive_statistics,02_data_distribution}`를 둔다. 다섯 단계의 중간 expression matrix는 중복 저장하지 않고, 공통 단계별 removed/retained gene 목록과 최종 56개 matrix를 저장한다. 공통 schema와 순서가 보장되므로 이 목록으로 모든 중간 단계를 재현할 수 있다.

현재 56개 네트워크 모두 최종 matrix와 10개 QC 산출물을 생성했다. 55개는 PCA가 완료됐고, 표본이 2개인 `threshold_10_cells / cell_separate / CSF / Non_classical_Monocytes` 1개만 명시적 `not_run` 상태다.

## 13. Transcriptomics CPM/variance feature-selection 비교

`experiments/04_variance_filtered/transcriptomics`는 `03_gene_filtered`의 18,045-gene matrix를 입력으로 사용한다. CIS 제외와 seed-0 terminal replicate 선택은 이미 입력에 반영되어 있다. threshold 5/10, cell-separated/integrated, CSF/PB의 동일한 56개 네트워크에 다음 세 전략을 각각 독립 적용한다.

1. `01_cpm_detection`: 원본 CPM workbook에서 `CPM >= 1`인 row가 `min(현재 row 수, max(4, ceil(현재 row 수의 10%)))` 이상인 gene을 유지한다. corrected value를 CPM으로 다시 계산하지 않는다. cell-integrated 분석은 cell type별로 이 조건을 평가하고 한 cell type에서라도 통과한 gene을 유지한다.
2. `02_hvg_top9000`: 각 네트워크의 저장된 batch-corrected value에서 population variance(`ddof=0`)가 큰 gene 9,000개를 유지한다. 동률은 기존 gene 저장 순서를 유지한다.
3. `03_hvg_top5000`: 동일한 variance 규칙으로 5,000개를 유지한다. proteomics network와 feature 규모를 비슷하게 맞춘 민감도 분석이다.

세 전략 모두 feature 선택 전에 추가 log, normalization, scaling 또는 결측 대체를 하지 않는다. CPM은 source CPM에서 detection support만 판정하고, 실제 QC와 PyWGCNA 입력값은 선택된 batch-corrected matrix의 stored value다.

### 13.1 Preprocess와 자동 sample cut

- 선택된 matrix로 before-preprocess descriptive 5개와 distribution 5개를 생성한다.
- PyWGCNA 2.2.1의 `preprocess()`를 먼저 `cut=inf`, 기본 `TPMcutoff=1`, 기본 `goodSamplesGenes()`로 실행한다.
- sample tree는 Euclidean distance와 average linkage로 계산한다.
- sample-distance robust z가 5 이상인 row를 cut 후보로 flag한다.
- 인접 linkage height gap의 절대 percentile과 상대 percentile이 모두 0.75 이상이고, 상대 gap이 0.05 이상인 지점만 검토한다.
- cut 아래에 유일한 최대 cluster가 있어야 하며, 제거율이 5% 이하여야 하고, 제거되는 모든 row가 distance flag여야 하며, metadata level을 소진하지 않아야 한다.
- 조건을 만족하는 후보 중 `min(absolute gap percentile, relative gap percentile)`이 가장 큰 지점의 gap midpoint를 selected cut으로 사용한다. 후보가 없으면 `cut=inf`를 유지한다.
- finite cut은 최대 cluster가 PyWGCNA 내부 label 0이 되도록 label만 교환한 뒤 새 객체에서 `preprocess()`를 다시 실행한다. expression value는 바꾸지 않는다.
- selected-cut matrix로 after-preprocess descriptive 5개와 distribution 5개를 다시 계산한다.
- 표본이 4개 미만이면 QC는 보존하고 PyWGCNA preprocess는 `not_run_insufficient_samples`로 명시한다.

### 13.2 저장 구조와 완료 상태

```text
04_variance_filtered/transcriptomics/
├── 00_shared/cpm_detection_cache/
├── threshold_5_cells/
│   ├── 01_cpm_detection/
│   ├── 02_hvg_top9000/
│   └── 03_hvg_top5000/
├── threshold_10_cells/
│   ├── 01_cpm_detection/
│   ├── 02_hvg_top9000/
│   └── 03_hvg_top5000/
├── cpm_source_audit.csv
├── variance_filter_preprocess_summary.csv
├── variance_filter_preprocess_comparison.png
└── run_manifest.json
```

각 network 아래에는 `01_feature_selection`, `02_before_preprocess`, `03_preprocess/01_cut_diagnostic`, `03_preprocess/02_selected_cut`, `04_after_preprocess`, `run_manifest.json`을 둔다. 큰 matrix는 gzip, 이후 생성한 PyWGCNA 객체는 zstd로 압축하며 이전에 완료된 객체 일부는 gzip 형식을 유지한다. Plotly HTML은 저장 공간을 줄이기 위해 CDN loader를 사용한다.

총 168개 전략-network 조합을 기록했다. 165개가 selected-cut preprocess까지 완료됐고, 동일한 2-sample network가 세 전략에서 각각 한 번씩 명시적으로 skip됐다. finite cut은 CPM 18개, top 9,000 19개, top 5,000 17개이며 각각 18, 20, 19개 sample row를 제거했다. CPM 선택 gene 수는 network에 따라 4,504–12,223개다.

### 13.3 기본 findModules와 module dendrogram gallery

- 입력은 각 network의 `03_preprocess/02_selected_cut` PyWGCNA 객체다. selected-cut matrix에 추가 feature filtering, log, normalization, scaling 또는 결측 대체를 적용하지 않는다.
- `findModules()`는 PyWGCNA 2.2.1 기본 power 후보 `[1,2,3,4,5,6,7,8,9,10,11,13,15,17,19]`에서 network별 soft threshold를 자동 선택한다.
- module parameter는 `networkType='signed hybrid'`, `TOMType='signed'`, `deepSplit=2`, `minModuleSize=50`, `MEDissThres=0.20`, `pamRespectsDendro=False`를 사용한다.
- PyWGCNA 2.2.1 호환을 위해 TOM 입력 adjacency만 쓰기 가능한 C-order NumPy copy로 전달한다. 값이나 계산식은 바꾸지 않는다.
- PyWGCNA가 단일 module 병합에서 `SystemExit`하는 경우 동일한 `MEDissThres=0.20`으로 ME 병합을 재현하고, 최종 proper module이 하나이면 표준화한 gene matrix의 첫 PC를 eigengene으로 저장한다. 이 fallback 여부는 `analysis_summary.json`에 기록한다.
- raw dynamic-tree label 0을 unassigned로 기록한다. 색 이름만으로 grey 계열을 임의 제외하지 않는다.
- dense adjacency, TOM과 전체 module 객체는 저장하지 않는다. 대신 module assignment, module summary, eigengene, soft-threshold 표, gene-tree linkage, parameter/provenance JSON과 figure를 저장한다. 동일 결과는 selected-cut preprocess 객체와 manifest로 재계산할 수 있다.
- 각 network의 출력 위치는 `05_find_modules/01_default_ds2_mms50_med0p20/`이다. `gene_dendrogram_modules.png`에는 전체 gene tree 아래에 dynamic cut과 ME 병합 후 module color strip을 함께 표시하며 이는 보조 진단 그림으로 유지한다.
- 주 ME dendrogram은 dynamic tree cut 직후, ME 병합 전 proper dynamic module마다 signed first PC를 계산한 뒤 `1 - Pearson correlation`을 average linkage로 군집화한다. raw dynamic label 0은 unassigned이므로 제외한다. `module_eigengene_dendrogram.png`에는 `MEDissThres=0.20` 가로선을 표시하며, 모든 y축은 이 기준선을 포함한다.
- dynamic module eigengene 값은 `dynamic_module_eigengenes.csv`, ME linkage는 `module_eigengene_linkage.csv`로 함께 저장한다.
- 전역 `module_eigengene_dendrogram_gallery.html`은 모든 ME dendrogram을 상대경로로 lazy-load하는 주 offline HTML이다. threshold, feature-selection 전략, cell mode, specimen, cell type, status로 필터링할 수 있고 dynamic/final module 수를 함께 표시한다. 기존 `module_dendrogram_gallery.html`은 gene dendrogram 보조 gallery로 유지한다.
- 전역 `find_modules_summary.csv`는 입력 차원, 선택 power, dynamic/final module 수, unassigned gene 수, 단일-module collapse 여부와 실행 상태를 한 행씩 기록한다.
- `validate_transcriptomics_variance_find_modules.py`는 168행 상태, 완료 결과별 14개 필수 파일, gene/ME linkage와 final/dynamic eigengene 차원, tested/selected power, 두 gallery의 record와 모든 상대 이미지 링크를 다시 검증하고 `find_modules_validation.json`을 저장한다.

완료 상태는 165개 `complete_findModules`, 표본이 2개인 동일 network의 전략별 3개 `not_run_preprocess_incomplete`, 0개 failure다. 완료 network의 final proper module 수는 1–51개이며, 9개 network가 ME 병합 후 하나의 proper module로 collapse했다. headless Edge로 `file://` HTML을 렌더링해 168개 카드와 로컬 이미지 표시도 확인했다.

### 13.4 ME 병합 후보 진단과 예측 파라미터 실제 재계산

- ME dendrogram으로 직접 추정할 수 있는 parameter는 module eigengene 병합 기준인 `MEDissThres`다. gene tree를 다시 자르는 `deepSplit`은 ME dendrogram만으로 예측하지 않는다.
- 후보 조건은 `deepSplit=2`, `minModuleSize=50`, `pamRespectsDendro=False`, 선택 power, adjacency, TOM, gene tree와 dynamic module을 기본 결과와 동일하게 유지하고 `MEDissThres`만 `0.20`에서 `0.10`으로 바꾼다.
- 이 조건에서는 adjacency/TOM과 dynamic tree cut을 다시 계산할 필요가 없다. 저장된 dynamic module에서 PyWGCNA 공식 `mergeCloseModules(cutHeight=0.10, iterate=True)`를 실행해 변경된 ME 병합 단계만 정확히 재계산한다.
- 재계산 전에 같은 dynamic module을 `MEDissThres=0.20`으로 병합해 기존 baseline partition과 color label을 모두 재현해야 한다. 어느 하나라도 다르면 해당 network는 후보 계산을 중단한다.
- 단순 dendrogram 절단으로 얻은 `predicted_diagnostic_module_count`는 진단값이다. 실제 `mergeCloseModules()`는 병합 후 ME를 다시 계산하며 반복하므로 `actual_final_module_count`가 예측보다 작을 수 있다. 최종 판단에는 실제 재계산값을 사용한다.
- 출력은 각 network의 `05_find_modules/02_predicted_ds2_mms50_med0p10/`에 baseline과 같은 compact artifact 구조로 저장한다. dense adjacency/TOM과 전체 model object는 중복 저장하지 않는다.
- 전역 결과는 `07_module_merge_diagnostics/actual_rerun_summary.csv`, `actual_rerun_analysis_summary.json`, `actual_rerun_errors.json`, `actual_rerun_manifest.json`에 저장하고 `07_module_merge_diagnostics/index.html`에서 baseline, 진단 예측, 실제 결과를 함께 비교한다.
- 완료된 165개 baseline 중 merge 후보로 선정된 120개를 실제 재계산했다. 120개 모두 baseline partition 및 exact color-label 재현을 통과했고 failure는 0개였다. 총 proper module 수는 baseline 1,774개, 단순 진단 예측 2,363개, 실제 iterative merge 2,339개였으며 baseline에서 collapse한 후보 4개는 실제 `MEDissThres=0.10` 결과에서 모두 해소됐다.

## 14. 실험 inspection dashboard

`dashboard/index.html`은 여러 실험 결과를 사람이 직접 찾고 비교하며 최종 판단을 기록하기 위한 offline inspection 화면이다. 결과 디렉토리를 HTML에 직접 하드코딩하면 실험이 추가될 때마다 화면 코드를 수정해야 하므로, 다음 세 층을 분리한다.

1. **Source adapter**: 실험 family 고유의 CSV·JSON·figure 경로를 읽는다.
2. **Normalized registry**: 각 실험을 공통 `record / metrics / artifacts / files / provenance` schema로 변환해 `dashboard-data.js`에 기록한다.
3. **Stable UI and review state**: `index.html`은 registry만 읽고 필터·검색·비교·검토 기능을 제공한다. 판단과 메모는 `localStorage`에 임시 보존하고 CSV로 export/import한다.

현재 예시는 `04_variance_filtered/transcriptomics`의 168개 feature-selection/findModules 조합을 연결한다. 실험 record에는 threshold, feature-selection 전략, cell mode, specimen, cell type, 상태, sample/gene 수, power, module 수, unassigned 비율, kME와 runtime을 포함한다. artifact에는 ME dendrogram, gene dendrogram, soft-threshold, module size, feature-selection, before/after heatmap·PCA, sample-cut figure를 연결한다. 원본 summary·assignment·eigengene·linkage·manifest 파일도 같은 record에서 열 수 있다.

기본 화면은 **Cell type comparison**이다. source와 cell type을 하나 고르면 해당 cell의 모든 threshold × feature-selection × specimen 조건을 유지한 채, 산출물 종류를 탭으로 전환한다. 한 탭에는 같은 산출물을 조건별 2열 panel로 모두 표시하고 `threshold → feature strategy → specimen` 순으로 자연 정렬해 같은 조건의 CSF/PB가 연속해서 나타나게 한다. 현재 transcriptomics record에서는 cell type당 최대 12개 조건과 16개 artifact 탭을 제공한다. 각 panel은 sample/gene 수, dynamic→final module 수와 review 상태를 표시하며 원본 artifact 또는 기존 단일-record inspection으로 이동할 수 있다.

검토자는 각 record에 `Use / Hold / Do not use`를 지정하고 soft threshold, module 병합·크기, sample outlier, downstream 준비 여부를 체크하며 근거 메모를 남긴다. 여러 실험은 최대 4개까지 comparison workspace에 올려 ME dendrogram과 핵심 수치를 함께 본다. 판단 완료 후 export한 review CSV를 의사결정 원본으로 보존하며, 브라우저 `localStorage` 자체는 공식 산출물로 간주하지 않는다.

새 실험 family를 추가할 때는 `code/src/build_experiment_inspection_dashboard.py`에 동일 schema를 반환하는 adapter만 추가하고 builder를 다시 실행한다. `index.html`의 화면 로직은 바꾸지 않는다. dashboard 자체는 외부 network 요청 없이 `file://`에서 동작해야 하며, 연결한 artifact의 상대경로 존재 여부와 headless browser 렌더링을 생성 후 확인한다.

## 15. 공식 PyWGCNA kME와 module-diagnosis 분석

`findModules`가 완료된 proteomics 8개 parameter condition과 transcriptomics 165개 variance-filter network에 다음 downstream 연산을 적용한다.

1. selected-cut PyWGCNA 객체의 저장 expression matrix를 읽는다. 추가 log, scaling, filtering 또는 결측 대체를 하지 않는다.
2. 저장된 final `module_color`와 final module eigengene을 feature ID와 sample ID로 정렬해 객체에 다시 연결한다. raw dynamic label 0으로 확인된 fully-unassigned eigengene은 정상적인 `MEgrey` 제외 의도와 동일하게 분석 입력에서 제외한다.
3. 공개 메서드 `PyWGCNA.WGCNA.CalculateSignedKME()`로 모든 feature × 모든 proper final ME의 signed Pearson correlation을 계산한다.
4. 공개 메서드 `PyWGCNA.WGCNA.module_trait_relationships_heatmap()`에 diagnosis trait와 `alternative='two-sided'`를 전달해 공식 module-trait heatmap을 생성한다.

PyWGCNA 2.2.1의 `CalculateSignedKME()`는 전달한 모든 feature 사이의 사용하지 않는 correlation block까지 함께 생성한다. 5,000–12,000 feature에서 O(p²) 메모리를 만들지 않도록 feature를 1,000개씩 독립적으로 전달해 공개 메서드를 반복 호출하고 feature-to-ME block만 결합한다. 각 feature의 ME 상관은 다른 feature의 포함 여부와 무관하므로 결과는 전체 호출과 동일하다. 기존에 저장된 assigned kME와의 최대 절대 차이는 `2.45e-15` 미만이어야 한다.

범주형 diagnosis는 PyWGCNA `getDatTraits()`가 생성한 one-hot 값을 사용한다. PyWGCNA 2.2.1과 pandas 3 조합에서 이 값이 object dtype으로 남아 SciPy `pearsonr`가 실패하므로 값은 바꾸지 않고 float dtype만 확정한 뒤 공식 heatmap 메서드의 numeric-trait 경로에 전달한다. heatmap cell에는 PyWGCNA 공식 형식대로 Pearson `r`과 보정 전 양측 p-value를 표시한다. BH-FDR은 network 내부의 module × diagnosis 검정 전체에 대해 별도 CSV에서만 제공하며 공식 heatmap annotation을 대체하지 않는다.

각 결과 디렉토리에는 다음을 저장한다.

- `feature_module_kme.csv.gz`: 모든 module에 대한 공식 kME, stored module kME, proper module assigned kME와 module 내 hub rank
- `module_kme_summary.csv`
- `module_diagnosis_correlations.csv`
- `module_diagnosis_pvalues.csv`
- `module_diagnosis_statistics.csv`: long form r, p, within-network BH q
- `diagnosis_trait_encoding.csv`
- `figures/module_diagnosis_heatmap.png`
- `analysis_summary.json`, `run_manifest.json`

Proteomics는 `01_batch_corrected/proteomics/06_kme_module_diagnosis/<specimen>/<condition>/`, transcriptomics는 각 network의 `06_kme_module_diagnosis/01_default_ds2_mms50_med0p20/`에 저장한다. 전역 summary는 각각 `summary.csv`와 `module_diagnosis_summary.csv`다.

검증 대상은 총 173개 network, 16,716,991개 공식 kME 값, 8,111개 module-diagnosis 검정과 173개 heatmap이다. raw-label-0 fully-unassigned eigengene을 제외한 proper module에서는 추정 불가능한 kME column이나 module-diagnosis 검정이 0개였다. non-estimable 항목이 생기는 향후 실행에 대비해 `experiments/module_diagnosis_non_estimable.csv`를 유지하고, 전체 검증 결과는 `experiments/module_diagnosis_validation.json`에 저장한다.

Inspection dashboard에는 proteomics 8개 condition을 새 source로 추가하고 transcriptomics를 포함한 각 record의 `Module-diagnosis` figure, 공식 kME, r/p/q table과 encoding/manifest를 연결한다. 전체 dashboard record는 176개이며 module-diagnosis figure를 기본 visual evidence로 표시한다.

## 16. Transcriptomics-proteomics kME-weighted Jaccard

`experiments/05_cross_omics_module_similarity/`에서는 완료된 transcriptomics network마다 같은 specimen의 proteomics module 전체를 비교한다. CSF transcriptomics는 CSF proteomics와, PB transcriptomics는 PB proteomics와만 비교한다. 실험 간 기준이 달라지는 것을 막기 위해 proteomics reference는 동료 결과 재현 조건인 `04_ds3_mms30_med0p10`으로 고정한다.

비교 가중치는 proper final module에 배정된 feature의 `absolute official assigned own-module kME`이다. raw dynamic label 0에서 유래한 fully-unassigned group은 양쪽 모두 제외한다. transcript module A와 proteomics module B의 similarity는 다음 fuzzy weighted Jaccard로 계산한다.

```text
w_A(g) = |kME_A(g)| if gene g is assigned to A, otherwise 0
w_B(g) = |kME_B(g)| if gene g is assigned to B, otherwise 0

J_w(A,B) = sum_g min(w_A(g), w_B(g)) / sum_g max(w_A(g), w_B(g))
```

gene universe는 각 transcriptomics 실험에서 측정 가능한 gene과 proteomics에서 측정 가능한 simple gene symbol의 교집합이다. 따라서 한 assay에서 애초에 측정되지 않은 gene 때문에 similarity가 낮아지지 않는다. gene symbol은 trim 후 uppercase로 정규화한다. `|`, `;`, `,`, `/`를 포함하는 다중 또는 모호한 proteomics symbol은 단일 gene 매핑에서 제외한다. 동일 gene에 proteomics probe가 여러 개이면 WGCNA `collapseRowsUsingKME`의 목적을 따라 proper module probe 중 모든 proper ME에 대한 최대 signed kME가 가장 큰 probe 하나를 대표로 선택하고, 동률은 `PROBEID`로 결정한다. 대표 probe의 기존 final module assignment는 바꾸지 않는다.

공통 universe에서 gene이 0개가 된 module도 사용자의 “모든 module” 비교 요구를 지키기 위해 heatmap 축과 pair table에 남긴다. 이 module의 non-empty module 대비 similarity는 0으로 계산하지만 best match나 reciprocal best로 표시하지 않는다. 각 heatmap은 transcriptomics module을 행, proteomics module을 열로 두며 모든 165개 figure에서 동일한 0–0.3133996733854167 color scale을 사용한다. `•`는 각 RNA module의 best protein match, `★`와 주황색 테두리는 양방향 reciprocal best를 뜻한다.

각 transcriptomics network 디렉토리에는 다음 산출물을 저장한다.

- `kme_weighted_jaccard_heatmap.png`
- `module_pair_similarity.csv`: 모든 RNA module × protein module 조합
- `transcript_module_best_matches.csv`
- `reciprocal_best_gene_contributions.csv.gz`
- `gene_mapping_audit.json`
- `analysis_summary.json`, `run_manifest.json`

전역 산출물은 `cross_omics_similarity_summary.csv`, `method.json`, `validation_summary.json`, `kme_weighted_jaccard_gallery.html`이다. proteomics representative gene table과 collapse audit는 `00_shared/proteomics_gene_collapse/`에 저장한다. 기존 `dashboard/index.html`에도 각 transcriptomics record의 기본 visual evidence와 관련 CSV 링크로 연결한다.

완료 상태는 transcriptomics network 165개, module pair 30,270개, heatmap 165개이다. validator는 모든 proper module set과 Cartesian pair 수, 0–1 값 범위, 공통 gene universe, gallery link, 전역 color scale을 검사하고, 165개 실험 각각의 최대 similarity pair를 원본 kME에서 독립적으로 다시 계산해 일치함을 확인했다. 공통 gene이 0개인 module 축 2개도 누락 없이 유지되었다.

이 값은 module membership profile의 기술적(descriptive) overlap이다. p-value나 방향성 있는 expression correlation이 아니며, WGCNA의 공식 `modulePreservation()` 통계량을 대신하지 않는다. 조건 간 차이는 후보 패턴을 찾는 용도로 해석하고, 정식 보존성 또는 생물학적 연관성 결론에는 별도 검정이 필요하다.

## 17. Direct disease-versus-HD module–diagnosis heatmap

module–diagnosis 분석의 주 결과는 PyWGCNA의 다범주 one-hot trait가 아니라 `MS vs HD`, `MOGAD vs HD`, `NMOSD vs HD` 직접 대조로 고정한다. PyWGCNA의 다범주 `getDatTraits()`는 각 diagnosis를 나머지 모든 diagnosis와 비교하므로 `MS` 행에 MOGAD와 NMOSD가 대조군으로 들어간다. 이는 이 프로젝트의 질환별 HD 대조 질문과 다르다.

각 직접 대조는 다음 순서로 계산한다.

1. 해당 target disease와 HD sample만 남기고 다른 두 질환 sample은 제외한다.
2. HD를 0, target disease를 1로 encoding한다.
3. 저장된 proper final module eigengene과 binary indicator의 `scipy.stats.pointbiserialr`를 계산한다. 이는 binary trait에 대한 Pearson correlation과 수학적으로 같다.
4. 양의 `r`은 target disease의 평균 ME가 HD보다 높음을, 음의 `r`은 낮음을 뜻한다.
5. HD 또는 target disease가 없거나 pairwise subset의 ME variance가 0이면 값을 0으로 대체하지 않고 `NA`로 유지한다.
6. 각 group이 3개 이상인 cell만 inference-ready로 두고, 해당 network 안의 inference-ready module × disease 검정에 BH-FDR을 적용한다.

heatmap은 3개 직접 대조를 행, final proper module을 열로 사용한다. cell에는 signed `r`와 two-sided unadjusted `p`를 표시하고, within-network BH `q<0.05`에는 `*`, group n<3에는 `†`를 붙인다. 열 label에는 module feature 수를 같이 표시한다.

기존 PyWGCNA one-hot 결과는 구현 감사와 비교를 위해 삭제하지 않고 `One-hot diagnosis (audit)`로 유지한다. 새 직접 대조 결과는 같은 `06_kme_module_diagnosis` 디렉토리에 다음 이름으로 추가한다.

- `figures/module_diagnosis_pairwise_hd_heatmap.png`
- `module_diagnosis_pairwise_hd_correlations.csv`
- `module_diagnosis_pairwise_hd_pvalues.csv`
- `module_diagnosis_pairwise_hd_qvalues.csv`
- `module_diagnosis_pairwise_hd_statistics.csv`
- `module_diagnosis_pairwise_hd_sample_counts.csv`
- `pairwise_hd_analysis_summary.json`
- `pairwise_hd_run_manifest.json`

재현 명령은 `code/src/run_pairwise_hd_module_diagnosis.py --omics all --overwrite`이다. 전역 summary는 proteomics의 `06_kme_module_diagnosis/pairwise_hd_summary.csv`와 transcriptomics의 `pairwise_hd_module_diagnosis_summary.csv`다. Inspection dashboard에서는 `Module-diagnosis (vs HD)`를 기본 visual evidence로 표시하고 기존 one-hot figure는 보조 audit tab으로 연결한다. Windows `file://`에서 260자를 넘는 원본 transcriptomics 경로가 깨지는 것을 막기 위해 dashboard builder는 heatmap을 `dashboard/assets/pairwise_hd/<path-hash>.png`로 mirror한다. 원본 분석 산출물과 manifest 경로는 바꾸지 않는다.

## 18. Pairwise-HD offline GitHub release

전체 `results_20261119`는 약 41GB이며 100MB 초과 원본 matrix와 PyWGCNA 객체를 포함하므로 GitHub 배포 대상으로 사용하지 않는다. `code/src/build_pairwise_hd_offline_release.py`가 `results_20261119_pairwise_hd_offline/`에 별도의 privacy-minimized 배포본을 생성한다.

배포본에는 완료된 173개 network의 `Module-diagnosis (vs HD)` heatmap, module-level r/p/q와 sample count, analysis summary, sanitised manifest, dashboard, 생성 코드와 문서만 포함한다. raw/formatted expression matrix, 개인 sample ID가 있는 sample-level table, `.p/.p.gz/.p.zst` 객체와 전체 feature-level kME table은 포함하지 않는다.

dashboard는 모든 asset과 table을 저장소 내부 상대경로로 연결하며 CDN, `fetch()`, web server가 필요하지 않다. 각 파일의 SHA-256은 `SHA256SUMS.txt`에 기록하고 `python verify_release.py`로 clone 이후 무결성을 확인한다. `.gitattributes`의 `* -text` 규칙으로 Windows clone 시 CRLF 자동 변환이 hash를 바꾸지 않게 한다.

private 저장소는 `https://github.com/sumincho0830/gene-analysis-pairwise-hd-offline`이다. 다른 컴퓨터에서는 clone 후 repository root의 `index.html`을 직접 열어 사용한다.

## 19. Inspection dashboard visual system

`dashboard/index.html`의 visual layer는 `results_20261015/design-system.md`를 기준으로 한다. 데이터 adapter, normalized record schema, filter, search, comparison, review, CSV import/export와 artifact 연결 방식은 바꾸지 않고 다음 표시 규칙만 적용한다.

- Arial/Helvetica/맑은 고딕 계열, 흰 배경, `#262626` 본문, 중립 회색 보조 텍스트를 사용한다.
- 카드형 배경, 큰 radius, 강한 그림자는 제거하고 `#ddd`, `#eee`의 가는 선으로 영역을 나눈다.
- 현재 record와 현재 artifact는 배경색보다 2px 선과 밑줄을 우선해 표시한다.
- KPI는 장식 카드가 아니라 한 줄 요약 표처럼 표시하고 숫자는 tabular numeral로 정렬한다.
- filter control은 흰 배경, 가는 테두리, inset chevron을 사용한다.
- metric, review, figure tab, module table, provenance accordion은 같은 간격과 선 체계를 공유한다.
- 일반 본문은 16px을 기준으로 하고, 자주 읽는 control과 label은 13–14px, 부차 metadata만 12px대로 제한한다.
- 최대 콘텐츠 폭은 1,560px이다. 820px 이하에서는 record list와 detail을 세로로 쌓고, 520px 이하에서는 filter와 metric을 1–2열로 전환한다.
- 인쇄 시 filter, record list, review control과 comparison workspace를 숨기고 선택한 결과 본문만 남긴다.

이 visual layer는 `code/src/build_experiment_inspection_dashboard.py`의 `DESIGN_SYSTEM_CSS`에 보관하고 builder가 생성 HTML의 기존 stylesheet 뒤에 한 번만 주입한다. 따라서 dashboard를 재생성해도 디자인과 사용자 기능이 함께 유지된다. 브라우저 100%와 200% 배율에서 2열/단일열 전환, cross-omics 기본 figure, module table과 provenance 목록을 확인했다.

## 20. 기존 PCA outlier candidate 표시와 전역 집계

PCA outlier 규칙과 저장된 `pca_outlier` 또는 `pca_outlier_candidate` 값은 변경하지 않는다. PCA figure는 이 PCA 전용 후보만 빨간색 빈 원으로 덧그리며 metadata color와 marker를 그대로 보존한다. descriptive-only 후보는 PCA figure에 추가하지 않는다. 저장된 PCA 후보의 source sample ID를 표시하고 각 metadata panel과 3D figure에는 `PCA outlier candidate` legend를 추가한다. Interactive PCA 재생성 시에도 동일 후보를 별도 trace로 표시한다. 별도의 top-N 라벨 선택은 적용하지 않는다.

`code/src/refresh_pca_outlier_annotations.py`는 기존 `pca_scores.csv`만 읽으며 PCA 좌표나 candidate flag를 다시 계산하지 않는다. Proteomics와 transcriptomics의 `pca_metadata.png`, `pca_3d.png`를 모두 갱신하고, `--interactive`를 지정한 경우에만 `pca_3d_interactive.html`도 갱신한다. PCA가 실행되지 않아 좌표가 non-finite인 scope는 기존 not-run 산출물을 보존하고 skip한다.

전역 audit 결과는 `experiments/outlier_candidate_summary/`에 저장한다.

- `outlier_candidate_occurrences.csv`: experiment 경로, PCA candidate sample ID, metadata, 실제 PCA distance·threshold와 실험 내 distance rank
- `outlier_candidate_experiment_counts.csv`: 실험별 후보 수와 최대 distance
- `outlier_candidate_sample_counts.csv`: source sample별 총 등장 횟수, 서로 다른 실험 수와 analysis-row 수
- `summary.json`: 후보 규칙과 전체 집계

동일 source sample이 threshold, cell type, feature-selection strategy, before/after-preprocess에서 반복 검출되면 각각 별도 occurrence로 센다. `distinct_experiment_count`는 서로 다른 `pca_scores.csv`에서 검출된 횟수이고, `distinct_analysis_row_count`는 cell-integrated ID를 포함한 서로 다른 분석 row 수다.

## 21. Transcriptomics 확정 outlier 제거 후 no-cut downstream

Variance-filtered transcriptomics의 민감도 분석에서는 공통적으로 심각한 outlier로 확인된 subject `Ng_MOG1`, `Ng_HD3`만 명시적으로 제거한다. 다른 sample은 자동 제거하지 않는다. 각 network에서 두 subject 제거 후 PyWGCNA `preprocess()`를 `cut=inf`로 실행하며, 생성된 `outlier_removal/03_preprocess/01_cut_diagnostic/matrix_before_cut.csv.gz`가 downstream의 유일한 expression 입력이다.

이 branch에는 `03_preprocess/02_selected_cut`을 만들지 않는다. dendrogram과 largest-gap 후보는 진단 근거로만 보존하며 finite sample-tree cut을 적용하지 않는다. 따라서 `findModules` 직전에도 추가 feature filtering, log 변환, normalization, scaling 또는 결측 대체를 하지 않는다.

실행 순서는 다음과 같다.

1. `run_transcriptomics_outlier_removed_find_modules.py`: 자동 soft-threshold와 기본 module 조건 `deepSplit=2`, `minModuleSize=50`, `MEDissThres=0.20`, `pamRespectsDendro=False`
2. `run_transcriptomics_outlier_removed_module_diagnosis.py`: 공식 `CalculateSignedKME()`, 공식 one-hot module-trait heatmap, disease-vs-HD pairwise heatmap
3. `run_cross_omics_kme_weighted_jaccard.py --branch after_outlier_removal`: 동일 specimen의 proteomics module과 kME-weighted Jaccard 비교

각 network의 결과는 각각 `outlier_removal/05_find_modules/01_default_ds2_mms50_med0p20/`와 `outlier_removal/06_kme_module_diagnosis/01_default_ds2_mms50_med0p20/`에 저장한다. cross-omics 결과는 `experiments/05_cross_omics_module_similarity/after_outlier_removal/`에 저장해 기존 결과를 덮어쓰지 않는다. Dashboard에서는 Before/After Outlier Removal을 같은 stage 안에서 전환하되, After Outlier Removal의 Preprocess에는 `Cut diagnostic`만 표시하고 `Selected cut`은 표시하지 않는다.

## 22. Transcriptomics batch-correction versioning

새로 전달받은 보정 CSV는 기존 보정값에 다시 ComBat을 적용한 순차적 2차 보정이 아니다. 동일한 upstream CPM 행렬에서 design matrix만 바꾸어 다시 적합한 대체 보정본이다. 따라서 다음 세 버전을 서로 덮어쓰지 않고 병렬 비교한다.

1. `uncorrected_cpm`: 전달받은 CPM workbook의 저장값
2. `combat_v1_cell_diagnosis`: cell type과 diagnosis를 보존한 기존 ComBat 결과
3. `combat_v2_cell_diagnosis_tissue`: cell type, diagnosis, tissue/specimen을 보존한 새 ComBat 결과

새 보정본의 재현 실행 루트는 `runs/combat_v2_cell_diagnosis_tissue/`이다. 입력 보정 CSV는 연구팀이 제공한 결과이므로 이 저장소에서 ComBat 자체를 다시 적합했다고 기록하지 않는다. 대신 받은 파일의 위치, 보존 변수, downstream 파라미터와 산출물 경로를 `RUN_CONFIG.md`와 각 단계 manifest에 기록한다.

Batch-correction QC는 동일한 `threshold × cell type × specimen` 범위에서 세 버전을 비교한다. heatmap과 PCA에는 저장된 expression 값을 그대로 사용하며 추가 log 변환, feature scaling, 결측치 대체, top-N feature 선택을 하지 않는다. PCA는 현재 분석 범위에서 결측값이 하나라도 있는 gene 전체와 zero-variance gene을 제외한 뒤 sklearn의 평균 중심화만 사용한다.

`Before Preprocess`는 raw/uncorrected라는 뜻이 아니다. 각 correction version에 대해 PyWGCNA `preprocess()`를 호출하기 전의 행렬이라는 stage 의미를 유지한다. Dashboard 최상위 `Batch correction` 탭이 correction version을 선택하고, 그 아래 기존 Before Preprocess / Preprocess / After Preprocess / Find Modules / Downstream 구조를 유지한다.

## 23. Final Jaccard의 보정 버전과 순열 검정

- Final Heatmap의 행, 순위, CSV와 모듈 상세는 선택한 transcriptomics 보정 버전의 record ID를 공유한다. CSF cell-separate RNA에서 `other`를 제외하고 같은 단백질 참조 조건 `04_ds3_mms30_med0p10`의 CSF 및 PB 네트워크를 비교한다.
- 메인 dashboard 인덱스를 만든 뒤 `build_final_jaccard_protein_scope.py`를 보정 버전별로 실행한다. v2에는 `--batch-correction combat_v2_cell_diagnosis_tissue`, v2의 `--run-root`, 공통 단백질 결과 루트인 `--proteomics-run-root`, 메인 `--dashboard-root`를 함께 지정한다. `build_cross_omics_view.py --dashboard-root`는 같은 대시보드에 diagnosis forest plot 데이터를 생성한다.
- `build_final_jaccard_bubble_figure.py --batch-correction`은 CSF RNA × CSF protein, HVG 9,000, Ng_MOG1·Ng_HD3 제거 후의 5-cell과 10-cell 조건을 각각 검정한다. 공통 측정 유전자 집합 내 RNA 라벨을 모듈 배정과 kME 가중치의 짝을 유지하며 네트워크별로 10,000회 섞는다. 기본 seed는 0으로 고정하고, 정해진 cell type 및 threshold 순서에서 네트워크의 0-based 순번을 더한다. p는 `(b+1)/(B+1)`이며, BH-FDR은 보정 버전과 threshold 하나에 포함된 모든 모듈 쌍에 적용한다. q ≤ 0.25는 탐색적 후보 기준이다.
- v2의 통계 원본, 후보 목록, 공유 gene 및 diagnosis 근거, PNG/PDF/SVG는 v2 루트의 `experiments/05_cross_omics_module_similarity/final_jaccard_permutation/`에 있다. 뷰어 경로는 `dashboard/figures/final_jaccard_bubble/combat_v2/`이다. 기존 v1 통계는 `--from-saved-statistics`로 뷰어 데이터를 생성할 수 있다. `build_full_dashboard_offline_release.py --final-jaccard-only`는 이 산출물과 overlap/diagnosis 데이터를 기존 오프라인 번들에 반영하고 checksum을 갱신한다.
