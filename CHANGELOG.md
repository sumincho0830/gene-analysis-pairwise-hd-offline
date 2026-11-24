# Pipeline change log

## 2026-11-24 - Compare CSF RNA with CSF and PB protein networks

- Preserved the established specimen-matched `CSF RNA x CSF protein` results.
- Added a sensitivity branch that recomputes the same absolute-own-kME
  weighted Jaccard statistic for `CSF RNA x PB protein`; no CSF result is
  relabelled or reused as PB.
- Added `All protein networks`, `CSF protein`, and `PB protein` controls to the
  Final Jaccard page. In the combined view, one RNA row is retained while CSF
  and PB protein-module columns are labelled separately.
- Added 138 PB-protein comparison records: 12 cell types per threshold-5
  condition and 11 saved cell types per threshold-10 condition, for three
  feature strategies and both outlier-removal states.
- Extended module-pair inspection so PB selections load PB protein membership
  and diagnosis evidence while reusing the exact selected CSF RNA network.

## 2026-11-24 - Restore post-preprocess outlier-removal comparison

- Generated the standard ten descriptive/distribution outputs from each saved
  `outlier_removal/03_preprocess/01_cut_diagnostic/matrix_before_cut.csv.gz`.
- These matrices are the outputs of PyWGCNA `preprocess(cut=inf)` after removing
  `Ng_MOG1` and `Ng_HD3`; no finite sample cut, imputation, or additional value
  transform was introduced.
- Added `Before Outlier Removal` and `After Outlier Removal` controls to the
  `After Preprocess` stage while preserving the selected output type when the
  user switches removal state.
- Packaged complete post-preprocess heatmap, PCA, interactive PCA 3D, and PCA
  score artifacts for 165 eligible transcriptomics networks.

## 2026-11-24 - Preserve Final Jaccard view during pair inspection

- Changed embedded Final Jaccard heatmap cells and ranking links to open the
  existing RNA-protein module comparison modal in the parent dashboard.
- The Final Jaccard page, selected condition, and iframe scroll position now
  remain unchanged while inspecting or closing a module pair.
- Kept the original deep-link navigation as the fallback when
  `final-jaccard.html` is opened as a standalone page.
- Restricted parent-window messages to the embedded Final Jaccard iframe and
  known cross-omics records.

## 2026-11-24 - Permutation-tested final CSF Jaccard bubble figures

- Fixed the comparison to CSF, cell-separated transcriptomics, HVG top 9,000,
  removal of `Ng_MOG1` and `Ng_HD3`, and exclusion of `cell_type=other`.
- Generated separate threshold-5 and threshold-10 bubble figures in PNG, PDF,
  and SVG. Dot area encodes the absolute-own-kME-weighted Jaccard similarity;
  fill color encodes permutation BH-FDR q.
- Used 10,000 deterministic gene-label permutations per cell network. The null
  preserves both module sizes and both omics' absolute own-module kME weights,
  while breaking RNA-protein gene identity correspondence.
- Computed one-sided empirical p-values as `(b+1)/(B+1)` and applied
  Benjamini-Hochberg correction across every module pair within each threshold
  figure.
- Threshold 5 had no pair at q<=0.25. Threshold 10 had five exploratory pairs
  at q<=0.25, but neither threshold had a pair at q<=0.10.
- Exported all pair statistics, q<=0.25 candidate rows, ranked shared genes, and
  RNA/protein module diagnosis evidence for MS, MOGAD, and NMOSD versus HD.

## 2026-11-23 - After-outlier-removal cut recommendation review

- Changed the review suggestion to the exact midpoint of the global largest
  consecutive sample-linkage height gap.
- Added warnings when the largest-cluster interpretation is ambiguous or would
  remove more than 5% of rows; warnings do not suppress the recommendation.
- Added a one-network-at-a-time review page with specimen, threshold, feature
  strategy, mode, cell type, and decision filters.
- Added editable cut values, `Use selected cut`, `Hold`, `Do not use`, browser
  persistence, and CSV export/import.
- Kept every finite cut unapplied pending manual review.

## 2026-11-23 - Confirmed outlier-removal comparison branch

- Preserved existing variance-filtered results as `Before Outlier Removal`.
- Excluded `Ng_MOG1` and `Ng_HD3` across CSF/PB and all cell-type rows.
- Regenerated descriptive/distribution QC with no feature or value transform.
  The preprocess input feature set is unchanged; default `goodSamplesGenes`
  removals after the sample change are recorded as preprocess results.
- Re-ran PyWGCNA `preprocess(cut=inf)` only.
- Saved a review-only global largest-gap midpoint; a warning marks suggestions
  whose largest-cluster interpretation would remove more than 5% of rows.
  No selected cut was applied.
- Added nested dashboard navigation for before/after outlier removal.

## 2026-11-23 — 기존 PCA outlier candidate 라벨과 전역 집계 추가

- PCA 좌표와 candidate 판정은 다시 계산하지 않고 기존 `pca_scores.csv`의
  `pca_outlier` 또는 `pca_outlier_candidate`를 그대로 사용했다.
- 저장된 PCA scope 554개를 순회해 좌표가 유효한 545개의 2D metadata figure와
  3D 정적 figure를 다시 생성했다. PCA 후보가 있는 373개 scope에는 빨간 빈 원,
  SampleID와 `PCA outlier candidate` legend를 추가했고, PCA가 실행되지 않은 9개는
  기존 not-run 산출물을 보존했다.
- descriptive-only 후보와 임의의 top-N 선택은 PCA figure에 추가하지 않으며,
  저장된 PCA 후보만 표시한다.
- experiment별 후보와 실제 distance를 저장한 occurrence table, experiment별
  후보 수, source sample별 반복 검출 횟수 표를
  `experiments/outlier_candidate_summary/`에 추가했다.
- PCA 후보는 총 1,638회, source sample 131개, experiment 373개에서 기록됐다.
- 재현용 `code/src/refresh_pca_outlier_annotations.py`를 추가했다.

## 2026-11-21 — 모든 cross-omics module pair의 공통 유전자 탐색 추가

- reciprocal-best pair에만 저장하던 gene contribution을 모든 30,270개 module pair로
  확장해 `all_pair_gene_contributions.csv.gz`로 저장했다.
- 165개 network에서 245,662개 공통 gene 행을 생성했다.
- pair별 shared-gene 수와 weighted-intersection 합을 원래 similarity table과 직접
  비교하는 validator와 단위 테스트를 추가했다.
- Dashboard에 RNA/protein module selector, gene search, kME 정렬, pagination,
  CSV export를 갖춘 `Cross-omics Gene Overlap` explorer를 추가했다.
- offline release builder가 overlap index, 165개 lazy-load shard, 계산·검증 코드를
  함께 포함하도록 변경했다.

## 2026-11-19 — Dashboard 산출물 탭 단계화

- 한 줄에 혼재하던 산출물 탭을 `Before Preprocess`, `Preprocess`,
  `After Preprocess`, `Find Modules`, `Downstream`의 상위 stage로 분리했다.
- 선택한 stage의 세부 산출물 탭만 아래에 나타나도록 변경했다.
- Cell type comparison과 Individual run에 동일한 분류 및 순서를 적용했다.

## 2026-11-19 — Module Explorer 추가

- transcriptomics와 proteomics의 공식 kME 결과를 공통 schema로 변환하는
  offline Module Explorer를 dashboard에 추가했다.
- 173개 network, 2,177개 module, 1,339,762개 feature를 module별 lazy-load shard로
  연결했다.
- module 선택, identifier 검색, hub rank/|kME|/signed kME/feature ID 정렬,
  pagination, CSV export를 추가했다.
- transcriptomics는 gene ID를, proteomics는 PROBEID/SYMBOL/UniProt을 함께 표시한다.
- `is_raw_unassigned=True` feature는 목록에 남기되 kME와 rank를 `NA`로 표시한다.
- Cell type comparison, Individual run, module table, comparison card에서 동일한
  Explorer를 열도록 연결했다.
- full offline release builder가 `module-members/` index와 shard도 함께 복사하도록
  갱신했다.

## 2026-11-19 — Cell type 중심 조건별 산출물 비교 화면 추가

- 기존의 단일 실험 record 중심 dashboard에 `Cell type comparison` 기본 화면을 추가했다.
- cell type을 한 번 선택하면 모든 threshold, feature-selection 전략과 CSF/PB 조건을 유지하고, 산출물 종류만 탭으로 전환하도록 화면 축을 바꿨다.
- 각 탭에는 같은 산출물을 모든 조건에 대해 2열로 배치하고, 동일 threshold·전략의 CSF/PB가 연속해서 보이도록 자연 정렬했다.
- transcriptomics B cell 기준 12개 조건과 16개 artifact 탭이 한 페이지에서 동작하며, 각 panel에서 원본 artifact 또는 기존 단일-record inspection으로 이동할 수 있다.
- 기존 review `localStorage`, CSV import/export, 개별 실험 상세, module table, provenance와 최대 4개 임의 비교 기능은 보존했다.
- Python compile, 176-record link validation과 headless Edge 실제 렌더링을 통과했고 preview PNG를 최신 화면으로 갱신했다.

## 2026-11-19 — 예측 `MEDissThres=0.10` 실제 재계산 추가

- ME dendrogram 후보 120개에 대해 `deepSplit=2`, `minModuleSize=50`, power, adjacency, TOM, gene tree와 dynamic module은 유지하고 `MEDissThres`만 0.20에서 0.10으로 변경했다.
- 전체 network를 처음부터 중복 계산하지 않고, 변경된 단계에 해당하는 PyWGCNA 공식 `mergeCloseModules(cutHeight=0.10, iterate=True)`를 저장된 dynamic module에 실행했다.
- 후보 계산 전 `MEDissThres=0.20`으로 기존 baseline partition과 exact color label을 재현하는 gate를 추가했고 120/120개가 통과했다.
- 120개 실제 계산은 failure 0개로 완료됐다. 총 module 수는 baseline 1,774개, 단순 dendrogram 예측 2,363개, 실제 iterative merge 2,339개였고 baseline collapse 후보 4개는 모두 해소됐다.
- 각 network의 `05_find_modules/02_predicted_ds2_mms50_med0p10/`에 baseline과 동일한 compact CSV/JSON/figure를 저장했다.
- 전역 비교 CSV/JSON과 dashboard를 갱신해 baseline, 예측값, 실제값 및 예측 오차를 한 화면에서 확인할 수 있게 했다.

## 2026-11-19 — Pairwise-HD offline GitHub release 추가

- 41GB 전체 실행 결과 대신 pairwise-HD 집계 결과만 포함한 약 41MB private 배포 저장소를 만들었다.
- 173개 heatmap, 1,211개 module-level CSV/JSON, offline dashboard와 재현 코드를 포함했다.
- raw/formatted matrix, 개인 sample ID, PyWGCNA 객체와 feature-level kME는 배포에서 제외했다.
- `SHA256SUMS.txt`와 `verify_release.py`를 추가하고 Windows line-ending 변환을 비활성화했다.
- GitHub에서 새로 clone한 복사본의 hash와 `file://` 렌더링을 확인했다.

## 2026-11-19 — Module–diagnosis를 direct disease-vs-HD 대조로 수정

- 기존 PyWGCNA 다범주 one-hot heatmap이 `MS vs HD`가 아니라 `MS vs all other diagnoses`를 계산한다는 차이를 확인했다.
- 각 network에서 다른 질환을 제외하고 `HD=0`, target disease=1로 encoding한 `MS vs HD`, `MOGAD vs HD`, `NMOSD vs HD` point-biserial correlation을 새 주 결과로 생성했다.
- signed `r`, two-sided `p`, inference-ready cell에 대한 within-network BH `q`, 질환별 sample 수와 target-minus-HD ME 평균 차이를 CSV에 보존했다.
- 기존 one-hot figure와 표는 삭제하지 않고 audit 결과로 유지했다.
- dashboard의 기본 module–diagnosis figure를 새 `vs HD` heatmap으로 바꾸고 기존 결과는 `One-hot diagnosis (audit)` tab으로 분리했다.
- Windows 260자 경로 제한으로 긴 transcriptomics heatmap이 `file://`에서 깨지는 문제를 짧은 dashboard asset mirror로 해결했다.
- 생성 규칙, 단위 테스트와 network별 manifest를 추가했다.

이 파일은 생성 규칙의 변경 이력을 누적해서 기록한다. 현재 유효한 규칙 자체는 `PIPELINE_SPEC.md`를 기준으로 한다.

## 2026-11-19 — Inspection dashboard visual redesign

- `results_20261015/design-system.md`의 흰 배경, 중립 회색, Arial, 가는 구분선과 밑줄형 활성 상태를 현재 experiment dashboard에 적용했다.
- 기존 카드·배지·둥근 버튼 중심 표현을 논문형 experiment atlas 레이아웃으로 정리했다.
- KPI, filter toolbar, experiment list, metric grid, review control, artifact tabs, module table과 provenance accordion의 typography와 spacing을 통일했다.
- 데이터 schema, 176개 record, 165개 cross-omics artifact, filter/search, review localStorage, CSV import/export와 comparison 기능은 변경하지 않았다.
- responsive breakpoint를 추가해 좁은 화면과 200% 확대에서 list/detail이 단일열로 전환되도록 했다.
- visual 규칙을 builder의 `DESIGN_SYSTEM_CSS`에 보관해 향후 dashboard 재생성에도 유지되도록 했다.
- HTML parse, 176개 record 링크 검증과 Chrome 100%/200% 실제 렌더링을 확인했다.

## 2026-11-19 — Cross-omics kME-weighted Jaccard 추가

- 완료된 transcriptomics network 165개 각각에 대해 specimen-matched proteomics reference `04_ds3_mms30_med0p10`의 모든 proper module과 비교했다.
- 공통 측정 gene universe에서 `absolute own-module kME`를 가중치로 한 fuzzy Jaccard를 계산했다. raw-label-0 fully-unassigned group은 제외했다.
- proteomics symbol은 uppercase simple symbol로 정규화하고 모호한 다중 symbol을 제외했다. 중복 probe는 최대 signed kME 기반 대표 probe 하나로 collapse했다.
- 모든 module 조합 30,270개를 CSV로 저장하고, 실험 간 동일한 color scale을 적용한 heatmap 165개와 filterable gallery를 생성했다.
- 공통 gene이 0개인 module도 축에는 유지하되 best/reciprocal-best로 표시하지 않는 경계 규칙을 추가했다.
- 단위 테스트 4개와 전수 validator를 추가했다. validator는 165개 각 실험의 최대 pair를 원본 kME에서 재계산했고 오류 0건으로 통과했다.
- 기존 experiment dashboard에 cross-omics heatmap, 공통 gene 수, 최대/중앙 best similarity, reciprocal pair 수, 세부 CSV·manifest 링크를 연결하고 브라우저에서 최신 화면과 이미지 로딩을 확인했다.
- weighted Jaccard는 기술적 overlap이며 공식 WGCNA module preservation 또는 유의성 검정이 아니라는 해석 한계를 명시했다.

## 2026-11-19 — 공식 PyWGCNA kME 및 module-diagnosis heatmap 추가

- Proteomics 8개 condition과 transcriptomics 완료 network 165개에서 공개 `CalculateSignedKME()`를 실행해 모든 feature × final ME signed kME를 저장했다.
- PyWGCNA `getDatTraits()`의 diagnosis one-hot encoding과 공개 `module_trait_relationships_heatmap()`을 사용해 공식 `r (p)` heatmap 173개를 생성했다.
- pandas 3에서 범주형 encoding이 object dtype으로 남는 호환 문제는 one-hot 값을 변경하지 않고 float dtype만 확정해 공식 numeric-trait 경로로 처리했다.
- 공식 heatmap은 보정 전 양측 Pearson p-value를 유지하고, within-network BH-FDR은 별도 long-form CSV에 기록했다.
- 1,000-feature chunk마다 공개 kME 메서드를 호출해 사용하지 않는 O(p²) feature-feature correlation 할당을 제한했다. 기존 assigned kME와 최대 절대 차이는 `2.4424906541753444e-15`였다.
- raw dynamic label 0 fully-unassigned eigengene을 proper module 분석에서 제외했다. 이는 PyWGCNA가 정상적인 grey module에서 의도하는 ME 제외 규칙을 색 이름이 아닌 저장된 raw label audit로 복원한 것이다.
- 173개 결과에서 proper module 대상 공식 kME 16,716,991개와 module-diagnosis 검정 8,111개를 검증했다. 추정 불가능한 proper-module kME column과 diagnosis 검정은 0개였다.
- Inspection dashboard에 proteomics adapter와 module-diagnosis artifact/kME/statistics 링크를 추가해 총 176개 record를 `file://`에서 검토할 수 있게 했다.

## 2026-11-19 — 확장 가능한 offline experiment inspection dashboard 예시 추가

- 실험별 경로를 화면에 하드코딩하지 않고 source adapter가 공통 record schema를 만드는 구조로 분리했다.
- `04_variance_filtered/transcriptomics`의 168개 실험을 threshold, feature 전략, cell mode, specimen, cell type과 실행 상태로 검색·필터링할 수 있게 했다.
- sample/gene 수, soft threshold, module 수, unassigned 비율, kME, runtime과 ME/gene dendrogram, QC figure, CSV·JSON provenance를 한 record에서 확인하도록 구성했다.
- `Use / Hold / Do not use`, 네 가지 inspection checklist, 근거 메모를 브라우저에 보존하고 CSV로 export/import할 수 있게 했다.
- 최대 4개 실험의 ME dendrogram과 핵심 수치를 comparison workspace에서 나란히 볼 수 있게 했다.
- 외부 library나 server 없이 `file://`에서 동작하도록 만들고, 168개 record ID와 모든 artifact/file 상대경로, 165개 ME dendrogram 연결, headless Edge 렌더링을 검증했다.

## 2026-11-19 — Transcriptomics dendrogram gallery를 ME dendrogram 기준으로 보완

- 최초 gallery가 gene dendrogram을 표시한 것은 사용자 의도와 달랐다. proteomics `module_eigengene_dendrogram.png`와 같은 pre-merge dynamic-module eigengene dendrogram을 별도로 생성했다.
- selected-cut matrix와 저장된 raw dynamic assignment로 proper dynamic module별 signed first PC를 계산했다. raw label 0 unassigned group은 제외했다.
- ME 거리는 `1 - Pearson correlation`, clustering은 average linkage이며, `MEDissThres=0.20` 가로선을 표시했다. module 간 거리가 모두 0.20보다 작아도 기준선이 보이도록 y축을 보정했다.
- 165개 완료 network에 `module_eigengene_dendrogram.png`, `dynamic_module_eigengenes.csv`, `module_eigengene_linkage.csv`와 전용 manifest를 추가했다.
- 주 gallery `module_eigengene_dendrogram_gallery.html`에는 dynamic/final module 수를 함께 표시했다. 기존 gene dendrogram gallery는 보조 진단용으로 보존했다.
- validator를 확장해 165개 ME 결과와 두 gallery의 168개 record 및 상대 이미지 링크를 모두 검증했다.

## 2026-11-19 — Transcriptomics variance-filter 기본 findModules와 dendrogram gallery 추가

- CPM detection, corrected-value variance 상위 9,000개, 상위 5,000개 전략의 selected-cut preprocess 객체를 입력으로 기본 `findModules()`를 실행하도록 재사용 가능한 runner를 추가했다.
- 추가 feature filtering이나 값 변환 없이 PyWGCNA 기본 power 후보를 사용해 network별 soft threshold를 자동 선택했다.
- 기본 module parameter(`deepSplit=2`, `minModuleSize=50`, `MEDissThres=0.20`, `pamRespectsDendro=False`)를 고정하고 parameter tuning 전 baseline으로 명시했다.
- dense adjacency/TOM과 전체 module 객체 대신 module assignment·summary·eigengene, soft-threshold 표, gene-tree linkage, figure와 manifest를 저장하도록 했다.
- dynamic cut과 ME 병합 후 module color를 함께 표시한 전체 gene dendrogram을 network별로 생성하고, 모든 dendrogram을 필터링해 볼 수 있는 offline `module_dendrogram_gallery.html`을 추가했다.
- PyWGCNA 2.2.1의 read-only TOM 입력과 단일-module `SystemExit` 경로는 계산 의미를 보존하는 scoped compatibility 처리로 대응하고 그 사용 여부를 결과 JSON에 기록했다.
- 장시간 dense TOM 계산은 완료 결과를 재사용하는 resume 규칙과 network-disjoint worker 분할을 사용하며, 전역 CSV·HTML은 모든 worker 종료 후 한 번만 합친다.
- 168개 조합 중 165개에서 `findModules()`를 완료했고, 표본이 2개인 동일 network 3개는 이전 preprocess 상태에 따라 명시적으로 skip했다. failure는 0개다.
- 완료 network의 final proper module 수는 1–51개였고 9개 network가 하나의 proper module로 collapse했다. 선택 power는 자동 탐색 결과 1–19 범위였다.
- 재실행 가능한 전용 validator를 추가해 165개 완료 output set의 필수 파일과 내부 차원, 168개 gallery record와 모든 상대 이미지 링크를 검증했다. headless Edge `file://` 렌더링에서도 168/168 카드와 로컬 dendrogram 이미지 표시를 확인했다.

## 2026-11-19 — Transcriptomics CPM/variance feature-selection preprocess 비교 생성

- `03_gene_filtered`의 18,045-gene batch-corrected matrix를 기준으로 CPM detection, corrected-value variance 상위 9,000개, corrected-value variance 상위 5,000개의 세 전략을 추가했다.
- CPM을 corrected matrix에서 재계산하지 않고 source CPM workbook 12개를 직접 읽어 `CPM >= 1` detection support를 계산했다. 모든 workbook에서 gene-filtered 18,045개 gene이 매칭됨을 audit에 기록했다.
- CPM 통과 기준은 현재 분석 row의 `max(4, ceil(10%))`이며, row가 4개 미만이면 사용 가능한 전체 row로 제한한다. cell-integrated에서는 cell type별로 평가한 뒤 어느 한 cell type에서라도 통과한 gene을 유지했다.
- HVG는 별도 log·scaling 없이 각 network의 저장된 corrected value에 대한 population variance(`ddof=0`)로 계산했다.
- 각 전략-network에서 before-preprocess QC, `cut=inf` diagnostic, 보수적인 자동 selected cut, selected-cut PyWGCNA preprocess, after-preprocess QC를 생성했다.
- 자동 cut은 절대·상대 linkage gap percentile, 최소 상대 gap, robust distance flag, 최대 5% 제거, metadata 보존 조건을 동시에 만족할 때만 적용했다.
- 168개 조합 중 165개를 완료하고, 표본이 2개인 동일 network 3개를 명시적으로 skip했다. finite cut은 전략별 18/19/17개였고 제거 sample row는 18/20/19개였다.
- 큰 matrix는 gzip과 NTFS hard link로 중복을 줄였고, 신규 PyWGCNA 객체는 zstd로 압축했다. Plotly 3D HTML은 CDN loader를 사용했다.
- 재개 실행에서 전역 summary의 입력 차원과 상태가 비는 문제를 수정하고 168행 전체에서 필수 필드 결측이 없음을 검증했다.
- 이 단계에서는 transcriptomics `findModules()`를 실행하지 않았다.

## 2026-11-19 — Transcriptomics seed-0 replicate filtering QC 생성

- batch-corrected CIS 제외 checkpoint에서 terminal specimen replicate 선택을 실행했다. 기존 selected-cut 결과는 입력으로 사용하지 않았다.
- sample ID 끝의 specimen 숫자만 제거해 group을 만들고 local `random.Random(0)`으로 하나를 유지했다.
- 두 threshold의 source sample 집합이 동일함을 확인하고, 하나의 전역 선택을 cell-separated와 cell-integrated에 공통 적용했다.
- 131개 source sample 중 반복 그룹 16개에서 20개를 제외하고 111개를 유지했다. 제외 ID는 2026-10-22의 `2_random_sampled` 결과와 정확히 일치했다.
- threshold 5/10, cell-separated/integrated, CSF/PB의 56개 네트워크에 filtered matrix와 descriptive statistics 5개, data distribution 5개를 생성했다.
- 55개 PCA가 완료됐다. threshold 10 / CSF / Non-classical Monocytes는 표본 2개라 PCA placeholder와 명시적 상태를 저장했다.
- 사용자의 중단 요청에 따라 PyWGCNA preprocess, dendrogram, finite cut, `.p` 객체와 `findModules`는 실행하지 않았다.

## 2026-11-19 — Transcriptomics reviewed selected-cut preprocess 생성

- `selected_cuts_review.csv`의 확정 결과 44개(`use` 20, `reject` 24)를 단일 의사결정 입력으로 사용했다.
- flagged=0이라 review에서 제외된 11개 network는 자동 no-cut, 표본 수 2개인 1개 network는 기존 기본값을 완화하지 않고 skip으로 기록했다.
- 각 network를 CIS 제외 checkpoint의 새 PyWGCNA 객체에서 다시 preprocess하고, cut 직전 행렬이 기존 `01_cut_diagnostic/matrix_before_cut.csv`와 동일한지 검증했다.
- finite cut은 선택한 높이 아래의 유일한 최대 sample cluster를 유지하도록 확정했다. 20개 finite-cut network에서 총 21개 row를 제거했다.
- PyWGCNA 2.2.1이 임의의 cluster label 0만 유지하는 문제를 막기 위해, 최대 군집이 label 0이 아닌 경우 두 label을 `preprocess()` 내부에서만 교환했다. 이 처리는 threshold 5/10의 PB CD4 T cell 두 network에서 발생했다.
- `02_selected_cut` 56개 디렉토리, 완료 객체·행렬 55개, 명시적 skip 1개를 생성했다. no-cut 35개 행렬은 diagnostic matrix와 SHA-256이 일치함을 검증했다.
- after-preprocess QC와 transcriptomics `findModules`는 실행하지 않았다.

## 2026-11-18 — Transcriptomics raw QC와 corrected cut diagnostic 생성

- 13개 cell type의 raw `cell_separate` QC와 `SampleID + specimen + cell type` 기반 `cell_integrated` QC를 생성했다. 각 분석에 ALL과 CSF/PB별 distribution을 저장했다.
- raw/corrected cell-integrated ALL matrix를 metadata의 specimen 값으로 분할한 CSF/PB formatted matrix를 추가하고 formatted manifest를 갱신했다.
- `threshold_5_cells`, `threshold_10_cells` 각각에서 `cell_separate`와 `cell_integrated`의 feature provenance, CIS 제외, before-preprocess QC, `cut=inf` diagnostic을 생성했다.
- stored value를 그대로 사용했으며 추가 log, scaling, 결측 대체, MAD/top-N feature 선택을 적용하지 않았다. PCA는 complete-case nonconstant gene 전체를 mean-center하여 계산했다.
- PyWGCNA 2.2.1의 default `goodSamplesGenes()`가 pandas 2.2 이상에서 bool Series에 integer count를 대입하며 실패하는 부분을 dtype 호환 shim으로 수정했다. 기본 threshold와 반복 판정식은 변경하지 않았다.
- threshold-5 결과는 `goodSamplesGenes()` 제거 gene이 0개여서 원본 구현으로 완료됐고, 반복 제거가 발생한 threshold-10에 호환 shim을 적용했다.
- transcriptomics heatmap에도 specimen/cohort group label을 추가했으며 cell-integrated heatmap에는 cell type까지 표시했다. outlier band는 expression panel 밖에 유지했다.
- `refresh_transcriptomics_heatmaps.py`로 heatmap 98개만 다시 생성했고, QC manifest에 heatmap-only refresh임을 기록했다.
- 전체 56개 diagnostic 중 55개를 완료했다. `threshold_10_cells / cell_separate / CSF / Non_classical_Monocytes`는 CIS 제외 후 sample 2개여서 PCA와 preprocess를 계산하지 않고 명시적 상태 파일/placeholder를 저장했다.
- finite cut, selected-cut preprocess, after-preprocess QC와 transcriptomics `findModules`는 실행하지 않았다.

## 2026-11-18 — Proteomics findModules parameter sensitivity study

- 기존 CSF/PB 기본값 결과를 `01_default_ds2_mms50_med0p20` 하위 디렉토리로 이동했다.
- 동일 specimen에서 power, adjacency, TOM, gene tree를 고정하고 `deepSplit`, `minModuleSize`, `MEDissThres`만 한 단계씩 바꾸는 네 조건을 실행했다.
- condition 디렉토리명에 비교 parameter를 기록하고, 실제 실행값은 각 `run_manifest.json`에 다시 저장했다.
- color 이름 기반 `dimgrey=unassigned` 규칙을 제거하고, `cutreeHybrid()`가 color 변환 전에 반환하는 raw `Name == 0`을 미배정 기준으로 사용했다.
- CSF color group은 7→10→15→19개, PB는 4→7→11→13개로 증가했다.
- `04_ds3_mms30_med0p10`은 과거 동료 결과의 module color별 protein 수를 CSF와 PB 모두 정확히 재현했다.
- 최종 condition은 자동 선택하지 않았고 `selected_condition.json`을 `pending_review`로 저장했다.
- module-trait association과 enrichment는 실행하지 않았다.

## 2026-11-18 — 시험 module-trait 단계 제거 및 module label audit

- 시험 실행했던 `06_module_trait_association` 결과와 전용 실행 코드를 제거했다. 파이프라인 gate를 baseline `findModules()` 이후로 되돌렸다.
- PyWGCNA 2.2.1의 `cutreeHybrid()` 내부 label 재정렬을 추적해 color 이름만으로 unassigned를 판정할 수 없음을 확인했다.
- 기존 export의 `dimgrey=unassigned` 가정은 잘못됐으며, 원래 dynamic label 0은 현재 CSF에서 `white` 37개, PB에서 `lightgrey` 158개 feature에 대응했다.
- module assignment 자체는 바꾸지 않았지만 현재 baseline의 unassigned 관련 summary 필드는 수정·재생성 전까지 해석하지 않는다.
- 동료 parameter 결과를 재현하기 전에 raw dynamic label을 별도 저장하는 규칙과 parameter sensitivity study 기준을 먼저 확정한다.

## 2026-11-18 — Proteomics findModules 완료

- selected-cut 이후 CSF 254개 sample과 PB 113개 sample에 대해 specimen별로 독립적인 `findModules()`를 실행했다.
- 두 분석 모두 5,284개 protein 전체를 사용했으며 추가 결측 대체, feature filtering, MAD/variance ranking, top-N 선택 또는 값 변환을 하지 않았다.
- PyWGCNA 2.2.1 기본 자동 soft-threshold 탐색을 사용해 CSF power 4, PB power 3을 선택했다. PB power 2는 scale-free fit은 통과했지만 mean connectivity 154.388로 `MeanCut=100`을 통과하지 못했다.
- CSF는 assigned module 6개와 unassigned 375개 feature, PB는 assigned module 3개와 unassigned 770개 feature를 얻었으며 두 network 모두 하나의 assigned module로 collapse하지 않았다.
- `findModules()` 후 `CalculateSignedKME()`를 실행하고 module color, eigengene, 전체 kME, assigned kME, soft-threshold 표 및 진단 figure를 저장했다.
- PyWGCNA 2.2.1의 read-only adjacency 호환성 오류에는 값이 동일한 writable copy만 전달했다. adjacency 값, 수식 및 분석 parameter는 변경하지 않았다.
- `analyseWGCNA()`, module-trait association, enrichment 및 기타 downstream 해석은 실행하지 않았다.

## 2026-11-18 — Selected cut preprocess와 after-preprocess QC 완료

- 동료의 기존 분석과 동일하게 CSF `cut=50`, PB `cut=80`을 적용했다.
- 현재 Euclidean/average-linkage dendrogram에서 CSF 257→254, PB 115→113 sample이 되어 기존 결과를 정확히 재현했다.
- CSF에서 3개, PB에서 2개 제거 sample과 cluster assignment를 CSV로 기록했다.
- diagnostic pass와 selected-cut pass의 cut 직전 matrix가 sample, feature, 값 및 hash에서 동일한지 검증했다.
- selected-cut `.p`, `matrix_after_cut.csv`, dendrogram, 제거 목록, summary 및 manifest를 저장했다.
- cut 이후 CSF/PB 각각에 descriptive statistics 5개와 distribution/PCA 5개를 생성했다.
- `findModules`는 실행하지 않았다.

## 2026-11-18 — Before-preprocess ALL QC 추가

- CIS 제외 후 specimen으로 분리하기 전 ALL 372개 sample에 대한 descriptive statistics와 distribution 산출물을 추가했다.
- ALL PCA는 CSF/PB PCA를 합치지 않고, CIS 제외 ALL matrix에서 독립적으로 적합하도록 했다.
- 중복 ALL checkpoint matrix는 저장하지 않고 formatted ALL matrix와 `sample_selection.csv`를 ALL manifest 입력으로 기록했다.
- PyWGCNA preprocess 대상은 기존과 동일하게 CSF/PB로 유지했다.

## 2026-11-18 — 재현성 계약과 저장 구조 통합

- 지금까지 확정된 raw QC, batch-corrected QC, CIS 제외, before-preprocess QC, `cut=inf` diagnostic 순서를 하나의 기준 문서로 통합했다.
- CIS 제외 전 전체 baseline QC와 CIS 제외 후 canonical before-preprocess QC의 용도를 구분했다.
- raw proteomics에서는 preprocess를 실행하지 않는다는 결정을 명시했다.
- `02_selected_cut`, `04_after_preprocess`, `findModules`가 아직 실행되지 않았고 finite cut 결정 전 gate에서 멈춘 상태임을 명시했다.
- state-changing matrix만 저장하는 checkpoint 규칙과 future selected-cut 재실행 조건을 기록했다.
- 문서 정리만 수행했으며 분석 결과와 기존 manifest는 변경하지 않았다.

## 2026-11-18 — Batch-corrected proteomics cut=inf diagnostic 완료

- Metadata diagnosis가 CIS인 13개 sample을 분석 설계 기준으로 제외하고 selection 목록과 CSF/PB checkpoint matrix를 저장했다.
- CIS 제외 후 CSF 257개, PB 115개 sample에 대한 before-preprocess descriptive statistics와 distribution을 생성했다.
- CSF/PB 각각 PyWGCNA 기본 `TPMcutoff=1`, patch하지 않은 `goodSamplesGenes()`, `cut=inf`로 diagnostic preprocess를 실행했다.
- 두 specimen 모두 TPMcutoff/goodSamplesGenes/dendrogram cut에서 실제 제거된 feature와 sample은 0개였다.
- `matrix_before_cut.csv`, PyWGCNA `.p`, linkage, leaf order, dendrogram, distance metrics, flagged-only 후보 및 manifest를 저장했다.
- CSF 6개와 PB 5개 sample을 cohort 내 median-distance robust z 기준 diagnostic 후보로만 기록했으며 제거하지 않았다.
- finite cut, selected-cut preprocess, after-preprocess QC와 `findModules`는 실행하지 않았다.

## 2026-11-18 — Batch-corrected proteomics 결측 제외 feature 기록

- raw union 11,083개 feature 중 결측이 하나라도 있는 feature와 corrected 출력에서 빠진 feature가 정확히 일치하는지 검증하는 코드를 추가했다.
- 검증된 제외 feature 목록을 `01_batch_corrected/proteomics/00_feature_provenance/excluded_features_due_to_missingness.csv`에 저장하도록 했다.
- 목록에는 feature identity와 panel 부재 결측, panel 내부 결측, 전체 결측 수 및 비율을 기록한다.
- 전달받은 corrected matrix를 비교한 provenance 기록이며 batch correction 자체를 재실행하지 않는다.

## 2026-11-18 — Specimen별 proteomics heatmap 추가

- raw와 batch-corrected 단계에 CSF/PB 전용 heatmap을 각각 추가했다.
- 기존 전체 heatmap은 `02_data_distribution/heatmap.png`에 유지한다.
- specimen 전용 heatmap은 `02_data_distribution/by_specimen/<specimen>/heatmap.png`에 저장한다.
- 전용 heatmap의 outlier annotation은 descriptive outlier와 해당 specimen에서 독립 계산한 PCA outlier의 합집합을 사용한다.

## 2026-11-18 — Proteomics PCA 결측 처리 수정 및 specimen별 PCA 추가

- PCA의 feature별 median 결측 대체를 제거했다.
- 각 분석 범위에서 결측값이 하나라도 있는 common `PROBEID`는 feature 전체를 제외하도록 수정했다.
- complete-case feature 중 variance가 0인 feature만 추가로 제외하며, MAD/variance ranking과 top-N 선택은 사용하지 않는다.
- raw와 batch-corrected 단계 모두 `ALL`, `CSF`, `PB` PCA를 독립적으로 계산한다.
- 기존 10개 기본 산출물은 같은 위치에 유지하고, CSF/PB 전용 PCA는 `02_data_distribution/by_specimen/<specimen>/`에 분리했다.
- PCA outlier에 영향을 받는 sample statistics, outlier candidates, intensity figure와 heatmap을 포함해 두 단계 전체 산출물을 다시 생성하도록 했다.

## 2026-11-18 — Batch-corrected proteomics QC 생성

- raw와 batch-corrected proteomics가 같은 QC engine과 동일한 10개 산출물 규칙을 사용하도록 코드를 공용화했다.
- batch-corrected 입력에는 저장된 ComBat-corrected 값을 추가 log 변환 없이 사용했다.
- heatmap은 specimen/cohort 정렬과 라벨, 별도 outlier annotation 규칙을 그대로 적용했다.
- descriptive outlier와 PCA를 corrected-value 기준으로 계산했다.
- PCA는 전체 common nonconstant `PROBEID`, 평균 중심화, variance scaling 없음 규칙을 적용했다.
- `experiments/01_batch_corrected/proteomics`에 결과와 독립된 `run_manifest.json`을 생성했다.

## 2026-11-18 — Raw proteomics QC 계산 기준 통일

- PCA 입력을 `log2(raw)`에서 raw intensity로 변경했다.
- PCA 결측 cell 대체를 log2 scale feature median에서 raw scale feature median으로 변경했다.
- PCA의 feature별 z-score scaling을 제거하고 평균 중심화만 적용하도록 변경했다.
- sample descriptive outlier의 median과 IQR 지표를 log2 값에서 raw 값으로 변경했다.
- sample/protein statistics에서 log2 파생 열을 제거하고 raw 통계만 저장하도록 변경했다.
- intensity distribution figure를 sample median raw intensity 기준으로 변경했다.
- heatmap, descriptive statistics, outlier candidate 및 PCA가 동일한 raw intensity 기준을 사용하도록 통일했다.

## 2026-11-18 — PyWGCNA formatted data 이관

- `results_20261102/outputs/1_formatted_files`의 기존 formatted 파일을 `results_20261119/formatted`의 raw/batch-corrected 구조로 복사했다.
- 이전 날짜 결과는 삭제하지 않고 그대로 보존했다.
- 기존 `ALL/CSF/PB` 및 cell-type integrated 구조를 유지했다.
- 변환 코드 `format_files.py`와 `format_corrected_transcriptomics_all_fast.py`를 새 코드 디렉토리로 가져와 `results_20261119` 경로에 맞게 수정했다.
- cell-type integrated matrix 재생성 코드와 formatted 파일 manifest 생성 코드를 추가했다.

## 2026-11-18 — Raw proteomics heatmap 규칙 확정

- intensity panel을 `log2(raw)`에서 raw expression value로 변경했다.
- sample 정렬을 `specimen → cohort → sample_id`로 확정했다.
- missingness와 intensity 두 패널 모두 `specimen / cohort` 라벨을 표시하도록 했다.
- CSF가 위, PB가 아래에 오도록 specimen 순서를 확정했다.
- heatmap 위의 빨간 outlier 행 테두리를 제거했다.
- outlier candidate를 오른쪽 별도 annotation 영역의 빨간 사각형과 sample ID로 표시하도록 변경했다.
- raw intensity의 색상 범위는 전체 관측값의 1st–99th percentile로 유지했다.

## 2026-11-18 — 재현성 기록 시작

- 현재 생성 규칙의 기준 문서 `PIPELINE_SPEC.md`를 추가했다.
- 변경 이력을 위한 `CHANGELOG.md`를 추가했다.
- 전체 실행 시 입력, 출력, 코드의 SHA-256과 parameter 및 package version을 `run_manifest.json`에 자동 기록하도록 했다.

## 2026-11-19 — Transcriptomics gene-filtered QC 생성

- `02_sample_filtered`의 seed-0 matrix를 입력으로 `03_gene_filtered/transcriptomics` 단계를 생성했다.
- 2026-10-22 실행과 동일하게 housekeeping, mitochondrial, ribosomal, non-coding, human immunoglobulin의 다섯 필터를 순차 적용했다.
- reference 파일을 `data/reference/gene_filtering`으로 복사해 새 결과 디렉토리만으로 필터를 다시 실행할 수 있게 했다.
- 36,617 gene 중 18,572개를 제거하고 18,045개를 유지했다.
- 단계별 removed/retained 목록과 최종 gene 순서가 2026-10-22 저장 결과와 정확히 일치하는 assertion 및 alignment report를 추가했다.
- 56개 네트워크 각각에 최종 filtered matrix와 descriptive 5개, distribution/PCA 5개 산출물을 생성했다.
- gene filtering은 sample row와 stored expression value를 변경하지 않으며 log, scaling, normalization, 결측 대체, top-N 선택을 하지 않는다.
- 중간 expression matrix 다섯 벌은 중복 저장하지 않고 단계별 gene 목록과 최종 matrix만 저장하도록 했다.
- 55개 네트워크에서 PCA를 완료했고, 표본이 2개인 1개 네트워크는 명시적 `not_run` placeholder를 유지했다.
- 이 단계에서는 PyWGCNA preprocess, dendrogram cut, `.p`, `findModules`를 실행하지 않았다.
