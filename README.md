# Gene-analysis complete offline inspection dashboard

이 저장소는 `results_20261119`에서 로컬로 보던 experiment inspection dashboard를
다른 컴퓨터에서도 같은 구조로 검토하기 위한 **private offline bundle**입니다.

## 다른 컴퓨터에서 보기

```powershell
git clone https://github.com/sumincho0830/gene-analysis-pairwise-hd-offline.git
cd gene-analysis-pairwise-hd-offline
start index.html
```

별도 Python 환경이나 web server가 필요하지 않습니다. `index.html`을 직접 열면
기본 **Cell type comparison** 화면에서 동일 cell type의 모든 조건을 산출물 탭별로
비교할 수 있습니다. **Individual run inspection**과 review CSV import/export도
로컬 dashboard와 동일하게 유지됩니다.

## 포함된 범위

- experiment records: 176
- visual/interactive artifacts: 4565
- downloadable CSV/JSON/provenance files: 10253
- feature selection, before/after heatmap·PCA, interactive PCA 3D
- sample-cut diagnostic, soft threshold, gene/ME dendrogram, module sizes
- `Ng_MOG1`·`Ng_HD3` 제거 전/후 결과 비교; 제거 후 branch에는 finite sample cut 없음
- module–diagnosis 및 RNA–protein module similarity
- dashboard와 결과를 재생성하는 핵심 Python 코드와 pipeline 기록

PNG figure는 decoded pixel이 원본과 완전히 같은 lossless WebP로만 재포장했습니다.
Plotly HTML은 저장소에 포함된 `assets/vendor/plotly.min.js`를 사용하므로 인터넷 없이
회전·확대할 수 있습니다.

## 데이터 취급

이 배포본은 현재 로컬 dashboard와 같은 inspection 경험을 보존하므로 interactive
PCA hover와 일부 QC 파일에 source sample ID가 포함될 수 있습니다. 저장소를 private로
유지하고 승인된 사용자에게만 접근 권한을 부여해야 합니다. raw/formatted expression
matrix와 PyWGCNA pickle object는 포함하지 않습니다.

## 무결성 확인

```powershell
python verify_release.py
```

모든 파일이 `SHA256SUMS.txt`와 일치하면 `OK`를 출력합니다.
