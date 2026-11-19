# Pairwise disease-vs-HD module diagnosis offline release

이 저장소는 `results_20261119`의 **집계된 module–diagnosis 결과만** 오프라인에서 검토하기 위한 배포본입니다.

## 다른 컴퓨터에서 보기

```powershell
git clone https://github.com/sumincho0830/gene-analysis-pairwise-hd-offline.git
cd gene-analysis-pairwise-hd-offline
start index.html
```

별도 Python 환경이나 web server가 필요하지 않습니다. 인터넷 연결을 끊은 상태에서도 `index.html`을 직접 열 수 있습니다.

## 포함된 결과

- 실험 record: 173개
- `MS vs HD`, `MOGAD vs HD`, `NMOSD vs HD` heatmap
- signed point-biserial `r`, two-sided `p`, within-network BH `q`
- 질환별 집계 sample 수, module feature 수, 분석 summary와 provenance
- 결과를 생성한 Python 코드와 단위 테스트

## 의도적으로 제외한 데이터

- raw/formatted expression matrix
- 개인 sample ID가 포함된 sample-level table
- PyWGCNA `.p`, `.p.gz`, `.p.zst` 객체
- 전체 feature-level kME table

따라서 이 저장소는 **결과 열람과 방법 감사용**입니다. 분석을 처음부터 재실행하려면 별도로 보관된 원본 데이터와 전체 `results_20261119` 실행 디렉토리가 필요합니다.

## 무결성 확인

```powershell
python verify_release.py
```

모든 파일이 `SHA256SUMS.txt`와 일치하면 `OK`를 출력합니다.
