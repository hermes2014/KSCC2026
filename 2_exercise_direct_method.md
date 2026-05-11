## Q. 참고구간 설정을 위한 파이썬 프로그램 요청

방금 요약해준 CLSI EP28-A3c의 내용에 따라서 참고집단(reference population)을 대상으로 검사를 수행하고 결과를 수집했어. 이제 통계분석을 수행해야 하는데, 가이드라인에 따라 통계분석을 수행해서 이 검사의 참고구간을 제시하는 파이썬 프로그램을 작성해줄래? 검사결과는 엑셀파일로 준비했어.

## A. LLM의 응답

가능합니다. 엑셀의 컬럼명이 검사실마다 다를 수 있어서, 실행 시 결과값 컬럼, 선택적으로 성별/연령군 같은 partition 컬럼을 인자로 지정하는 방식의 재사용 가능한 스크립트로 만들겠습니다. 이상치는 자동 삭제하지 않고 “flag 후 검토”가 기본값이 되도록 구성하겠습니다.



파이썬 스크립트를 작성해두었습니다. 핵심 기능은 다음과 같습니다.

* 엑셀 파일에서 검사결과 컬럼을 읽어 참고구간 계산
* 성별, 연령군 등 `partition`별 분석 가능
* 기본값: **중앙 95% 참고구간**, 즉 2.5–97.5 percentile
* 각 reference limit의 **90% confidence interval** 계산
* `n ≥ 120` 여부를 partition별로 판정
* Tukey 및 Reed/Dixon 방식의 이상치 후보 flag
* 이상치는 기본적으로 **자동 제외하지 않음**
* 결과를 `RI_summary`, `Data_with_flags`, `Audit` 시트가 포함된 엑셀 보고서로 저장

실행 예시는 아래와 같습니다.

```bash
pip install pandas numpy scipy openpyxl

python ep28_reference_interval_analysis.py data.xlsx \
  --sheet Sheet1 \
  --value-col Result \
  --partition-cols Sex AgeGroup \
  --units U/L \
  --out ep28_RI_report.xlsx
```

전체를 하나의 참고집단으로 분석하려면:

```bash
python ep28_reference_interval_analysis.py data.xlsx \
  --value-col Result \
  --units U/L
```

수동 제외 컬럼이 있다면:

```bash
python ep28_reference_interval_analysis.py data.xlsx \
  --value-col Result \
  --manual-exclude-col Exclude
```

엑셀 파일을 업로드해주면, 컬럼 구조에 맞춰 이 스크립트를 바로 실행 가능한 형태로 조정하거나 실제 참고구간 보고서까지 만들어드릴 수 있습니다.
