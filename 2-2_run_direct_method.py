####################################### 준비

# !는 Shell Commands
# %는 Magic Commands

%pwd # 현재 작업 디렉토리 경로 출력

!git clone https://github.com/hermes2014/KSCC2026.git # 파이썬 주피터노트북 안에서 셸 커맨드를 실행

%ls ./KSCC2026

%cd ./KSCC2026

%pwd

%ls .

# !pip install pandas numpy scipy openpyxl


####################################### 예제 실행
# import ep28_reference_interval_analysis as ep28
import ep28_reference_interval_analysis_harris_boyd as ep28

input_path = 'CLSI_EP28_Example_ALT.xlsx'
output_path = 'ep28_RI_report_ALT.xlsx'
# input_path = 'CLSI_EP28_Example_Calcium.xlsx'
# output_path = 'ep28_RI_report_Calcium.xlsx'
sheet_name = 'Sheet1'
value_col = 'Result'
partition_cols = ['Sex']
units = 'U/L'

cfg = ep28.AnalysisConfig(
    value_col=value_col, # 검사결과값 컬럼명
    partition_cols=partition_cols, # 파티션에 사용할 컬럼명 리스트
    manual_exclude_col=None,
    central_fraction=0.95, # 참고구간의 중앙 범위; 0.95가 기본
    limit_ci=0.9, # 참고구간의 CI; 0.9를 위해서는 각 120명 이상; 0.95를 위해서는 각 146명 이상; 0.99를 위해서는 각 210명 이상
    min_n_nonparametric=120, # 최소 검체 수
    outlier_method="both", # 아웃라이어 검출법; "none", "tukey", "reed", "both"
    exclude_flagged_outliers=False, # 지정된 아웃라이어 값 제거 여부
    units=units, # 검사 단위
    decimal_places=1, # 참고치 소수점 단위

    harris_boyd_enabled=True, # 파티셔닝 필요성 검사여부
    harris_boyd_transform="boxcox", # 파티셔닝 계산을 위한 변환방법; "none", "log", "boxcox"
    harris_boyd_min_n=3, # 각 파티션의 최소 갯수
    normality_alpha=0.05, # 정규성 검토시 사용할 Shapiro-Wilk alpha값; 0.05가 기본
)

result = ep28.run_analysis(input_path=input_path, sheet_name=sheet_name, cfg=cfg, output_path=output_path)
