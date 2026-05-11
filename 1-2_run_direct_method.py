####################################### 준비

# !는 Shell Commands
# %는 Magic Commands

%pwd # 현재 작업 디렉토리 경로 출력

!git clone https://github.com/hermes2014/KSCC2026.git

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
    value_col=value_col,
    partition_cols=partition_cols,
    manual_exclude_col=None,
    central_fraction=0.95,
    limit_ci=0.9,
    min_n_nonparametric=120,
    outlier_method="both",
    exclude_flagged_outliers=False,
    units=units,
    decimal_places=1,

    harris_boyd_enabled=True,
    harris_boyd_transform="boxcox", # "none", "log", "boxcox"
    harris_boyd_min_n=3,
    normality_alpha=0.05,    
)

result = ep28.run_analysis(input_path=input_path, sheet_name=sheet_name, cfg=cfg, output_path=output_path)
