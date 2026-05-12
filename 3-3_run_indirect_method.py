####################################### 준비
import pandas as pd
# from pyrefine_indirect import find_ri
from pyrefine_indirect_with_colors import find_ri

####################################### 예제 데이터 선택
df = pd.read_csv("testcase1.csv") # refineR 예제 중 testcase1 사용
# df = pd.read_csv("testcase2.csv")
# df = pd.read_csv("testcase3.csv")
# df = pd.read_csv("testcase4.csv")
# df = pd.read_csv("testcase5.csv")

####################################### 간접법으로 참고범위 추정
x = df["result"].to_numpy()

fit = find_ri(
    x,
    model="BoxCox",        # "BoxCox", "modBoxCoxFast", "modBoxCox"; 기본은 "BoxCox"
    n_bootstrap=0,         # 부트스트랩 반복 횟수; 너무 오래 걸려서 일단 pass, 검토 단계에서는 0~30, 최종 분석에서는 충분히 크게
    seed=123,              # 임의추출시 사용할 seed
    n_jobs=1               # 부트스트랩에 사용할 thread 수 (multi-cpu 지원)
)

print(fit.summary())
print(fit.get_ri(percentiles=(0.025, 0.975)))

####################################### 전체 데이터 중 추정한 정상 집단/비정상 집단 표시
import matplotlib.pyplot as plt

# ax = fit.plot()
ax = fit.plot(
    show_pathological_difference = True,
    healthy_color="green",
    pathological_color="red",
)
plt.show()
