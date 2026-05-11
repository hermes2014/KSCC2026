import pandas as pd
# from pyrefine_indirect import find_ri
from pyrefine_indirect_with_colors import find_ri

df = pd.read_csv("testcase1.csv")
# df = pd.read_csv("testcase2.csv")
# df = pd.read_csv("testcase3.csv")
# df = pd.read_csv("testcase4.csv")
# df = pd.read_csv("testcase5.csv")

x = df["result"].to_numpy()

fit = find_ri(
    x,
    model="BoxCox",        # "BoxCox", "modBoxCoxFast", "modBoxCox"
    n_bootstrap=0,         # 검토 단계에서는 0~30, 최종 분석에서는 충분히 크게
    seed=123,
    n_jobs=1
)

print(fit.summary())
print(fit.get_ri(percentiles=(0.025, 0.975)))


import matplotlib.pyplot as plt

# ax = fit.plot()
ax = fit.plot(
    show_pathological_difference = True,
    healthy_color="green",
    # healthy_color="#00CD00",
    pathological_color="red",
    # pathological_color="#CD0000",
)
plt.show()
