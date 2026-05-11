# Q2. refineR에 대한 설명 요청

좋아, 그러면 refineR을 이용한 indirect method로 이 검사의 참고구간(reference interval)을 수립하는게 좋겠어. refineR이 참고구간을 구하는 알고리즘에 대해 자세히 설명해 줄래? refineR의 자세한 정보는 아래 링크를 참고해. https://cran.r-project.org/web/packages/refineR/index.html

# A2. LLM의 응답

refineR 자체 설명과 원 논문/매뉴얼을 기준으로 정리하겠습니다. 알고리즘은 “mixed distribution에서 병적 결과의 영향을 줄이고 중심 기준분포를 재구성하는 절차”로 이해하면 됩니다.

확인해 보니 CRAN 현재 문서 기준으로 refineR 2.0.0은 findRI()가 중심 함수이고, 2021년 알고리즘 논문에 더해 bootstrap CI, RI 검증 관련 기능이 확장되어 있습니다. 아래는 “수식보다 작동 원리” 중심으로 풀어 설명드리겠습니다.

Thought for 32s > 

네. **refineR의 핵심은 “환자값이 섞인 전체 분포에서, 비병적(non-pathological) 성분만을 가장 잘 설명하는 Box-Cox 변환 정규분포를 역방향으로 찾아내는 알고리즘”**이라고 이해하면 됩니다.

CRAN 현재 문서 기준 refineR는 real-world data, 즉 병적·비병적 결과가 섞인 routine 검사값을 입력으로 받아 비병적 분포를 모델링하고, 그 분포에서 참고구간을 산출하는 R 패키지입니다. 주 함수는 `findRI()`이고, 현재 CRAN에는 version 2.0.0이 올라와 있습니다. ([CRAN][1])

## 1. refineR이 가정하는 것

refineR는 다음 세 가지를 전제로 합니다.

첫째, 입력 데이터의 **대부분은 비병적 결과**라고 가정합니다.
둘째, 비병적 결과의 분포는 **Box-Cox transformed normal distribution**으로 표현 가능하다고 봅니다. 즉 원자료에서는 정규분포가 아니어도, 적절한 Box-Cox 변환 후에는 정규분포에 가까워진다고 가정합니다.
셋째, 전체 결과값 범위 중 어딘가에는 **병적 결과의 비율이 매우 낮은 구간**이 존재한다고 가정합니다. 병적 결과의 위치나 모양은 별도로 가정하지 않습니다. ([Nature][2])

이 세 번째 가정이 중요합니다. 예를 들어 CRP, troponin처럼 “정상군 중심 peak”가 뚜렷하고 병적 값이 주로 한쪽 꼬리로 붙는 검사는 비교적 적합할 수 있습니다. 반대로 정상군과 환자군이 넓게 겹치거나, 입원환자 위주라서 병적 결과가 다수이면 불안정해질 수 있습니다.

## 2. 모델의 형태: Box-Cox 변환 정규분포

refineR는 비병적 결과를 다음과 같은 형태로 봅니다.

원자료 검사값을 (x)라고 하면, Box-Cox 변환 후 값 (y)가

[
y \sim N(\mu, \sigma^2)
]

를 따른다고 가정합니다.

Box-Cox 변환은 보통 다음과 같습니다.

[
y =
\begin{cases}
\frac{x^\lambda - 1}{\lambda}, & \lambda \neq 0 \
\log(x), & \lambda = 0
\end{cases}
]

여기서 refineR가 추정하는 주요 모수는 다음입니다.

| 모수        | 의미                                  |
| --------- | ----------------------------------- |
| (\lambda) | 분포의 왜도, 즉 Box-Cox 변환 정도             |
| (\mu)     | 변환 공간에서 비병적 분포의 평균                  |
| (\sigma)  | 변환 공간에서 비병적 분포의 표준편차                |
| (P)       | 전체 데이터 중 비병적 결과의 추정 비율, NP fraction |
| shift     | modified Box-Cox 모델에서 쓰는 이동 모수      |

기본 `model = "BoxCox"`는 1-parameter Box-Cox 모델이고, `model = "modBoxCox"` 또는 `"modBoxCoxFast"`는 shift parameter를 추가한 2-parameter modified Box-Cox 모델입니다. modified Box-Cox는 0에서 떨어져 있는 심한 왜도 분포를 더 잘 다루기 위한 옵션입니다. ([CRAN][3])

## 3. “inverse modeling”이 refineR의 특징

기존 indirect method 중 일부는 먼저 데이터를 Box-Cox 변환하고, 변환된 공간에서 모델을 맞춘 뒤, 다시 원래 단위로 되돌립니다. refineR 논문은 이를 **forward modeling**이라고 설명합니다.

refineR는 반대로, 여러 후보 (\lambda, \mu, \sigma)로부터 **원래 검사값 단위에서 예상되는 분포를 만들어 놓고**, 실제 histogram과 얼마나 잘 맞는지를 비교합니다. 이것을 논문에서는 **inverse modeling approach**라고 설명합니다. 장점은 변환 공간에서 작은 오차가 원래 단위에서 큰 오차로 증폭되는 문제를 줄이려는 데 있습니다. ([Nature][2])

즉 refineR는 이렇게 묻습니다.

> “이 (\lambda, \mu, \sigma) 조합이 진짜 비병적 분포라면, 원래 검사값 축에서 어떤 histogram이 나와야 하는가? 그 예상 histogram이 실제 데이터의 중심부를 얼마나 잘 설명하는가?”

## 4. 알고리즘 단계

### 4-1. 입력 데이터 전처리

refineR 자체는 LIS에서 가져온 데이터를 자동으로 임상적으로 정제하지 않습니다. CRAN vignette도 prefiltering, cleaning, partitioning은 refineR 패키지 밖에서 미리 수행해야 한다고 설명합니다. ([CRAN][3])

따라서 실제 적용 전에는 적어도 다음을 해야 합니다.

검사법·장비·시약 lot 변화가 없는 기간으로 제한하고, 단위가 섞이지 않았는지 확인하며, 검출한계 이하 값, 명백한 입력 오류, 반복측정 편향, 극단값 처리 원칙을 정해야 합니다. 2023년 practical guide도 장비 간 결과 비교 가능성, unit consistency, revalidation 이후/이전 자료 제외 등을 실무 조건으로 제시합니다. ([OUP Academic][4])

### 4-2. main peak 찾기

refineR는 먼저 전체 데이터의 밀도분포를 보고 **main peak**를 찾습니다. 여기서 main peak는 단순히 가장 높은 봉우리가 아니라, **area under the curve가 가장 큰 peak**로 정의됩니다. 이는 “대부분은 비병적 결과”라는 가정에 기반합니다. ([Nature][2])

이 과정에서 여러 (\lambda) 후보에 대해 Box-Cox 변환을 적용하고, 변환 공간에서 분포 밀도를 추정합니다. 논문에서는 density estimation에 `ash` 알고리즘을 사용한다고 설명합니다. ([Nature][2])

### 4-3. (\lambda, \mu, \sigma) 탐색 범위 설정

refineR는 모든 가능한 모수를 무작정 찾지 않고, 먼저 합리적인 search region을 좁힙니다.

(\lambda)는 분포의 왜도를 조절합니다. (\lambda = 0)이면 log-normal 형태, (\lambda = 1)이면 Gaussian 형태에 해당합니다. 논문에서는 (\lambda)의 기본 탐색 범위를 0에서 1.5까지로 설명합니다. 이는 대부분의 검사값 분포가 log-normal, right-skewed, normal 범위에 들어가며, 드문 left-skewed 분포를 일부 허용하기 위한 것입니다. ([Nature][2])

(\mu)와 (\sigma)는 main peak의 폭을 여러 높이에서 측정해 추정합니다. 예를 들어 peak 높이의 50%, 55%, …, 95% 지점에서 분포 폭을 계산하고, 그 중심점들로 (\mu) 탐색 범위를 잡으며, 폭을 표준편차로 환산해 (\sigma) 탐색 범위를 잡습니다. ([Nature][2])

### 4-4. 관심 histogram 구간 만들기

그다음 refineR는 main peak 주변의 관심 구간을 정합니다. 논문에서는 선택된 (\lambda)에서 평균적인 (\mu)와 (\sigma)를 이용해 변환 공간에서 대략

[
\mu \pm 3\sigma
]

범위를 만들고, 이를 inverse Box-Cox 변환으로 원래 검사값 단위로 되돌린다고 설명합니다. 그 범위 안에서 overlapping bins를 이용한 histogram representation을 만듭니다. ([Nature][2])

overlapping bins를 쓰는 이유는 너무 거친 histogram의 noise를 줄이면서도 분포 모양을 부드럽게 표현하기 위해서입니다. ([Nature][2])

### 4-5. 후보 비병적 분포를 만들고 실제 histogram과 비교

각 후보 (\lambda, \mu, \sigma) 조합에 대해 refineR는 Box-Cox transformed normal distribution을 원래 단위로 되돌려 예상 분포 (M)을 만듭니다. 이 예상 분포는 전체 데이터가 아니라 비병적 데이터만 설명해야 하므로, 높이를 (P)라는 scaling factor로 조정합니다. 여기서 (P \leq 1)이고, 최종 출력에서 **NP fraction**으로 보이는 값입니다. ([Nature][2])

예를 들어 출력에서 `NP fraction: 0.909`가 나오면, refineR가 전체 데이터의 약 90.9%를 비병적 분포 성분으로 추정했다는 뜻입니다. CRAN vignette의 예시 출력에서도 `lambda`, `mu`, `sigma`, `shift`, `cost`, `NP fraction`이 함께 제시됩니다. ([CRAN][3])

### 4-6. asymmetric confidence band로 “믿을 만한 bin”만 선택

이 부분이 refineR의 중요한 아이디어입니다.

실제 histogram에는 비병적 결과와 병적 결과가 섞여 있습니다. 병적 결과는 특정 bin의 count를 “늘리는 방향”으로 작용합니다. 어떤 값이 낮은 쪽 병적이든 높은 쪽 병적이든, histogram count 관점에서는 특정 위치의 count가 추가되는 것입니다.

그래서 refineR는 후보 모델 (M) 주변에 **asymmetric confidence band**를 만듭니다. 이 confidence band는 Poisson 통계에 기반하며, 모델이 맞을 때 비병적 histogram count가 포함될 법한 범위를 나타냅니다. 논문은 이 band가 99% 사례에서 비병적 count를 덮도록 구성되며, 병적 결과로 인한 positive bias를 줄이기 위해 비대칭적으로 설계된다고 설명합니다. ([Nature][2])

그다음 실제 histogram bin 중 이 band 안에 들어오는 bin만 “비병적 분포를 잘 반영하는 bin”으로 보고 cost function 계산에 사용합니다. band 밖에 있는 bin은 병적 결과가 섞였을 가능성이 커서 모델 fitting에 덜 반영됩니다. ([Nature][2])

### 4-7. cost function 최소화

각 후보 모델에 대해 refineR는 선택된 bin들에서 **negative log-likelihood with regularization**을 계산합니다. histogram bin count는 Poisson 분포를 따른다고 보고, 관측 count가 후보 모델의 예상 count로부터 나왔을 가능도를 평가합니다. ([Nature][2])

regularization은 높은 count를 가진 bin, 즉 비병적 중심분포일 가능성이 큰 bin에는 보상을 주고, 낮은 count의 bin이나 너무 많은 bin이 선택되는 경우에는 불리하게 작용하도록 설계되어 있습니다. ([Nature][2])

결국 refineR는 다음을 찾습니다.

[
(\lambda^*, \mu^*, \sigma^*, P^*) =
\arg\min \text{cost}
]

이 최적 모수가 “비병적 분포를 가장 잘 설명하는 모델”이 됩니다.

### 4-8. multi-level grid search

최적화는 multi-level grid search로 수행됩니다. 처음에는 비교적 성긴 grid에서 (\lambda, \mu, \sigma) 조합을 탐색하고, cost가 낮은 영역을 찾은 뒤, 그 주변을 더 촘촘한 grid로 다시 탐색합니다. CRAN 문서도 `findRI()`가 multi-level grid search로 negative log-likelihood with regularization을 최소화한다고 설명합니다. ([CRAN][5])

## 5. 참고구간 산출

최적화가 끝나면 refineR는 추정된 비병적 분포에서 원하는 percentile을 계산합니다. 기본은 2.5th와 97.5th percentile입니다. CRAN vignette의 예시에서도 `print(fit)`은 기본적으로 lower limit [2.5%], upper limit [97.5%]를 출력합니다. ([CRAN][3])

개념적으로는 다음과 같습니다.

[
RI_{lower} = \text{invBoxCox}(\mu^* + z_{0.025}\sigma^*)
]

[
RI_{upper} = \text{invBoxCox}(\mu^* + z_{0.975}\sigma^*)
]

즉 “전체 1만 건의 2.5–97.5 percentile”이 아니라, **전체 자료에서 추정한 비병적 분포의 2.5–97.5 percentile**입니다. 이 차이가 indirect method의 핵심입니다.

## 6. bootstrap confidence interval

refineR는 `NBootstrap` 옵션을 통해 참고한계의 신뢰구간을 계산할 수 있습니다. bootstrap은 원자료에서 같은 크기로 복원추출한 데이터를 여러 번 만들고, 각 bootstrap sample에 대해 refineR를 반복 적용한 뒤, 산출된 lower/upper reference limit들의 중앙 구간을 신뢰구간으로 삼는 방식입니다. 논문에서는 예를 들어 200회 반복을 설명하고, CRAN vignette는 실제 real-world analysis에서는 `NBootstrap >= 200`을 권장한다고 설명합니다. ([Nature][2])

실무에서는 다음처럼 해석합니다.

```r
fit <- findRI(Data = x, NBootstrap = 200)
getRI(fit, RIperc = c(0.025, 0.975), pointEst = "medianBS")
```

`pointEst = "medianBS"`는 bootstrap 결과들의 median model을 point estimate로 쓰는 방식이고, CRAN vignette도 bootstrap을 사용한 경우 medianBS 사용을 권장한다고 설명합니다. ([CRAN][3])

## 7. 출력값 해석

`print(fit)` 또는 `getRI(fit)`에서 주로 보는 항목은 다음입니다.

| 항목                  | 해석                          |
| ------------------- | --------------------------- |
| lower limit [2.5%]  | 추정 비병적 분포의 2.5 percentile   |
| upper limit [97.5%] | 추정 비병적 분포의 97.5 percentile  |
| lambda              | Box-Cox 변환 모수               |
| mu                  | 변환 공간에서 비병적 분포 평균           |
| sigma               | 변환 공간에서 비병적 분포 표준편차         |
| shift               | modified Box-Cox에서 쓰는 이동값   |
| cost                | 최적화 목적함수 값                  |
| NP fraction         | 전체 데이터 중 비병적 성분의 추정 비율      |
| CILow/CIHigh        | bootstrap 사용 시 각 참고한계의 신뢰구간 |

CRAN vignette 예시에서도 `getRI()`는 percentile별 point estimate, CI, uncertainty margin을 출력하고, `plot()`은 원자료 histogram과 추정된 비병적 분포 곡선, 참고한계를 함께 보여줍니다. ([CRAN][3])

## 8. refineR가 “환자값을 제거”하는 방식은 아니다

중요한 점은 refineR가 개별 결과값을 하나씩 “정상/비정상”으로 분류하는 알고리즘은 아니라는 것입니다. 전체 분포를 보고 **비병적 성분의 통계적 모델**을 추정합니다.

plot에서 `showPathol = TRUE`를 쓰면 병적 분포처럼 보이는 빨간 곡선을 표시할 수 있지만, CRAN vignette는 이것이 raw histogram과 추정 비병적 모델의 차이를 표시한 것이므로 해석에 주의해야 한다고 설명합니다. ([CRAN][3])

따라서 refineR 결과를 보고 “이 1만 건 중 어느 환자가 건강인인가”를 말할 수는 없습니다. 말할 수 있는 것은 “이 mixed distribution에서 비병적 성분은 대략 이런 분포로 추정된다”입니다.

## 9. 성능과 한계

원 논문에서 refineR는 simulated datasets에서 kosmic, direct method N=120, direct method N=400과 비교되었습니다. 평균 percentage error는 refineR가 가장 낮았고, ±1 total error deviation 기준에서는 direct method N=400보다는 낮았지만 kosmic 및 direct method N=120보다 좋은 성능을 보였다고 보고되었습니다. ([Nature][2])

다만 성능은 병적 분포의 위치와 비율에 따라 달라집니다. 논문은 병적 분포와 비병적 분포의 겹침이 크거나, 병적 결과 비율이 매우 높은 경우에는 ground truth에서 점점 벗어났다고 설명합니다. 특히 다수의 검체가 비병적이라는 가정과 병적 결과가 거의 없는 구간이 존재한다는 가정이 깨지면 결과가 불안정해질 수 있습니다. ([Nature][2])

## 10. 선생님 상황에 적용하면

검사값 1만 건만 있고 건강인/환자 구분이 없다면, refineR는 매우 적절한 후보입니다. 다만 “검사값만” 있다는 조건 때문에 다음과 같이 접근하는 것이 좋겠습니다.

1. 먼저 전체 분포를 그려서 단봉성 main peak가 있는지 확인합니다.
2. 기본 `model = "BoxCox"`로 `findRI()`를 실행합니다.
3. 왜도가 크거나 하한 쪽 모양이 부자연스러우면 `model = "modBoxCox"`도 함께 실행합니다.
4. `NBootstrap = 200` 이상으로 참고한계의 CI를 산출합니다.
5. `plot(fit)`으로 raw histogram과 추정 비병적 분포가 중심부를 잘 설명하는지 확인합니다.
6. 기존 제조사 RI, 문헌 RI, 임상적 flag rate와 비교합니다.
7. 나이·성별 정보가 없으므로 partitioned RI가 불가능하다는 한계를 보고서에 명시합니다.

실제 R 코드의 골격은 이 정도입니다.

```r
library(refineR)

# x: numeric vector, 검사 결과값 10,000건
x <- x[is.finite(x)]
x <- x[x > 0]  # Box-Cox 적용상 0 이하 값이 있으면 별도 검토 필요

fit.bc <- findRI(Data = x, model = "BoxCox", NBootstrap = 200, seed = 123)
print(fit.bc)
getRI(fit.bc, RIperc = c(0.025, 0.975), pointEst = "medianBS")
plot(fit.bc, pointEst = "medianBS", showCI = TRUE)

fit.mbc <- findRI(Data = x, model = "modBoxCox", NBootstrap = 200, seed = 123)
print(fit.mbc)
getRI(fit.mbc, RIperc = c(0.025, 0.975), pointEst = "medianBS")
plot(fit.mbc, pointEst = "medianBS", showCI = TRUE)
```

정리하면, refineR는 **“전체 결과값의 percentile을 구하는 방법”이 아니라, mixed routine data 안에서 비병적 분포를 가장 잘 설명하는 Box-Cox 기반 모델을 찾고, 그 모델의 percentile을 참고구간으로 쓰는 방법**입니다. 선생님처럼 건강인/환자 구분 없이 검사 결과값만 1만 건 있는 상황에서 사용할 수 있는 대표적인 indirect RI 알고리즘입니다.

[1]: https://cran.r-project.org/web/packages/refineR/index.html "CRAN: Package refineR"
[2]: https://www.nature.com/articles/s41598-021-95301-2 "refineR: A Novel Algorithm for Reference Interval Estimation from Real-World Data | Scientific Reports"
[3]: https://cran.r-project.org/web/packages/refineR/vignettes/refineR_package.html "refineR: Reference Interval Estimation using Real-World Data (RWD)"
[4]: https://academic.oup.com/jalm/article/8/1/84/6847917?utm_source=chatgpt.com "Estimation of Reference Intervals from Routine Data Using the ..."
[5]: https://cran.r-project.org/web/packages/refineR/refman/refineR.html "Help for package refineR"


