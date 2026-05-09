# Q. LLM에게 가이드라인 요약을 요청하기
나는 임상검사실을 운영하는 진단검사의학과 전문의야. 이번에 새로 세팅한 검사의 참고구간(reference interval)을 설정해야 하는데, 관련 가이드라인(CLSI EP28-A3C)을 따라서 이 작업을 진행하려고 해. 저 가이드라인을 자세하게 정리해 줄 수 있니?

# A. LLM의 응답

(가이드라인 자체는 저작권 문서라 원문을 길게 인용하지는 않고, 공개적으로 확인 가능한 자료와 표준적인 통계 절차를 바탕으로 실무에 바로 쓸 수 있게 구조화해서 정리하겠습니다. 먼저 현재 문서 상태와 핵심 원칙을 확인해 보겠습니다.)

Thought for 1m 23s>

네. **CLSI EP28-A3c**는 임상검사실에서 **정량검사의 참고값(reference values), 참고한계(reference limits), 참고구간(reference intervals, RI)을 설정·이전·검증**하는 방법을 다루는 핵심 가이드라인입니다. 현재 CLSI 페이지 기준으로 EP28은 **2010년 10월 19일 발간된 제3판**이고, **2020년 4월 내용 개정 없이 재확인**된 문서입니다. CLSI-IFCC 공동 프로젝트이며, FDA가 규제 요구사항 충족에 사용할 수 있는 consensus standard로 인정한 문서로 안내되어 있습니다. ([CLSI](https://clsi.org/shop/standards/ep28/ "EP28 | Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory"))

## 1. EP28-A3c의 핵심 개념

**Reference interval은 “정상범위”가 아니라, 정의된 참고 모집단에서 얻은 검사값 분포의 일정 구간**입니다. 일반적으로 양측 참고구간은 건강 관련 참고집단의 **중앙 95%**, 즉 하위 2.5백분위수부터 상위 97.5백분위수까지로 설정합니다. 단측 해석이 임상적으로 타당한 검사는 하한 또는 상한 한쪽만 설정할 수 있습니다. ([De Gruyter Brill](https://www.degruyterbrill.com/document/doi/10.1515/tjb-2017-0299/html?srsltid=AfmBOoqGQPuSVHQrngUDn5uNvVNZ4AI8PCcVOzC_K7togvNefn-dKmbR "Establishing and using reference intervals"))

중요한 점은 **reference interval과 clinical decision limit은 다르다**는 것입니다. 예를 들어 HbA1c, 콜레스테롤처럼 국제적·국가적 합의나 outcome 기반으로 정해진 절단값이 있는 검사는 새로 참고구간을 설정하는 것이 아니라, 검사법의 정확도·traceability·QC·PT를 통해 그 decision limit을 신뢰성 있게 적용하는 데 초점을 둬야 합니다. ([Bfldr](https://cdn.bfldr.com/YLD4EVFU/at/scjm7t2858p8thq2gkbwb8/ep28a3ce_sample.pdf "EP28: Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory, 3rd Edition"))

EP28-A3c가 다루는 큰 흐름은 다음 세 가지입니다.

| 상황  | 권장 접근 |
| --- | --- |
| 새 분석물 또는 새 검사법이고 신뢰할 외부 RI가 없음 | 직접 설정, establishment |
| 제조사, 문헌, 다기관 연구 RI가 있고 내 검사실에 적용 가능해 보임 | 이전 및 검증, transference/verification |
| 이미 합의된 clinical decision limit이 있음 | RI 설정보다 정확도·표준화·검사 성능 관리 |

CLSI 공식 설명에서도 EP28은 참고 대상자 선택, 전처리·분석 요인, 계산법, 참고구간 이전, 결과 제시와 활용을 포괄한다고 설명합니다. ([CLSI](https://clsi.org/shop/standards/ep28/ "EP28 | Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory"))

---

## 2. 직접 참고구간을 설정하는 절차

직접 설정은 가장 원칙적인 방법입니다. EP28-A3c는 **각 partition, 예: 성별·연령군별로 최소 120명 이상의 적격 참고개체를 확보하여 비모수법으로 계산하는 것**을 전통적 최선의 방법으로 제시합니다. ([Bfldr](https://cdn.bfldr.com/YLD4EVFU/at/scjm7t2858p8thq2gkbwb8/ep28a3ce_sample.pdf "EP28: Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory, 3rd Edition"))

### 2.1 검사와 intended use 정의

먼저 아래 항목을 고정해야 합니다.

| 항목  | 정해야 할 내용 |
| --- | --- |
| Measurand | 무엇을 측정하는가 |
| 검체  | serum, plasma, whole blood, urine 등 |
| 검사법 | analyzer, reagent, calibrator, lot, traceability |
| 단위  | SI 또는 통상 단위 |
| 결과 보고 형식 | 소수점 자리, flag 기준 |
| 임상 사용 목적 | 선별, 진단 보조, 모니터링, 위험도 평가 등 |

검사법이 아직 안정화되지 않았거나 calibration, precision, bias, analytical measuring range, carryover, interference 검증이 끝나지 않았다면 참고구간 설정을 먼저 진행하면 안 됩니다. RI는 분석시스템이 안정적이라는 전제 위에서만 의미가 있습니다.

### 2.2 참고 모집단과 partition 정의

“건강한 사람”을 막연히 모집하는 방식은 좋지 않습니다. EP28-A3c의 핵심은 **reference population을 사전에 정의**하는 것입니다.

예를 들어 성인 혈청 검사라면 다음을 명확히 해야 합니다.

| 구분  | 예시  |
| --- | --- |
| 포함 기준 | 20–65세, 특정 지역 거주, 자발적 동의, 일반 건강상태 양호 |
| 제외 기준 | 급성질환, 만성질환, 임신, 특정 약제, 음주·흡연 기준, 최근 수술·입원, 비정상 선별검사 |
| 생리적 요인 | 성별, 연령, 월경주기, 폐경, 임신, BMI |
| 생활 요인 | 식이, 운동, 음주, 흡연 |
| 채혈 조건 | 공복 여부, 채혈 시간, 자세, tourniquet 시간 |

Partition은 생물학적 타당성과 데이터 양쪽으로 판단합니다. 대표적으로 성별, 연령, 소아 성장 단계, 임신, 폐경, 신기능 상태, 검체 종류가 partition 후보가 됩니다. 다만 partition을 너무 많이 나누면 각 partition별 120명 확보가 어려워지고, 통계적 불확실성이 커집니다.

### 2.3 참고개체 모집과 설문

EP28-A3c는 참고 대상자 선택, 제외 기준, sample questionnaire를 별도 절차로 다룹니다. 목차상으로도 reference individuals의 selection, exclusion/partitioning, questionnaire가 핵심 장으로 배치되어 있습니다. ([Bfldr](https://cdn.bfldr.com/YLD4EVFU/at/scjm7t2858p8thq2gkbwb8/ep28a3ce_sample.pdf "EP28: Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory, 3rd Edition"))

실무적으로는 다음을 권합니다.

| 단계  | 실무 포인트 |
| --- | --- |
| 사전 설문 | 병력, 약물, 보충제, 음주, 흡연, 임신, 최근 감염, 운동 |
| 기본 선별 | 필요 시 CBC, CRP, 간·신장 기능 등 |
| 제외 기준 적용 | 사전에 정한 기준으로만 제외 |
| 기록 보존 | 제외 사유, 채혈일, 검체상태, 분석일, lot 기록 |

검사 특성상 “진정한 건강인” 확보가 어렵다면, 선별검사 또는 문진 기준을 더 정교하게 설계해야 합니다. 소아·노인·임산부에서는 윤리적·실무적 제약이 크므로 다기관 RI 또는 문헌 RI의 verification 접근이 현실적일 수 있습니다.

### 2.4 전처리 및 분석 조건 통제

Reference interval은 biological variation뿐 아니라 preanalytical variation의 영향을 크게 받습니다. 따라서 실제 환자검체 처리와 동일한 조건으로 reference samples를 다뤄야 합니다.

확인해야 할 전처리 항목은 다음과 같습니다.

| 영역  | 예시  |
| --- | --- |
| 대상자 준비 | 공복, 채혈 시간, 자세, 운동 제한 |
| 채혈  | tube 종류, 항응고제, tourniquet 시간 |
| 운송  | 온도, 시간, 원심 전 대기시간 |
| 처리  | 원심 조건, 분리 시간, hemolysis/lipemia/icterus 기준 |
| 보관  | 분석 전 보관 온도, freeze-thaw 허용 여부 |
| 분석  | 같은 장비·시약·calibrator·QC 체계 사용 |

검체 안정성이 낮은 analyte는 모집보다 logistics가 더 중요합니다. 채혈 후 분석까지의 허용시간과 보관 조건을 SOP로 먼저 고정해야 합니다.

---

## 3. 통계 분석 원칙

### 3.1 기본 계산

가장 표준적인 방식은 **비모수 percentile method**입니다. 양측 95% RI는 정렬된 reference values에서 2.5th percentile과 97.5th percentile을 구합니다. CLSI 방식에서는 각 reference limit에 대해 보통 **90% confidence interval**을 함께 제시하는 것이 권장됩니다. ([MedCalc](https://www.medcalc.org/en/manual/referenceinterval.php "Reference Interval (Reference Range) Calculation | MedCalc"))

실무 보고 예시는 다음과 같습니다.

> 검사명: X  
> 검체: Serum  
> 분석법: Y analyzer / Z reagent  
> 참고집단: 성인 남성 20–64세, n=132  
> Reference interval: 4.1–8.9 U/L  
> 90% CI for lower limit: 3.8–4.4 U/L  
> 90% CI for upper limit: 8.4–9.5 U/L  
> 산출법: Nonparametric percentile method, CLSI EP28-A3c 기반

### 3.2 표본수

| 방법  | 표본수 |
| --- | --- |
| 비모수법 | 권장: partition당 ≥120 |
| Robust method | 120 미만일 때 대안 가능, 단 통계 역량 필요 |
| Verification | 보통 partition당 최소 20명 |

비모수 percentile method에서 120명은 단순 관행이 아니라, reference limit과 그 confidence interval을 계산할 수 있는 최소 수준으로 이해하는 것이 좋습니다. 120명보다 많을수록 reference limit의 confidence interval이 좁아져 신뢰도가 올라갑니다. ([MedCalc](https://www.medcalc.org/en/manual/referenceinterval.php "Reference Interval (Reference Range) Calculation | MedCalc"))

### 3.3 이상치 처리

이상치는 자동 삭제하면 안 됩니다. 먼저 hemolysis, 검체오류, 약물, 숨겨진 질환, 분석오류, 입력오류를 확인해야 합니다. EP28-A3c 기반 실무 문헌에서는 Reed/Dixon 또는 Tukey 방법을 사용해 outlier를 평가하고, 제거한 경우 replacement sample을 확보하는 절차를 설명합니다.

간단히 말하면:

| 상황  | 처리  |
| --- | --- |
| 명백한 검체·분석·입력 오류 | 제외 가능, 사유 기록 |
| 통계적 outlier이나 생물학적으로 설명 가능 | 무조건 제외하지 말고 검토 |
| 제외한 경우 | 해당 partition의 n을 유지하도록 대체 검체 확보 |
| 반복적 outlier 발생 | 모집 기준, hidden disease, partition 필요성 재검토 |

### 3.4 분포 확인과 변환

분석 전에 histogram, box plot, Q-Q plot을 확인합니다. 정규분포가 아니어도 n≥120이면 비모수법을 쓰면 됩니다. 다만 소표본에서 parametric 또는 robust approach를 고려할 때는 log transformation, Box-Cox transformation 등 분포 변환이 필요할 수 있습니다. Analyse-it의 EP28-A3c workflow 설명도 nonparametric, robust, bootstrap, transformation, partitioning, Tukey outlier screening, transfer/verification을 포함합니다. ([Analyse-it](https://analyse-it.com/products/method-validation/reference-intervals "Reference intervals per CLSI EP28-A3C | Analyse-it®"))

---

## 4. 외부 참고구간의 이전과 검증

현실적으로 많은 검사실은 직접 120명 이상을 모집하기 어렵습니다. EP28-A3c도 실제 검사실에서는 직접 RI 설정이 드물다는 점을 인정하고, **잘 설정된 외부 RI를 내 검사실에 적용 가능한지 검증하는 접근**을 중요하게 다룹니다. ([Bfldr](https://cdn.bfldr.com/YLD4EVFU/at/scjm7t2858p8thq2gkbwb8/ep28a3ce_sample.pdf "EP28: Defining, Establishing, and Verifying Reference Intervals in the Clinical Laboratory, 3rd Edition"))

### 4.1 먼저 외부 RI의 품질을 평가

외부 RI를 그대로 쓰기 전에 다음을 확인해야 합니다.

| 항목  | 질문  |
| --- | --- |
| 참고집단 | 내 검사실 환자군/지역/인구학과 유사한가 |
| 검사법 | 동일 장비·시약·calibrator인가 |
| traceability | 표준화 또는 commutability가 확보되었나 |
| 검체  | serum/plasma/whole blood 등 동일한가 |
| 전처리 | 공복, 채혈시간, 보관조건이 유사한가 |
| partition | 성별·연령 구분이 내 population에 맞는가 |
| 통계법 | n, outlier 처리, percentile/CI 산출법이 적절한가 |
| 최신성 | 오래된 방법·population에 근거한 RI는 아닌가 |

외부 RI의 method와 population이 내 검사실과 다르면, 단순 verification 전에 method comparison 또는 직접 설정이 필요할 수 있습니다. 문헌상 EP28-A3c 기반 transference는 분석법과 population의 비교 가능성을 핵심으로 보며, 방법 비교에서는 regression의 slope가 1에 가깝고 intercept가 0에 가까우며 상관성이 높아야 한다고 설명됩니다. ([De Gruyter Brill](https://www.degruyterbrill.com/document/doi/10.1515/tjb-2017-0299/html?srsltid=AfmBOoqGQPuSVHQrngUDn5uNvVNZ4AI8PCcVOzC_K7togvNefn-dKmbR "Establishing and using reference intervals"))

### 4.2 20명 verification 절차

가장 실용적인 EP28-A3c 방식은 **local healthy reference individuals 20명으로 검증**하는 것입니다. 검증 전에는 전처리, 분석법, local population이 원 RI와 충분히 유사한지 확인해야 합니다.

절차는 다음과 같습니다.

| 단계  | 내용  |
| --- | --- |
| 1   | 해당 partition에 맞는 local healthy individuals 20명 모집 |
| 2   | 원 RI 설정 때와 유사한 조건으로 검체 채취·분석 |
| 3   | Reed/Dixon 또는 Tukey 등으로 outlier 확인 |
| 4   | outlier를 제거했다면 replacement sample 확보 |
| 5   | 20개 결과 중 외부 RI 밖에 있는 값의 개수 계산 |

판정은 다음과 같이 적용합니다.

| 20명 결과 | 해석  |
| --- | --- |
| RI 밖 ≤2개 | 해당 RI를 수용 가능 |
| RI 밖 3–4개 | 새 20명으로 한 번 더 검증 |
| RI 밖 ≥5개 | 해당 RI는 부적합 가능성이 큼 |

두 번째 20명을 검사했을 때도 RI 밖 값이 3개 이상이면, 분석법 차이, 전처리 차이, population 차이, partition 필요성, 원 RI의 품질을 재검토하고 직접 설정을 고려해야 합니다. 이 흐름은 EP28-A3c verification 절차를 설명한 문헌과 도식에서 동일하게 제시됩니다.

---

## 5. 새로 세팅한 검사에 대한 실무 적용 전략

새 검사를 도입한 상황이라면, 저는 보통 아래 순서로 접근하는 것이 가장 안전하다고 봅니다.

### Step 1. 이 검사가 RI 대상인지 decision limit 대상인지 결정

먼저 이 검사의 해석이 건강인 분포 기반인지, outcome 또는 consensus 기반 cutoff인지 확인합니다.

| 유형  | 예   |
| --- | --- |
| Reference interval 대상 | 대부분의 chemistry, hematology 정량검사 |
| Clinical decision limit 대상 | HbA1c, lipids 일부, troponin 일부, 약물치료농도 등 |

Decision limit 대상이면 RI 설정보다 traceability, LoQ, imprecision, bias, lot-to-lot, QC 전략이 더 중요합니다.

### Step 2. 제조사 RI 또는 문헌 RI를 검토

제조사 insert, peer-reviewed multicenter study, 국내외 데이터, 유사 장비 RI를 수집합니다. 단, “제조사가 제시했으므로 사용 가능”이 아니라, 내 검사실 조건과 비교하여 **adopt 가능한 RI 후보**인지 평가해야 합니다.

### Step 3. 가능하면 transference/verification으로 시작

검사법이 제조사와 동일하거나 traceable하고, population 차이가 크지 않다면 20명 verification이 가장 현실적입니다. 특히 새 검사 세팅 초기에는 직접 120명 이상을 모집하기보다, 외부 RI의 품질을 평가하고 local verification을 수행하는 방식이 효율적입니다.

### Step 4. verification 실패 시 원인 분석

실패했다고 곧바로 “우리 population이 다르다”고 결론 내리면 안 됩니다. 다음 순서로 확인합니다.

| 원인 후보 | 확인 사항 |
| --- | --- |
| 분석 문제 | calibration, reagent lot, QC, bias, carryover |
| 전처리 차이 | fasting, tube, centrifugation, storage |
| 검체 matrix | serum vs plasma, anticoagulant 영향 |
| 모집 문제 | hidden disease, 약물, 연령 분포 |
| partition 문제 | 성별·연령 구분 필요성 |
| 외부 RI 문제 | n 부족, 오래된 방법, 부적절한 통계 |

### Step 5. 직접 설정이 필요한 경우

다음이면 de novo establishment를 고려해야 합니다.

| 상황  | 이유  |
| --- | --- |
| 새 분석물이고 신뢰할 RI가 없음 | 외부 근거 부족 |
| 검사법이 기존법과 matrix/bias가 크게 다름 | transfer 부적절 |
| local verification 반복 실패 | population 또는 method 차이 가능 |
| 특정 특수집단 검사 | 소아, 임산부, 노인, 특정 질환 배제집단 |
| 규제·인증상 자체 근거 필요 | 문서화 필요성 |

---

## 6. 결과 보고와 문서화

RI는 계산보다 문서화가 중요합니다. 검사실 SOP와 validation file에 다음을 남기는 것이 좋습니다.

| 문서 항목 | 포함 내용 |
| --- | --- |
| RI 출처 | 직접 설정, 제조사, 문헌, 다기관 연구 |
| 대상 검사법 | 장비, 시약, calibrator, lot, traceability |
| 검체 조건 | 검체 종류, 채혈·보관·분석 조건 |
| 참고집단 | 포함·제외 기준, 인구학적 특성 |
| partition | 성별, 연령, 기타 구분 근거 |
| 통계법 | 비모수법, robust method, outlier 처리법 |
| 표본수 | 전체 및 partition별 n |
| 결과  | reference limits, 90% CI, 단위 |
| verification 결과 | 20명 결과, RI 밖 개수, outlier 처리 |
| 승인  | 책임 전문의 검토 및 승인일 |
| 재검토 기준 | 장비·시약 변경, method comparison 결과, population 변화 |

환자 보고서에는 최소한 **검사명, 결과, 단위, 해당 환자에게 적용되는 RI, flag, 필요 시 성별·연령 기준**이 명확해야 합니다. RI가 검사법 의존적인 경우에는 장비·검사법 변경 시 비교 가능성도 관리해야 합니다.

---

## 7. 주의해야 할 함정

가장 흔한 문제는 다음입니다.

| 함정  | 왜 문제인가 |
| --- | --- |
| 제조사 RI를 무검증 사용 | local population·preanalytics 차이 반영 안 됨 |
| “정상범위”라는 표현 | 건강/질병을 이분법적으로 오해시킴 |
| decision limit과 RI 혼동 | HbA1c, lipid 등에서 잘못된 기준 적용 가능 |
| partition 과다 | n 부족과 불안정한 limit 발생 |
| partition 부족 | 성별·연령 차이를 숨김 |
| outlier 자동 삭제 | 실제 생물학적 변이를 제거할 위험 |
| 오래된 문헌 RI 사용 | 검사법·population 변화 반영 안 됨 |
| 혈청/혈장 혼용 | matrix effect 가능 |
| verification 실패 후 원인분석 생략 | 잘못된 RI 채택 또는 불필요한 자체 연구로 이어짐 |

---

## 8. 새 검사 RI 설정용 간단 체크리스트

| 확인  | 항목  |
| --- | --- |
| □   | 이 검사는 RI 대상인가, decision limit 대상인가 |
| □   | 분석 성능 검증은 완료되었는가 |
| □   | 검체 종류와 전처리 조건이 고정되었는가 |
| □   | 외부 RI 후보의 출처와 품질을 검토했는가 |
| □   | 내 검사실 method와 population에 적용 가능한가 |
| □   | 필요한 partition을 사전에 정의했는가 |
| □   | 각 partition별 20명 verification 또는 120명 establishment 계획이 있는가 |
| □   | outlier 처리 기준을 사전에 정했는가 |
| □   | reference limits와 90% CI를 산출했는가 |
| □   | 책임 전문의가 검토·승인했는가 |
| □   | 장비·시약·calibrator 변경 시 재검토 기준을 정했는가 |

요약하면, EP28-A3c의 실무 메시지는 **"직접 설정은 partition당 ≥120명의 건강 참고개체로 비모수법을 적용하는 것이 원칙이지만, 현실적인 검사실에서는 신뢰할 수 있는 외부 RI를 엄격히 평가한 뒤 local population에서 최소 20명으로 검증하는 접근이 매우 중요하다**입니다.
