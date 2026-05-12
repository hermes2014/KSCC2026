getwd() # 현재 작업 디렉토리 위치 확인

# git clone https://github.com/hermes2014/KSCC2026.git # GitHub에서 코드와 자료 내려받기 (colab 터미널 창에서 실행)

setwd("/content/KSCC2026") # 작업 디렉토리 변경

getwd() # 작업 디렉토리가 잘 변경되었는지 확인

# refineR 패키지 설치 (설치되지 않은 경우)
pkg <- "refineR"
if (!requireNamespace(pkg, quietly = TRUE)) {
  install.packages(pkg)
}

# refineR 패키지 로딩
library(refineR) 


# refineR 도움말 및 문서 페이지 보기 (패키지 설치 확인용)
?refineR

# refineR 예제 로딩
data(testcase1)
data(testcase2)
data(testcase3)
data(testcase4)
data(testcase5)

# 검사 결과의 예로 refineR 예제 입력
data <- testcase1

# 참고구간 추정
fit <- findRI(Data=data)

# 참고구간 상/하한의 신뢰구간을 추정하기 위해 bootstrapping 수행
fit.bs <- findRI(Data=data, NBootstrap=0) # 실해에 시간이 너무 오래 걸려서 일단 패스, 실제는 NBootstrap=200 권장


# 추정한 참고구간의 정보 출력
print(fit)
 
# bootstrapping을 포함해서 추정한 참고구간의 정보 출력
print(fit.bs)

# RIperc로 주어진 범위로 참고구간 만들기 (default: 0.025 ~ 0.975)
# CIprop: 참고범위 신뢰구간
# pointEst: "fullDataEst" (Default), "medianBS" (권장), "meanBS"
getRI(fit.bs, RIperc=c(0.025, 0.975), CIprop=0.95, pointEst="medianBS")  

# 그래프 그리기
plot(fit) # 기본형
plot(fit, showPathol = TRUE) # pathologic population도 추가
 



# (Optional) 원데이터의 치우침이 심할 경우는 아래와 같이 modified BoxCox (2-p BoxCox) 모델링 사용
fit.mbc <- findRI(Data=data, model="modBoxCox")
print(fit.mbc)
plot(fit.mbc)
