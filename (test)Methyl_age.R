
# DNA 연령 관련 분석을 위한 작업 디렉토리 설정
setwd("C:/Users/ko911/OneDrive/바탕 화면/UnsortedStudy/DNAage")

# 필요한 패키지 로드
library(WGCNA)  # Weighted gene co-expression network analysis
library(sqldf)  # SQL 데이터베이스 쿼리를 사용하여 데이터 처리
library(dplyr)
library(tidyverse)
## Bioconductor 설치
# install.packages("BiocManager")

## "impute" 패키지 설치 (Bioconductor 패키지)
#if (!requireNamespace("BiocManager", quietly = TRUE)) {
#  install.packages("BiocManager")
#}
#BiocManager::install("impute")

## devtools 패키지 설치 (만약 설치되지 않은 경우)
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}

# install.packages("RPMM")
# library(RPMM)
source("C:/Users/ko911/OneDrive/바탕 화면/UnsortedStudy/DNAage/AdditionalFile24NROMALIZATION.R.txt")

# 나이 변환과 프로브 주석 함수 정의
trafo= function(x, adult.age=20) { x = (x + 1) / (1 + adult.age); y = ifelse(x <= 1, log(x), x - 1); y }
anti.trafo= function(x, adult.age=20) { ifelse(x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age) }

# 21k DNA 메틸화 데이터와 프로브 주석 데이터 읽기
probeAnnotation21kdatMethUsed = read.csv("AdditionalFile22probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k = read.csv("AdditionalFile21datMiniAnnotation27k.csv")
datClock = read.csv("AdditionalFile23predictor.csv")

# sqldf 패키지 로드
library(sqldf)

# DNA 메틸화 데이터 (beta 값) 읽기
dat0 = read.csv.sql("AdditionalFile26MethylationDataExample55.csv")
nSamples = dim(dat0)[[2]] - 1
nProbes = dim(dat0)[[1]]
dat0[,1] = gsub(x = dat0[,1], pattern = "\"", replacement = "")

# 로그 파일 생성 및 오류 체크
file.remove("LogFile.txt")
file.create("LogFile.txt")
DoNotProceed = FALSE
cat(paste("메틸화 데이터에는", nSamples, "개의 샘플 (배열)과", nProbes, "개의 프로브가 있습니다."), file = "LogFile.txt")
if (nSamples == 0) {
  DoNotProceed = TRUE
  cat("\n ERROR: 샘플이 없는 것 같습니다. 데이터 파일을 올바르게 입력했는지 확인하세요.", file = "LogFile.txt", append = TRUE)
}
if (nProbes == 0) {
  DoNotProceed = TRUE
  cat("\n ERROR: 프로브가 없는 것 같습니다. 데이터 파일을 올바르게 입력했는지 확인하세요.", file = "LogFile.txt", append = TRUE)
}
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
  nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
  for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
  if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="LogFile.txt",append=TRUE)  } 
  XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
  selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
  selectXchromosome[is.na(selectXchromosome)]=FALSE
  meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
  if (   sum(selectXchromosome) >=500 )  {
    meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
  if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="LogFile.txt",append=TRUE)  } 
}


match1 = match(probeAnnotation21kdatMethUsed$Name, dat0[, 1])
if (sum(is.na(match1)) > 0) { 
  missingProbes = probeAnnotation21kdatMethUsed$Name[!is.element(probeAnnotation21kdatMethUsed$Name, dat0[, 1])] 
  DoNotProceed = TRUE 
  cat(paste("\n \n Input error: You forgot to include the following", length(missingProbes), "CpG probes (or probe names):\n", paste(missingProbes, sep = "", collapse = ", ")), file = "LogFile.txt", append = TRUE) 
}

# 단계 2: 21k 프로브로 데이터 제한 및 숫자 형식 확인
match1 = match(probeAnnotation21kdatMethUsed$Name, dat0[, 1])
if (sum(is.na(match1)) > 0) stop(paste(sum(is.na(match1)), "CpG probes cannot be matched"))
dat1 = dat0[match1, ]
asnumeric1 = function(x) { as.numeric(as.character(x)) }
dat1[, -1] = apply(as.matrix(dat1[, -1]), 2, asnumeric1)

# 단계 3: 결과 파일 datout 생성
set.seed(1)
# 데이터 정규화 수행 여부 (권장)
normalizeData = TRUE
source("AdditionalFile25StepwiseAnalysis.txt")


# 단계 4: 결과 출력
if (sum(datout$Comment != "") == 0) { 
  cat("\n 개별 샘플은 정상적으로 처리되었습니다.", file = "LogFile.txt", append = TRUE) 
} 
if (sum(datout$Comment != "") > 0) { 
  cat(paste("\n 경고: 다음 샘플에 대해 경고가 생성되었습니다.\n", datout[, 1][datout$Comment != ""], "\n 자세한 내용은 로그 파일을 확인하세요."), file = "LogFile.txt", append = TRUE) 
}

# 결과를 디렉토리에 출력
write.table(datout, "Output.csv", row.names = FALSE, sep = ",")

######test###############

# 이 작업을 수행하기 위해 연령 데이터가 포함된 샘플 주석 데이터를 읽습니다.

######Here!!!! We can change the test data!
# 여기에서 테스트 데이터를 변경할 수 있습니다.
datSample = read.csv("AdditionalFile27SampleAnnotationExample55.csv")



## ------------------------------------- 새로 만든 코드

library(dplyr)
library(tidyverse)

#datClock %>% colnames()
datMethUsedNormalized %>% head()
datSample %>% head()

# 데이터 재생성
#pre23 <- datClock
meth26 <- datMethUsedNormalized %>% as.data.frame()

# 데이터 전처리
#pre23 <- pre23[pre23$CpGmarker!="(Intercept)",]

#join(CpGmarker로)하기 위해 데이터 맞춤
meth26 <- meth26 %>% rownames() %>% cbind(., meth26) %>% as.data.frame()
colnames(meth26)[1]<- 'id'

datAll <- datSample %>% right_join(meth26, by = 'id')
datSe <- datAll %>% select(Age, colnames(datAll)[grepl("^cg", colnames(datAll))])

# datSe <- datSe %>% t() %>% 

x_matrix <- as.matrix(datSe[, -which(names(datSe) == 'Age')])
y_vector <- datSe$Age


library(caret)
# 데이터를 훈련 데이터와 테스트 데이터로 나눔 (70% 훈련, 30% 테스트)
set.seed(123)  # 재현성을 위한 시드 설정
index <- createDataPartition(y_vector, p = 0.7, list = FALSE)
x_train <- x_matrix[index, ]
y_train <- y_vector[index]
x_test <- x_matrix[-index, ]
y_test <- y_vector[-index]


## Elastic net, Ridge, Lasso
library(glmnet)
# use 10 fold cross validation to estimate the lambda parameter 
# in the training data
alpha = 0.5 # α=0은 ridge,α=1은 lasso
#alpha=0.5는 Elastic Net 모델에서 L1 penalty와 L2 penalty를 반반씩 적용
glmnet.Training.CV = cv.glmnet(x = x_train, y = y_train, nfolds = 10, alpha = alpha, family = "gaussian")# The definition of the lambda parameter:
plot(glmnet.Training.CV)
lambda.glmnet.Training = glmnet.Training.CV$lambda.min # 오차가 제일 작은 람다
#x축은 로그 스케일로 람다(lambda) 값의 경로가 나타납니다. 람다 값은 규제의 강도를 조절하는 파라미터로, 높은 람다 값은 모델의 복잡성을 줄이고 더 많은 특성을 선택하지 않게 합니다.
#y축은 모델 성능을 측정하는 평가 지표입니다. 이 경우에는 평균 제곱 오차(Mean-Squared Error)가 사용되었습니다. 이 값은 작을수록 모델의 성능이 더 좋다는 의미입니다.

# alpha값에 따라 !경우의 수를 나눠보자!
alphas <- seq(0, 1, by = 0.1)  # Create a sequence of alpha values from 0 to 1
cv_results  <- list()
lambda.glmnet.Training_list <- list()

for (alpha in alphas) {
  glmnet.Training.CV_list <- cv.glmnet(x = x_train, y = y_train, nfolds = 10, alpha = alpha, family = "gaussian")
  cv_results[[as.character(alpha)]] <- glmnet.Training.CV_list
  lambda.glmnet.Training_list[[as.character(alpha)]] <- glmnet.Training.CV_list$lambda.min
}
lambda.glmnet.Training_list # alpha값에 따른 최적의 lambda

par(mfrow = c(2, 3))  # Adjust the layout according to your preference

for (i in 1:length(alphas)) {
  alpha <- as.character(alphas[i])
  plot(cv_results[[alpha]], main = paste("Alpha =", alpha))
}


# Fit the elastic net predictor to the training data
glmnet.Training = glmnet(x = x_train, y = y_train, family="gaussian", alpha=0.5, nlambda=100)
# Arrive at an estimate of of DNAmAge
DNAmAgeBasedOnTraining=inverse.F(predict(glmnet.Training,datout,type="response",s=lambda.glmnet.Training))




## Randomforest
library(randomForest)
library(ggplot2)
# RandomForest 모델 학습
#RandomForest 모델은 트리의 개수(ntree)와 트리의 최대 깊이(mtry) 등 여러 매개변수를 갖음
# 경우의 수에 따른 그래프 그리기
ntree_values <- seq(100, 1000, by = 100)  # 트리 개수의 경우의 수
mtry_values <- seq(1, ncol(x_train), by = 1)  # mtry의 경우의 수

results <- matrix(NA, nrow = length(ntree_values), ncol = length(mtry_values))

# ! 엄청 오래걸림 주의 !
for (i in 1:length(ntree_values)) {
  for (j in 1:length(mtry_values)) {
    ntree <- ntree_values[i]
    mtry <- mtry_values[j]
    rf_model <- randomForest(x = x_train, y = y_train, ntree = ntree, mtry = mtry)
    rf_predictions <- predict(rf_model, newdata = x_test)
    mse <- mean((rf_predictions - y_test)^2)
    results[i, j] <- mse
    print(paste(i, '/', length(ntree_values), j, '/', length(mtry_values)))
  }
}

# 결과를 그래프로 나타내기
results_df <- as.data.frame(results)
colnames(results_df) <- mtry_values
rownames(results_df) <- ntree_values

# 그래프 그리기
ggplot(results_df, aes(x = as.factor(colnames(results_df)), y = as.factor(rownames(results_df)), fill = as.vector(results_df))) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "RandomForest Mean Squared Error by ntree and mtry",
       x = "mtry", y = "ntree") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## Xgboost
library(xgboost)
#ntree와 max_depth 매개변수를 변경
# 경우의 수에 따른 그래프 그리기
ntree_values2 <- seq(100, 1000, by = 100)  # 트리 개수의 경우의 수
max_depth_values <- seq(1, 10, by = 1)  # max_depth의 경우의 수

results2 <- expand.grid(ntree_values = ntree_values, max_depth_values = max_depth_values)
results2$MSE <- NA

for (i in 1:nrow(results2)) {
  ntree2 <- results2$ntree_values[i]
  max_depth <- results2$max_depth_values[i]
  xgb_model <- xgboost(data = x_train, label = y_train, nrounds = ntree2, max_depth = max_depth, objective = "reg:squarederror")
  xgb_predictions <- predict(xgb_model, newdata = x_test)
  mse2 <- mean((xgb_predictions - y_test)^2)
  results2$MSE[i] <- mse2
  print(paste(i, '/', nrow(results2)))
}

# factor화
results2$ntree_values2 <- as.factor(results2$ntree_values2)
results2$max_depth_values <- as.factor(results2$max_depth_values)

# 결과를 그래프로 나타내기
ggplot(results2, aes(x = max_depth_values, y = ntree_values2, fill = MSE)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "XGBoost Mean Squared Error by ntree and max_depth",
       x = "max_depth", y = "ntree") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# x 축 (max_depth): 트리의 최대 깊이(max_depth)에 해당합니다. 이 값은 트리의 깊이가 얼마나 깊어지는지를 나타내며, 너무 낮으면 모델이 복잡성을 캡처하지 못하고, 너무 높으면 과적합의 위험이 있습니다.
# 
# y 축 (ntree): 트리의 개수(ntree)에 해당합니다. 이 값은 앙상블 내에서 생성되는 개별 의사결정 트리의 개수입니다.
# 
# 색상: 색상은 평균 제곱 오차(MSE)를 나타내며, 색상이 진할수록 낮은 MSE 값을 가짐을 의미합니다. 즉, 색상이 진할수록 모델의 예측 성능이 좋은 것을 나타냅니다.
# 그래프를 해석하는 과정에서 주목해야 할 점은:
#   
#   최적값 탐색: 그래프에서 MSE가 가장 낮은 영역을 찾아봅니다. 보통 색상이 진한 부분이 그 중요한 영역일 수 있습니다.
# 과적합 확인: max_depth가 너무 높아질수록 MSE가 증가하는 부분이 있을 수 있습니다. 이는 모델이 훈련 데이터에 과적합되는 경향을 보이는 것을 나타낼 수 있습니다.
# 트리의 개수와 성능: ntree 값에 따라서도 성능이 변하는 경향을 확인할 수 있습니다. 초기에 트리 개수가 적을 때는 성능이 낮을 수 있으나, 일정 수준 이상에서는 큰 성능 향상이 없거나 증가할 수 있습니다.
# 상호작용: max_depth와 ntree가 서로 어떻게 상호작용하는지 관찰합니다. 예를 들어, 더 큰 max_depth 값이 주어졌을 때 ntree 값이 어떤 영향을 미치는지 확인합니다.
#--------------------------

DNAmAge = datout$DNAmAge
medianAbsDev = function(x, y) median(abs(x - y), na.rm = TRUE)
medianAbsDev1 = signif(medianAbsDev(DNAmAge, datSample$Age), 2)
par(mfrow = c(1, 1))
verboseScatterplot(DNAmAge, datSample$Age, xlab = "DNAm Age", ylab = "Chronological Age", main = paste("All, err=", medianAbsDev1))
abline(0, 1) 