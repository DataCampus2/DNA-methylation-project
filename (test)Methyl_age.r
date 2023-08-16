
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


# library(caret)
# # 데이터를 훈련 데이터와 테스트 데이터로 나눔 (70% 훈련, 30% 테스트)
# set.seed(123)  # 재현성을 위한 시드 설정
# index <- createDataPartition(y_vector, p = 0.7, list = FALSE)
# x_train <- x_matrix[index, ]
# y_train <- y_vector[index]
# x_test <- x_matrix[-index, ]
# y_test <- y_vector[-index]

library(glmnet)
# use 10 fold cross validation to estimate the lambda parameter 
# in the training data
alpha = 0.5 # α=0은 ridge,α=1은 lasso
#alpha=0.5는 Elastic Net 모델에서 L1 penalty와 L2 penalty를 반반씩 적용
glmnet.Training.CV = cv.glmnet(x = x_matrix, y = y_vector, nfolds = 10, alpha = alpha, family = "gaussian")# The definition of the lambda parameter:
plot(glmnet.Training.CV)
lambda.glmnet.Training = glmnet.Training.CV$lambda.min # 오차가 제일 작은 람다
#x축은 로그 스케일로 람다(lambda) 값의 경로가 나타납니다. 람다 값은 규제의 강도를 조절하는 파라미터로, 높은 람다 값은 모델의 복잡성을 줄이고 더 많은 특성을 선택하지 않게 합니다.
#y축은 모델 성능을 측정하는 평가 지표입니다. 이 경우에는 평균 제곱 오차(Mean-Squared Error)가 사용되었습니다. 이 값은 작을수록 모델의 성능이 더 좋다는 의미입니다.

# alpha값에 따라 !경우의 수를 나눠보자!
alphas <- seq(0, 1, by = 0.1)  # Create a sequence of alpha values from 0 to 1
cv_results  <- list()
lambda.glmnet.Training_list <- list()

for (alpha in alphas) {
  glmnet.Training.CV_list <- cv.glmnet(x = x_matrix, y = y_vector, nfolds = 10, alpha = alpha, family = "gaussian")
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
glmnet.Training = glmnet(x = x_matrix, y = y_vector, family="gaussian", alpha=0.5, nlambda=100)
# Arrive at an estimate of of DNAmAge
DNAmAgeBasedOnTraining=inverse.F(predict(glmnet.Training,datout,type="response",s=lambda.glmnet.Training))


#--------------------------

DNAmAge = datout$DNAmAge
medianAbsDev = function(x, y) median(abs(x - y), na.rm = TRUE)
medianAbsDev1 = signif(medianAbsDev(DNAmAge, datSample$Age), 2)
par(mfrow = c(1, 1))
verboseScatterplot(DNAmAge, datSample$Age, xlab = "DNAm Age", ylab = "Chronological Age", main = paste("All, err=", medianAbsDev1))
abline(0, 1) 