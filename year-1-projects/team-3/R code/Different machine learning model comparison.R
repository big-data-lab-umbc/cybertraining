
library(caret)
library(randomForest)
library(pROC)
library(caret)
library(ROCR)
library(e1071)
library(nnet)
## read the MODIS and CALIPSO data file;
raw715<-read.table(file="~/umbc/is698/CalMod_07152007_include_geometry.csv",sep=",",header=TRUE)
raw622<-read.table(file="~/umbc/is698/CalMod_22062009_angles.csv",sep=",",header=TRUE)
raw715<-na.omit(raw715)
raw622<-na.omit(raw622)
# create model formula using all band variables;
PredictorVariables <- paste("Band", 1:38, sep="")
Formula <- formula(paste("dust ~ ", paste(PredictorVariables, collapse=" + ")))

ds_prep <- function(fname) {
  p1<-as.data.frame(fname)
  p1$id<-NULL
  p1$X<-NULL
  p1$dust<-as.integer(p1$dust)
  typeof(p1$dust)
  p1$dust[which(p1$dust ==2)]<-0
  p1$dust[which(p1$dust ==3)]<-1
  p1$v1 = p1$Band2/p1$Band31 
  p1$v2 = p1$Band3/p1$Band2 
  p1$v3 = p1$Band2/p1$Band32 
  return(p1)
}

ds715<-ds_prep(raw715)
ds622<-ds_prep(raw622)
ds622<-na.omit(ds622)
ds715<-na.omit(ds715)
set.seed(1)
index<-sample(1:nrow(ds715), size=0.3*nrow(ds715))
#ds715.test<-ds715[index,]
ds715.test<-ds622
ds715.train<-ds715

## the folloing function is to get predicted value and calcuate AUC and accuracy;
pred <- function(fname,mname) {
  if(grepl("network",fname)==FALSE){
  prob<-predict(mname,ds715.test,type="response")}
  else {
    #prob<-compute(fit.nn, ds715.test.nn)$net.result[,1]
    prob<-predict(fit.nn,ds715.test.nn,type="raw")[,1]
  }
  
  ds715.test$pred2<-(prob>0.5)*1
  auc<-auc(ds715.test$dust,prob)
  test<-confusionMatrix(as.factor(ds715.test$pred2),as.factor(ds715.test$dust))
  accuracy<-test$overall['Accuracy']
  sesp<-t(test$byClass)
  se<-sesp[1]
  sp<-sesp[2] 
  youden<-se+sp-1
  df<-rbind(df,data.frame(fname,auc,accuracy,youden))
  return(df)
}

#logistic regression model
fit.lr<-glm(Formula, data=ds715.train, family="binomial")

# random forest model
fit.rf<-randomForest(Formula, data=ds715.train,importance=TRUE)

# scale the data set for ANN model;
ds715.nn<-as.data.frame(scale(ds715.train[c(PredictorVariables)]))
head(ds715.nn)
ds715.nn$dust<-ds715.train$dust
# scale the testing data for ANN model;
ds715.test.nn<-as.data.frame(scale(ds715.test[c(PredictorVariables)]))
#support vector machine method;
fit.svm<-svm(Formula,data=ds715.train)
#ANN model;
fit.nn<- nnet(Formula, data=ds715.nn, size=10, decay=1.0e-5, maxit=50)

df<-data.frame(df_name=character(),auc=integer(),accuracy=double(),youden=double())
df<-pred("Logistic regression",fit.lr)
df<-pred("Random forests",fit.rf)
df<-pred("Neural network",fit.nn)
df<-pred("Support vector machine",fit.svm)
## the follwing is to get predicted dust probability from each model;
ds715.test$prob1<-predict(fit.lr,ds715.test,type="response")
ds715.test$prob2<-predict(fit.rf,ds715.test,type="response")
ds715.test$prob3<-predict(fit.svm,ds715.test,type="response")
ds715.test$prob4<-predict(fit.nn, ds715.test.nn,type="raw")[,1]
## stacking classifier by averaging the predicted values from the above 4 models;
ds715.test$prob<-rowMeans(ds715.test[,c("prob1", "prob2","prob3","prob4")], na.rm=TRUE)
ds715.test$pred2<-(ds715.test$prob>0.5)*1
auc<-auc(ds715.test$dust,ds715.test$prob)
test<-confusionMatrix(as.factor(ds715.test$pred2),as.factor(ds715.test$dust))
accuracy<-test$overall['Accuracy']
sesp<-t(test$byClass)
se<-sesp[1]
sp<-sesp[2] 
youden<-se+sp-1
fname<-"Ensemble"
df<-rbind(df,data.frame(fname,auc,accuracy,youden))
df
write.csv(df,file="~//is698//is698 table 1.csv")

