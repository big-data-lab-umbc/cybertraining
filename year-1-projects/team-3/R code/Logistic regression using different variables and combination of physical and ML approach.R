library(randomForest)
library(pROC)
library(caret)
library(ROCR)
## read the MODIS and CALIPSO data file;
raw715<-read.table(file="~/umbc/is698/CalMod_07152007_include_geometry.csv",sep=",",header=TRUE)
raw622<-read.table(file="~/umbc/is698/CalMod_22062009_angles.csv",sep=",",header=TRUE)
raw1105<-read.table(file="~/umbc/is698/CalMod_11052009.csv",sep=",",header=TRUE)

ds_prep <- function(fname) {
  p1<-as.data.frame(fname)
  p1$id<-NULL
  p1$X<-NULL
  p1$dust<-as.integer(p1$dust)
  typeof(p1$dust)
  # regroup dust variable, dust or dust with cloud is consider 1, otherwise is 0;
  p1$dust[which(p1$dust ==2)]<-0
  p1$dust[which(p1$dust ==3)]<-1
  
  # create physical variables 
  p1$v1 = p1$Band2/p1$Band31 # Band 860 nm vs 11 um
  p1$v2 = p1$Band3/p1$Band2 # Band 460 nm / 860 nm
  p1$v3 = p1$Band2/p1$Band32 # Band 860 nm vs 11 um
  p1$v4 = p1$Band32-p1$Band33 # Band 11um - 12 um
  return(p1)
}

ds715<-ds_prep(raw715)
ds622<-ds_prep(raw622)
ds1105<-ds_prep(raw1105)

ds7152<-rbind(ds715,ds1105)
set.seed(1)
# create training and testing datset for the day 715 ;
index<-sample(1:nrow(ds715), size=0.3*nrow(ds715))
ds715.test<-ds715[index,]
ds715.train<-ds715[-index,]
# create model formula using all Band variables;
PredictorVariables <- paste("Band", 1:38, sep="")
Formula <- formula(paste("dust ~ ", paste(PredictorVariables, collapse=" + ")))


nothing <- glm(dust ~ 1,data=ds715.train,family=binomial)
fullmod   <-glm(Formula, data=ds715.train, family="binomial")
## stepwise selection method;
bothways <-step(nothing, list(lower=formula(nothing),upper=formula(fullmod)),data=ds715,direction="both",trace=0)

# after the above step, we can get 16 Band variables for formula1;
formula1<-formula(bothways)

# formula with physical varibles by adding physical variables;
formula.physical<-as.formula((paste("dust~",paste(as.character(formula1[3]),"v1+v2+v3+v4",sep="+"))))

# formula with angle varibles;
formula.angle<-as.formula((paste("dust~",paste(as.character(formula1[3]),"Band39+Band40+Band41+Band42",sep="+"))))



## the folloing function is to get predicted value and calcuate AUC and accuracy;
pred <- function(fname,ds.test,form, dstrain) {
  
mylogit<-glm(as.formula(form), data=as.data.frame(dstrain), family="binomial")
p1.test<-data.frame(ds.test)
prob<-predict(mylogit,p1.test,type="response")
p1.test$pred2<-(prob>0.5)*1
auc<-auc(p1.test$pred2,p1.test$dust)
test<-confusionMatrix(as.factor(p1.test$pred2),as.factor(p1.test$dust))
accuracy<-test$overall['Accuracy']
sesp<-t(test$byClass)
se<-sesp[1]
sp<-sesp[2] 
youden<-se+sp-1
df<-rbind(df,data.frame(fname,auc,accuracy,youden))
return(df)
}

## combine the output from different models;
df<-data.frame(df_name=character(),auc=integer(),accuracy=double(),youden=double())
df<-pred("logit with all data 715",ds715.test,"dust~.",ds715.train)
df<-pred("logit with ml 715",ds715.test,formula1,ds715.train)
df<-pred("logit with physical 715",ds715.test,formula.physical,ds715.train)
df<-pred("logit with angle 715",ds715.test,formula.angle,ds715.train)

## prediction using a different day data set;
df<-pred("logit with all data 622",ds622,"dust~.",ds715.train)
df<-pred("logit with ml 622",ds622,formula1,ds715.train)
df<-pred("logit with physical 622",ds622,formula.physical,ds715.train)
df

## Physical algorithm BTD using data on 715, threshold is set 0.8

ds715.test$physical<-((ds715.test$Band32-ds715.test$Band33)>0.8)*1
accuracy <- table(ds715.test$physical, ds715.test$dust)
sum(diag(accuracy))/sum(accuracy)


## Physical algorithm BTD using data on 622, threshold is set 0.8

ds622$physical<-((ds622$Band32-ds622$Band33)>0.8)*1
accuracy <- table(ds622$physical, ds622$dust)
sum(diag(accuracy))/sum(accuracy)
df


df<-pred("logit with angle 715",ds622,"dust ~ Band16+Band35+Band29+Band26+Band9+Band23+v1+v2+v3")

df

write.csv(df,file="~//UMBC//is698//performance comparsion combined.csv")

## the following is to predict dust for MODIS regrion on 622 day

mylogit<-glm(formula.physical, data=ds715.train, family="binomial")

raw622_gran<-read.table(file="~/umbc/is698/ModGranule_06222009.csv",sep=",",header=TRUE)
raw622_gran[is.na(raw622_gran)] <- 0
raw622_gran$v1 = raw622_gran$Band2/raw622_gran$Band31 # R0.47um
raw622_gran$v2 = raw622_gran$Band3/raw622_gran$Band2 # band 460 num / 860 num
raw622_gran$v3 = raw622_gran$Band2/raw622_gran$Band32 
raw622_gran$v3 = raw622_gran$Band32 - raw622_gran$Band33 

prob<-predict(mylogit,raw622_gran,type="response")


mylogit<-glm("dust ~Band16+Band35+Band29+Band26+Band9+Band23+v1+v2+v3+v4", ds715.train, family="binomial")
prob<-predict(mylogit,raw622_gran,type="response")
write.table(as.data.frame(prob),file="~//is698//622predict2.txt",sep="\t",col.names = F,row.names=FALSE)


