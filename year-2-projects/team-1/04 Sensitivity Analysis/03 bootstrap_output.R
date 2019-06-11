## This file evaluates the mean and SE of the precipitation spatially and temporally based on inflated variance
library(ggplot2)
library(tidyr)
library(gridExtra)

#Set WD to where yor output folder is
#setwd("D:/Google Drive/UMBC/2019 02 Spring/Cybertraining/Project/data/output")

bootcomp = function(size) #Size is the number of bootstrap iterations
{
results <- vector("list", size) 
a = list.dirs(getwd())
a=paste0(a,'/')
for(j in 1:size)
{
  file_names = list.files(a[j+1])
  file_names2 = paste0(a[j+1],file_names)
  files=lapply( file_names2 , read.table ,header = F,sep = '\t',numerals = 'no.loss',stringsAsFactors = F) # gets files
  
  lat = substr(file_names,8,14) # extract lat and long from the filenames vector
  lat = as.numeric(lat)
  long = substr(file_names,16,23)
  long = as.numeric(long)
  
  # empty dataframes
  precip.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
  evapo.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
  runoff.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
  balance.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
  baseflow.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
  rb.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
  
  for(i in 1:387)
  {
    
    files[[i]] = data.frame(files[[i]][457:639,c(2,4:7)]) # Keep only Apr-Sep2017 in row, drop date and year variable in col
    names(files[[i]])=c('month','precip','evapo','runoff','baseflow') # name the variable
    
    precip.df$lat[i] = lat[i]
    precip.df$long[i]=long[i]
    precip.df$apr[i] = sum(files[[i]]$precip[files[[i]]$month==4])
    precip.df$may[i] = sum(files[[i]]$precip[files[[i]]$month==5])
    precip.df$jun[i] = sum(files[[i]]$precip[files[[i]]$month==6])
    precip.df$jul[i] = sum(files[[i]]$precip[files[[i]]$month==7])
    precip.df$aug[i] = sum(files[[i]]$precip[files[[i]]$month==8])
    precip.df$sep[i] = sum(files[[i]]$precip[files[[i]]$month==9])
    precip.df$tot[i] = sum(precip.df[i,3:8])
    
    evapo.df$lat[i] = lat[i]
    evapo.df$long[i]=long[i]
    evapo.df$apr[i] = sum(files[[i]]$evapo[files[[i]]$month==4])
    evapo.df$may[i] = sum(files[[i]]$evapo[files[[i]]$month==5])
    evapo.df$jun[i] = sum(files[[i]]$evapo[files[[i]]$month==6])
    evapo.df$jul[i] = sum(files[[i]]$evapo[files[[i]]$month==7])
    evapo.df$aug[i] = sum(files[[i]]$evapo[files[[i]]$month==8])
    evapo.df$sep[i] = sum(files[[i]]$evapo[files[[i]]$month==9])
    evapo.df$tot[i] = sum(evapo.df[i,3:8])
    evapo.df = replace(evapo.df,is.na(evapo.df),0) #there were some NAs, I replaced them by 0
    
    runoff.df$lat[i] = lat[i]
    runoff.df$long[i]=long[i]
    runoff.df$apr[i] = sum(files[[i]]$runoff[files[[i]]$month==4])
    runoff.df$may[i] = sum(files[[i]]$runoff[files[[i]]$month==5])
    runoff.df$jun[i] = sum(files[[i]]$runoff[files[[i]]$month==6])
    runoff.df$jul[i] = sum(files[[i]]$runoff[files[[i]]$month==7])
    runoff.df$aug[i] = sum(files[[i]]$runoff[files[[i]]$month==8])
    runoff.df$sep[i] = sum(files[[i]]$runoff[files[[i]]$month==9])
    runoff.df$tot[i] = sum(runoff.df[i,3:8])
    
    baseflow.df$lat[i] = lat[i]
    baseflow.df$long[i]=long[i]
    baseflow.df$apr[i] = sum(files[[i]]$baseflow[files[[i]]$month==4])
    baseflow.df$may[i] = sum(files[[i]]$baseflow[files[[i]]$month==5])
    baseflow.df$jun[i] = sum(files[[i]]$baseflow[files[[i]]$month==6])
    baseflow.df$jul[i] = sum(files[[i]]$baseflow[files[[i]]$month==7])
    baseflow.df$aug[i] = sum(files[[i]]$baseflow[files[[i]]$month==8])
    baseflow.df$sep[i] = sum(files[[i]]$baseflow[files[[i]]$month==9])
    baseflow.df$tot[i] = sum(baseflow.df[i,3:8])
  
    
  }
  rb.df[,3:9] = runoff.df[,3:9] + baseflow.df[,3:9]
  rb.df$lat = lat
  rb.df$long = long
  balance.df[,3:9] = precip.df[,3:9] - evapo.df[,3:9] - rb.df[,3:9]
  balance.df$lat = lat
  balance.df$long = long
  
results[[j]] = balance.df #dataframe with each month for the j-th bootstrap
}
return(results)
}

setwd("D:/VIC/output_sn2/cvx15/")
results15 = bootcomp(66) #For some reason 67 onwards was not being read properly for cv x 1.5
setwd("D:/VIC/output_sn2/cvx1/")
results1 = bootcomp(100)
setwd("D:/VIC/output_sn2/cvx2/")
results2 = bootcomp(100)


totalwb = function(results)
{
  size = length(results)
r = data.frame(lat = numeric(387),long=numeric(387),wb = numeric(387))
for(i in 1:size)
  r$wb = r$wb + results[[i]]$tot**2 #standard errors of the simulations

r$wb = sqrt(r$wb/size)
r$lat = results[[1]]$lat
r$long = results[[1]]$long
return(r)
}
r15 = totalwb(results15)
r2 = totalwb(results2)
r1 = totalwb(results1)


b1 = ggplot(r1, aes(x=long,y=lat,fill=wb))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'SE for original CV') + theme(legend.position = 'bottom')

b2 = ggplot(r2, aes(x=long,y=lat,fill=wb))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'SE for 2xCV') + theme(legend.position = 'bottom')

b15 = ggplot(r15, aes(x=long,y=lat,fill=wb))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'SE for 1.5xCV') + theme(legend.position = 'bottom')

b1
b15
b2

cv = data.frame(lat=lat,long=long,cv1 = r1$wb, cv15=r15$wb, cv2 = r2$wb)
head(cv)

write.csv(cv,file='perturb.csv',row.names = F)
varb1 = as.data.frame(matrix(NA,387,100))
varb2 = as.data.frame(matrix(NA,387,100))
varb15 = as.data.frame(matrix(NA,387,66))

for(i in 1:100)
{
  varb1[i] = results1[[i]]$tot
  varb2[i] = results2[[i]]$tot
}

for(i in 1:66)
{
  varb15[i] = results15[[i]]$tot
}
View(varb1)
iqr1 = apply(varb1, 1, iqr)
iqr15 = apply(varb15, 1, iqr)
iqr2 = apply(varb2, 1, iqr)
iqr = data.frame(lat=lat,long=long,iqr1,iqr15,iqr2)

c1 = ggplot(iqr, aes(x=long,y=lat,fill=iqr1))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'IQR for original CV') + theme(legend.position = 'bottom')

c2 = ggplot(iqr, aes(x=long,y=lat,fill=iqr2))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'IQR for 2xCV') + theme(legend.position = 'bottom')

c15 = ggplot(iqr, aes(x=long,y=lat,fill=iqr15))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'IQR for 1.5xCV') + theme(legend.position = 'bottom')
write.csv(iqr,file='perturb2.csv',row.names = F)

grid.arrange(b1,b15,b2,c1,c15,c2,nrow=2)

#####################################################
### This gives the plots for mean water balance each month for the 3 different variability levels
wbdist = function(results)
{
  size = length(results)
  r = data.frame(month = month.abb[4:9],wb=NA,p5=NA,p95=NA)
  month = data.frame(apr=NA,may=NA,jun=NA,jul=NA,aug=NA,sep=NA)
  for(i in 1:size)
  {
    temp = apply(results[[i]], 2, mean,na.rm=T)
    temp = temp[3:8]
    month[i,] = temp
  }
  month=month[-1,]
  for(j in 1:6)
  {
    r[j,2] = mean(month[,j])
    r[j,3] = quantile(month[,j],.025)
    r[j,4] = quantile(month[,j],.975)
  }
  return(r)
}

mcv1 = wbdist(results1)
mcv1$cv = 'original CV'
mcv15 = wbdist(results15)
mcv15$cv = 'CV x 1.5'
mcv2 = wbdist(results2)
mcv2$cv = 'CV x 2'

cvdist = rbind(mcv1,mcv15,mcv2)
cvdist$month = factor(cvdist$month,levels(cvdist$month)[c(1,5,4,3,2,6)])
cvdist$cv = as.factor(cvdist$cv)
cvdist$cv = factor(cvdist$cv,levels(cvdist$cv)[c(3,1,2)])
ggplot(cvdist,aes(x=month,y=wb,fill=cv)) + geom_bar(stat = 'identity',position = position_dodge()) +
  geom_errorbar(aes(ymin=p5,ymax=p95),width=.2,position = position_dodge(.9)) + labs(y = 'Estimated mean water balance for the Potomac basin')
write.csv(cvdist,'cvdist.csv')
