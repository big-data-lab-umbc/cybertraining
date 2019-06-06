library(ggplot2)
library(tidyr)
library(gridExtra)

## Setwd to the folder where your flux output folder is
#setwd("D:/Google Drive/UMBC/2019 02 Spring/Cybertraining/Project/data/output")


file_names=list.files(pattern="fluxes_*") #gets file names from working directory; can specify path

files=lapply( file_names , read.table ,header = F,sep = '\t',numerals = 'no.loss',stringsAsFactors = F) # gets files

#View(files[[1]])

lat = substr(file_names,8,14) # extract lat and long from the filenames vector
lat = as.numeric(lat)
long = substr(file_names,16,23)
long = as.numeric(long)

# empty dataframes
precip.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
evapo.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
runoff.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
balance.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
#soil.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
baseflow.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))
rb.df = data.frame(lat=numeric(387),long=numeric(387),apr = numeric(387),may = numeric(387),jun = numeric(387),jul = numeric(387),aug = numeric(387),sep = numeric(387),tot = numeric(387))

for(i in 1:387)
{

  files[[i]] = data.frame(files[[i]][457:639,c(2,4:11)]) # Keep only Apr-Sep2017 in row, drop date and year variable in col
  names(files[[i]])=c('month','precip','evapo','runoff','baseflow','swe','soil1','soil2','soil3') # name the variable
  
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
  
  #soil.df$lat[i] = lat[i]
  #soil.df$long[i]=long[i]
  #soil.df$apr[i] = sum(files[[i]]$soil1[files[[i]]$month==4] + files[[i]]$soil2[files[[i]]$month==4] + files[[i]]$soil3[files[[i]]$month==4])
  #soil.df$may[i] = sum(files[[i]]$soil1[files[[i]]$month==5] + files[[i]]$soil2[files[[i]]$month==5] + files[[i]]$soil3[files[[i]]$month==5])
  #soil.df$jun[i] = sum(files[[i]]$soil1[files[[i]]$month==6] + files[[i]]$soil2[files[[i]]$month==6] + files[[i]]$soil3[files[[i]]$month==6])
  #soil.df$jul[i] = sum(files[[i]]$soil1[files[[i]]$month==7] + files[[i]]$soil2[files[[i]]$month==7] + files[[i]]$soil3[files[[i]]$month==7])
  #soil.df$aug[i] = sum(files[[i]]$soil1[files[[i]]$month==8] + files[[i]]$soil2[files[[i]]$month==8] + files[[i]]$soil3[files[[i]]$month==8])
  #soil.df$sep[i] = sum(files[[i]]$soil1[files[[i]]$month==9] + files[[i]]$soil2[files[[i]]$month==9] + files[[i]]$soil3[files[[i]]$month==9])
  #soil.df$tot[i] = sum(soil.df[i,3:8])
  
}
rb.df[,3:9] = runoff.df[,3:9] + baseflow.df[,3:9] #variable for runoff + baseflow
rb.df$lat = lat
rb.df$long = long
balance.df[,3:9] = precip.df[,3:9] - evapo.df[,3:9] - rb.df[,3:9] #Water balance, monthly
balance.df$lat = lat
balance.df$long = long

write.csv(precip.df,file = 'precip.csv',row.names = F)
write.csv(evapo.df,file = 'evapo.csv',row.names = F)
write.csv(rb.df,file = 'rb.csv',row.names = F)
write.csv(balance.df,file = 'balance.csv',row.names = F)
#write.csv(soil.df,file = 'soil.csv',row.names = F)
#write.csv(baseflow.df,file = 'baseflow.csv',row.names = F)
#summary(soil.df)

########## TOtal for the whole basic across all months
total.df = data.frame(lat,long,precip = precip.df$tot,evapo=evapo.df$tot,balance=balance.df$tot,rb=rb.df$tot)
precip.df = precip.df[,-9]
runoff.df = runoff.df[,-9]
evapo.df = evapo.df[,-9]
balance.df = balance.df[,-9]
soil.df = soil.df[,-9]
baseflow.df = baseflow.df[,-9]


######plots
p1 = ggplot(precip.df, aes(x=long,y=lat,fill=apr))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'April precipitation') +theme(legend.position = 'bottom')
p2 = ggplot(precip.df, aes(x=long,y=lat,fill=may))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'May precipitation') + theme(legend.position = 'bottom')
p3 = ggplot(precip.df, aes(x=long,y=lat,fill=jun))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'June precipitation') + theme(legend.position = 'bottom')
p4 = ggplot(precip.df, aes(x=long,y=lat,fill=jul))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'July precipitation') + theme(legend.position = 'bottom')
p5 = ggplot(precip.df, aes(x=long,y=lat,fill=aug))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Aug precipitation') + theme(legend.position = 'bottom')
p6 = ggplot(precip.df, aes(x=long,y=lat,fill=sep))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'September Precipitation') + theme(legend.position = 'bottom')

grid.arrange(p1,p2,p3,p4,p5,p6,nrow=2)

e1 = ggplot(evapo.df, aes(x=long,y=lat,fill=apr))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'April evapotranspiration') +theme(legend.position = 'bottom')
e2 = ggplot(evapo.df, aes(x=long,y=lat,fill=may))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'May evapotranspiration') + theme(legend.position = 'bottom')
e3 = ggplot(evapo.df, aes(x=long,y=lat,fill=jun))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'June evapotranspiration') + theme(legend.position = 'bottom')
e4 = ggplot(evapo.df, aes(x=long,y=lat,fill=jul))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'July evapotranspiration') + theme(legend.position = 'bottom')
e5 = ggplot(evapo.df, aes(x=long,y=lat,fill=aug))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Aug evapotranspiration') + theme(legend.position = 'bottom')
e6 = ggplot(evapo.df, aes(x=long,y=lat,fill=sep))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'September evapotranspiration') + theme(legend.position = 'bottom')

grid.arrange(e1,e2,e3,e4,e5,e6,nrow=2)

r1 = ggplot(runoff.df, aes(x=long,y=lat,fill=apr))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'April runoff') +theme(legend.position = 'bottom')
r2 = ggplot(runoff.df, aes(x=long,y=lat,fill=may))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'May runoff') + theme(legend.position = 'bottom')
r3 = ggplot(runoff.df, aes(x=long,y=lat,fill=jun))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'June runoff') + theme(legend.position = 'bottom')
r4 = ggplot(runoff.df, aes(x=long,y=lat,fill=jul))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'July runoff') + theme(legend.position = 'bottom')
r5 = ggplot(runoff.df, aes(x=long,y=lat,fill=aug))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Aug runoff') + theme(legend.position = 'bottom')
r6 = ggplot(runoff.df, aes(x=long,y=lat,fill=sep))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'September runoff') + theme(legend.position = 'bottom')

grid.arrange(r1,r2,r3,r4,r5,r6,nrow=2)

sw1 = ggplot(baseflow.df, aes(x=long,y=lat,fill=apr))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'April baseflow') +theme(legend.position = 'bottom')
sw2 = ggplot(baseflow.df, aes(x=long,y=lat,fill=may))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'May baseflow') + theme(legend.position = 'bottom')
sw3 = ggplot(baseflow.df, aes(x=long,y=lat,fill=jun))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'June baseflow') + theme(legend.position = 'bottom')
sw4 = ggplot(baseflow.df, aes(x=long,y=lat,fill=jul))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'July baseflow') + theme(legend.position = 'bottom')
sw5 = ggplot(baseflow.df, aes(x=long,y=lat,fill=aug))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Aug baseflow') + theme(legend.position = 'bottom')
sw6 = ggplot(baseflow.df, aes(x=long,y=lat,fill=sep))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'September baseflow') + theme(legend.position = 'bottom')

grid.arrange(sw1,sw2,sw3,sw4,sw5,sw6,nrow=2)

#s1 = ggplot(soil.df, aes(x=long,y=lat,fill=apr))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  #labs(fill = 'mm', x='Longitude', y='Latitude',title = 'April soil moisture') +theme(legend.position = 'bottom')
#s2 = ggplot(soil.df, aes(x=long,y=lat,fill=may))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  #labs(fill = 'mm', x='Longitude', y='Latitude',title = 'May soil moisture') + theme(legend.position = 'bottom')
#s3 = ggplot(soil.df, aes(x=long,y=lat,fill=jun))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  #labs(fill = 'mm', x='Longitude', y='Latitude',title = 'June soil moisture') + theme(legend.position = 'bottom')
#s4 = ggplot(soil.df, aes(x=long,y=lat,fill=jul))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  #labs(fill = 'mm', x='Longitude', y='Latitude',title = 'July soil moisture') + theme(legend.position = 'bottom')
#s5 = ggplot(soil.df, aes(x=long,y=lat,fill=aug))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  #labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Aug soil moisture') + theme(legend.position = 'bottom')
#s6 = ggplot(soil.df, aes(x=long,y=lat,fill=sep))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  #labs(fill = 'mm', x='Longitude', y='Latitude',title = 'September soil moisture') + theme(legend.position = 'bottom')

#grid.arrange(s1,s2,s3,s4,s5,s6,nrow=2)

b1 = ggplot(balance.df, aes(x=long,y=lat,fill=apr))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'April water balance') +theme(legend.position = 'bottom')
b2 = ggplot(balance.df, aes(x=long,y=lat,fill=may))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'May water balance') + theme(legend.position = 'bottom')
b3 = ggplot(balance.df, aes(x=long,y=lat,fill=jun))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'June water balance') + theme(legend.position = 'bottom')
b4 = ggplot(balance.df, aes(x=long,y=lat,fill=jul))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'July water balance') + theme(legend.position = 'bottom')
b5 = ggplot(balance.df, aes(x=long,y=lat,fill=aug))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Aug water balance') + theme(legend.position = 'bottom')
b6 = ggplot(balance.df, aes(x=long,y=lat,fill=sep))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'September water balance') + theme(legend.position = 'bottom')

grid.arrange(b1,b2,b3,b4,b5,b6,nrow=2)

t1 = ggplot(total.df, aes(x=long,y=lat,fill=precip))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Total Precipitation') +theme(legend.position = 'bottom')
t2 = ggplot(total.df, aes(x=long,y=lat,fill=evapo))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Total evapotranspiration') + theme(legend.position = 'bottom')
t3 = ggplot(total.df, aes(x=long,y=lat,fill=rb))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Total runoff and baseflow') + theme(legend.position = 'bottom')
t4 = ggplot(total.df, aes(x=long,y=lat,fill=balance))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'mm', x='Longitude', y='Latitude',title = 'Total water balance') + theme(legend.position = 'bottom')


grid.arrange(t1,t2,t3,t4,nrow=2)

#april water budget
grid.arrange(p1,e1,r1,b1,nrow=2)
#may water budget
grid.arrange(p2,e2,r2,b2,nrow=2)
#june water budget
grid.arrange(p3,e3,r3,b3,nrow=2)
#july water budget
grid.arrange(p4,e4,r4,b4,nrow=2)
#aug water budget
grid.arrange(p5,e5,r5,b5,nrow=2)
#sept water budget
grid.arrange(p6,e6,r6,b6,nrow=2)

### Multi-bar bar chart for monthly water budget components
combined = data.frame(month=numeric(6),precip=numeric(6),evapo=numeric(6),rb=numeric(6))
combined$month = paste(as.character(4:9),month.name[4:9],sep = '-') #Had to do 04-Apr and so on to keep month order. Look in sensitivity analysis code to see how it can be done using relevel
#Means across all locations, one value for each month
combined$precip = apply(precip.df[,3:8],2,mean)
combined$evapo = apply(evapo.df[,3:8],2,mean)
combined$rb = apply(rb.df[,3:8],2,mean,na.rm=T)
names(combined) = c('Month','Precipitation','Evapotranspiration','Runoff and Baseflow')
#Turn it into a long dataset so that we can plot it nicely
combined_long = gather(combined,component, measurement,2:4,factor_key = T)

ggplot(combined_long,aes(fill=component,y=measurement,x=Month)) + geom_bar(position = 'dodge',stat = 'identity') + 
  labs(x = 'Month',y='Value (mm)',title = 'Water Budget components between April and September 2017')