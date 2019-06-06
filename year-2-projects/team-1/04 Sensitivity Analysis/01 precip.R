## This file does some basic spatio temporal analysis of the precipitation data

# install.packages("lubridate")
# install.packages("TSA")
# install.packages("forecast")
# install.packages("RColorBrewer")
# install.packages("ggplot2")
# install.packages("gstat")
# install.packages("sp")
# install.packages("spatial")
library(lubridate)
library(TSA)
library(forecast)
library(RColorBrewer)
library(ggplot2)
library(gstat)
library(sp)
library(spatial)
library(ncf)

###################################
##################################

#setwd("D:/Google Drive/UMBC/2019 02 Spring/Cybertraining/Project/meteo_forcing")
#setwd("C:/Users/reetam/Google Drive/2019 02 Spring/Cybertraining/Project/meteo_forcing")


file_names=list.files(pattern="forcing_*") #gets file names from working directory; can specify path

files=lapply( file_names , read.table ,header = F,sep = ' ',numerals = 'no.loss',stringsAsFactors = F) # gets files

for(i in 1:387)
{
  #files[[i]] = files[[i]][1:456,]
  files[[i]] = data.frame(files[[i]][,1]) # Just keep precip
  names(files[[i]])='precip' # name the variable
}
#names(files[[1]])
#dim(files[[1]])
lat = substr(file_names,9,15) # extract lat and long from the filenames vector
long = substr(file_names,17,24)

for(i in 1:387) # add lat and long to each element in list. might not be necessary
{
  files[[i]]$lat = rep(as.numeric(as.character(lat[i])),639) #makes lat and long numeric.
  files[[i]]$long = rep(as.numeric(as.character(long[i])),639)
}
#head(files[[2]])


#### Temporal Analysis #######
##############################

#create temporal data frame; each row is a grid point

mean=numeric(387);mean0=numeric(387);zero=numeric(387)
sd = numeric(387); sd0=numeric(387)
for(i in 1:387)
{
  mean0[i]=mean(files[[i]]$precip)
  mean[i]=mean(files[[i]]$precip[files[[i]]$precip>0])
  sd0[i]=sd(files[[i]]$precip)
  sd[i]=sd(files[[i]]$precip[files[[i]]$precip>0])
  zero[i]=sum(files[[i]]$precip==0)/639
}

params = data.frame(lat,long,mean,sd,mean0,sd0,zero)
params$lat = as.numeric(as.character(params$lat))
params$long = as.numeric(as.character(params$long))
head(params)
#write.csv(params,"heatmap.csv")
ggplot(params, aes(x=long,y=lat,fill=mean0))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'Mean', x='Longitude', y='Latitude',title = 'Heat map of annual mean precipitation')

ggplot(params, aes(x=long,y=lat,fill=sd0))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
  labs(fill = 'Mean', x='Longitude', y='Latitude',title = 'Heat map of annual SD of precipitation')

ggplot(params, aes(x=long,y=lat,fill=zero))+geom_raster() + scale_fill_gradientn(colours=rainbow(7)) +
    labs(fill = 'Mean', x='Longitude', y='Latitude',title = 'Heat map of proportion of days with no rainfall')


for(i in 1:387)
{
  params$ts.model[i]=as.character(auto.arima(ts(files[[i]]$precip))) #extract auto.arima model for each location
}
xtabs(~ts.model,data=params)

params$ts.model = as.factor(substr(params$ts.model,6,12)) # keep only the values of p,d,q
ggplot(params,aes(x=long,y=lat,label=ts.model,col=ts.model))+ geom_raster() +
       geom_text() + scale_color_brewer(palette = "Spectral") + scale_fill_gradientn(colors=terrain.colors(10)) +
        labs(col='(p,d,q)',x='Longitude',y='Latitude',title='ARIMA(p,d,q) models for annual precipitation' )


##### SPatial Analysis ##########
#################################

## Correlogram
dates  =seq(as.Date("2016-01-01"),as.Date("2017-03-31"),by='days')

sp.temp = data.frame(matrix(NA,nrow = 387,ncol = 458))
names(sp.temp) = c("lat","long",as.character(dates))
sp.temp$lat = params$lat
sp.temp$long = params$long
for(i in 1:387)
{
  sp.temp[i,3:641]=files[[i]][,1]
}
rain=numeric(639)
for(i in 1:639)
  rain[i] = mean(sp.temp[,i+2]>0)
summary(rain)
sum(rain==0)
sum(rain<.3)
plot.ecdf(rain)


# using spatial - individual correlograms
cor1 = data.frame(x=params$long,y=params$lat,z=sp.temp[,31])
head(cor1)
topo.cor1 = surf.ls(2,cor1)
cor1 = correlogram(topo.cor1,25)
  head(cor1)

#Get the first lag of correlogram at each data point
max.cor = numeric(456)
for(i in 1:456)
{
  temp = data.frame(x=params$long,y=params$lat,z=sp.temp[,i+2])
  topo.temp = surf.ls(2,temp)
  correl = correlogram(topo.temp,25)
  max.cor[i] = correl$y[2]
}

ggplot(as.data.frame(max.cor),aes(x=max.cor))+geom_histogram(bins=15) #hist of correlations from correlograms
boxplot(max.cor)
hist(max.cor)
loc8 = as.numeric(sp.temp[8,3:458])
hist(loc8[loc8>0],probability = T)
summary(loc8[0<loc8 & loc8<63])
# Using package ncf
cor.mult=correlog(sp.temp$long,sp.temp$lat,as.matrix(sp.temp[,3:458]),increment = 1,latlon = T,resamp = 100)
plot(cor.mult)

#Variogram analysis (partial)
params2=params[,c(1,2,5)]
coordinates(params2) = ~ long+lat
g = gstat(id='mean0',formula = mean0~1,data = params2)
expvar = variogram(g)
head(expvar)
plot(expvar)
ggplot(expvar,aes(x=dist,y=gamma,size=np)) + geom_point()

expvar2 <- variogram(g,width=3,cutoff=5,map=TRUE)
plot(expvar2)