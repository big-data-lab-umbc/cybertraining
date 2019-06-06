## This file creates simultion ensembles using parametric bootstrap

setwd("D:/VIC/forcing") #Point to orignal folder with the forcing files
library(EnvStats)

file_names=list.files(pattern="forcing_*")                        #gets file names from working directory; can specify path
files=lapply(file_names , read.table ,header = F,sep = ' ',numerals = 'no.loss',stringsAsFactors = F) #read in all file names
mu0 = as.numeric(387);cv0 = as.numeric(387) #gamma parameters for each grid point
for(i in 1:387)
{
  precip = files[[i]]$V1                                        #extract precip into a different variable
  param = egammaAlt(precip[precip>0])                           #gamme parameters for positive part
  mu0[i] = as.vector(param$parameters[1])                       #mean estimate for i-th grid point
  cv0[i] = as.vector(param$parameters[2])                       #cv estimate for i-th grid point
  
}
path = as.character(1:100)                                      #create folder names for 100 new folders for bootstrap
for(j in 1:100)
{
  dir.create(path[j])                                           #create numbered directories
  file_names2 = paste0(path[j],'/',file_names)                  #path for writing out the files
  
  for(i in 1:387)
  {
    precip[precip>0] = rgammaAlt(sum(precip>0),mu0[i],cv0[i]*2) 
    files[[i]]$V1 = precip                                       #rewrite it back to the file
    
    write.table(files[[i]],file_names2[i],row.names = F,col.names = F,quote = F)  #write files out in folder
  }
  
}

###############################################
params = data.frame(p0,mu0,cv0)
write.csv(params,"params.csv")
getwd()


##write out global param files
##original global.params need to be in working directory. Creates 300 new files
cvlist = c('cvx1','cvx15','cvx2')
for(i in 1:3)
{
  for(j in 1:100)
  {
    x=readLines("global.params")
    replace=paste0('meteo_sn/',cvlist[i],'/',path[j]) #replace original meteo_sn path with something like cvx2/20
    filename = paste0('global.params.',cvlist[i],'.',path[j]) #Save it as a global param corresponding to the foldername from above
    y=gsub("meteo_forcing",replace,x)
    outloc = paste0("output_sn/",cvlist[i],'/',path[j],'/') #output in a similar directory
    z=gsub('output_sn/',outloc,y)
    cat(z,file=filename,sep = '\n')
  }
}


## create blank directories in taki for output because for some reason VIC output would not create new dirs, only write in existing dirs
for(c in cvlist)
{
  dir.create(c)
}
for(i in 1:100)
  dir.create(as.character(i))
