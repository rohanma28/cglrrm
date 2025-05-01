
## Install/load packages
library(matrixStats)
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)

##  Initiate starting and ending years
ystart = 1980
yend = 2021

## Specify location where CGLRRM output files should be stored
fileoutloc = "G:/My Drive/Winter 2023/NOAA-USACE-BIL/Data"


## read in LL observations (LLobs)
setwd("G:/My Drive/Winter 2023/NOAA-USACE-BIL/Data/LL/Monthly Mean")
LLobs = read.csv('GLHYD.csv', skip = 12, header = TRUE)
LLobs = LLobs[which(LLobs$year >= ystart),]
LLobs = LLobs[which(LLobs$year <= yend+1),]


## initialize variables for stepping through code without running entire loop
y = ystart # y can be thought of as the modeled time (year the forecast starts)
m = 1 # month of the simulation start

## "i" is used as a index on the number of times the loop has completed
i=0

## This loop starts at the first year of simulation and runs until the last year
## while reading in the monthly mean water levels from Lake Superior and storing
## them in SUP

for (y in ystart:yend){
  
  ## this loop runs through all months of the year
  for (m in 1:12){
    
    
    ##  Set up SUP comparison dataframe starting with 1900 (SUPERIOR)
    setwd(file.path("C:/Users/cwegener/Desktop/CGLRRM output/",as.character(y),
                    "/",as.character(month.abb[m]),"/1900/CGLRRM/output"))
    SUP = read.table("superior.test",header=F,skip= 20,fill=TRUE,row.names=NULL,
                     colClasses = c(rep("character", 4),
                                    rep("numeric", 3),
                                    rep("NULL", 20)))
    colnames(SUP) = c("day","month","year","period","BOM","mean","EOP")
    SUP = unite(SUP,date,day, month, year)
    SUP = select(SUP,-c(period, EOP, BOM))
    SUP$date = lubridate::dmy(SUP$date)
    colnames(SUP) = c("date","1900")
    SUP = SUP[,1:2]
    
    ##  Then go into the loop for 1901 until the end of the simulation period
    ##  stopping one year from yend because there is no model output for the
    ##  current simulation year
    
    j = 3  #current column counter
    sy = 1901 #simulation year (the year the NBS actually occurred)
    
    
    ## This gets the output that is associated with each historical realization
    ## of the net basin supply for a specific month and year in the simulation
    ## period
    
    for (sy in 1901:(y-1)){  
      setwd(file.path("C:/Users/cwegener/Desktop/CGLRRM output/",as.character(y),
                      "/",as.character(month.abb[m]),"/",as.character(as.integer(sy)),
                      "/CGLRRM/output"))
      
      SUP_next = read.table("superior.test",header=F,skip= 20,fill=TRUE,row.names=NULL,
                            colClasses = c(rep("character", 4),
                                           rep("numeric", 3),
                                           rep("NULL", 20)))
      colnames(SUP_next) = c("day","month","year","period","BOM","mean","EOP")
      SUP_next = unite(SUP_next,date,day,month,year)
      SUP_next = select(SUP_next,-c(period,BOM,EOP))
      bottom = matrix(data=NA,ncol=ncol(SUP_next),nrow = (nrow(SUP)-nrow(SUP_next)))
      colnames(bottom) = c('date','Mean')
      SUP_next = rbind(SUP_next,bottom)
      colnames(SUP_next) = c('date','mean')
      SUP = cbind(SUP,SUP_next$mean)
      colnames(SUP)[j]=as.character(sy)
      j = j+1
    }
    
    ## Aggregate by month (take average)
    SUP$month = month(ymd(SUP$date))
    SUP$year = year(ymd(SUP$date))
    # SUP = aggregate (. ~ month + year, 
                     # data = SUP, mean)
    SUP = subset(SUP, select = -c(date))
    SUP$horizon=c(1:12)
    SUP = cbind(SUP[,(ncol(SUP)-2):ncol(SUP)],SUP[,1:(ncol(SUP)-3)])
    SUP$startmo = paste(as.character(month.abb[m]),as.character(y))
    SUP = cbind(SUP[,1:3],SUP[,ncol(SUP)],SUP[,4:(ncol(SUP)-1)])
    colnames(SUP)[1:4] = c("month","year","horizon","startmo")     
    SUP$obsLL = LLobs$Superior[(i+1):(i+12)]
    SUP = cbind(SUP[,1:4],SUP$obsLL,SUP[,5:(ncol(SUP)-1)])
    colnames(SUP)[5] = "obsLL"
    
    ## if it's the first year, store SUP into SUP_main            
    if(i == 0){ 
      SUP_main = SUP
      
    ## if it's not the first year, append the SUP dataframe to the bottom of SUP_main
    }else{
      
      ## after each simulation year, there is another year of historical NBS
      ## that can be used so there will be another column in SUP so SUP_main
      ## is given a column of NAs so that the number of columns match up for
      ## using the rbind command
      
      if(ncol(SUP_main)==ncol(SUP)){
        SUP_main = rbind(SUP_main,SUP)
      }else{
        newcol = rep(NA,nrow(SUP_main))
        SUP_main = cbind(SUP_main,newcol)
        colnames(SUP_main)[ncol(SUP_main)] = colnames(SUP[ncol(SUP)])
        SUP_main = rbind(SUP_main,SUP)
      }
      
    }
    
    i=i+1
  }
}


## Write output to .csv file

setwd(fileoutloc)
write.csv(SUP_main,file=paste("CGLRRMout_SUP_",ystart,"_to_",yend,"_fewcolumns.csv",sep=""))




###########################################################################
############################################################################

                    #     Michigan Huron (MHU)

###########################################################################
############################################################################

y = ystart
m = 1
i=0

for (y in ystart:yend){
  for (m in 1:12){
    
    
    ##  Set up MHU comparison dataframe starting with 1900 (MHU)
    setwd(file.path("C:/Users/cwegener/Desktop/CGLRRM output/",as.character(y),
                    "/",as.character(month.abb[m]),"/1900/CGLRRM/output"))
    MHU = read.table("micherie.test",header=F,skip= 26,fill=TRUE,row.names=NULL,
                     colClasses = c(rep("character", 6),
                                    rep("numeric", 2),
                                    rep("NULL", 20)))
    colnames(MHU) = c("yrday","day","day2","month","year","period","USslip","mean")
    MHU = unite(MHU,date,day, month, year)
    MHU = select(MHU,-c(yrday,day2,period,USslip))
    MHU$date = lubridate::dmy(MHU$date)
    colnames(MHU) = c("date","1900")
    MHU = MHU[,1:2]
    
    ##  Then go into the loop for 1901:2017  
    
    j = 3
    sy = 1901
    
    for (sy in 1901:(y-1)){  
      setwd(file.path("C:/Users/cwegener/Desktop/CGLRRM output/",as.character(y),
                      "/",as.character(month.abb[m]),"/",as.character(as.integer(sy)),
                      "/CGLRRM/output"))
      
      MHU_next = read.table("micherie.test",header=F,skip= 26,fill=TRUE,row.names=NULL,
                            colClasses = c(rep("character", 6),
                                           rep("numeric", 2),
                                           rep("NULL", 20)))
      colnames(MHU_next) = c("yrday","day","day2","month","year","period","USslip","mean")
      MHU_next = unite(MHU_next,date,day, month, year)
      MHU_next = select(MHU_next,-c(yrday,day2,period,USslip))
      bottom = matrix(data=NA,ncol=ncol(MHU_next),nrow = (nrow(MHU)-nrow(MHU_next)))
      colnames(bottom) = c('date','Mean')
      MHU_next = rbind(MHU_next,bottom)
      colnames(MHU_next) = c('date','mean')
      MHU = cbind(MHU,MHU_next$mean)
      colnames(MHU)[j]=as.character(sy)
      j = j+1
    }
    
    ## Aggregate by month (take average)
    MHU$month = month(ymd(MHU$date))
    MHU$year = year(ymd(MHU$date))
    MHU = aggregate (. ~ month + year, 
                     data = MHU, mean)
    MHU = subset(MHU, select = -c(date))
    
    MHU$horizon=c(1:12)
    MHU = cbind(MHU[,1],MHU$horizon,MHU[,2:(ncol(MHU)-1)])
    MHU$startmo = paste(as.character(month.abb[m]),as.character(y))
    MHU = cbind(MHU[,1],MHU[,3],MHU$startmo,MHU[,2],MHU[,4:(ncol(MHU)-1)])
    colnames(MHU)[1:4] = c("month","year","startmo","horizon") 
    
    MHU$obsLL = LLobs$Michigan.Huron[(i+1):(i+12)]
    MHU = cbind(MHU[,1:4],MHU$obsLL,MHU[,5:(ncol(MHU)-1)])
    colnames(MHU)[5] = "obsLL"
    
    if(i == 0){
      MHU_main = MHU
    }
    else{
      
      if(ncol(MHU_main)==ncol(MHU)){
        MHU_main = rbind(MHU_main,MHU)
      }
      else{
        newcol = rep(NA,nrow(MHU_main))
        MHU_main = cbind(MHU_main,newcol)
        colnames(MHU_main)[ncol(MHU_main)] = colnames(MHU[ncol(MHU)])
        MHU_main = rbind(MHU_main,MHU)
      }
      
    }
    
    i = i+1
  }
}

## Write output to .csv file

setwd(fileoutloc)
write.csv(MHU_main,file=paste("CGLRRMout_MHU_",ystart,"_to_",yend,"_fewcolumns.csv",sep=""))


###########################################################################
############################################################################

#     Erie (ERI)

###########################################################################
############################################################################


y = ystart
m = 1
i=0

for (y in ystart:yend){
  for (m in 1:12){
    
    
    ##  Set up ERI comparison dataframe starting with 1900 (ERI)
    setwd(file.path("C:/Users/cwegener/Desktop/CGLRRM output/",as.character(y),
                    "/",as.character(month.abb[m]),"/1900/CGLRRM/output"))
    ERI = read.table("micherie.test",header=F,skip= 26,fill=TRUE,row.names=NULL,
                     colClasses = c(rep("character", 6),
                                    rep("numeric", 8),
                                    rep("NULL", 20)))
    colnames(ERI) = c("yrday","day","day2","month","year","period","USslip","MHUmean"
                      ,"MHUEOP","STCriv","STCmean","STCEOP","DETriv","mean")
    ERI = unite(ERI,date,day, month, year)
    ERI = select(ERI,-c(yrday,day2,period,USslip,MHUmean,MHUEOP,STCriv,STCmean,
                        STCEOP,DETriv))
    ERI$date = lubridate::dmy(ERI$date)
    colnames(ERI) = c("date","1900")
    
    ##  Then go into the loop for 1901:2017  
    
    j = 3
    sy = 1901
    
    for (sy in 1901:(y-1)){  
      setwd(file.path("C:/Users/cwegener/Desktop/CGLRRM output/",as.character(y),
                      "/",as.character(month.abb[m]),"/",as.character(as.integer(sy)),
                      "/CGLRRM/output"))
      
      ERI_next = read.table("micherie.test",header=F,skip= 26,fill=TRUE,row.names=NULL,
                            colClasses = c(rep("character", 6),
                                           rep("numeric", 8),
                                           rep("NULL", 20)))
      colnames(ERI_next) = c("yrday","day","day2","month","year","period","USslip","MHUmean"
                             ,"MHUEOP","STCriv","STCmean","STCEOP","DETriv","mean")
      ERI_next = unite(ERI_next,date,day, month, year)
      ERI_next = select(ERI_next,-c(yrday,day2,period,USslip,MHUmean,MHUEOP,STCriv,STCmean,
                                    STCEOP,DETriv))
      bottom = matrix(data=NA,ncol=ncol(ERI_next),nrow = (nrow(ERI)-nrow(ERI_next)))
      colnames(bottom) = c('date','Mean')
      ERI_next = rbind(ERI_next,bottom)
      colnames(ERI_next) = c('date','mean')
      ERI = cbind(ERI,ERI_next$mean)
      colnames(ERI)[j]=as.character(sy)
      j = j+1
    }
    
    ## Aggregate by month (take average)
    ERI$month = month(ymd(ERI$date))
    ERI$year = year(ymd(ERI$date))
    ERI = aggregate (. ~ month + year, 
                     data = ERI, mean)
    ERI = subset(ERI, select = -c(date))
    ERI$horizon=c(1:12)
    ERI = cbind(ERI[,1],ERI$horizon,ERI[,2:(ncol(ERI)-1)])
    ERI$startmo = paste(as.character(month.abb[m]),as.character(y))
    ERI = cbind(ERI[,1],ERI[,3],ERI$startmo,ERI[,2],ERI[,4:(ncol(ERI)-1)])
    colnames(ERI)[1:4] = c("month","year","startmo","horizon") 
    
    ERI$obsLL = LLobs$Erie[(i+1):(i+12)]
    ERI = cbind(ERI[,1:4],ERI$obsLL,ERI[,5:(ncol(ERI)-1)])
    colnames(ERI)[5] = "obsLL"
    
    
    if(i == 0){
      ERI_main = ERI
    }
    else{
      
      if(ncol(ERI_main)==ncol(ERI)){
        ERI_main = rbind(ERI_main,ERI)
      }
      else{
        newcol = rep(NA,nrow(ERI_main))
        ERI_main = cbind(ERI_main,newcol)
        colnames(ERI_main)[ncol(ERI_main)] = colnames(ERI[ncol(ERI)])
        ERI_main = rbind(ERI_main,ERI)
      }
      
    }
    
    
    i = i+1
  }
}

## Write output to .csv file

setwd(fileoutloc)
write.csv(ERI_main,file=paste("CGLRRMout_ERI_",ystart,"_to_",yend,"_fewcolumns.csv",sep=""))

