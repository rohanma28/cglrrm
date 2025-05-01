## The purpose of this script is to recreate the USACE Ensemble product output
## for 2018 using residual NBS sequences between 1900 and 2017


## For now, assume that the input file structure is set up and just work with
## opening the executable file and dealing with the output

##  Load necessary packages

library(tidyr)
library(dplyr)
library(lubridate)
library(stringr)
library(matrixStats)


## Read in starting LL

setwd("G:/My Drive/Winter 2023/NOAA-USACE-BIL/Data/LL/Beginning of Month")

SUP_SLL = read.csv("SUP_BOMLL.csv",skip=9,header=T)
MHU_SLL = read.csv("MHU_BOMLL.csv",skip=9,header=T)
ERI_SLL = read.csv("ERI_BOMLL.csv",skip=9,header=T)
STC_SLL = read.csv("STC_BOMLL.csv",skip=9,header=T)


##  Define start and end years
mo = 1980
mf = 2021


##  Start for loops that modify input parameter files and run CGLRRM


## Loop through 2000-2022 (m = simulation year)
for(m in mo:mf){
  
  setwd("C:/Users/Caleb/Desktop/CGLRRM_model_run/")
  dir.create(as.character(m))

  ## Loop through the 12 months of simulation (l = simulation month)
  for (l in 1:12){
    setwd(file.path("C:/Users/Caleb/Desktop/CGLRRM_model_run/",as.character(m)))
    dir.create(as.character(month.abb[l]))
    
    ##  Loop through years 1900 until 2017 (i = NBS year)
    for (i in 1900:(m-1)){
      
      ##  Set model run directory

      setwd(file.path("C:/Users/Caleb/Desktop/CGLRRM_model_run/",as.character(m),"/",as.character(month.abb[l])))
      dir.create(as.character(i))
      file.copy("C:/Users/Caleb/Desktop/CGLRRM_model_run/CGLRRM",
                file.path("C:/Users/Caleb/Desktop/CGLRRM_model_run/",
                as.character(m),"/",as.character(month.abb[l]),"/",as.character(i)),
                copy.mode = T,copy.date = T,recursive = T)
      setwd(file.path("C:/Users/Caleb/Desktop/CGLRRM_model_run/",
        as.character(m),"/",as.character(month.abb[l]),"/",as.character(i),"/CGLRRM"))
      
      
      ##  Change the parameter file to read in the next year of NBS sequences
      params = readLines("CGLRRM_params.2008")
      params[1] = str_c("Sup NBS Data: ./input/",as.character(i),"/",as.character(month.abb[l]),"/MNBS_2008_sp.txt")
      params[2] = str_c("MHu NBS data: ./input/",as.character(i),"/",as.character(month.abb[l]),"/MNBS_2008_mh.txt")
      params[3] = str_c("Eri NBS data: ./input/",as.character(i),"/",as.character(month.abb[l]),"/MNBS_2008_er.txt")
      params[4] = str_c("Stc NBS data: ./input/",as.character(i),"/",as.character(month.abb[l]),"/MNBS_2008_sc.txt")
      params[5] = str_c("Output Directory: ./output/")
      params[6] = str_c("Sup Start Level      : ",as.character(SUP_SLL[which(SUP_SLL$Year == m),(l+1)])," m")
      params[7] = str_c("MHu Start Level      : ",as.character(MHU_SLL[which(MHU_SLL$Year == m),(l+1)])," m")
      params[8] = str_c("Eri Start Level      : ",as.character(ERI_SLL[which(ERI_SLL$Year == m),(l+1)])," m")
      params[9] = str_c("St. C Start Level    : ",as.character(STC_SLL[which(STC_SLL$Year == m),(l+1)])," m")
      
      #Set start and end dates of simulation
      
      # January start months
      if(l==1){
        params[10] = str_c("Start Date: ", as.character(m), ",1,1")
        params[11] =str_c("End Date: ",as.character(m), ",12,31")
      }  
      
      # start months with corresponding end months with 31 days (Feb,Apr,Jun,Aug,Nob)
      if(l==2|l==4|l==6|l==8|l==9|l==11){
        params[10] = str_c("Start Date: ", as.character(m), ",", as.character(l),",1")
        params[11] =str_c("End Date: ",as.character(m+1),"," , as.character(l-1),",31")
      }
      
      # March starting month and February ending month.....
      # Assume that February is always 28 days long (ignore leap years)
      if(l==3){
        params[10] = str_c("Start Date: ", as.character(m), ",3,1")
        
        if(leap_year(m+1)){
        params[11] =str_c("End Date: ",as.character(m+1), ",2,29")
        }
        else{
          params[11] =str_c("End Date: ",as.character(m+1), ",2,28")
        }
      }  
      
      if(l==5|l==7|l==10|l==12){
        params[10] = str_c("Start Date: ", as.character(m), ",", as.character(l),",1")
        params[11] =str_c("End Date: ", as.character(m+1), ",", as.character(l-1),",30")
      } 
      
      ##  Rewrite the parameter file for the next round
      sink("CGLRRM_params.2008")
      for (k in 1:249){
        cat(as.character(params[k]))
        cat("\n")
      }
      
      ##  Close the file
      sink()
      
      ##  Run the model
      system("cglrrm")
      
      ##  Delete input files
      unlink('input', recursive = TRUE)
      unlink('CGLRRM.exe', recursive = TRUE)
       
      setwd(file.path("C:/Users/Caleb/Desktop/CGLRRM_model_run/",
                      as.character(m),"/",as.character(month.abb[l]),"/",
                      as.character(i),"/CGLRRM/output"))
      
      
    }
  }
}

