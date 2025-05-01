## The purpose of this code is to setup the input file structure for reading
## in NBS sequences into the CGLRRM for use in recreation of the USACE 
## historical NBS sequence Lake Level Ensemble Product

ystart = 2021
yend = 2090

#########################################

##    SUPERIOR

#########################################

## Read in the data
setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/Copula/output")
NBS = read.csv("Superior_2020_2090_NBS_Forecast.csv", header=T, row.names = 1)
superior_area = 82100000000

## Subset the data
#NBS = subset(NBS,Year=="2018")

#setwd(file.path("G:\\My Drive\\Winter 2023\\NOAA-USACE-BIL\\CGLRRM\\TEST\\1may2023_feb2018\\input\\",paste(as.character(j))))
## Open output file

setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/input")

sink("MNBS_2008_sp.txt")

## write the output
cat("# LAKE SUPERIOR MONTHLY NET BASIN SUPPLY (10m3s)")
cat("\n")
cat("# FORTRAN FORMAT: I4 12F8.0")
cat("\n")
cat("# 2021-2090")
cat("\n")
cat("UNITS:m3s")
cat("\n")
cat("INTERVAL:m")
cat("\n")
cat("####     Jan     Feb     Mar     Apr     May     Jun     Jul     Aug     Sep     Oct     Nov     Dec")
cat("\n")

for (k in ystart:yend){
  cat(as.character(k),
      format(round(NBS[k-ystart+1,1]*superior_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,2]*superior_area/1000/(86400*28)),width=7),
      format(round(NBS[k-ystart+1,3]*superior_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,4]*superior_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,5]*superior_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,6]*superior_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,7]*superior_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,8]*superior_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,9]*superior_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,10]*superior_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,11]*superior_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,12]*superior_area/1000/(86400*31)),width=7))
  cat("\n")
}

## Close the output file
sink()

#########################################

##    MICHIGAN.HURON

#########################################


## Read in the data
setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/Copula/output")
NBS = read.csv("Michigan_2020_2090_NBS_Forecast.csv", header=T, row.names = 1)
michigan_area = 117400000000

## Subset the data
#NBS = subset(NBS,Year=="2018")

setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/input")

## Open output file
sink("MNBS_2008_mh.txt")

## write the output
cat("# LAKE MICHIGAN-HURON MONTHLY NET BASIN SUPPLY (10m3s)")
cat("\n")
cat("# FORTRAN FORMAT: I4 12F8.0")
cat("\n")
cat("# 2021-2090")
cat("\n")
cat("UNITS:m3s")
cat("\n")
cat("INTERVAL:m")
cat("\n")
cat("####     Jan     Feb     Mar     Apr     May     Jun     Jul     Aug     Sep     Oct     Nov     Dec")
cat("\n")

for (k in ystart:yend){
  cat(as.character(k),
      format(round(NBS[k-ystart+1,1]*michigan_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,2]*michigan_area/1000/(86400*28)),width=7),
      format(round(NBS[k-ystart+1,3]*michigan_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,4]*michigan_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,5]*michigan_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,6]*michigan_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,7]*michigan_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,8]*michigan_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,9]*michigan_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,10]*michigan_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,11]*michigan_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,12]*michigan_area/1000/(86400*31)),width=7))
  cat("\n")
}
## Close the output file
sink()

#########################################

##    ERIE

#########################################


## Read in the data
setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/Copula/output")
NBS = read.csv("Erie_2020_2090_NBS_Forecast.csv", header=T, row.names = 1)
erie_area = 25700000000

## Subset the data
#NBS = subset(NBS,Year=="2018")

setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/input")

## Open output file
sink("MNBS_2008_er.txt")

## write the output
cat("# LAKE ERIE MONTHLY NET BASIN SUPPLY (10m3s)")
cat("\n")
cat("# FORTRAN FORMAT: I4 12F8.0")
cat("\n")
cat("# 2021-2090")
cat("\n")
cat("UNITS:m3s")
cat("\n")
cat("INTERVAL:m")
cat("\n")
cat("####     Jan     Feb     Mar     Apr     May     Jun     Jul     Aug     Sep     Oct     Nov     Dec")
cat("\n")

for (k in ystart:yend){
  cat(as.character(k),
      format(round(NBS[k-ystart+1,1]*erie_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,2]*erie_area/1000/(86400*28)),width=7),
      format(round(NBS[k-ystart+1,3]*erie_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,4]*erie_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,5]*erie_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,6]*erie_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,7]*erie_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,8]*erie_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,9]*erie_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,10]*erie_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,11]*erie_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,12]*erie_area/1000/(86400*31)),width=7))
  cat("\n")
}
## Close the output file
sink()

#########################################

##    ST.CLAIR

#########################################


## Read in the data
setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/Copula/output")
NBS = read.csv("StClair_2020_2090_NBS_Forecast.csv", header=T, row.names = 1)
stclair_area = 1114000000

## Subset the data
#NBS = subset(NBS,Year=="2018")

setwd("/Volumes/Hydro/projects/routing_models/CGLRRM_Rohan/input/")
    
## Open output file
sink("MNBS_2008_sc.txt")

## write the output
cat("# LAKE ST. CLAIR MONTHLY NET BASIN SUPPLY (10m3s)")
cat("\n")
cat("# FORTRAN FORMAT: I4 12F8.0")
cat("\n")
cat("# 2021-2090")
cat("\n")
cat("UNITS:m3s")
cat("\n")
cat("INTERVAL:m")
cat("\n")
cat("####     Jan     Feb     Mar     Apr     May     Jun     Jul     Aug     Sep     Oct     Nov     Dec")
cat("\n")

for (k in ystart:yend){
  cat(as.character(k),
      format(round(NBS[k-ystart+1,1]*stclair_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,2]*stclair_area/1000/(86400*28)),width=7),
      format(round(NBS[k-ystart+1,3]*stclair_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,4]*stclair_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,5]*stclair_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,6]*stclair_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,7]*stclair_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,8]*stclair_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,9]*stclair_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,10]*stclair_area/1000/(86400*31)),width=7),
      format(round(NBS[k-ystart+1,11]*stclair_area/1000/(86400*30)),width=7),
      format(round(NBS[k-ystart+1,12]*stclair_area/1000/(86400*31)),width=7))
  cat("\n")
}
## Close the output file
sink()
