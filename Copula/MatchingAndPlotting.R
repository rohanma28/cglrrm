setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\Copula")
library(copula)
library(VineCopula)
library(rvinecopulib)
library(FAdist)
library(brms)
library(png)
library(gridExtra)
library(grid)

# set number of samples as a global variable 
n <- 100
forecast_month <- 1

# Define a list of month and lake abbreviations
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul","Aug","Sep","Oct","Nov","Dec")
lake_names  <- c("Eri", "Ont", "MiH", "Sup", "StC") 

setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\Copula\\formatted_input")
cal_inps <- data.matrix(read.csv(paste(month.abb[forecast_month],"_12Forecast_3Ant_CopulaInput.csv",sep=""))[1:70,2:226])

# load the names of the vine copula model files 
model_name <- paste(month.abb[forecast_month],"_15month_1950_2020_RvineDist.rds",sep="")
samp_name  <- paste(month.abb[forecast_month],"_15month_1950_2020_samp.rds",sep="")

# load the L2S confidence intervals
setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\Copula\\R_objects")
CI95s <- readRDS("1950_2020_L2S_Avg_CI95_width.rds")

### Define a function that randomly samples n sequences
random_sample <- function(sample_num){
  n      <- sample_num
  
  # read in the proper vine distribution model
  vd     <- readRDS(model_name)
  
  samps  <- readRDS(samp_name)
  colnames(samps) <- colnames(cal_inps)
  
  # Create empty matrix of the proper dimensions and column titles to hold forecast values and antecedent matches
  samp_row     <- cal_inps[1,]
  samp_mat     <- t(replicate(n,samp_row,simplify="row"))
  
  # Grab random rows from the sample matrix 
  samp_mat <- samps[sample(1:nrow(samps),n,replace=FALSE),]
  
  return(samp_mat)
}



## Now, define a function that takes a conditionally sampled matrix and returns a vector of NBS 12 months into the future
forecast_lake_month_NBS <- function(matched_sample, lake_name, month_name){
  lake_month <- matched_sample[,intersect(which(substr(colnames(matched_sample),1,3)=="for"),
                                              intersect(which(substr(colnames(matched_sample),9,11)==lake_name),
                                                        which(substr(colnames(matched_sample),5,7)==month_name)))]
  NBS <- as.matrix(lake_month[,which(substr(colnames(lake_month),13,13)=="P")]-
                     lake_month[,which(substr(colnames(lake_month),13,13)=="E")]+
                     lake_month[,which(substr(colnames(lake_month),13,13)=="R")],ncol=1)
  colnames(NBS) <- paste("for.",month_name,".",lake_name,".NBS",sep="")
  return(NBS)
}

## Now define a function that creates a by-lake list of matrices of all month-lake NBS combinations
forecast_complete_NBS <- function(matched_sample, sample_num){
  # create matrix to store NBS 
  complete_NBS <- vector(mode="list",length=length(lake_names))
  # loop through lakes and months, but want months in the order that they appear in the forecast
  if (forecast_month==1){ordered_months=month_names}
  else{ordered_months = c(month_names[forecast_month:12],month_names[1:(forecast_month-1)])}
  for (lake in 1:length(lake_names)){
    complete_NBS[[lake]] = matrix(nrow=sample_num,ncol=0)
    for (month in 1:length(ordered_months)){
      complete_NBS[[lake]] <- cbind(complete_NBS[[lake]],forecast_lake_month_NBS(matched_sample,lake_name=lake_names[lake],month_name=ordered_months[month]))
    }
    complete_NBS[[lake]] = data.matrix(as.data.frame(complete_NBS[[lake]]))
  }
  names(complete_NBS) <- lake_names
  return(complete_NBS)
}

######################
# Generate Random Samples for each validation month
######################
setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\Copula\\R_objects")

rval_samps <- vector(mode="list", length=70)
count = 1
for (year in 1:70){
  print(count)
  names(rval_samps)[count] <- paste("year",year,sep="")
  rsamples  <- random_sample(sample_num=n)
  rval_samps[[count]]      <- forecast_complete_NBS(matched_sample = rsamples, sample_num=n)
  count = count+1
}

saveRDS(rval_samps,file="random_validation_samples_100k_100_2021_2090.rds")
rval_samps <- readRDS(file="random_validation_samples_100k_100_2021_2090.rds")

dir.create("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\output\\flatshiftrunoffsim2\\")

for (sim in 1:n) {

  setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\Copula\\output")
  
  eri_sample <- data.frame(matrix(0, 70, 12))
  mih_sample <- data.frame(matrix(0, 70, 12))
  sup_sample <- data.frame(matrix(0, 70, 12))
  stc_sample <- data.frame(matrix(0, 70, 12))
  
  for (year in 1:70) {
    for (month in 1:12) {
      eri_sample[year, month] <- rval_samps[[year]][[1]][sim,month]
      mih_sample[year, month] <- rval_samps[[year]][[3]][sim,month]
      sup_sample[year, month] <- rval_samps[[year]][[4]][sim,month]
      stc_sample[year, month] <- rval_samps[[year]][[5]][sim,month]
    }
  }
  
  colnames(eri_sample) <- colnames(rval_samps[[1]][[1]])
  colnames(mih_sample) <- colnames(rval_samps[[1]][[3]])
  colnames(sup_sample) <- colnames(rval_samps[[1]][[4]])
  colnames(stc_sample) <- colnames(rval_samps[[1]][[5]])
  
  write.csv(eri_sample, "Erie_2020_2090_NBS_Forecast.csv")
  write.csv(mih_sample, "Michigan_2020_2090_NBS_Forecast.csv")
  write.csv(sup_sample, "Superior_2020_2090_NBS_Forecast.csv")
  write.csv(stc_sample, "StClair_2020_2090_NBS_Forecast.csv")
  
  source("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\r_code\\input_file_structure.R")
  
  setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan")
  
  system("cglrrm_test.exe")
  
  new_dir <- paste("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\output\\flatshiftrunoffsim2\\", sim, sep = "")
  dir.create(new_dir)
  setwd(new_dir)

  write.csv(read.table("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\output\\spmmlv.test", skip = 18), "supforecast.csv", row.names = F)
  write.csv(read.table("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\output\\scmmlv.test", skip = 19), "stclairforecast.csv", row.names = F)
  write.csv(read.table("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\output\\mhmmlv.test", skip = 19), "mihurforecast.csv", row.names = F)
  write.csv(read.table("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\output\\ermmlv.test", skip = 19), "erieforecast.csv", row.names = F)

  source("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan\\r_code\\plotting.R")
  
}

setwd("C:\\Users\\rathreya\\Downloads\\CGLRRM_Rohan")
