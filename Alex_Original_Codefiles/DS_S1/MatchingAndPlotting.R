rm(list=ls())
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
library(copula)
library(VineCopula)
library(rvinecopulib)
library(FAdist)
library(brms)
library(png)
library(gridExtra)
library(grid)

# set number of samples as a global variable 
n <- 1000

# Define a list of month and lake abbreviations
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul","Aug","Sep","Oct","Nov","Dec")
lake_names  <- c("Eri", "Ont", "MiH", "Sup") 

# Load validation data for 2006-2020, separate out antecedents
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/formatted_input")
val_inps <- vector("list", length=12)
val_ants <- vector("list", length=12)
for (month in 1:12){
  val_inps[[month]] <- read.csv(paste(month.abb[month],"_12Forecast_3Ant_CopulaInput.csv",sep=""))[57:70,2:181]
  val_ants[[month]] <- val_inps[[month]][,which(substr(colnames(val_inps[[month]]),1,3)=="ant")]
}
names(val_inps) <- month.abb
names(val_ants) <- month.abb

cal_inps <- vector("list",length=12)
for (month in 1:12){cal_inps[[month]] <- data.matrix(read.csv(paste(month.abb[month],"_12Forecast_3Ant_CopulaInput.csv",sep=""))[1:56,2:181])}
names(cal_inps) = month.abb


# load the names of the vine copula model files 
model_names <- vector(length=12)
samp_names  <- vector(length=12)
for (month in 1:12){
  model_names[month] <- paste(month.abb[month],"_15month_1950_2006_RvineDist.rds",sep="")
  samp_names[month]  <- paste(month.abb[month],"_15month_1950_2006_samp.rds",sep="")
}

# load the L2S confidence intervals
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
CI95s <- readRDS("1950_2006_L2S_Avg_CI95_width.rds")

#################
### Function definitions
#################

# Matching function which selects n samples from the provided set of samples based on antecedent matching
# The matching is based on minimizing the difference between antecedent observations and sampled values, 
# normalized by the standard deviation (approximated as 95%CI/4) of the L2S samples for that variable. 
# This gives matching "scores" for each member where a lower score indicates a closer match. 
# To then weight each member, take the inverse of its matching score and divide it by the min of the inverse of all matching scores. 
# This will give a number >1 which (rounded) is the number of times to repeat that member. 
# Then, thin the new weighted ensemble down to n members. 
matching_sample <- function(antecedents,forecast_month,sample_num){
  n      <- sample_num
  ants   <- antecedents
  fmonth <- forecast_month
  
  # read in the proper vine distribution model
  vd     <- readRDS(model_names[fmonth])
  
  samps  <- readRDS(samp_names[fmonth])
  colnames(samps) <- colnames(val_inps[[fmonth]])

  # read in the proper L2S estimated average 95% intervals
  CIs <- t(data.matrix(as.data.frame(CI95s[[fmonth]])))
  
  # Create empty matrix of the proper dimensions and column titles to hold forecast values and antecedent matches
  samp_row     <- val_inps[[fmonth]][1,]
  samp_row[1,] <- NA
  samp_mat     <- t(replicate(n,samp_row,simplify="row"))
  
  # Create matrix of antecedent matching values (not summed across variables yet)
  match_mat <- abs(sweep(data.matrix(as.data.frame(samps[,colnames(ants)])),2,data.matrix(ants))/(CIs[,colnames(ants)]/4))
  
  # Take the inverse of the matching scores
  inv_match <- 1/rowSums(match_mat)
  
  # Normalize to the minimum of the inverse of the matching scores
  weights <- round(inv_match/min(inv_match))

  # Weight each member by repeating it the number of times as its weight
  weighted_samps <- samps[rep(1:nrow(samps), times=weights),]
  
  # Sample n members from the weighted ensemble
  samp_mat <- weighted_samps[sample(1:nrow(weighted_samps), size=n, replace=FALSE),]
  
  ### Superseded
  # # Select the n samples which minimize the matching formula
  # samp_row_nums <- order(rowSums(match_mat))[1:n]
  # 
  # # Grab the rows from the sample matrix corresponding to these matched sample indices 
  # samp_mat <- samps[samp_row_nums,]
  
  return(samp_mat)
  
}

### Define a function that randomly samples n sequences
random_sample <- function(forecast_month,sample_num){
  n      <- sample_num
  fmonth <- forecast_month
  
  # read in the proper vine distribution model
  vd     <- readRDS(model_names[fmonth])
  
  samps  <- readRDS(samp_names[fmonth])
  colnames(samps) <- colnames(val_inps[[fmonth]])
  
  # Create empty matrix of the proper dimensions and column titles to hold forecast values and antecedent matches
  samp_row     <- val_inps[[fmonth]][1,]
  samp_row[1,] <- NA
  samp_mat     <- t(replicate(n,samp_row,simplify="row"))
  
  # Grab random rows from the sample matrix 
  samp_mat <- samps[sample(1:100000,n,replace=FALSE),]
  
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
forecast_complete_NBS <- function(matched_sample, forecast_month, sample_num){
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

## Now define a function that converts a list of NBS matrices to cumulative NBS
cumulative <- function(NBS_list){
  cumulative_NBS <- vector(mode="list",length=length(lake_names))
  for (lake in 1:length(lake_names)){
    cumulative_NBS[[lake]] <- t(apply(NBS_list[[lake]],1,cumsum))
    colnames(cumulative_NBS[[lake]]) <- paste(colnames(NBS_list[[lake]]),".cum",sep="")
  }
  names(cumulative_NBS) <- lake_names
  return(cumulative_NBS)
}

#################

######################
# Calculate calibration NBS sequences
######################
cal_NBS        <- vector(mode="list",length=12)
cum_cal_NBS    <- vector(mode="list",length=12)
names(cal_NBS) <- month_names
names(cum_cal_NBS)    <- month_names
for (month in 1:12){
  cal_NBS[[month]]     <- forecast_complete_NBS(cal_inps[[month]],forecast_month=month,sample_num=length(cal_inps[[1]][,1]))
  cum_cal_NBS[[month]] <- cumulative(cal_NBS[[month]])
}

cum_cal_NBS_ts_1  <- vector(mode="list", length=4)
cum_cal_NBS_ts_3  <- vector(mode="list", length=4)
cum_cal_NBS_ts_6  <- vector(mode="list", length=4)
cum_cal_NBS_ts_12 <- vector(mode="list", length=4)
names(cum_cal_NBS_ts_1)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_cal_NBS_ts_3)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_cal_NBS_ts_6)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_cal_NBS_ts_12) <- c("Eri", "Ont", "MiH", "Sup")
for (lake in 1:4){
  cum_cal_NBS_ts_1[[lake]] <- data.frame(year=rep(seq(1951,2006), each=12),
                                         month=rep(1:12, times=56),
                                         forecast_month = rep(1:12, times = 56),
                                         cum_NBS = rep(NA, 672))
  rep_count = 1
  for (year in 1:56){
    for (month in 1:12){
      value <- cum_cal_NBS[[month]][[lake]][year,1]
      cum_cal_NBS_ts_1[[lake]]$cum_NBS[rep_count] <- value
      rep_count = rep_count+1
    }
  }
}
for (lake in 1:4){
  cum_cal_NBS_ts_3[[lake]] <- data.frame(year=rep(seq(1951,2006), each=12),
                                         month=rep(1:12, times=56),
                                         forecast_month = rep(NA, times = 672),
                                         cum_NBS = rep(NA, 672))
  cum_cal_NBS_ts_3[[lake]]$forecast_month <- cum_cal_NBS_ts_3[[lake]]$month-2
  cum_cal_NBS_ts_3[[lake]]$forecast_month[cum_cal_NBS_ts_3[[lake]]$forecast_month < 1] <- cum_cal_NBS_ts_3[[lake]]$forecast_month[cum_cal_NBS_ts_3[[lake]]$forecast_month < 1]+12
  cum_cal_NBS_ts_3[[lake]]$forecast_month[c(1,2)] <- NA
  rep_count <- 3
  exit_flag=F
  for (year in 1:56){
    for (month in 1:12){
      if(rep_count>672){exit_flag=T}
      if(exit_flag){break}
      value <- cum_cal_NBS[[month]][[lake]][year,3]
      cum_cal_NBS_ts_3[[lake]]$cum_NBS[rep_count] <- value
      rep_count <- rep_count + 1
    }
    if(exit_flag){break}
  }
}
for (lake in 1:4){
  cum_cal_NBS_ts_6[[lake]] <- data.frame(year=rep(seq(1951,2006), each=12),
                                         month=rep(1:12, times=56),
                                         forecast_month = rep(NA, times = 672),
                                         cum_NBS = rep(NA, 672))
  cum_cal_NBS_ts_6[[lake]]$forecast_month <- cum_cal_NBS_ts_6[[lake]]$month-5
  cum_cal_NBS_ts_6[[lake]]$forecast_month[cum_cal_NBS_ts_6[[lake]]$forecast_month < 1] <- cum_cal_NBS_ts_6[[lake]]$forecast_month[cum_cal_NBS_ts_6[[lake]]$forecast_month < 1]+12
  cum_cal_NBS_ts_6[[lake]]$forecast_month[c(1:5)] <- NA
  rep_count <- 6
  exit_flag=F
  for (year in 1:56){
    for (month in 1:12){
      if(rep_count>672){exit_flag=T}
      if(exit_flag){break}
      value <- cum_cal_NBS[[month]][[lake]][year,6]
      cum_cal_NBS_ts_6[[lake]]$cum_NBS[rep_count] <- value
      rep_count <- rep_count + 1
    }
    if(exit_flag){break}
  }
}
for (lake in 1:4){
  cum_cal_NBS_ts_12[[lake]] <- data.frame(year=rep(seq(1951,2006), each=12),
                                          month=rep(1:12, times=56),
                                          forecast_month = rep(NA, times = 672),
                                          cum_NBS = rep(NA, 672))
  cum_cal_NBS_ts_12[[lake]]$forecast_month <- cum_cal_NBS_ts_12[[lake]]$month-11
  cum_cal_NBS_ts_12[[lake]]$forecast_month[cum_cal_NBS_ts_12[[lake]]$forecast_month < 1] <- cum_cal_NBS_ts_12[[lake]]$forecast_month[cum_cal_NBS_ts_12[[lake]]$forecast_month < 1]+12
  cum_cal_NBS_ts_12[[lake]]$forecast_month[c(1:11)] <- NA
  rep_count <- 12
  exit_flag=F
  for (year in 1:56){
    for (month in 1:12){
      if(rep_count>672){exit_flag=T}
      if(exit_flag){break}
      value <- cum_cal_NBS[[month]][[lake]][year,12]
      cum_cal_NBS_ts_12[[lake]]$cum_NBS[rep_count] <- value
      rep_count <- rep_count + 1
    }
    if(exit_flag){break}
  }
}
######################

######################
# Calculate validation NBS sequences
######################
val_NBS        <- vector(mode="list",length=12)
cum_val_NBS    <- vector(mode="list",length=12)
names(val_NBS) <- month_names
names(cum_val_NBS)    <- month_names
for (month in 1:12){
  val_NBS[[month]]     <- forecast_complete_NBS(val_inps[[month]],forecast_month=month,sample_num=length(val_inps[[1]][,1]))
  cum_val_NBS[[month]] <- cumulative(val_NBS[[month]])
}

cum_val_NBS_ts_1  <- vector(mode="list", length=4)
cum_val_NBS_ts_3  <- vector(mode="list", length=4)
cum_val_NBS_ts_6  <- vector(mode="list", length=4)
cum_val_NBS_ts_12 <- vector(mode="list", length=4)
names(cum_val_NBS_ts_1)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_val_NBS_ts_3)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_val_NBS_ts_6)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_val_NBS_ts_12) <- c("Eri", "Ont", "MiH", "Sup")
for (lake in 1:4){
  cum_val_NBS_ts_1[[lake]] <- data.frame(year=rep(seq(2007,2021), each=12),
                                         month=rep(1:12, times=15),
                                         forecast_month = rep(1:12, times = 15),
                                         cum_NBS = rep(NA, 180))
  rep_count = 1
  for (year in 1:14){
    for (month in 1:12){
      value <- cum_val_NBS[[month]][[lake]][year,1]
      cum_val_NBS_ts_1[[lake]]$cum_NBS[rep_count] <- value
      rep_count = rep_count+1
    }
  }
}
for (lake in 1:4){
  cum_val_NBS_ts_3[[lake]] <- data.frame(year=rep(seq(2007,2021), each=12),
                                         month=rep(1:12, times=15),
                                         forecast_month = rep(NA, times = 180),
                                         cum_NBS = rep(NA, 180))
  cum_val_NBS_ts_3[[lake]]$forecast_month <- cum_val_NBS_ts_3[[lake]]$month-2
  cum_val_NBS_ts_3[[lake]]$forecast_month[cum_val_NBS_ts_3[[lake]]$forecast_month < 1] <- cum_val_NBS_ts_3[[lake]]$forecast_month[cum_val_NBS_ts_3[[lake]]$forecast_month < 1]+12
  cum_val_NBS_ts_3[[lake]]$forecast_month[c(1,2,171:180)] <- NA
  rep_count <- 3
  for (year in 1:14){
    for (month in 1:12){
        value <- cum_val_NBS[[month]][[lake]][year,3]
        cum_val_NBS_ts_3[[lake]]$cum_NBS[rep_count] <- value
        rep_count <- rep_count + 1
    }
  }
}
for (lake in 1:4){
  cum_val_NBS_ts_6[[lake]] <- data.frame(year=rep(seq(2007,2021), each=12),
                                         month=rep(1:12, times=15),
                                         forecast_month = rep(NA, times = 180),
                                         cum_NBS = rep(NA, 180))
  cum_val_NBS_ts_6[[lake]]$forecast_month <- cum_val_NBS_ts_6[[lake]]$month-5
  cum_val_NBS_ts_6[[lake]]$forecast_month[cum_val_NBS_ts_6[[lake]]$forecast_month < 1] <- cum_val_NBS_ts_6[[lake]]$forecast_month[cum_val_NBS_ts_6[[lake]]$forecast_month < 1]+12
  cum_val_NBS_ts_6[[lake]]$forecast_month[c(1:5,174:180)] <- NA
  rep_count <- 6
  for (year in 1:14){
    for (month in 1:12){
      value <- cum_val_NBS[[month]][[lake]][year,6]
      cum_val_NBS_ts_6[[lake]]$cum_NBS[rep_count] <- value
      rep_count <- rep_count + 1
    }
  }
}
for (lake in 1:4){
  cum_val_NBS_ts_12[[lake]] <- data.frame(year=rep(seq(2007,2021), each=12),
                                         month=rep(1:12, times=15),
                                         forecast_month = rep(NA, times = 180),
                                         cum_NBS = rep(NA, 180))
  cum_val_NBS_ts_12[[lake]]$forecast_month <- cum_val_NBS_ts_12[[lake]]$month-11
  cum_val_NBS_ts_12[[lake]]$forecast_month[cum_val_NBS_ts_12[[lake]]$forecast_month < 1] <- cum_val_NBS_ts_12[[lake]]$forecast_month[cum_val_NBS_ts_12[[lake]]$forecast_month < 1]+12
  cum_val_NBS_ts_12[[lake]]$forecast_month[c(1:11)] <- NA
  rep_count <- 12
  for (year in 1:14){
    for (month in 1:12){
      value <- cum_val_NBS[[month]][[lake]][year,12]
      cum_val_NBS_ts_12[[lake]]$cum_NBS[rep_count] <- value
      rep_count <- rep_count + 1
    }
  }
}
###########################

######################
# Generate matched samples for each validation month
######################
# val_samps <- vector(mode="list", length=168)
# count = 1
# for (year in 1:14){
#   for (month in 1:12){
#     print(count)
#     names(val_samps)[count] <- paste("year",year,"month",month,sep="")
#     msamples  <- matching_sample(antecedents = val_ants[[month]][year,],forecast_month = month, sample_num=n)
#     mNBS      <- forecast_complete_NBS(matched_sample = msamples,       forecast_month = month, sample_num=n)
#     val_samps[[count]]  <- cumulative(mNBS)
#     count = count+1
#   }
# }
# 
# saveRDS(val_samps,file="validation_samples_100k_100_2007_2021.rds")
mval_samps <- readRDS(file="validation_samples_100k_100_2007_2021.rds")

# convert to the time series format 
cum_samp_NBS_ts_1  <- vector(mode="list", length=4)
cum_samp_NBS_ts_3  <- vector(mode="list", length=4)
cum_samp_NBS_ts_6  <- vector(mode="list", length=4)
cum_samp_NBS_ts_12 <- vector(mode="list", length=4)
names(cum_samp_NBS_ts_1)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_samp_NBS_ts_3)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_samp_NBS_ts_6)  <- c("Eri", "Ont", "MiH", "Sup")
names(cum_samp_NBS_ts_12) <- c("Eri", "Ont", "MiH", "Sup")
for (lake in 1:4){
  cum_samp_NBS_ts_1[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_1[[lake]]$year,
                                                month          = cum_val_NBS_ts_1[[lake]]$month,
                                                forecast_month = cum_val_NBS_ts_1[[lake]]$forecast_month),
                                      matrix(NA,nrow=180,ncol=n))
  cum_samp_NBS_ts_3[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_3[[lake]]$year,
                                                month          = cum_val_NBS_ts_3[[lake]]$month,
                                                forecast_month = cum_val_NBS_ts_3[[lake]]$forecast_month),
                                     matrix(NA,nrow=180,ncol=n))
  cum_samp_NBS_ts_6[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_6[[lake]]$year,
                                                month          = cum_val_NBS_ts_6[[lake]]$month,
                                                forecast_month = cum_val_NBS_ts_6[[lake]]$forecast_month),
                                     matrix(NA,nrow=180,ncol=n))
  cum_samp_NBS_ts_12[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_12[[lake]]$year,
                                                 month          = cum_val_NBS_ts_12[[lake]]$month,
                                                 forecast_month = cum_val_NBS_ts_12[[lake]]$forecast_month),
                                      matrix(NA,nrow=180,ncol=n))
  rep_count1  <- 1
  rep_count3  <- 3
  rep_count6  <- 6
  rep_count12 <- 12
  for (year in 1:14){
    for (month in 1:12){
      values1 <- mval_samps[[rep_count1]][[lake]][,1]
      cum_samp_NBS_ts_1[[lake]][rep_count1,4:(n+3)] <- values1
      
      values3 <- mval_samps[[rep_count1]][[lake]][,3]
      cum_samp_NBS_ts_3[[lake]][rep_count3,4:(n+3)] <- values3
      
      values6 <- mval_samps[[rep_count1]][[lake]][,6]
      cum_samp_NBS_ts_6[[lake]][rep_count6,4:(n+3)] <- values6
      
      values12 <- mval_samps[[rep_count1]][[lake]][,12]
      cum_samp_NBS_ts_12[[lake]][rep_count12,4:(n+3)] <- values12
      
      rep_count1  <- rep_count1  + 1
      rep_count3  <- rep_count3  + 1
      rep_count6  <- rep_count6  + 1
      rep_count12 <- rep_count12 + 1
    }
  }
}
#####################

######################
# Generate Random Samples for each validation month
######################
# rval_samps <- vector(mode="list", length=168)
# count = 1
# for (year in 1:14){
#   for (month in 1:12){
#     print(count)
#     names(rval_samps)[count] <- paste("year",year,"month",month,sep="")
#     rsamples  <- random_sample(forecast_month = month, sample_num=n)
#     rNBS      <- forecast_complete_NBS(matched_sample = rsamples,       forecast_month = month, sample_num=n)
#     rval_samps[[count]]  <- cumulative(rNBS)
#     count = count+1
#   }
# }
# 
# saveRDS(rval_samps,file="random_validation_samples_100k_100_2007_2021.rds")
rval_samps <- readRDS(file="random_validation_samples_100k_100_2007_2021.rds")

# convert to the time series format 
rcum_samp_NBS_ts_1  <- vector(mode="list", length=4)
rcum_samp_NBS_ts_3  <- vector(mode="list", length=4)
rcum_samp_NBS_ts_6  <- vector(mode="list", length=4)
rcum_samp_NBS_ts_12 <- vector(mode="list", length=4)
names(rcum_samp_NBS_ts_1)  <- c("Eri", "Ont", "MiH", "Sup")
names(rcum_samp_NBS_ts_3)  <- c("Eri", "Ont", "MiH", "Sup")
names(rcum_samp_NBS_ts_6)  <- c("Eri", "Ont", "MiH", "Sup")
names(rcum_samp_NBS_ts_12) <- c("Eri", "Ont", "MiH", "Sup")
for (lake in 1:4){
  rcum_samp_NBS_ts_1[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_1[[lake]]$year,
                                                month           = cum_val_NBS_ts_1[[lake]]$month,
                                                forecast_month  = cum_val_NBS_ts_1[[lake]]$forecast_month),
                                     matrix(NA,nrow=180,ncol=n))
  rcum_samp_NBS_ts_3[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_3[[lake]]$year,
                                                month           = cum_val_NBS_ts_3[[lake]]$month,
                                                forecast_month  = cum_val_NBS_ts_3[[lake]]$forecast_month),
                                     matrix(NA,nrow=180,ncol=n))
  rcum_samp_NBS_ts_6[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_6[[lake]]$year,
                                                month           = cum_val_NBS_ts_6[[lake]]$month,
                                                forecast_month  = cum_val_NBS_ts_6[[lake]]$forecast_month),
                                     matrix(NA,nrow=180,ncol=n))
  rcum_samp_NBS_ts_12[[lake]] <- cbind(data.frame(year           = cum_val_NBS_ts_12[[lake]]$year,
                                                 month           = cum_val_NBS_ts_12[[lake]]$month,
                                                 forecast_month  = cum_val_NBS_ts_12[[lake]]$forecast_month),
                                      matrix(NA,nrow=180,ncol=n))
  rep_count1  <- 1
  rep_count3  <- 3
  rep_count6  <- 6
  rep_count12 <- 12
  for (year in 1:14){
    for (month in 1:12){
      values1 <- rval_samps[[rep_count1]][[lake]][,1]
      rcum_samp_NBS_ts_1[[lake]][rep_count1,4:(n+3)] <- values1
      
      values3 <- rval_samps[[rep_count1]][[lake]][,3]
      rcum_samp_NBS_ts_3[[lake]][rep_count3,4:(n+3)] <- values3
      
      values6 <- rval_samps[[rep_count1]][[lake]][,6]
      rcum_samp_NBS_ts_6[[lake]][rep_count6,4:(n+3)] <- values6
      
      values12 <- rval_samps[[rep_count1]][[lake]][,12]
      rcum_samp_NBS_ts_12[[lake]][rep_count12,4:(n+3)] <- values12
      
      rep_count1  <- rep_count1  + 1
      rep_count3  <- rep_count3  + 1
      rep_count6  <- rep_count6  + 1
      rep_count12 <- rep_count12 + 1
    }
  }
}
#####################

# Generate table of fraction of observations within the 90% prediction interval

# Lake Erie
Eri1_90s      <- apply(as.matrix(cum_samp_NBS_ts_1$Eri[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri1_90s     <- apply(as.matrix(rcum_samp_NBS_ts_1$Eri[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Eri1_90_frac  <- round(mean((cum_val_NBS_ts_1$Eri$cum_NBS[1:168]>Eri1_90s[1,])&(cum_val_NBS_ts_1$Eri$cum_NBS[1:168]<Eri1_90s[2,])),3)
rEri1_90_frac <- round(mean((cum_val_NBS_ts_1$Eri$cum_NBS[1:168]>rEri1_90s[1,])&(cum_val_NBS_ts_1$Eri$cum_NBS[1:168]<rEri1_90s[2,])),3)

Eri3_90s      <- apply(as.matrix(cum_samp_NBS_ts_3$Eri[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri3_90s     <- apply(as.matrix(rcum_samp_NBS_ts_3$Eri[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Eri3_90_frac  <- round(mean((cum_val_NBS_ts_3$Eri$cum_NBS[3:170]>Eri3_90s[1,])&(cum_val_NBS_ts_3$Eri$cum_NBS[3:170]<Eri3_90s[2,])),3)
rEri3_90_frac <- round(mean((cum_val_NBS_ts_3$Eri$cum_NBS[3:170]>rEri3_90s[1,])&(cum_val_NBS_ts_3$Eri$cum_NBS[3:170]<rEri3_90s[2,])),3)

Eri6_90s      <- apply(as.matrix(cum_samp_NBS_ts_6$Eri[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri6_90s     <- apply(as.matrix(rcum_samp_NBS_ts_6$Eri[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Eri6_90_frac  <- round(mean((cum_val_NBS_ts_6$Eri$cum_NBS[6:173]>Eri6_90s[1,])&(cum_val_NBS_ts_6$Eri$cum_NBS[6:173]<Eri6_90s[2,])),3)
rEri6_90_frac <- round(mean((cum_val_NBS_ts_6$Eri$cum_NBS[6:173]>rEri6_90s[1,])&(cum_val_NBS_ts_6$Eri$cum_NBS[6:173]<rEri6_90s[2,])),3)

Eri12_90s      <- apply(as.matrix(cum_samp_NBS_ts_12$Eri[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri12_90s     <- apply(as.matrix(rcum_samp_NBS_ts_12$Eri[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Eri12_90_frac  <- round(mean((cum_val_NBS_ts_12$Eri$cum_NBS[12:179]>Eri12_90s[1,])&(cum_val_NBS_ts_12$Eri$cum_NBS[12:179]<Eri12_90s[2,])),3)
rEri12_90_frac <- round(mean((cum_val_NBS_ts_12$Eri$cum_NBS[12:179]>rEri12_90s[1,])&(cum_val_NBS_ts_12$Eri$cum_NBS[12:179]<rEri12_90s[2,])),3)

# Lake Ontario
Ont1_90s      <- apply(as.matrix(cum_samp_NBS_ts_1$Ont[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt1_90s     <- apply(as.matrix(rcum_samp_NBS_ts_1$Ont[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Ont1_90_frac  <- round(mean((cum_val_NBS_ts_1$Ont$cum_NBS[1:168]>Ont1_90s[1,])&(cum_val_NBS_ts_1$Ont$cum_NBS[1:168]<Ont1_90s[2,])),3)
rOnt1_90_frac <- round(mean((cum_val_NBS_ts_1$Ont$cum_NBS[1:168]>rOnt1_90s[1,])&(cum_val_NBS_ts_1$Ont$cum_NBS[1:168]<rOnt1_90s[2,])),3)

Ont3_90s      <- apply(as.matrix(cum_samp_NBS_ts_3$Ont[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt3_90s     <- apply(as.matrix(rcum_samp_NBS_ts_3$Ont[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Ont3_90_frac  <- round(mean((cum_val_NBS_ts_3$Ont$cum_NBS[3:170]>Ont3_90s[1,])&(cum_val_NBS_ts_3$Ont$cum_NBS[3:170]<Ont3_90s[2,])),3)
rOnt3_90_frac <- round(mean((cum_val_NBS_ts_3$Ont$cum_NBS[3:170]>rOnt3_90s[1,])&(cum_val_NBS_ts_3$Ont$cum_NBS[3:170]<rOnt3_90s[2,])),3)

Ont6_90s      <- apply(as.matrix(cum_samp_NBS_ts_6$Ont[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt6_90s     <- apply(as.matrix(rcum_samp_NBS_ts_6$Ont[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Ont6_90_frac  <- round(mean((cum_val_NBS_ts_6$Ont$cum_NBS[6:173]>Ont6_90s[1,])&(cum_val_NBS_ts_6$Ont$cum_NBS[6:173]<Ont6_90s[2,])),3)
rOnt6_90_frac <- round(mean((cum_val_NBS_ts_6$Ont$cum_NBS[6:173]>rOnt6_90s[1,])&(cum_val_NBS_ts_6$Ont$cum_NBS[6:173]<rOnt6_90s[2,])),3)

Ont12_90s      <- apply(as.matrix(cum_samp_NBS_ts_12$Ont[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt12_90s     <- apply(as.matrix(rcum_samp_NBS_ts_12$Ont[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Ont12_90_frac  <- round(mean((cum_val_NBS_ts_12$Ont$cum_NBS[12:179]>Ont12_90s[1,])&(cum_val_NBS_ts_12$Ont$cum_NBS[12:179]<Ont12_90s[2,])),3)
rOnt12_90_frac <- round(mean((cum_val_NBS_ts_12$Ont$cum_NBS[12:179]>rOnt12_90s[1,])&(cum_val_NBS_ts_12$Ont$cum_NBS[12:179]<rOnt12_90s[2,])),3)

# Lake Michigan-Huron
MiH1_90s      <- apply(as.matrix(cum_samp_NBS_ts_1$MiH[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH1_90s     <- apply(as.matrix(rcum_samp_NBS_ts_1$MiH[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
MiH1_90_frac  <- round(mean((cum_val_NBS_ts_1$MiH$cum_NBS[1:168]>MiH1_90s[1,])&(cum_val_NBS_ts_1$MiH$cum_NBS[1:168]<MiH1_90s[2,])),3)
rMiH1_90_frac <- round(mean((cum_val_NBS_ts_1$MiH$cum_NBS[1:168]>rMiH1_90s[1,])&(cum_val_NBS_ts_1$MiH$cum_NBS[1:168]<rMiH1_90s[2,])),3)

MiH3_90s      <- apply(as.matrix(cum_samp_NBS_ts_3$MiH[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH3_90s     <- apply(as.matrix(rcum_samp_NBS_ts_3$MiH[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
MiH3_90_frac  <- round(mean((cum_val_NBS_ts_3$MiH$cum_NBS[3:170]>MiH3_90s[1,])&(cum_val_NBS_ts_3$MiH$cum_NBS[3:170]<MiH3_90s[2,])),3)
rMiH3_90_frac <- round(mean((cum_val_NBS_ts_3$MiH$cum_NBS[3:170]>rMiH3_90s[1,])&(cum_val_NBS_ts_3$MiH$cum_NBS[3:170]<rMiH3_90s[2,])),3)

MiH6_90s      <- apply(as.matrix(cum_samp_NBS_ts_6$MiH[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH6_90s     <- apply(as.matrix(rcum_samp_NBS_ts_6$MiH[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
MiH6_90_frac  <- round(mean((cum_val_NBS_ts_6$MiH$cum_NBS[6:173]>MiH6_90s[1,])&(cum_val_NBS_ts_6$MiH$cum_NBS[6:173]<MiH6_90s[2,])),3)
rMiH6_90_frac <- round(mean((cum_val_NBS_ts_6$MiH$cum_NBS[6:173]>rMiH6_90s[1,])&(cum_val_NBS_ts_6$MiH$cum_NBS[6:173]<rMiH6_90s[2,])),3)

MiH12_90s      <- apply(as.matrix(cum_samp_NBS_ts_12$MiH[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH12_90s     <- apply(as.matrix(rcum_samp_NBS_ts_12$MiH[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
MiH12_90_frac  <- round(mean((cum_val_NBS_ts_12$MiH$cum_NBS[12:179]>MiH12_90s[1,])&(cum_val_NBS_ts_12$MiH$cum_NBS[12:179]<MiH12_90s[2,])),3)
rMiH12_90_frac <- round(mean((cum_val_NBS_ts_12$MiH$cum_NBS[12:179]>rMiH12_90s[1,])&(cum_val_NBS_ts_12$MiH$cum_NBS[12:179]<rMiH12_90s[2,])),3)

# Lake Superior
Sup1_90s      <- apply(as.matrix(cum_samp_NBS_ts_1$Sup[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup1_90s     <- apply(as.matrix(rcum_samp_NBS_ts_1$Sup[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Sup1_90_frac  <- round(mean((cum_val_NBS_ts_1$Sup$cum_NBS[1:168]>Sup1_90s[1,])&(cum_val_NBS_ts_1$Sup$cum_NBS[1:168]<Sup1_90s[2,])),3)
rSup1_90_frac <- round(mean((cum_val_NBS_ts_1$Sup$cum_NBS[1:168]>rSup1_90s[1,])&(cum_val_NBS_ts_1$Sup$cum_NBS[1:168]<rSup1_90s[2,])),3)

Sup3_90s      <- apply(as.matrix(cum_samp_NBS_ts_3$Sup[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup3_90s     <- apply(as.matrix(rcum_samp_NBS_ts_3$Sup[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Sup3_90_frac  <- round(mean((cum_val_NBS_ts_3$Sup$cum_NBS[3:170]>Sup3_90s[1,])&(cum_val_NBS_ts_3$Sup$cum_NBS[3:170]<Sup3_90s[2,])),3)
rSup3_90_frac <- round(mean((cum_val_NBS_ts_3$Sup$cum_NBS[3:170]>rSup3_90s[1,])&(cum_val_NBS_ts_3$Sup$cum_NBS[3:170]<rSup3_90s[2,])),3)

Sup6_90s      <- apply(as.matrix(cum_samp_NBS_ts_6$Sup[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup6_90s     <- apply(as.matrix(rcum_samp_NBS_ts_6$Sup[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Sup6_90_frac  <- round(mean((cum_val_NBS_ts_6$Sup$cum_NBS[6:173]>Sup6_90s[1,])&(cum_val_NBS_ts_6$Sup$cum_NBS[6:173]<Sup6_90s[2,])),3)
rSup6_90_frac <- round(mean((cum_val_NBS_ts_6$Sup$cum_NBS[6:173]>rSup6_90s[1,])&(cum_val_NBS_ts_6$Sup$cum_NBS[6:173]<rSup6_90s[2,])),3)

Sup12_90s      <- apply(as.matrix(cum_samp_NBS_ts_12$Sup[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup12_90s     <- apply(as.matrix(rcum_samp_NBS_ts_12$Sup[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
Sup12_90_frac  <- round(mean((cum_val_NBS_ts_12$Sup$cum_NBS[12:179]>Sup12_90s[1,])&(cum_val_NBS_ts_12$Sup$cum_NBS[12:179]<Sup12_90s[2,])),3)
rSup12_90_frac <- round(mean((cum_val_NBS_ts_12$Sup$cum_NBS[12:179]>rSup12_90s[1,])&(cum_val_NBS_ts_12$Sup$cum_NBS[12:179]<rSup12_90s[2,])),3)

unmatched_90_frac_table <- data.frame(
  "1-month"  = c(rSup1_90_frac,rMiH1_90_frac,rEri1_90_frac,rOnt1_90_frac),
  "3-month"  = c(rSup3_90_frac,rMiH3_90_frac,rEri3_90_frac,rOnt3_90_frac),
  "6-month"  = c(rSup6_90_frac,rMiH6_90_frac,rEri6_90_frac,rOnt6_90_frac),
  "12-month" = c(rSup12_90_frac,rMiH12_90_frac,rEri12_90_frac,rOnt12_90_frac),
  row.names=c("Superior","Michigan-Huron","Erie","Ontario")
)

matched_90_frac_table <- data.frame(
  "1-month"  = c(Sup1_90_frac,MiH1_90_frac,Eri1_90_frac,Ont1_90_frac),
  "3-month"  = c(Sup3_90_frac,MiH3_90_frac,Eri3_90_frac,Ont3_90_frac),
  "6-month"  = c(Sup6_90_frac,MiH6_90_frac,Eri6_90_frac,Ont6_90_frac),
  "12-month" = c(Sup12_90_frac,MiH12_90_frac,Eri12_90_frac,Ont12_90_frac),
  row.names=c("Superior","Michigan-Huron","Erie","Ontario")
)



# Autocorrelation plots
#################
# Load original l2s outputs
# create index of where the primary variable in each plot occurs ex (13,25...)
# calculate correlation between that variable and other values that are +- 1:12 away from it S
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/L2SWBM_med_1950-2022")
cal_Sup_P <- read.csv("SupPrecip_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_Sup_E <- read.csv("SupEvap_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_Sup_R <- read.csv("SupRunoff_analysis19502022_prior19001969_1m.csv")[1:696,1:3]

cal_MiH_P <- read.csv("MiHurPrecip_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_MiH_E <- read.csv("MiHurEvap_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_MiH_R <- read.csv("MiHurRunoff_analysis19502022_prior19001969_1m.csv")[1:696,1:3]

cal_Eri_P <- read.csv("EriePrecip_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_Eri_E <- read.csv("ErieEvap_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_Eri_R <- read.csv("ErieRunoff_analysis19502022_prior19001969_1m.csv")[1:696,1:3]

cal_Ont_P <- read.csv("OntPrecip_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_Ont_E <- read.csv("OntEvap_analysis19502022_prior19001969_1m.csv")[1:696,1:3]
cal_Ont_R <- read.csv("OntRunoff_analysis19502022_prior19001969_1m.csv")[1:696,1:3]

Jan_rows  <- seq(13,673,by=12) 
Apr_rows  <- seq(16,676,by=12)
Jul_rows  <- seq(19,679,by=12)
Oct_rows  <- seq(22,682,by=12)

Month_rows <- list("Jan" = Jan_rows, "Apr" = Apr_rows, "Jul" = Jul_rows, "Oct" = Oct_rows)

obsCors = matrix(ncol=12,nrow=25)
corRows           <- c("-12month","-11month","-10month","-09month","-08month",
                      "-07month","-06month","-05month","-04month","-03month",
                      "-02month","-01month","+00month","+01month","+02month",
                      "+03month","+04month","+05month","+06month","+07month",
                      "+08month","+09month","+10month","+11month","+12month")
rownames(obsCors) <- corRows
colnames(obsCors) <- c("cal_Sup_P", "cal_Sup_E", "cal_Sup_R", "cal_MiH_P", "cal_MiH_E", "cal_MiH_R", 
                       "cal_Eri_P", "cal_Eri_E", "cal_Eri_R", "cal_Ont_P", "cal_Ont_E", "cal_Ont_R")
obsCors[13,] <- 1
obsCorsList <- list("Jan" = obsCors, "Apr" = obsCors, "Jul" = obsCors, "Oct" = obsCors)

for (rows in 1:4){
  for (lag in 1:12){
    obsCorsList[[rows]][lag,1] = cor(x=cal_Sup_P$Median[Month_rows[[rows]]], y=cal_Sup_P$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,2] = cor(x=cal_Sup_E$Median[Month_rows[[rows]]], y=cal_Sup_E$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,3] = cor(x=cal_Sup_R$Median[Month_rows[[rows]]], y=cal_Sup_R$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    
    obsCorsList[[rows]][lag,4] = cor(x=cal_MiH_P$Median[Month_rows[[rows]]], y=cal_MiH_P$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,5] = cor(x=cal_MiH_E$Median[Month_rows[[rows]]], y=cal_MiH_E$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,6] = cor(x=cal_MiH_R$Median[Month_rows[[rows]]], y=cal_MiH_R$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    
    obsCorsList[[rows]][lag,7] = cor(x=cal_Eri_P$Median[Month_rows[[rows]]], y=cal_Eri_P$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,8] = cor(x=cal_Eri_E$Median[Month_rows[[rows]]], y=cal_Eri_E$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,9] = cor(x=cal_Eri_R$Median[Month_rows[[rows]]], y=cal_Eri_R$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    
    obsCorsList[[rows]][lag,10] = cor(x=cal_Ont_P$Median[Month_rows[[rows]]], y=cal_Ont_P$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,11] = cor(x=cal_Ont_E$Median[Month_rows[[rows]]], y=cal_Ont_E$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
    obsCorsList[[rows]][lag,12] = cor(x=cal_Ont_R$Median[Month_rows[[rows]]], y=cal_Ont_R$Median[(Month_rows[[rows]]-(12-lag+1))], method="s")
  }
  for (look in 1:12){
    obsCorsList[[rows]][(look+13),1] = cor(x=cal_Sup_P$Median[Month_rows[[rows]]], y=cal_Sup_P$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),2] = cor(x=cal_Sup_E$Median[Month_rows[[rows]]], y=cal_Sup_E$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),3] = cor(x=cal_Sup_R$Median[Month_rows[[rows]]], y=cal_Sup_R$Median[(Month_rows[[rows]]+look)], method="s")
    
    obsCorsList[[rows]][(look+13),4] = cor(x=cal_MiH_P$Median[Month_rows[[rows]]], y=cal_MiH_P$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),5] = cor(x=cal_MiH_E$Median[Month_rows[[rows]]], y=cal_MiH_E$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),6] = cor(x=cal_MiH_R$Median[Month_rows[[rows]]], y=cal_MiH_R$Median[(Month_rows[[rows]]+look)], method="s")
    
    obsCorsList[[rows]][(look+13),7] = cor(x=cal_Eri_P$Median[Month_rows[[rows]]], y=cal_Eri_P$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),8] = cor(x=cal_Eri_E$Median[Month_rows[[rows]]], y=cal_Eri_E$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),9] = cor(x=cal_Eri_R$Median[Month_rows[[rows]]], y=cal_Eri_R$Median[(Month_rows[[rows]]+look)], method="s")
    
    obsCorsList[[rows]][(look+13),10] = cor(x=cal_Ont_P$Median[Month_rows[[rows]]], y=cal_Ont_P$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),11] = cor(x=cal_Ont_E$Median[Month_rows[[rows]]], y=cal_Ont_E$Median[(Month_rows[[rows]]+look)], method="s")
    obsCorsList[[rows]][(look+13),12] = cor(x=cal_Ont_R$Median[Month_rows[[rows]]], y=cal_Ont_R$Median[(Month_rows[[rows]]+look)], method="s")
  }
}


### simulated autocorrelation
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
JanSamps <- readRDS("Jan_15month_1950_2006_samp.rds")
colnames(JanSamps) = colnames(val_inps[["Jan"]])
AprSamps <- readRDS("Apr_15month_1950_2006_samp.rds")
colnames(AprSamps) = colnames(val_inps[["Apr"]])
JulSamps <- readRDS("Jul_15month_1950_2006_samp.rds")
colnames(JulSamps) = colnames(val_inps[["Jul"]])
OctSamps <- readRDS("Oct_15month_1950_2006_samp.rds")
colnames(OctSamps) = colnames(val_inps[["Oct"]])

seasonal_samps = list("Jan" = JanSamps, "Apr" = AprSamps, "Jul" = JulSamps, "Oct" = OctSamps)


simCors = matrix(ncol=12,nrow=15)
corRows           <- c("-03month",
                       "-02month","-01month","+00month","+01month","+02month",
                       "+03month","+04month","+05month","+06month","+07month",
                       "+08month","+09month","+10month","+11month")
rownames(simCors) <- corRows
colnames(simCors) <- c("cal_Sup_P", "cal_Sup_E", "cal_Sup_R", "cal_MiH_P", "cal_MiH_E", "cal_MiH_R", 
                       "cal_Eri_P", "cal_Eri_E", "cal_Eri_R", "cal_Ont_P", "cal_Ont_E", "cal_Ont_R")
simCors[4,] <- 1
simCorsList <- list("Jan" = simCors, "Apr" = simCors, "Jul" = simCors, "Oct" = simCors)

for (season in 1:4){
  for (fmonth in 1:15){
    simCorsList[[season]][fmonth,1] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Sup.P", colnames(seasonal_samps[[season]]))][4]],
                                         y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Sup.P", colnames(seasonal_samps[[season]]))][fmonth]],
                                         method="s")
    simCorsList[[season]][fmonth,2] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Sup.E", colnames(seasonal_samps[[season]]))][4]],
                                         y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Sup.E", colnames(seasonal_samps[[season]]))][fmonth]],
                                         method="s")
    simCorsList[[season]][fmonth,3] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Sup.R", colnames(seasonal_samps[[season]]))][4]],
                                         y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Sup.R", colnames(seasonal_samps[[season]]))][fmonth]],
                                         method="s")
    
    simCorsList[[season]][fmonth,4] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("MiH.P", colnames(seasonal_samps[[season]]))][4]],
                                          y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("MiH.P", colnames(seasonal_samps[[season]]))][fmonth]],
                                          method="s")
    simCorsList[[season]][fmonth,5] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("MiH.E", colnames(seasonal_samps[[season]]))][4]],
                                          y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("MiH.E", colnames(seasonal_samps[[season]]))][fmonth]],
                                          method="s")
    simCorsList[[season]][fmonth,6] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("MiH.R", colnames(seasonal_samps[[season]]))][4]],
                                          y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("MiH.R", colnames(seasonal_samps[[season]]))][fmonth]],
                                          method="s")
    
    simCorsList[[season]][fmonth,7] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Eri.P", colnames(seasonal_samps[[season]]))][4]],
                                          y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Eri.P", colnames(seasonal_samps[[season]]))][fmonth]],
                                          method="s")
    simCorsList[[season]][fmonth,8] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Eri.E", colnames(seasonal_samps[[season]]))][4]],
                                          y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Eri.E", colnames(seasonal_samps[[season]]))][fmonth]],
                                          method="s")
    simCorsList[[season]][fmonth,9] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Eri.R", colnames(seasonal_samps[[season]]))][4]],
                                          y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Eri.R", colnames(seasonal_samps[[season]]))][fmonth]],
                                          method="s")
    
    simCorsList[[season]][fmonth,10] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Ont.P", colnames(seasonal_samps[[season]]))][4]],
                                         y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Ont.P", colnames(seasonal_samps[[season]]))][fmonth]],
                                         method="s")
    simCorsList[[season]][fmonth,11] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Ont.E", colnames(seasonal_samps[[season]]))][4]],
                                         y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Ont.E", colnames(seasonal_samps[[season]]))][fmonth]],
                                         method="s")
    simCorsList[[season]][fmonth,12] = cor(x=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Ont.R", colnames(seasonal_samps[[season]]))][4]],
                                         y=seasonal_samps[[season]][,colnames(seasonal_samps[[season]])[grep("Ont.R", colnames(seasonal_samps[[season]]))][fmonth]],
                                         method="s")
  }
}


### generate individual plots with simulation values
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/Figures")
acfmonths = c("Jan","Apr","Jul","Oct")

png(filename = "SupSimACF.png", height=800,width=800)
par(mfrow=c(4,3), mar=c(1,1,1,1), oma=c(4,4,4,1.5))
for (month in 1:4){
  for (comp in 1:3){
    plot(c(-12:12),obsCorsList[[month]][,comp],type="n", xlim=c(-12,12), ylim=c(-.5,1),pch=19,
         cex=0.4, axes=F,xlab="",ylab="")
    rect(xleft=-3,ybottom=-0.5,xright=11,ytop=1.0,col="lightgray", border=NA)
    abline(h=0, cex=0.3,col="darkgray")
    for (j in 1:25){
      segments(x0=j-13,x1=j-13,y0=0,y1=obsCorsList[[month]][j,comp], col="black", lwd=1)
    }    
    points(x=seq(-12,12), y=obsCorsList[[month]][,comp], pch=19, cex=1)
    points(x=seq(-3,11), y=simCorsList[[month]][,comp], pch=19, cex=1, col="red")
    axis(2, label=F, tick=T, tck=-0.04)
    axis(1, label=F, tick=T, tck=-0.05, cex.lab=0.5, mgp=c(0.3,0,0), cex=0.8)
    box()
    
    if (comp==1){
      mtext(text=acfmonths[month], side=2, line=2.5, cex=0.8)
      axis(2,label=T,tick=F)
    }
    
    if (comp==1 & month==1){mtext(text="Precipitation",side=3,line=1.5,cex=0.8)}
    if (comp==2 & month==1){
      mtext(text="Evaporation",side=3,line=1.5,cex=0.8)
      mtext(text="Superior", side=3, line=3, cex=1.2)}
    if (comp==3 & month==1){mtext(text="Runoff",side=3,line=1.5,cex=0.8)}
    
    if (month==4){
      mtext(text="Lag [months]", side=1, line=2.5, cex=0.8)
      axis(1,label=T,tick=T)
    }
  }
}
dev.off()


png(filename = "MiHSimACF.png", height=800,width=800)
par(mfrow=c(4,3), mar=c(1,1,1,1), oma=c(4,4,4,1.5))
for (month in 1:4){
  for (comp in 1:3){
    plot(c(-12:12),obsCorsList[[month]][,comp],type="n", xlim=c(-12,12), ylim=c(-.5,1),pch=19,
         cex=0.4, axes=F,xlab="",ylab="")
    rect(xleft=-3,ybottom=-0.5,xright=11,ytop=1.0,col="lightgray", border=NA)
    abline(h=0, cex=0.3,col="darkgray")
    for (j in 1:25){
      segments(x0=j-13,x1=j-13,y0=0,y1=obsCorsList[[month]][j,comp+3], col="black", lwd=1)
    }    
    points(x=seq(-12,12), y=obsCorsList[[month]][,comp+3], pch=19, cex=1)
    points(x=seq(-3,11), y=simCorsList[[month]][,comp+3], pch=19, cex=1, col="red")
    axis(2, label=F, tick=T, tck=-0.04)
    axis(1, label=F, tick=T, tck=-0.05, cex.lab=0.5, mgp=c(0.3,0,0), cex=0.8)
    box()
    
    if (comp==1){
      mtext(text=acfmonths[month], side=2, line=2.5, cex=0.8)
      axis(2,label=T,tick=F)
    }
    
    if (comp==1 & month==1){mtext(text="Precipitation",side=3,line=1.5,cex=0.8)}
    if (comp==2 & month==1){
      mtext(text="Evaporation",side=3,line=1.5,cex=0.8)
      mtext(text="Michigan-Huron", side=3, line=3, cex=1.2)}
    if (comp==3 & month==1){mtext(text="Runoff",side=3,line=1.5,cex=0.8)}
    
    if (month==4){
      mtext(text="Lag [months]", side=1, line=2.5, cex=0.8)
      axis(1,label=T,tick=T)
    }
  }
}
dev.off()


png(filename = "EriSimACF.png", height=800,width=800)
par(mfrow=c(4,3), mar=c(1,1,1,1), oma=c(4,4,4,1.5))
for (month in 1:4){
  for (comp in 1:3){
    plot(c(-12:12),obsCorsList[[month]][,comp],type="n", xlim=c(-12,12), ylim=c(-.5,1),pch=19,
         cex=0.4, axes=F,xlab="",ylab="")
    rect(xleft=-3,ybottom=-0.5,xright=11,ytop=1.0,col="lightgray", border=NA)
    abline(h=0, cex=0.3,col="darkgray")
    for (j in 1:25){
      segments(x0=j-13,x1=j-13,y0=0,y1=obsCorsList[[month]][j,comp+6], col="black", lwd=1)
    }    
    points(x=seq(-12,12), y=obsCorsList[[month]][,comp+6], pch=19, cex=1)
    points(x=seq(-3,11), y=simCorsList[[month]][,comp+6], pch=19, cex=1, col="red")
    axis(2, label=F, tick=T, tck=-0.04)
    axis(1, label=F, tick=T, tck=-0.05, cex.lab=0.5, mgp=c(0.3,0,0), cex=0.8)
    box()
    
    if (comp==1){
      mtext(text=acfmonths[month], side=2, line=2.5, cex=0.8)
      axis(2,label=T,tick=F)
    }
    
    if (comp==1 & month==1){mtext(text="Precipitation",side=3,line=1.5,cex=0.8)}
    if (comp==2 & month==1){
      mtext(text="Evaporation",side=3,line=1.5,cex=0.8)
      mtext(text="Erie", side=3, line=3, cex=1.2)}
    if (comp==3 & month==1){mtext(text="Runoff",side=3,line=1.5,cex=0.8)}
    
    if (month==4){
      mtext(text="Lag [months]", side=1, line=2.5, cex=0.8)
      axis(1,label=T,tick=T)
    }
  }
}
dev.off()



png(filename = "OntSimACF.png", height=800,width=800)
par(mfrow=c(4,3), mar=c(1,1,1,1), oma=c(4,4,4,1.5))
for (month in 1:4){
  for (comp in 1:3){
    plot(c(-12:12),obsCorsList[[month]][,(comp+9)],type="n", xlim=c(-12,12), ylim=c(-.5,1),pch=19,
         cex=0.4, axes=F,xlab="",ylab="")
    rect(xleft=-3,ybottom=-0.5,xright=11,ytop=1.0,col="lightgray", border=NA)
    abline(h=0, cex=0.3,col="darkgray")  
    for (j in 1:25){
      segments(x0=j-13,x1=j-13,y0=0,y1=obsCorsList[[month]][j,(comp+9)], col="black", lwd=1)
    }
    points(x=seq(-12,12), y=obsCorsList[[month]][,(comp+9)], pch=19, cex=1)
    points(x=seq(-3,11), y=simCorsList[[month]][,comp+9], pch=19, cex=1, col="red")
    axis(2, label=F, tick=T, tck=-0.04)
    axis(1, label=F, tick=T, tck=-0.05, cex.lab=0.5, mgp=c(0.3,0,0), cex=0.8)
    box()
    
    if (comp==1){
      mtext(text=acfmonths[month], side=2, line=2.5, cex=0.8)
      axis(2,label=T,tick=F)
    }
    
    if (comp==1 & month==1){mtext(text="Precipitation",side=3,line=1.5,cex=0.8)}
    if (comp==2 & month==1){
      mtext(text="Evaporation",side=3,line=1.5,cex=0.8)
      mtext(text="Ontario", side=3, line=3, cex=1.2)}
    if (comp==3 & month==1){mtext(text="Runoff",side=3,line=1.5,cex=0.8)}
    
    if (month==4){
      mtext(text="Lag [months]", side=1, line=2.5, cex=0.8)
      axis(1,label=T,tick=T)
    }
  }
}
dev.off()

# combine the Sup and Ont Plotting arrays
SupSimACF <- readPNG("SupSimACF.png")  #
OntSimACF <- readPNG("OntSimACF.png")  

# Save the combined plot as a new PNG file
png("combinedSim_ACF.png", width = 2400, height = 1200)  #
grid.arrange(rasterGrob(SupSimACF), rasterGrob(OntSimACF), ncol = 2)
dev.off()



#################

# NBS Spatial correlation plots
##################
lakes=c("Sup", "MiH", "Eri", "Ont")
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/Figures")

png("ErieOntarioSpatialCorr.png", width=600, height=600)
par(oma=c(2.7,2.7,0,0))
x1=cum_cal_NBS_ts_1$Eri$cum_NBS
y1=cum_cal_NBS_ts_1$Ont$cum_NBS
x2=na.omit(data.matrix(rcum_samp_NBS_ts_1$Eri)[,4:7])
y2=na.omit(data.matrix(rcum_samp_NBS_ts_1$Ont)[,4:7])
zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
x1hist = hist(x1, plot=FALSE, breaks=seq(-250,460,length.out=15))
y1hist = hist(y1, plot=FALSE, breaks=seq(-120,610,length.out=15))
x2hist = hist(x2, plot=FALSE, breaks=seq(-250,460,length.out=15))
y2hist = hist(y2, plot=FALSE, breaks=seq(-120,610,length.out=15))
top = max(c(x1hist$counts, y1hist$counts, x2hist$counts, y2hist$counts))
par(mar=c(3,3,1,1))
plot(x1,y1, pch=19, col="gray", cex=1.2, xlim=c(-250,460),ylim=c(-120,610))
points(x=x2, y=y2, pch=19, col=rgb(1,0.4,0.4))
mtext("Lake Erie Monthly Water Supply [mm]", side=1, outer=F, line=3.5, font=2)
mtext("Lake Ontario Monthly Water Supply [mm]", side=2, outer=F, line=3.5, font=2)
par(mar=c(0,3,1,1))
barplot(x1hist$counts, axes=FALSE, ylim=c(0, top), space=0, col="gray")
par(new=TRUE)
barplot(x2hist$counts, axes=FALSE, ylim=c(0, top), space=0, col=rgb(1,0.4,0.4, alpha=0.5))
par(mar=c(3,0,1,1))
barplot(y1hist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col="gray")
par(new=TRUE)
barplot(y2hist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col=rgb(1,0.4,0.4, alpha=0.5))
dev.off()


png("ErieSuperiorSpatialCorr.png", width=600, height=600)
par(oma=c(2.7,2.7,0,0))
x1=cum_cal_NBS_ts_1$Eri$cum_NBS
y1=cum_cal_NBS_ts_1$Sup$cum_NBS
x2=na.omit(data.matrix(rcum_samp_NBS_ts_1$Eri)[,4:7])
y2=na.omit(data.matrix(rcum_samp_NBS_ts_1$Sup)[,4:7])
zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
x1hist = hist(x1, plot=FALSE, breaks=seq(-250,460,length.out=15))
y1hist = hist(y1, plot=FALSE, breaks=seq(-110,330,length.out=15))
x2hist = hist(x2, plot=FALSE, breaks=seq(-250,460,length.out=15))
y2hist = hist(y2, plot=FALSE, breaks=seq(-110,330,length.out=15))
top = max(c(x1hist$counts, y1hist$counts, x2hist$counts, y2hist$counts))
par(mar=c(3,3,1,1))
plot(x1,y1, pch=19, col="gray", cex=1.2, xlim=c(-250,460),ylim=c(-110,330))
points(x=x2, y=y2, pch=19, col=rgb(1,0.4,0.4))
mtext("Lake Erie Monthly Water Supply [mm]", side=1, outer=F, line=3.5, font=2)
mtext("Lake Superior Monthly Water Supply [mm]", side=2, outer=F, line=3.5, font=2)
par(mar=c(0,3,1,1))
barplot(x1hist$counts, axes=FALSE, ylim=c(0, top), space=0, col="gray")
par(new=TRUE)
barplot(x2hist$counts, axes=FALSE, ylim=c(0, top), space=0, col=rgb(1,0.4,0.4, alpha=0.5))
par(mar=c(3,0,1,1))
barplot(y1hist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col="gray")
par(new=TRUE)
barplot(y2hist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE, col=rgb(1,0.4,0.4, alpha=0.5))
dev.off()


# combine the Sup and Ont plots
SupErispatial <- readPNG("ErieSuperiorSpatialCorr.png")  
OntErispatial <- readPNG("ErieontarioSpatialCorr.png")  

# Save the combined plot as a new PNG file
png("SpatialCorr.png", width = 1600, height = 800)  #
grid.arrange(rasterGrob(SupErispatial), rasterGrob(OntErispatial), ncol = 2)
dev.off()

# 1-month NBS correlation across all lake pairs (for Supplementary Information)
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/Figures")
png("1MonthSpatialCorrAll.png", width = 8, height = 10, units = "in", res = 300)
par(mfrow = c(3, 2), oma = c(1, 0.5, 0.5, 0), mar = c(4, 1, 1, 0), pty="s", mgp = c(2.3, 1, 0))
# Sup-MiH
plot(x=cum_cal_NBS_ts_1$MiH$cum_NBS, y=cum_cal_NBS_ts_1$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-150,380),ylim=c(-110,330), 
     xlab= "Lake Michigan-Huron Monthly Water Supply [mm]", ylab = "Lake Superior Monthly Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_1$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_1$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

#MiH-Eri
plot(x=cum_cal_NBS_ts_1$MiH$cum_NBS, y=cum_cal_NBS_ts_1$Eri$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-150,380),ylim=c(-250,460), 
     xlab= "Lake Michigan-Huron Monthly Water Supply [mm]", ylab = "Lake Erie Monthly Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_1$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_1$Eri)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

#Sup-Eri
plot(x=cum_cal_NBS_ts_1$Eri$cum_NBS, y=cum_cal_NBS_ts_1$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-250,460),ylim=c(-110,330), 
     xlab= "Lake Erie Monthly Water Supply [mm]", ylab = "Lake Superior Monthly Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_1$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_1$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

#MiH-Ont
plot(x=cum_cal_NBS_ts_1$MiH$cum_NBS, y=cum_cal_NBS_ts_1$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-150,380),ylim=c(-120,610), 
     xlab= "Lake Michigan-Huron Monthly Water Supply [mm]", ylab = "Lake Ontario Monthly Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_1$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_1$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

#Sup-Ont
plot(x=cum_cal_NBS_ts_1$Ont$cum_NBS, y=cum_cal_NBS_ts_1$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-120,610),ylim=c(-110,330), 
     xlab= "Lake Ontario Monthly Water Supply [mm]", ylab = "Lake Superior Monthly Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_1$Ont)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_1$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Eri-Ont
plot(x=cum_cal_NBS_ts_1$Eri$cum_NBS, y=cum_cal_NBS_ts_1$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-250,460),ylim=c(-120,610), 
     xlab= "Lake Erie Monthly Water Supply [mm]", ylab = "Lake Ontario Monthly Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_1$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_1$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)
dev.off()


# 3-month NBS correlation across all lake pairs (for Supplementary Information)  
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/Figures")
png("3MonthSpatialCorrAll.png", width = 8, height = 10, units = "in", res = 300)
par(mfrow = c(3, 2), oma = c(1, 0.5, 0.5, 0), mar = c(4, 1, 1, 0), pty="s", mgp = c(2.3, 1, 0))
# Sup-MiH
plot(x=cum_cal_NBS_ts_3$MiH$cum_NBS, y=cum_cal_NBS_ts_3$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-200,750), ylim=c(-220,680),
     xlab= "Lake Michigan-Huron 3-Month Water Supply [mm]", ylab = "Lake Superior 3-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_3$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_3$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# MiH-Eri
plot(x=cum_cal_NBS_ts_3$MiH$cum_NBS, y=cum_cal_NBS_ts_3$Eri$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-200,750), ylim=c(-400,850),
     xlab= "Lake Michigan-Huron 3-Month Water Supply [mm]", ylab = "Lake Erie 3-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_3$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_3$Eri)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Sup-Eri
plot(x=cum_cal_NBS_ts_3$Eri$cum_NBS, y=cum_cal_NBS_ts_3$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-400,850), ylim=c(-220,680),
     xlab= "Lake Erie 3-Month Water Supply [mm]", ylab = "Lake Superior 3-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_3$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_3$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# MiH-Ont
plot(x=cum_cal_NBS_ts_3$MiH$cum_NBS, y=cum_cal_NBS_ts_3$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-200, 750), ylim=c(-100,1300),
     xlab= "Lake Michigan-Huron 3-Month Water Supply [mm]", ylab = "Lake Ontario 3-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_3$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_3$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Sup-Ont
plot(x=cum_cal_NBS_ts_3$Ont$cum_NBS, y=cum_cal_NBS_ts_3$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-100, 1300), ylim=c(-220, 680),
     xlab= "Lake Ontario 3-Month Water Supply [mm]", ylab = "Lake Superior 3-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_3$Ont)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_3$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Eri-Ont
plot(x=cum_cal_NBS_ts_3$Eri$cum_NBS, y=cum_cal_NBS_ts_3$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-400, 850), ylim=c(-100, 1300),
     xlab= "Lake Erie 3-Month Water Supply [mm]", ylab = "Lake Ontario 3-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_3$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_3$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)
dev.off()

# 6-month NBS correlation across all lake pairs (for Supplementary Information)  
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/Figures")
png("6MonthSpatialCorrAll.png", width = 8, height = 10, units = "in", res = 300)
par(mfrow = c(3, 2), oma = c(1, 0.5, 0.5, 0), mar = c(4, 1, 1, 0), pty="s", mgp = c(2.3, 1, 0))
# Sup-MiH
plot(x=cum_cal_NBS_ts_6$MiH$cum_NBS, y=cum_cal_NBS_ts_6$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-250,1050), ylim=c(-270,1050),
     xlab= "Lake Michigan-Huron 6-Month Water Supply [mm]", ylab = "Lake Superior 6-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_6$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_6$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# MiH-Eri
plot(x=cum_cal_NBS_ts_6$MiH$cum_NBS, y=cum_cal_NBS_ts_6$Eri$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-250,1050), ylim=c(-550, 1350), 
     xlab= "Lake Michigan-Huron 6-Month Water Supply [mm]", ylab = "Lake Erie 6-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_6$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_6$Eri)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Sup-Eri
plot(x=cum_cal_NBS_ts_6$Eri$cum_NBS, y=cum_cal_NBS_ts_6$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-550,1350), ylim=c(-270,1050),
     xlab= "Lake Erie 6-Month Water Supply [mm]", ylab = "Lake Superior 6-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_6$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_6$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# MiH-Ont
plot(x=cum_cal_NBS_ts_6$MiH$cum_NBS, y=cum_cal_NBS_ts_6$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-250,1050), ylim=c(-50,1900),
     xlab= "Lake Michigan-Huron 6-Month Water Supply [mm]", ylab = "Lake Ontario 6-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_6$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_6$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Sup-Ont
plot(x=cum_cal_NBS_ts_6$Ont$cum_NBS, y=cum_cal_NBS_ts_6$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-50,1900), ylim=c(-270, 1050),
     xlab= "Lake Ontario 6-Month Water Supply [mm]", ylab = "Lake Superior 6-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_6$Ont)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_6$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Eri-Ont
plot(x=cum_cal_NBS_ts_6$Eri$cum_NBS, y=cum_cal_NBS_ts_6$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-550,1350), ylim=c(-50,1900),
     xlab= "Lake Erie 6-Month Water Supply [mm]", ylab = "Lake Ontario 6-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_6$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_6$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)
dev.off()

# 12-month NBS correlation across all lake pairs (for Supplementary Information)  
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/Figures")
png("12MonthSpatialCorrAll.png", width = 8, height = 10, units = "in", res = 300)
par(mfrow = c(3, 2), oma = c(1, 0.5, 0.5, 0), mar = c(4, 1, 1, 0), pty="s", mgp = c(2.3, 1, 0))
# Sup-MiH
plot(x=cum_cal_NBS_ts_12$MiH$cum_NBS, y=cum_cal_NBS_ts_12$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(350,1700), ylim=c(200,1250),
     xlab= "Lake Michigan-Huron 12-Month Water Supply [mm]", ylab = "Lake Superior 12-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_12$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_12$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# MiH-Eri
plot(x=cum_cal_NBS_ts_12$MiH$cum_NBS, y=cum_cal_NBS_ts_12$Eri$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(350,1700), ylim=c(-50,2100),
     xlab= "Lake Michigan-Huron 12-Month Water Supply [mm]", ylab = "Lake Erie 12-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_12$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_12$Eri)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Sup-Eri
plot(x=cum_cal_NBS_ts_12$Eri$cum_NBS, y=cum_cal_NBS_ts_12$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-50,2100), ylim=c(200,1250),
     xlab= "Lake Erie 12-Month Water Supply [mm]", ylab = "Lake Superior 12-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_12$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_12$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# MiH-Ont
plot(x=cum_cal_NBS_ts_12$MiH$cum_NBS, y=cum_cal_NBS_ts_12$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(350,1700), ylim=c(800,3500),
     xlab= "Lake Michigan-Huron 12-Month Water Supply [mm]", ylab = "Lake Ontario 12-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_12$MiH)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_12$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Sup-Ont
plot(x=cum_cal_NBS_ts_12$Ont$cum_NBS, y=cum_cal_NBS_ts_12$Sup$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(800,3500), ylim=c(200,1250),
     xlab= "Lake Ontario 12-Month Water Supply [mm]", ylab = "Lake Superior 12-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_12$Ont)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_12$Sup)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)

# Eri-Ont
plot(x=cum_cal_NBS_ts_12$Eri$cum_NBS, y=cum_cal_NBS_ts_12$Ont$cum_NBS, pch=19, col="gray", cex=0.8, xlim=c(-50,2100), ylim=c(800,3500),
     xlab= "Lake Erie 12-Month Water Supply [mm]", ylab = "Lake Ontario 12-Month Water Supply [mm]")
points(x=na.omit(data.matrix(rcum_samp_NBS_ts_12$Eri)[,4:7]), y=na.omit(data.matrix(rcum_samp_NBS_ts_12$Ont)[,4:7]), pch=19, col=rgb(1,0.4,0.4), cex=0.8)
dev.off()



##################
# 1 month Forecast Time Series Validation
##################
png(file="1MonthValidationTS.png", width=1150, height=750)
par(mfrow=c(4,1),mar=c(0,0,0,0), oma=c(4,7.5,2,5), xaxs="i", yaxs="i")

plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-250,700), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(-200,0,200,400, 600), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Superior", outer=FALSE, line=2, cex=1.2, font=1)
Sup1_90s  <- apply(as.matrix(cum_samp_NBS_ts_1$Sup[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup1_90s <- apply(as.matrix(rcum_samp_NBS_ts_1$Sup[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 1:168){
  rect(xleft=val_month-0.4, ybottom=rSup1_90s[1, val_month], xright=val_month+0.4, ytop=rSup1_90s[2, val_month], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Sup1_90s[1, val_month],Sup1_90s[2, val_month]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}
points(cum_val_NBS_ts_1$Sup$cum_NBS, pch=19)
box()

# legend(50,820, legend = c("Observed monthly water supply", "Random 90% prediction interval", "Matched 90% prediction interval"), col = c("black", "gray", rgb(1,0.4,0.4)), 
#        lty = c(0, 0, 0), pch = c(19, 15, 15), bty = "n", cex = 1.2, inset=-0.1, xpd=NA, ncol=3,
#        pt.bg=c("black", "gray", rgb(1,0.4,0.4)), pt.cex=c(1,2,2))

plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-250,700), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(-200,0,200,400, 600), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Michigan-Huron", outer=FALSE, line=2, cex=1.2, font=1)
MiH1_90s <- apply(as.matrix(cum_samp_NBS_ts_1$MiH[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH1_90s <- apply(as.matrix(rcum_samp_NBS_ts_1$MiH[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 1:168){
  rect(xleft=val_month-0.4, ybottom=rMiH1_90s[1, val_month], xright=val_month+0.4, ytop=rMiH1_90s[2, val_month], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(MiH1_90s[1, val_month],MiH1_90s[2, val_month]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_1$MiH$cum_NBS, pch=19)
box()

plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-250,700), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(-200,0,200,400, 600), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Erie", outer=FALSE, line=2, cex=1.2, font=1)
Eri1_90s <- apply(as.matrix(cum_samp_NBS_ts_1$Eri[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri1_90s <- apply(as.matrix(rcum_samp_NBS_ts_1$Eri[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 1:168){
  rect(xleft=val_month-0.4, ybottom=rEri1_90s[1, val_month], xright=val_month+0.4, ytop=rEri1_90s[2, val_month], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Eri1_90s[1, val_month],Eri1_90s[2, val_month]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_1$Eri$cum_NBS, pch=19)
box()

plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-250,700), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(-200,0,200,400,600), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Ontario", outer=FALSE, line=2, cex=1.2, font=1)
Ont1_90s <- apply(as.matrix(cum_samp_NBS_ts_1$Ont[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt1_90s <- apply(as.matrix(rcum_samp_NBS_ts_1$Ont[1:168, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 1:168){
  rect(xleft=val_month-0.4, ybottom=rOnt1_90s[1, val_month], xright=val_month+0.4, ytop=rOnt1_90s[2, val_month], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Ont1_90s[1, val_month],Ont1_90s[2, val_month]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_1$Ont$cum_NBS, pch=19)
axis(1,at=seq(1,181, by=1),labels=F,tck=-0.02)
axis(1,at=seq(1,181,by=12),labels=F,tck=-0.05)
axis(1,at=seq(7,177,by=12),tck=0,labels=seq(2007,2021), cex.axis=1.5, padj=0.5)
box()

mtext(side=2, text="1-month Cumulative Water Supply (mm)", outer=TRUE, line=5.3, font=2)

dev.off()


##################
# 3 month Forecast Time Series Validation
###################
png(file="3MonthValidationTS.png", width=1150, height=750)
par(mfrow=c(4,1),mar=c(0,0,0,0), oma=c(4,7.5,2,5), xaxs="i", yaxs="i")

plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-500,1500), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(-500,0,500,1000,1500), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Superior", outer=FALSE, line=2, cex=1.2, font=1)
Sup3_90s <- apply(as.matrix(cum_samp_NBS_ts_3$Sup[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup3_90s <- apply(as.matrix(rcum_samp_NBS_ts_3$Sup[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 3:170){
  rect(xleft=val_month-0.4, ybottom=rSup3_90s[1, val_month-2], xright=val_month+0.4, ytop=rSup3_90s[2, val_month-2], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Sup3_90s[1, val_month-2],Sup3_90s[2, val_month-2]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}
points(cum_val_NBS_ts_3$Sup$cum_NBS, pch=19)
box()

# legend(50,1800, legend = c("Observed monthly water supply", "Random 90% prediction interval", "Matched 90% prediction interval"), col = c("black", "gray", rgb(1,0.4,0.4)), 
#        lty = c(0, 0, 0), pch = c(19, 15, 15), bty = "n", cex = 1.2, inset=-0.1, xpd=NA, ncol=3,
#        pt.bg=c("black", "gray", rgb(1,0.4,0.4)), pt.cex=c(1,2,2))


plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-500,1500), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(-500,0,500,1000,1500), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Michigan-Huron", outer=FALSE, line=2, cex=1.2, font=1)
MiH3_90s <- apply(as.matrix(cum_samp_NBS_ts_3$MiH[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH3_90s <- apply(as.matrix(rcum_samp_NBS_ts_3$MiH[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 3:170){
  rect(xleft=val_month-0.4, ybottom=rMiH3_90s[1, val_month-2], xright=val_month+0.4, ytop=rMiH3_90s[2, val_month-2], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(MiH3_90s[1, val_month-2],MiH3_90s[2, val_month-2]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_3$MiH$cum_NBS, pch=19)
box()


plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-500,1500), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(-500,0,500,1000,1500), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Erie", outer=FALSE, line=2, cex=1.2, font=1)
Eri3_90s <- apply(as.matrix(cum_samp_NBS_ts_3$Eri[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri3_90s <- apply(as.matrix(rcum_samp_NBS_ts_3$Eri[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 3:170){
  rect(xleft=val_month-0.4, ybottom=rEri3_90s[1, val_month-2], xright=val_month+0.4, ytop=rEri3_90s[2, val_month-2], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Eri3_90s[1, val_month-2],Eri3_90s[2, val_month-2]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_3$Eri$cum_NBS, pch=19)
box()


plot(1, type="n", xlab="",ylab="", xlim=c(0,180),ylim=c(-500,1500), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(-500,0,500,1000,1500), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Ontario", outer=FALSE, line=2, cex=1.2, font=1)
Ont3_90s <- apply(as.matrix(cum_samp_NBS_ts_3$Ont[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt3_90s <- apply(as.matrix(rcum_samp_NBS_ts_3$Ont[3:170, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 3:170){
  rect(xleft=val_month-0.4, ybottom=rOnt3_90s[1, val_month-2], xright=val_month+0.4, ytop=rOnt3_90s[2, val_month-2], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Ont3_90s[1, val_month-2],Ont3_90s[2, val_month-2]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_3$Ont$cum_NBS, pch=19)

axis(1,at=seq(1,181, by=1),labels=F,tck=-0.02)
axis(1,at=seq(1,181,by=12),labels=F,tck=-0.05)
axis(1,at=seq(7,177,by=12),tck=0,labels=seq(2007,2021), cex.axis=1.5, padj=0.5)
box()

mtext(side=2, text="3-month Cumulative Water Supply (mm)", outer=TRUE, line=5.3, font=2)

dev.off()


##################
# 6 month Forecast Time Series Validation
###################
png(file="6MonthValidationTS.png", width=1150, height=750)
par(mfrow=c(4,1),mar=c(0,0,0,0), oma=c(4,7.5,2,5), xaxs="i", yaxs="i")

plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(-250,1000), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(0,500,1000), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Superior", outer=FALSE, line=2, cex=1.2, font=1)
Sup6_90s <- apply(as.matrix(cum_samp_NBS_ts_6$Sup[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup6_90s <- apply(as.matrix(rcum_samp_NBS_ts_6$Sup[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 6:173){
  rect(xleft=val_month-0.4, ybottom=rSup6_90s[1, val_month-5], xright=val_month+0.4, ytop=rSup6_90s[2, val_month-5], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Sup6_90s[1, val_month-5],Sup6_90s[2, val_month-5]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}
points(cum_val_NBS_ts_6$Sup$cum_NBS, pch=19)
box()

# legend(50,1200, legend = c("Observed monthly water supply", "Random 90% prediction interval", "Matched 90% prediction interval"), col = c("black", "gray", rgb(1,0.4,0.4)), 
#        lty = c(0, 0, 0), pch = c(19, 15, 15), bty = "n", cex = 1.2, inset=-0.1, xpd=NA, ncol=3,
#        pt.bg=c("black", "gray", rgb(1,0.4,0.4)), pt.cex=c(1,2,2))


plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(-250,1250), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(0,500,1000), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Michigan-Huron", outer=FALSE, line=2, cex=1.2, font=1)
MiH6_90s <- apply(as.matrix(cum_samp_NBS_ts_6$MiH[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH6_90s <- apply(as.matrix(rcum_samp_NBS_ts_6$MiH[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 6:173){
  rect(xleft=val_month-0.4, ybottom=rMiH6_90s[1, val_month-5], xright=val_month+0.4, ytop=rMiH6_90s[2, val_month-5], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(MiH6_90s[1, val_month-5],MiH6_90s[2, val_month-5]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_6$MiH$cum_NBS, pch=19)
box()


plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(-500,1500), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(-500,0,500,1000,1500), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Erie", outer=FALSE, line=2, cex=1.2, font=1)
Eri6_90s <- apply(as.matrix(cum_samp_NBS_ts_6$Eri[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri6_90s <- apply(as.matrix(rcum_samp_NBS_ts_6$Eri[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 6:173){
  rect(xleft=val_month-0.4, ybottom=rEri6_90s[1, val_month-5], xright=val_month+0.4, ytop=rEri6_90s[2, val_month-5], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Eri6_90s[1, val_month-5],Eri6_90s[2, val_month-5]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_6$Eri$cum_NBS, pch=19)
box()


plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(0,2000), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(0,500,1000,1500,2000), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Ontario", outer=FALSE, line=2, cex=1.2, font=1)
Ont6_90s <- apply(as.matrix(cum_samp_NBS_ts_6$Ont[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt6_90s <- apply(as.matrix(rcum_samp_NBS_ts_6$Ont[6:173, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 6:173){
  rect(xleft=val_month-0.4, ybottom=rOnt6_90s[1, val_month-5], xright=val_month+0.4, ytop=rOnt6_90s[2, val_month-5], 
       col="gray", border=NA)
  lines(x=c(val_month,val_month),y=c(Ont6_90s[1, val_month-5],Ont6_90s[2, val_month-5]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_6$Ont$cum_NBS, pch=19)
axis(1,at=seq(1,181, by=1),labels=F,tck=-0.02)
axis(1,at=seq(1,181,by=12),labels=F,tck=-0.05)
axis(1,at=seq(7,177,by=12),tck=0,labels=seq(2007,2021), cex.axis=1.5, padj=0.5)
box()

mtext(side=2, text="6-month Cumulative Water Supply (mm)", outer=TRUE, line=5.3, font=2)

dev.off()
###################

# 12 month Forecast Time Series Validation
###################
png(file="12MonthValidationTS.png", width=1150, height=750)
par(mfrow=c(4,1),mar=c(0,0,0,0), oma=c(4,7.5,2,5), xaxs="i", yaxs="i")

plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(150,1250), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(0,500,1000), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Superior", outer=FALSE, line=2, cex=1.2, font=1)
Sup12_90s <- apply(as.matrix(cum_samp_NBS_ts_12$Sup[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rSup12_90s <- apply(as.matrix(rcum_samp_NBS_ts_12$Sup[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 12:179){
  rect(xleft=val_month-0.4, ybottom=rSup12_90s[1, val_month-11], xright=val_month+0.4, ytop=rSup12_90s[2, val_month-11], 
       col="gray", border=NA)
  # rect(xleft=val_month-0.25, ybottom=Sup12_90s[1, val_month-11], xright=val_month+0.25, ytop=Sup12_90s[2, val_month-11], 
  #      col=rgb(1,0.4,0.4), border=NA)
  lines(x=c(val_month,val_month),y=c(Sup12_90s[1, val_month-11],Sup12_90s[2, val_month-11]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}
points(cum_val_NBS_ts_12$Sup$cum_NBS, pch=19)
box()

# legend(50,1400, legend = c("Observed monthly water supply", "Random 90% prediction interval", "Matched 90% prediction interval"), col = c("black", "gray", rgb(1,0.4,0.4)), 
#        lty = c(0, 0, 0), pch = c(19, 15, 15), bty = "n", cex = 1.2, inset=-0.1, xpd=NA, ncol=3,
#        pt.bg=c("black", "gray", rgb(1,0.4,0.4)), pt.cex=c(1,2,2))


plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(400,1700), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(0,500,1000,1500), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Michigan-Huron", outer=FALSE, line=2, cex=1.2, font=1)
MiH12_90s <- apply(as.matrix(cum_samp_NBS_ts_12$MiH[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rMiH12_90s <- apply(as.matrix(rcum_samp_NBS_ts_12$MiH[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 12:179){
  rect(xleft=val_month-0.4, ybottom=rMiH12_90s[1, val_month-11], xright=val_month+0.4, ytop=rMiH12_90s[2, val_month-11], 
       col="gray", border=NA)
  # rect(xleft=val_month-0.25, ybottom=MiH12_90s[1, val_month-11], xright=val_month+0.25, ytop=MiH12_90s[2, val_month-11], 
  #      col=rgb(1,0.4,0.4), border=NA)
  lines(x=c(val_month,val_month),y=c(MiH12_90s[1, val_month-11],MiH12_90s[2, val_month-11]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_12$MiH$cum_NBS, pch=19)
box()


plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(0,2000), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(2, at=c(-500,0,500,1000,1500,2000), labels=T, cex.axis=1.5, las=1)
mtext(side=4,text="Erie", outer=FALSE, line=2, cex=1.2, font=1)
Eri12_90s <- apply(as.matrix(cum_samp_NBS_ts_12$Eri[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rEri12_90s <- apply(as.matrix(rcum_samp_NBS_ts_12$Eri[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 12:179){
  rect(xleft=val_month-0.4, ybottom=rEri12_90s[1, val_month-11], xright=val_month+0.4, ytop=rEri12_90s[2, val_month-11], 
       col="gray", border=NA)
  # rect(xleft=val_month-0.25, ybottom=Eri12_90s[1, val_month-11], xright=val_month+0.25, ytop=Eri12_90s[2, val_month-11], 
  #      col=rgb(1,0.4,0.4), border=NA)
  lines(x=c(val_month,val_month),y=c(Eri12_90s[1, val_month-11],Eri12_90s[2, val_month-11]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_12$Eri$cum_NBS, pch=19)
box()


plot(1, type="n", xlab="",ylab="", xlim=c(1,180),ylim=c(750,3000), axes=F, main="")
abline(v=seq(1,181,by=12), col="lightgray")
axis(4, at=c(0,500,1000,1500,2000,2500,3000), labels=T, cex.axis=1.5, las=1)
mtext(side=2,text="Ontario", outer=FALSE, line=2, cex=1.2, font=1)
Ont12_90s <- apply(as.matrix(cum_samp_NBS_ts_12$Ont[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
rOnt12_90s <- apply(as.matrix(rcum_samp_NBS_ts_12$Ont[12:179, 4:59]), 1, quantile, probs = c(0.05, 0.95))
for (val_month in 12:179){
  rect(xleft=val_month-0.4, ybottom=rOnt12_90s[1, val_month-11], xright=val_month+0.4, ytop=rOnt12_90s[2, val_month-11], 
       col="gray", border=NA)
  # rect(xleft=val_month-0.25, ybottom=Ont12_90s[1, val_month-11], xright=val_month+0.25, ytop=Ont12_90s[2, val_month-11], 
  #      col=rgb(1,0.4,0.4), border=NA)
  lines(x=c(val_month,val_month),y=c(Ont12_90s[1, val_month-11],Ont12_90s[2, val_month-11]), col=rgb(1,0.4,0.4), lwd=2, lend="butt")
}

points(cum_val_NBS_ts_12$Ont$cum_NBS, pch=19)
axis(1,at=seq(1,181, by=1),labels=F,tck=-0.02)
axis(1,at=seq(1,181,by=12),labels=F,tck=-0.05)
axis(1,at=seq(7,177,by=12),tck=0,labels=seq(2007,2021), cex.axis=1.5, padj=0.5)
box()

mtext(side=2, text="12-month Cumulative Water Supply (mm)", outer=TRUE, line=5.3, font=2)

dev.off()
###################



