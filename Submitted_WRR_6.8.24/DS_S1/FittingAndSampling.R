rm(list=ls())
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")

library(copula)
library(VineCopula)
library(rvinecopulib)
library(FAdist)
library(brms)

######
## Load formatted input P, E, and R values from 1950-2006 (80% of data) for copula fitting
######
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/formatted_input")
cals <- vector("list",length=12)
for (month in 1:12){cals[[month]] <- data.matrix(read.csv(paste(month.abb[month],"_12Forecast_3Ant_CopulaInput.csv",sep=""))[1:56,2:181])}
names(cals) = month.abb
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
######

######
## Define file names to store/retrieve the copula models, pseudo samples, and samples
######
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
struct_names <- c("Jan_15month_1950_2006_RvineCop.rds",
                  "Feb_15month_1950_2006_RvineCop.rds",
                  "Mar_15month_1950_2006_RvineCop.rds",
                  "Apr_15month_1950_2006_RvineCop.rds",
                  "May_15month_1950_2006_RvineCop.rds",
                  "Jun_15month_1950_2006_RvineCop.rds",
                  "Jul_15month_1950_2006_RvineCop.rds",
                  "Aug_15month_1950_2006_RvineCop.rds",
                  "Sep_15month_1950_2006_RvineCop.rds",
                  "Oct_15month_1950_2006_RvineCop.rds",
                  "Nov_15month_1950_2006_RvineCop.rds",
                  "Dec_15month_1950_2006_RvineCop.rds")

dist_names  <-  c("Jan_15month_1950_2006_RvineDist.rds",
                  "Feb_15month_1950_2006_RvineDist.rds",
                  "Mar_15month_1950_2006_RvineDist.rds",
                  "Apr_15month_1950_2006_RvineDist.rds",
                  "May_15month_1950_2006_RvineDist.rds",
                  "Jun_15month_1950_2006_RvineDist.rds",
                  "Jul_15month_1950_2006_RvineDist.rds",
                  "Aug_15month_1950_2006_RvineDist.rds",
                  "Sep_15month_1950_2006_RvineDist.rds",
                  "Oct_15month_1950_2006_RvineDist.rds",
                  "Nov_15month_1950_2006_RvineDist.rds",
                  "Dec_15month_1950_2006_RvineDist.rds")

psamp_names <- c("Jan_15month_1950_2006_psamp.rds",
                 "Feb_15month_1950_2006_psamp.rds",
                 "Mar_15month_1950_2006_psamp.rds",
                 "Apr_15month_1950_2006_psamp.rds",
                 "May_15month_1950_2006_psamp.rds",
                 "Jun_15month_1950_2006_psamp.rds",
                 "Jul_15month_1950_2006_psamp.rds",
                 "Aug_15month_1950_2006_psamp.rds",
                 "Sep_15month_1950_2006_psamp.rds",
                 "Oct_15month_1950_2006_psamp.rds",
                 "Nov_15month_1950_2006_psamp.rds",
                 "Dec_15month_1950_2006_psamp.rds")

samp_names <- c("Jan_15month_1950_2006_samp.rds",
                "Feb_15month_1950_2006_samp.rds",
                "Mar_15month_1950_2006_samp.rds",
                "Apr_15month_1950_2006_samp.rds",
                "May_15month_1950_2006_samp.rds",
                "Jun_15month_1950_2006_samp.rds",
                "Jul_15month_1950_2006_samp.rds",
                "Aug_15month_1950_2006_samp.rds",
                "Sep_15month_1950_2006_samp.rds",
                "Oct_15month_1950_2006_samp.rds",
                "Nov_15month_1950_2006_samp.rds",
                "Dec_15month_1950_2006_samp.rds")

setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
######

######
## Generate vine structures to fit the copula model (commented out by default, takes ~1 hr)
###### 
# setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
# for (month in 1:12){
#   pseudo      <- pseudo_obs(cals[[month]])
#   struct_vine <- vinecop(pseudo,cores=5)
#   saveRDS(struct_vine,struct_names[month])
# }
# setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
######

######
## Generate pseudo samples from the copula models (commented out by default, superceded, takes ~15 min)
###### 
#setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
# for (month in 1:12){
#   cop   <- readRDS(struct_names[month])
#   psamp <- rvinecop(10000,cop)
#   saveRDS(psamp,psamp_names[month])
# }
#setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
######



######
## Fit marginal distributions to observed values 
######
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
marginal_pars <- vector(mode="list",length=12)
names(marginal_pars) <- month.abb

for (month in 1:12){
  obs_mat <- cals[[month]]
  marginal_pars[[month]] <- vector(mode="list",length=180)
  names(marginal_pars[[month]]) <- colnames(obs_mat)
  
  for (col in 1:180){
    comp <- substr(colnames(obs_mat)[col],start=nchar(colnames(obs_mat)[col]),stop=nchar(colnames(obs_mat)[col]))
    vals <- obs_mat[,col]
    
    if (comp=="P"){
      # Define gamma parameters
      marginal_pars[[month]][[col]]$x.bar_s <- mean(vals, na.rm=T)
      marginal_pars[[month]][[col]]$theta_s <- log(marginal_pars[[month]][[col]]$x.bar_s)-mean(log(vals),na.rm = TRUE)
      marginal_pars[[month]][[col]]$shape   <- 1/(4*marginal_pars[[month]][[col]]$theta_s)*(1+sqrt(1+4*marginal_pars[[month]][[col]]$theta_s/3))
      marginal_pars[[month]][[col]]$rate    <- marginal_pars[[month]][[col]]$shape/marginal_pars[[month]][[col]]$x.bar_s
      marginal_pars[[month]][[col]]$scale   <- 1/marginal_pars[[month]][[col]]$rate
      marginal_pars[[month]][[col]]$thres   <- 0
    }
    
    else if(comp=="E"){
      # Define normal parameters
      marginal_pars[[month]][[col]]$mean    <- mean(vals,na.rm=T)
      marginal_pars[[month]][[col]]$var     <- var(vals,na.rm=T)
      marginal_pars[[month]][[col]]$std     <- sqrt(marginal_pars[[month]][[col]]$var)
      marginal_pars[[month]][[col]]$prec    <- 1/marginal_pars[[month]][[col]]$var
    }
    
    else if(comp=="R"){
      # Define log-normal parameters
      marginal_pars[[month]][[col]]$logmean <- mean(log(vals))
      marginal_pars[[month]][[col]]$logvar  <- var(log(vals))
      marginal_pars[[month]][[col]]$logstd  <- sqrt(marginal_pars[[month]][[col]]$logvar)
      marginal_pars[[month]][[col]]$logprec <- 1/marginal_pars[[month]][[col]]$logvar
      marginal_pars[[month]][[col]]$shift   <- 0
    }
    
    else{
      print("Error in observation column titles")
    }
  }
}
saveRDS(marginal_pars,"1950_2006_Marginal_CDF_full_pars.rds")
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
######

######
## Create lists of the marginal distributions
######
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
marginal_lists <- marginal_pars
names(marginal_pars) <- month.abb

for(month in 1:12){
  for (col in 1:180){
    comp <- substr(names(marginal_pars[[month]][col]),13,13)
    
    if(comp=="P"){
      marginal_lists[[month]][[col]] <- list(distr="gamma", 
                                             shape = marginal_pars[[month]][[col]]$shape,
                                             rate  = marginal_pars[[month]][[col]]$rate)
    }
    
    else if(comp=="E"){
      marginal_lists[[month]][[col]] <- list(distr="norm",
                                             mean = marginal_pars[[month]][[col]]$mean,
                                             sd   = marginal_pars[[month]][[col]]$std)
    }
    
    else if(comp=="R"){
      marginal_lists[[month]][[col]] <- list(distr="lnorm",
                                             meanlog = marginal_pars[[month]][[col]]$logmean,
                                             sdlog   = marginal_pars[[month]][[col]]$logstd)
    }
    
  }
    
}
saveRDS(marginal_lists,"1950_2006_Marginal_CDF_lists.rds")
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")

######
## Create objects representing the copula distributions using the structure and bivariate copula found earlier and the marginal distributions
######
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
for (month in 1:12){
  struct <- readRDS(struct_names[month])
  dist   <- vine_dist(margins=marginal_lists[[month]], pair_copulas = struct$pair_copulas, structure = struct$structure)
  saveRDS(dist,dist_names[month])
}
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")


######
## Generate samples from the copula distributions (
###### 
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
for (month in 1:12){
  dist   <- readRDS(dist_names[month])
  samp  <- rvine(100000,dist, cores=5)
  colnames(samp) <- colnames(cals[[month]])
  saveRDS(samp,samp_names[month])
}


######


######
## Convert pseudo samples to component values using the marginal distributions (commented out by default, superceded takes ~2 mins)
######
# setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
# for (month in 1:12){
#   psamp_mat <- readRDS(psamp_names[[month]])
#   samp_mat  <- psamp_mat   # #Placeholder, real values will replace [0,1] below
#   
#   for (col in 1:180){
#     comp  <- substr(colnames(psamp_mat)[col],start=nchar(colnames(psamp_mat)[col]),stop=nchar(colnames(psamp_mat)[col]))
#     pvals <- psamp_mat[,col]
#     
#     if(comp=="P"){
#       samp_mat[,col] <- qgamma3(pvals,shape=marginal_pars[[month]][[col]]$shape, 
#                                       scale=marginal_pars[[month]][[col]]$scale,
#                                       thres=marginal_pars[[month]][[col]]$thres)
#     }
#     
#     else if(comp=="E"){
#       samp_mat[,col] <- qnorm(pvals,mean=marginal_pars[[month]][[col]]$mean,
#                                     sd=marginal_pars[[month]][[col]]$std)
#     }
#     
#     else if(comp=="R"){
#       samp_mat[,col] <- qshifted_lnorm(pvals,meanlog=marginal_pars[[month]][[col]]$logmean,
#                                              sdlog=marginal_pars[[month]][[col]]$logstd,
#                                              shift=marginal_pars[[month]][[col]]$shift)
#     }
#     
#     else{
#       print("Error in pseudo-sample column titles")
#     }
#     
#   }
#   saveRDS(samp_mat,samp_names[month])
# }
# setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
######

