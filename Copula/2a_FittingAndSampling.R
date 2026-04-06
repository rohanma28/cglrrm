rm(list=ls())
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")

library(copula)
library(VineCopula)
library(rvinecopulib)
library(FAdist)
library(brms)


forecast_month <- 1

######
## Load formatted input P, E, and R values from 1950-2020 for copula fitting
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/1_out_formatted_input")
cal <- data.matrix(read.csv(paste(month.abb[forecast_month],"_12Forecast_3Ant_CopulaInput.csv",sep=""))[1:70,2:226])
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")
######

######
## Define file names to store/retrieve the copula models, pseudo samples, and samples
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
struct_name <- paste(month.abb[forecast_month], "_15month_1950_2020_RvineCop.rds", sep = "")
dist_name <- paste(month.abb[forecast_month], "_15month_1950_2020_RvineDist.rds", sep = "")
psamp_name <- paste(month.abb[forecast_month], "_15month_1950_2020_psamp.rds", sep = "")
samp_name <- paste(month.abb[forecast_month], "_15month_1950_2020_samp.rds", sep = "")
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")
######

######
## Generate vine structures to fit the copula model (commented out by default, takes ~1 hr)
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
pseudo      <- pseudo_obs(cal)
struct_vine <- vinecop(pseudo,cores=5)
saveRDS(struct_vine,struct_name)
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")
######

######
## Generate pseudo samples from the copula models (commented out by default, superceded, takes ~15 min)
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
cop   <- readRDS(struct_name)
psamp <- rvinecop(10000,cop)
saveRDS(psamp,psamp_name)
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")
######

######
## Fit marginal distributions to observed values 
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")

obs_mat <- cal
marginal_pars <- vector(mode="list",length=225)
names(marginal_pars) <- colnames(obs_mat)

for (col in 1:225){
  comp <- substr(colnames(obs_mat)[col],start=nchar(colnames(obs_mat)[col]),stop=nchar(colnames(obs_mat)[col]))
  vals <- obs_mat[,col]
  
  if (comp=="P"){
    # Define gamma parameters
    marginal_pars[[col]]$x.bar_s <- mean(vals, na.rm=T)
    marginal_pars[[col]]$theta_s <- log(marginal_pars[[col]]$x.bar_s)-mean(log(vals),na.rm = TRUE)
    marginal_pars[[col]]$shape   <- 1/(4*marginal_pars[[col]]$theta_s)*(1+sqrt(1+4*marginal_pars[[col]]$theta_s/3))
    marginal_pars[[col]]$rate    <- marginal_pars[[col]]$shape/marginal_pars[[col]]$x.bar_s
    marginal_pars[[col]]$scale   <- 1/marginal_pars[[col]]$rate
    marginal_pars[[col]]$thres   <- 0
  }
  
  else if(comp=="E"){
    # Define normal parameters
    marginal_pars[[col]]$mean    <- mean(vals,na.rm=T)
    marginal_pars[[col]]$var     <- var(vals,na.rm=T)
    marginal_pars[[col]]$std     <- sqrt(marginal_pars[[col]]$var)
    marginal_pars[[col]]$prec    <- 1/marginal_pars[[col]]$var
  }
  
  else if(comp=="R"){
    # Define log-normal parameters
    marginal_pars[[col]]$logmean <- mean(log(vals))
    marginal_pars[[col]]$logvar  <- var(log(vals))
    marginal_pars[[col]]$logstd  <- sqrt(marginal_pars[[col]]$logvar)
    marginal_pars[[col]]$logprec <- 1/marginal_pars[[col]]$logvar
    marginal_pars[[col]]$shift   <- 0
  }
  
  else{
    print("Error in observation column titles")
  }
}
saveRDS(marginal_pars,"1950_2020_Marginal_CDF_full_pars.rds")
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")
######

######
## Create lists of the marginal distributions
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
marginal_lists <- marginal_pars

for (col in 1:225){
  comp <- substr(names(marginal_pars[col]),13,13)
  
  if(comp=="P"){
    marginal_lists[[col]] <- list(distr="gamma", 
                                           shape = marginal_pars[[col]]$shape,
                                           rate  = marginal_pars[[col]]$rate)
  }
  
  else if(comp=="E"){
    marginal_lists[[col]] <- list(distr="norm",
                                           mean = marginal_pars[[col]]$mean,
                                           sd   = marginal_pars[[col]]$std)
  }
  
  else if(comp=="R"){
    marginal_lists[[col]] <- list(distr="lnorm",
                                           meanlog = marginal_pars[[col]]$logmean,
                                           sdlog   = marginal_pars[[col]]$logstd)
  }
  
}
saveRDS(marginal_lists,"1950_2020_Marginal_CDF_lists.rds")
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")

######
## Create objects representing the copula distributions using the structure and bivariate copula found earlier and the marginal distributions
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
struct <- readRDS(struct_name)
dist   <- vine_dist(margins=marginal_lists, pair_copulas = struct$pair_copulas, structure = struct$structure)
saveRDS(dist,dist_name)
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")

######
## Generate samples from the copula distributions (
###### 
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
dist   <- readRDS(dist_name)
samp  <- rvine(100000,dist, cores=5)
colnames(samp) <- colnames(cal)
saveRDS(samp,samp_name)
######

######
## Convert pseudo samples to component values using the marginal distributions (commented out by default, superceded takes ~2 mins)
######
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
psamp_mat <- readRDS(psamp_name)
samp_mat  <- psamp_mat   # #Placeholder, real values will replace [0,1] below

for (col in 1:225){
  comp  <- substr(colnames(psamp_mat)[col],start=nchar(colnames(psamp_mat)[col]),stop=nchar(colnames(psamp_mat)[col]))
  pvals <- psamp_mat[,col]

  if(comp=="P"){
    samp_mat[,col] <- qgamma3(pvals,shape=marginal_pars[[col]]$shape,
                                    scale=marginal_pars[[col]]$scale,
                                    thres=marginal_pars[[col]]$thres)
  }

  else if(comp=="E"){
    samp_mat[,col] <- qnorm(pvals,mean=marginal_pars[[col]]$mean,
                                  sd=marginal_pars[[col]]$std)
  }

  else if(comp=="R"){
    samp_mat[,col] <- qshifted_lnorm(pvals,meanlog=marginal_pars[[col]]$logmean,
                                           sdlog=marginal_pars[[col]]$logstd,
                                           shift=marginal_pars[[col]]$shift)
  }

  else{
    print("Error in pseudo-sample column titles")
  }

}
saveRDS(samp_mat,samp_name)
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula")
######
