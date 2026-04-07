library(copula)
library(VineCopula)
library(rvinecopulib)
library(FAdist)
library(brms)


forecast_month <- 1

######
## Load formatted input P, E, and R values from 1950-2020 for copula fitting
######
cal <- data.matrix(read.csv(paste("1_out_formatted_input/", month.abb[forecast_month],"_12Forecast_3Ant_CopulaInput.csv",sep=""))[1:70,2:226])
######

######
## Define file names to store/retrieve the copula models, pseudo samples, and samples
######
struct_name <- paste("R_objects/", month.abb[forecast_month], "_15month_1950_2020_RvineCop.rds", sep = "")
dist_name <- paste("R_objects/", month.abb[forecast_month], "_15month_1950_2020_RvineDist.rds", sep = "")
psamp_name <- paste("R_objects/", month.abb[forecast_month], "_15month_1950_2020_psamp.rds", sep = "")
samp_name <- paste("R_objects/", month.abb[forecast_month], "_15month_1950_2020_samp.rds", sep = "")
######

######
## Generate vine structures to fit the copula model (commented out by default, takes ~1 hr)
######
pseudo      <- pseudo_obs(cal)
struct_vine <- vinecop(pseudo,cores=5)
saveRDS(struct_vine,struct_name)
######

######
## Generate pseudo samples from the copula models (commented out by default, superceded, takes ~15 min)
######
cop   <- readRDS(struct_name)
psamp <- rvinecop(10000,cop)
saveRDS(psamp,psamp_name)
######

######
## Fit marginal distributions to observed values 
######

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
saveRDS(marginal_pars,"R_objects/1950_2020_Marginal_CDF_full_pars.rds")
######

######
## Create lists of the marginal distributions
######
# marginal_lists <- marginal_pars
# 
# for (col in 1:225){
#   comp <- substr(names(marginal_pars[col]),13,13)
#   
#   if(comp=="P"){
#     marginal_lists[[col]] <- list(distr="gamma", 
#                                            shape = marginal_pars[[col]]$shape,
#                                            rate  = marginal_pars[[col]]$rate)
#   }
#   
#   else if(comp=="E"){
#     marginal_lists[[col]] <- list(distr="norm",
#                                            mean = marginal_pars[[col]]$mean,
#                                            sd   = marginal_pars[[col]]$std)
#   }
#   
#   else if(comp=="R"){
#     marginal_lists[[col]] <- list(distr="lnorm",
#                                            meanlog = marginal_pars[[col]]$logmean,
#                                            sdlog   = marginal_pars[[col]]$logstd)
#   }
#   
# }
# saveRDS(marginal_lists,"R_objects/1950_2020_Marginal_CDF_lists.rds")

######
## Create objects representing the copula distributions using the structure and bivariate copula found earlier and the marginal distributions
######
# struct <- readRDS(struct_name)
# dist   <- vine_dist(margins=marginal_lists, pair_copulas = struct$pair_copulas, structure = struct$structure)
# saveRDS(dist,dist_name)

######
## Generate samples from the copula distributions (
###### 
# dist   <- readRDS(dist_name)
# samp  <- rvine(100000,dist, cores=5)
# colnames(samp) <- colnames(cal)
######

######
## Convert pseudo samples to component values using the marginal distributions (commented out by default, superceded takes ~2 mins)
######

### PRECIPITATION

psamp_mat <- readRDS(psamp_name)
samp_mat  <- psamp_mat   # #Placeholder, real values will replace [0,1] below

for (col in 1:225){
  comp  <- substr(colnames(psamp_mat)[col],start=nchar(colnames(psamp_mat)[col]),stop=nchar(colnames(psamp_mat)[col]))
  pvals <- psamp_mat[,col]
  
  if(comp=="E"){
    samp_mat[,col] <- qnorm(pvals,mean=marginal_pars[[col]]$mean,
                            sd=marginal_pars[[col]]$std)
  }
  
  else if(comp=="R"){
    samp_mat[,col] <- qshifted_lnorm(pvals,meanlog=marginal_pars[[col]]$logmean,
                                     sdlog=marginal_pars[[col]]$logstd,
                                     shift=marginal_pars[[col]]$shift)
  }
  
  else if (comp != "P"){
    print("Error in pseudo-sample column titles")
  }
  
}

pweights <- c(1.0309, 1.0576, 1.0764, 1.0465, 1.0131, 1.0484, 1.0933, 1.1766, 
              1.059, 1.057, 1.0536, 1.0684, 1.0309, 1.0576, 1.0764)

for (col in c(1, 46, 91, 136, 181)) {
  meanvec <- as.numeric(unlist(sapply(marginal_pars, "[", 1)[col:(col+14)]))
  smoothmeans <- loess.smooth(c(1:15), meanvec * pweights, evaluation = 15)[[2]]
  shapevec <- as.numeric(unlist(sapply(marginal_pars, "[", 3)[col:(col+14)]))
  newscales <- smoothmeans / shapevec
  
  for (subcol in seq(col, col + 14, by = 1)) {
    comp  <- substr(colnames(psamp_mat)[subcol],start=nchar(colnames(psamp_mat)[subcol]),stop=nchar(colnames(psamp_mat)[subcol]))
    pvals <- psamp_mat[,subcol]
    
    samp_mat[,subcol] <- qgamma3(pvals,shape=marginal_pars[[subcol]]$shape,
                                 scale=newscales[subcol %% 45],
                                 thres=marginal_pars[[subcol]]$thres)
  }
}

saveRDS(samp_mat,samp_name)

### EVAPORATION
evapchange <- read.csv("2b_input_changefiles/evapchange.csv", row.names = 1)[,c(3,4,2,1)]
evapchange <- rbind(evapchange[10:12,], evapchange)
evapchange <- cbind(evapchange, rep(0, 15))

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

  else if(comp=="R"){
    samp_mat[,col] <- qshifted_lnorm(pvals,meanlog=marginal_pars[[col]]$logmean,
                                     sdlog=marginal_pars[[col]]$logstd,
                                     shift=marginal_pars[[col]]$shift)
  }
  
  else if (comp != "E"){
    print("Error in pseudo-sample column titles")
  }

}

for (col in c(16, 61, 106, 151, 196)) {
  for (subcol in seq(col, col + 14, by = 1))
  {
    comp  <- substr(colnames(psamp_mat)[subcol],start=nchar(colnames(psamp_mat)[subcol]),stop=nchar(colnames(psamp_mat)[subcol]))
    pvals <- psamp_mat[,subcol]
    
    samp_mat[,subcol] <- qnorm(pvals,
                            mean=(marginal_pars[[subcol]]$mean + evapchange[(subcol %% 45) - 15, (col %/% 45) + 1]),
                            sd=marginal_pars[[subcol]]$std)
  }
}

saveRDS(samp_mat,samp_name)

######

### RUNOFF
runoffchange <- ((read.csv("2b_input_changefiles/runoffchange.csv", row.names = 1)/100)+1)[,c(3,4,2,1)]
runoffchange <- rbind(runoffchange[10:12,], runoffchange)
runoffchange <- cbind(runoffchange, rep(1, 15))

psamp_mat <- readRDS(psamp_name)
samp_mat  <- psamp_mat   # #Placeholder, real values will replace [0,1] below
basesamp_mat <- psamp_mat

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
  
  else if (comp != "R"){
    print("Error in pseudo-sample column titles")
  }
  
}

for (col in c(31, 76, 121, 166, 211)) {
  for (subcol in seq(col, col + 14, by = 1))
  {
    comp  <- substr(colnames(psamp_mat)[subcol],start=nchar(colnames(psamp_mat)[subcol]),stop=nchar(colnames(psamp_mat)[subcol]))
    pvals <- psamp_mat[,subcol]
    
    samp_mat[,subcol] <- qshifted_lnorm(pvals,meanlog=(marginal_pars[[subcol]]$logmean + log(runoffchange[(subcol %% 45) - 30, (col %/% 45) + 1])),
                                                sdlog=marginal_pars[[subcol]]$logstd,
                                                shift=marginal_pars[[subcol]]$shift)
  }
}

saveRDS(samp_mat,samp_name)

######

### THREE-DIMENSIONAL P vs E PLOT

  # supwlmat <- matrix(nrow = 15, ncol = 15)
  # stcwlmat <- matrix(nrow = 15, ncol = 15)
  # mihwlmat <- matrix(nrow = 15, ncol = 15)
  # eriwlmat <- matrix(nrow = 15, ncol = 15)

  # for (precipweight in seq(0.86, 1.14, 0.02)) {
  #     for (evapweight in seq(0.86, 1.14, 0.02)) {
  #         psamp_mat <- readRDS(psamp_name)
  #         samp_mat  <- psamp_mat   # #Placeholder, real values will replace [0,1] below
          
  #         for (col in 1:225){
  #             comp  <- substr(colnames(psamp_mat)[col],start=nchar(colnames(psamp_mat)[col]),stop=nchar(colnames(psamp_mat)[col]))
  #             pvals <- psamp_mat[,col]
              
              
  #             if(comp=="R"){
  #                 samp_mat[,col] <- qshifted_lnorm(pvals,meanlog=marginal_pars[[col]]$logmean,
  #                                                 sdlog=marginal_pars[[col]]$logstd,
  #                                                 shift=marginal_pars[[col]]$shift)
  #             }
              
  #             else if (comp != "P" & comp != "E"){
  #                 print("Error in pseudo-sample column titles")
  #             }
          
  #         }
          
  #         pweights <- rep(precipweight, times=15)
          
  #         for (col in c(1, 46, 91, 136, 181)) {
  #             meanvec <- as.numeric(unlist(sapply(marginal_pars, "[", 1)[col:(col+14)]))
  #             smoothmeans <- meanvec * pweights
  #             shapevec <- as.numeric(unlist(sapply(marginal_pars, "[", 3)[col:(col+14)]))
  #             newscales <- smoothmeans / shapevec
              
  #             for (subcol in seq(col, col + 14, by = 1)) {
  #                 comp  <- substr(colnames(psamp_mat)[subcol],start=nchar(colnames(psamp_mat)[subcol]),stop=nchar(colnames(psamp_mat)[subcol]))
  #                 pvals <- psamp_mat[,subcol]
                  
  #                 samp_mat[,subcol] <- qgamma3(pvals,shape=marginal_pars[[subcol]]$shape,
  #                                             scale=newscales[subcol %% 45],
  #                                             thres=marginal_pars[[subcol]]$thres)
  #             }
  #         }
          
  #         eweights <- rep(evapweight, times=15)
          
  #         for (col in c(16, 61, 106, 151, 196)) {
  #             for (subcol in seq(col, col + 14, by = 1)) {
  #                 comp  <- substr(colnames(psamp_mat)[subcol],start=nchar(colnames(psamp_mat)[subcol]),stop=nchar(colnames(psamp_mat)[subcol]))
  #                 pvals <- psamp_mat[,subcol]
                  
  #                 samp_mat[,subcol] <- qnorm(pvals,
  #                                         mean=(marginal_pars[[subcol]]$mean*eweights[(subcol - 15) %% 45]),
  #                                         sd=marginal_pars[[subcol]]$std)
  #             }
  #         }
              
  #         saveRDS(samp_mat,samp_name)
          
  #         source("3_MatchingAndPlotting.R")
          
  #         supwlmat[round(((precipweight-0.86)*50)+1),round(((evapweight-0.86)*50)+1)] <- mean(supmeans, na.rm = TRUE)
  #         stcwlmat[round(((precipweight-0.86)*50)+1),round(((evapweight-0.86)*50)+1)] <- mean(stcmeans, na.rm = TRUE)
  #         mihwlmat[round(((precipweight-0.86)*50)+1),round(((evapweight-0.86)*50)+1)] <- mean(mihmeans, na.rm = TRUE)
  #         eriwlmat[round(((precipweight-0.86)*50)+1),round(((evapweight-0.86)*50)+1)] <- mean(erimeans, na.rm = TRUE)
  #     }
  # }

  # dir.create("/output/LongTermWLs")

  # write.csv(supwlmat, "/output/LongTermWLs/Superior_LongTerm_PEChange_WL.csv")
  # write.csv(stcwlmat, "/output/LongTermWLs/StClair_LongTerm_PEChange_WL.csv")
  # write.csv(mihwlmat, "/output/LongTermWLs/MichiganHuron_LongTerm_PEChange_WL.csv")
  # write.csv(eriwlmat, "/output/LongTermWLs/Erie_LongTerm_PEChange_WL.csv")
