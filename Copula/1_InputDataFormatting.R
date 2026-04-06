rm(list=ls())
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/L2SWBM_med_1950-2022")

L_forecast <- 12
L_ant      <- 3
forecast_month <- 1

EriP <- read.csv("eriePrecip_analysis19502022_prior19001969_100k.csv")[10:864,]
EriE <- read.csv("erieEvap_analysis19502022_prior19001969_100k.csv")[10:864,]
EriR <- read.csv("erieRunoff_analysis19502022_prior19001969_100k.csv")[10:864,]
OntP <- read.csv("ontarioPrecip_analysis19502022_prior19001969_100k.csv")[10:864,]
OntE <- read.csv("ontarioEvap_analysis19502022_prior19001969_100k.csv")[10:864,]
OntR <- read.csv("ontarioRunoff_analysis19502022_prior19001969_100k.csv")[10:864,]
MiHP <- read.csv("miHuronPrecip_analysis19502022_prior19001969_100k.csv")[10:864,]
MiHE <- read.csv("miHuronEvap_analysis19502022_prior19001969_100k.csv")[10:864,]
MiHR <- read.csv("miHuronRunoff_analysis19502022_prior19001969_100k.csv")[10:864,]
SupP <- read.csv("superiorPrecip_analysis19502022_prior19001969_100k.csv")[10:864,]
SupE <- read.csv("superiorEvap_analysis19502022_prior19001969_100k.csv")[10:864,]
SupR <- read.csv("superiorRunoff_analysis19502022_prior19001969_100k.csv")[10:864,]
StCP <- read.csv("clairPrecip_analysis19502022_prior19001969_100k.csv")[10:864,]
StCE <- read.csv("clairEvap_analysis19502022_prior19001969_100k.csv")[10:864,]
StCR <- read.csv("clairRunoff_analysis19502022_prior19001969_100k.csv")[10:864,]

comps_inp = lapply(list(EriP, EriE, EriR, OntP, OntE, OntR, MiHP, MiHE, MiHR, SupP, SupE, SupR, StCP, StCE, StCR), function(x) {rownames(x) <- c(1:855); x})
names(comps_inp) <- c("EriP", "EriE", "EriR", "OntP", "OntE", "OntR", "MiHP", "MiHE", "MiHR", "SupP", "SupE", "SupR", "StCP", "StCE", "StCR")

for (comp in 1:15){
  comps_inp[[comp]]$CI.95.W <- abs(comps_inp[[comp]]$X97.5.Percentile - comps_inp[[comp]]$X2.5.Percentile)
}

start_month_indices <- which(comps_inp$EriP$Month[4:843]==forecast_month)+3
cal_df <- data.frame()

for(lake_comp in 1:15){
  comp_cal <- data.frame()
  cal_in   <- comps_inp[[lake_comp]]$Median
  
  for (obs_year in 1:length(start_month_indices)){
    cal_row  <- cal_in[(start_month_indices[obs_year]-L_ant):(start_month_indices[obs_year]+L_forecast-1)]
    comp_cal <- rbind(comp_cal,cal_row)
  }
  
  if(ncol(cal_df)==0){cal_df <- comp_cal}
  else               {cal_df <- cbind(cal_df,comp_cal)}
}

start_month_indices <- which(comps_inp$EriP$Month[4:843]==forecast_month)+3
CI95_df <- data.frame()

for(lake_comp in 1:15){
  comp_CI <- data.frame()
  CI_in   <- comps_inp[[lake_comp]]$CI.95.W
  
  for (obs_year in 1:length(start_month_indices)){
    CI_row  <- CI_in[(start_month_indices[obs_year]-L_ant):(start_month_indices[obs_year]+L_forecast-1)]
    comp_CI <- rbind(comp_CI,CI_row)
  }
  
  if(ncol(CI95_df)==0){CI95_df <- comp_CI}
  else                {CI95_df <- cbind(CI95_df,comp_CI)}
}


############## Column Naming
############################
clas <- rep(c(replicate(3,"ant"),replicate(12,"for")),15)
lake  <- c(rep("Eri",45),rep("Ont",45),rep("MiH",45),rep("Sup",45),rep("StC",45))
comp  <- c(rep(c(rep("P",15),rep("E",15),rep("R",15)),5))

months <- vector(length = 15)
for (i in (forecast_month-3):(forecast_month+11)) {
  months[i - forecast_month + 4] = month.abb[((i - 1) %% 12) + 1]
}
months <- rep(months, 12)

colnames(cal_df) <- paste(clas, months, lake, comp, sep = ".")
colnames(CI95_df) <- paste(clas, months, lake, comp, sep = ".")

#################################

#################################
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/formatted_input")
filename <- paste(month.abb[forecast_month], "_12Forecast_3Ant_CopulaInput.csv", sep = "")
write.csv(cal_df, filename)

##################################


##################################
setwd("~/INSERT_WORKING_DIRECTORY_HERE/Copula/R_objects")
avg_CI95 <- apply(CI95_df[1:70,],2,mean)
saveRDS(avg_CI95,"1950_2020_L2S_Avg_CI95_width.rds")
