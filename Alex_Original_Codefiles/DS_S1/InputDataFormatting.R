rm(list=ls())
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code")
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/L2SWBM_med_1950-2022")

months = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
L_forecast <- 12
L_ant      <- 3

EriP <- read.csv("EriePrecip_analysis19502022_prior19001969_1m.csv")[10:864,]
EriE <- read.csv("ErieEvap_analysis19502022_prior19001969_1m.csv")[10:864,]
EriR <- read.csv("ErieRunoff_analysis19502022_prior19001969_1m.csv")[10:864,]
OntP <- read.csv("OntPrecip_analysis19502022_prior19001969_1m.csv")[10:864,]
OntE <- read.csv("OntEvap_analysis19502022_prior19001969_1m.csv")[10:864,]
OntR <- read.csv("OntRunoff_analysis19502022_prior19001969_1m.csv")[10:864,]
MiHP <- read.csv("MiHurPrecip_analysis19502022_prior19001969_1m.csv")[10:864,]
MiHE <- read.csv("MiHurEvap_analysis19502022_prior19001969_1m.csv")[10:864,]
MiHR <- read.csv("MiHurRunoff_analysis19502022_prior19001969_1m.csv")[10:864,]
SupP <- read.csv("SupPrecip_analysis19502022_prior19001969_1m.csv")[10:864,]
SupE <- read.csv("SupEvap_analysis19502022_prior19001969_1m.csv")[10:864,]
SupR <- read.csv("SupRunoff_analysis19502022_prior19001969_1m.csv")[10:864,]

comps_inp = lapply(list(EriP, EriE, EriR, OntP, OntE, OntR, MiHP, MiHE, MiHR, SupP, SupE, SupR), function(x) {rownames(x) <- c(1:855); x})
names(comps_inp) <- c("EriP", "EriE", "EriR", "OntP", "OntE", "OntR", "MiHP", "MiHE", "MiHR", "SupP", "SupE", "SupR")

for (comp in 1:12){
  comps_inp[[comp]]$CI.95.W <- abs(comps_inp[[comp]]$X97.5.Percentile - comps_inp[[comp]]$X2.5.Percentile)
}

cal_dfs  <- vector("list",length=12)
CI95_dfs <- vector("list",length=12)
names(cal_dfs)  <- month.abb
names(CI95_dfs) <- month.abb

for (forecast_month in 1:12){
  start_month_indices <- which(comps_inp$EriP$Month[4:843]==forecast_month)+3
  cal_dfs[[forecast_month]] <- data.frame()
  
  for(lake_comp in 1:12){
    comp_cal <- data.frame()
    cal_in   <- comps_inp[[lake_comp]]$Median
    
    for (obs_year in 1:length(start_month_indices)){
      cal_row  <- cal_in[(start_month_indices[obs_year]-L_ant):(start_month_indices[obs_year]+L_forecast-1)]
      comp_cal <- rbind(comp_cal,cal_row)
    }
    
    if(ncol(cal_dfs[[forecast_month]])==0){cal_dfs[[forecast_month]] <- comp_cal}
    else                                  {cal_dfs[[forecast_month]] <- cbind(cal_dfs[[forecast_month]],comp_cal)}
  }
}

for (forecast_month in 1:12){
  start_month_indices <- which(comps_inp$EriP$Month[4:843]==forecast_month)+3
  CI95_dfs[[forecast_month]] <- data.frame()
  
  for(lake_comp in 1:12){
    comp_CI <- data.frame()
    CI_in   <- comps_inp[[lake_comp]]$CI.95.W
    
    for (obs_year in 1:length(start_month_indices)){
      CI_row  <- CI_in[(start_month_indices[obs_year]-L_ant):(start_month_indices[obs_year]+L_forecast-1)]
      comp_CI <- rbind(comp_CI,CI_row)
    }
    
    if(ncol(CI95_dfs[[forecast_month]])==0){CI95_dfs[[forecast_month]] <- comp_CI}
    else                                   {CI95_dfs[[forecast_month]] <- cbind(CI95_dfs[[forecast_month]],comp_CI)}
  }
}


############## Column Naming
############################
clas <- rep(c(replicate(3,"ant"),replicate(12,"for")),12)
lake  <- c(rep("Eri",45),rep("Ont",45),rep("MiH",45),rep("Sup",45))
comp  <- rep(c(rep("P",15),rep("E",15),rep("R",15)),4)

Jan_month <- rep(c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),12)
Feb_month <- rep(c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan"),12)
Mar_month <- rep(c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb"),12)
Apr_month <- rep(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar"),12)
May_month <- rep(c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr"),12)
Jun_month <- rep(c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May"),12)
Jul_month <- rep(c("Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"),12)
Aug_month <- rep(c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul"),12)
Sep_month <- rep(c("Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug"),12)
Oct_month <- rep(c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"),12)
Nov_month <- rep(c("Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"),12)
Dec_month <- rep(c("Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"),12)


colnames(cal_dfs$Jan) <- paste(clas, Jan_month, lake, comp, sep = ".")
colnames(cal_dfs$Feb) <- paste(clas, Feb_month, lake, comp, sep = ".")
colnames(cal_dfs$Mar) <- paste(clas, Mar_month, lake, comp, sep = ".")
colnames(cal_dfs$Apr) <- paste(clas, Apr_month, lake, comp, sep = ".")
colnames(cal_dfs$May) <- paste(clas, May_month, lake, comp, sep = ".")
colnames(cal_dfs$Jun) <- paste(clas, Jun_month, lake, comp, sep = ".")
colnames(cal_dfs$Jul) <- paste(clas, Jul_month, lake, comp, sep = ".")
colnames(cal_dfs$Aug) <- paste(clas, Aug_month, lake, comp, sep = ".")
colnames(cal_dfs$Sep) <- paste(clas, Sep_month, lake, comp, sep = ".")
colnames(cal_dfs$Oct) <- paste(clas, Oct_month, lake, comp, sep = ".")
colnames(cal_dfs$Nov) <- paste(clas, Nov_month, lake, comp, sep = ".")
colnames(cal_dfs$Dec) <- paste(clas, Dec_month, lake, comp, sep = ".")

colnames(CI95_dfs$Jan) <- paste(clas, Jan_month, lake, comp, sep = ".")
colnames(CI95_dfs$Feb) <- paste(clas, Feb_month, lake, comp, sep = ".")
colnames(CI95_dfs$Mar) <- paste(clas, Mar_month, lake, comp, sep = ".")
colnames(CI95_dfs$Apr) <- paste(clas, Apr_month, lake, comp, sep = ".")
colnames(CI95_dfs$May) <- paste(clas, May_month, lake, comp, sep = ".")
colnames(CI95_dfs$Jun) <- paste(clas, Jun_month, lake, comp, sep = ".")
colnames(CI95_dfs$Jul) <- paste(clas, Jul_month, lake, comp, sep = ".")
colnames(CI95_dfs$Aug) <- paste(clas, Aug_month, lake, comp, sep = ".")
colnames(CI95_dfs$Sep) <- paste(clas, Sep_month, lake, comp, sep = ".")
colnames(CI95_dfs$Oct) <- paste(clas, Oct_month, lake, comp, sep = ".")
colnames(CI95_dfs$Nov) <- paste(clas, Nov_month, lake, comp, sep = ".")
colnames(CI95_dfs$Dec) <- paste(clas, Dec_month, lake, comp, sep = ".")
#################################


#################################
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/formatted_input")
write.csv(cal_dfs$cal_df_Jan, "Jan_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Feb, "Feb_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Mar, "Mar_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Apr, "Apr_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_May, "May_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Jun, "Jun_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Jul, "Jul_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Aug, "Aug_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Sep, "Sep_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Oct, "Oct_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Nov, "Nov_12Forecast_3Ant_CopulaInput.csv")
write.csv(cal_dfs$cal_df_Dec, "Dec_12Forecast_3Ant_CopulaInput.csv")
##################################


##################################
setwd("Q:/Hydro/writing/in_prep/alex_copula/Code/R_objects")
avg_CI95s <- CI95_dfs
for (forecast_month in 1:12){
  avg_CI95s[[forecast_month]] <- apply(CI95_dfs[[forecast_month]][1:56,],2,mean)
}
saveRDS(avg_CI95s,"1950_2006_L2S_Avg_CI95_width.rds")
s