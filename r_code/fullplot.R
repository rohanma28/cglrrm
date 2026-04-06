library("tidyverse")
library("zoo")

# WATER LEVEL PLOT #############################################################

setwd("~/INSERT_WORKING_DIRECTORY_HERE/output")

old.data <- read.csv("GLHYD_data_metric.csv")[997:1248,3:7]

setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")

old.sup.linearMonthly <- as.numeric(old.data[,1])
old.mihur.linearMonthly <- as.numeric(old.data[,2])
old.stclair.linearMonthly <- as.numeric(old.data[,3])
old.erie.linearMonthly <- as.numeric(old.data[,4])

old.sup.linearAnnual <- rep(NA, length(old.sup.linearMonthly))
old.mihur.linearAnnual <- rep(NA, length(old.mihur.linearMonthly))
old.stclair.linearAnnual <- rep(NA, length(old.stclair.linearMonthly))
old.erie.linearAnnual <- rep(NA, length(old.erie.linearMonthly))

for (i in seq(1, length(old.sup.linearMonthly), 12)) {
  old.sup.linearAnnual[i+5] = mean(old.sup.linearMonthly[i:(i+11)])
  old.mihur.linearAnnual[i+5] = mean(old.mihur.linearMonthly[i:(i+11)])
  old.stclair.linearAnnual[i+5] = mean(old.stclair.linearMonthly[i:(i+11)])
  old.erie.linearAnnual[i+5] = mean(old.erie.linearMonthly[i:(i+11)])
}

old.mihur.linearMonthly <- old.mihur.linearMonthly+4
old.mihur.linearAnnual <- old.mihur.linearAnnual+4

old.stclair.linearMonthly <- old.stclair.linearMonthly+2
old.stclair.linearAnnual <- old.stclair.linearAnnual+2

oldmonths <- length(old.sup.linearMonthly)
oldyears <- length(old.sup.linearAnnual)

months <- seq.Date(from = as.Date("2000-01-01"), to = as.Date("2090-12-01"), by = "1 month")

png("waterLevels.png", width = 7, height = 7, units = "in", res = 300)
par(mar = c(5, 5, 1, 1))
par(xpd=FALSE)

plot.new()

plot(months, c(old.sup.linearMonthly, rep(NA, 840)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.sup.linearAnnual, rep(NA, 840)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(old.mihur.linearMonthly, rep(NA, 840)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.mihur.linearAnnual, rep(NA, 840)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(old.stclair.linearMonthly, rep(NA, 840)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.stclair.linearAnnual, rep(NA, 840)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(old.erie.linearMonthly, rep(NA, 840)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.erie.linearAnnual, rep(NA, 840)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)

sampsims <- sample.int(1000, 100)
sampbold <- sampsims[100]

for (sim in sampsims) {
  
  print(sim)
  
  setwd(paste("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE/", sim, sep = ""))
  
  supdata <- data.frame(lapply(read.csv("supforecast.csv", header = F)[2:71,], as.numeric))
  mihurdata <- data.frame(lapply(read.csv("mihurforecast.csv", header = F)[2:71,], as.numeric))
  stclairdata <- data.frame(lapply(read.csv("stclairforecast.csv", header = F)[2:71,], as.numeric))
  eriedata <- data.frame(lapply(read.csv("erieforecast.csv", header = F)[2:71,], as.numeric))
  
  setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")
  
  sup.linearMonthly <- data.frame(pivot_longer(supdata, colnames(supdata)[2:13], values_to = "water_level"))[,3]
  mihur.linearMonthly <- data.frame(pivot_longer(mihurdata, colnames(mihurdata)[2:13], values_to = "water_level"))[,3]
  stclair.linearMonthly <- data.frame(pivot_longer(stclairdata, colnames(stclairdata)[2:13], values_to = "water_level"))[,3]
  erie.linearMonthly <- data.frame(pivot_longer(eriedata, colnames(eriedata)[2:13], values_to = "water_level"))[,3]
  
  sup.linearAnnual <- rep(NA, length(sup.linearMonthly))
  mihur.linearAnnual <- rep(NA, length(mihur.linearMonthly))
  stclair.linearAnnual <- rep(NA, length(stclair.linearMonthly))
  erie.linearAnnual <- rep(NA, length(erie.linearMonthly))
  
  for (i in seq(1, length(sup.linearMonthly), 12)) {
    sup.linearAnnual[i+5] = mean(sup.linearMonthly[i:(i+11)])
    mihur.linearAnnual[i+5] = mean(mihur.linearMonthly[i:(i+11)])
    stclair.linearAnnual[i+5] = mean(stclair.linearMonthly[i:(i+11)])
    erie.linearAnnual[i+5] = mean(erie.linearMonthly[i:(i+11)])
  }
  
  mihur.linearMonthly <- mihur.linearMonthly+4
  mihur.linearAnnual <- mihur.linearAnnual+4
  
  stclair.linearMonthly <- stclair.linearMonthly+2
  stclair.linearAnnual <- stclair.linearAnnual+2
  
  plot(months, c(rep(NA, oldmonths), sup.linearMonthly), type="l", xlab="Year", ylab="Water Surface Elevation (meters)", 
       col=ifelse(sim == sampbold, "darkblue", "skyblue"), cex=0.3, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
  par(new=TRUE)
  
  plot(months, c(rep(NA, oldmonths), mihur.linearMonthly), type="l", xlab="Year", ylab="Water Surface Elevation (meters)", 
       col=ifelse(sim == sampbold, "darkblue", "skyblue"), cex=0.3, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
  par(new=TRUE)
  
  plot(months, c(rep(NA, oldmonths), stclair.linearMonthly), type="l", xlab="Year", ylab="Water Surface Elevation (meters)", 
       col=ifelse(sim == sampbold, "darkblue", "skyblue"), cex=0.3, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
  par(new=TRUE)
  
  plot(months, c(rep(NA, oldmonths), erie.linearMonthly), type="l", xlab="Year", ylab="Water Surface Elevation (meters)", 
       col=ifelse(sim == sampbold, "darkblue", "skyblue"), cex=0.3, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
  par(new=TRUE)
  
}

axis.Date(1, at = seq.Date(as.Date("2000-06-01"), as.Date("2090-06-01"), by = "10 years"),
          format = "%Y", cex.axis = 0.75)
axis.Date(1, at = seq.Date(as.Date("2000-06-01"), as.Date("2090-06-01"), by = "2 years"), 
          labels = FALSE, tck = -0.01)

axis.Date(3, at = seq.Date(as.Date("2000-06-01"), as.Date("2090-06-01"), by = "10 years"), 
          labels = FALSE, cex.axis = 0.75)
axis.Date(3, at = seq.Date(as.Date("2000-06-01"), as.Date("2090-06-01"), by = "2 years"), 
          labels = FALSE, tck = -0.01)

axis(2, at = seq(183, 184, 1.0), las = 1, cex.axis = 0.75)
axis(2, at = seq(182.8, 184.0, 0.2), labels = FALSE, tck = -0.01)

axis(2, at = seq(180, 181, 1.0), labels=176:177,las = 1, cex.axis = 0.75)
axis(2, at = seq(179.6, 181.6, 0.2), labels = FALSE, tck = -0.01)

axis(2, at = seq(177, 178, 1.0), labels=175:176, las = 1, cex.axis = 0.75)
axis(2, at = seq(176.4, 178.0, 0.2), labels = FALSE, tck = -0.01)

axis(2, at = seq(174, 175, 1.0), las = 1, cex.axis = 0.75)
axis(2, at = seq(173.6, 175.2, 0.2), labels = FALSE, tck = -0.01)

text(as.Date("1999-01-01"), 184.4, "Lake Superior", pos = 4, cex = 0.85)
text(as.Date("1999-01-01"), 182, "Lake Michigan and Huron", pos = 4, cex = 0.85)
text(as.Date("1999-01-01"), 178.8, "Lake St. Clair", pos = 4, cex = 0.85)
text(as.Date("1999-01-01"), 176, "Lake Erie", pos = 4, cex = 0.85)

box()

dev.off()

# COMPARISON PLOTS #############################################################

## HISTORICAL DATA EXTRACTION ##################################################

setwd("~/INSERT_WORKING_DIRECTORY_HERE/output")

old.data <- read.csv("GLHYD_data_metric.csv")[397:1248,3:7]

old.sup.linearMonthly <- as.numeric(old.data[,1])
old.mihur.linearMonthly <- as.numeric(old.data[,2])
old.stclair.linearMonthly <- as.numeric(old.data[,3])
old.erie.linearMonthly <- as.numeric(old.data[,4])

old.supfull <- c(old.sup.linearMonthly, rep(NA, 840))
old.mihurfull <- c(old.mihur.linearMonthly, rep(NA, 840))
old.stclairfull <- c(old.stclair.linearMonthly, rep(NA, 840))
old.eriefull <- c(old.erie.linearMonthly, rep(NA, 840))

zsupreal <- zoo(as.matrix(old.sup.linearMonthly, ncol = 1))
supsdreal <- rollapply(zsupreal, width = 12, FUN = sd)
supsdreal <- as.numeric(supsdreal)
supsdreal <- c(rep(NA, 6), supsdreal, rep(NA, 845))

zmihurreal <- zoo(as.matrix(old.mihur.linearMonthly, ncol = 1))
mihursdreal <- rollapply(zmihurreal, width = 12, FUN = sd)
mihursdreal <- as.numeric(mihursdreal)
mihursdreal <- c(rep(NA, 6), mihursdreal, rep(NA, 845))

zstclairreal <- zoo(as.matrix(old.stclair.linearMonthly, ncol = 1))
stclairsdreal <- rollapply(zstclairreal, width = 12, FUN = sd)
stclairsdreal <- as.numeric(stclairsdreal)
stclairsdreal <- c(rep(NA, 6), stclairsdreal, rep(NA, 845))

zeriereal <- zoo(as.matrix(old.erie.linearMonthly, ncol = 1))
eriesdreal <- rollapply(zeriereal, width = 12, FUN = sd)
eriesdreal <- as.numeric(eriesdreal)
eriesdreal <- c(rep(NA, 6), eriesdreal, rep(NA, 845))

supmaxreal <- max(old.sup.linearMonthly)
mihurmaxreal <- max(old.mihur.linearMonthly)
stclairmaxreal <- max(old.stclair.linearMonthly)
eriemaxreal <- max(old.erie.linearMonthly)

supminreal <- min(old.sup.linearMonthly)
mihurminreal <- min(old.mihur.linearMonthly)
stclairminreal <- min(old.stclair.linearMonthly)
erieminreal <- min(old.erie.linearMonthly)

supavgreal <- mean(old.sup.linearMonthly)
mihuravgreal <- mean(old.mihur.linearMonthly)
stclairavgreal <- mean(old.stclair.linearMonthly)
erieavgreal <- mean(old.erie.linearMonthly)

## BASELINE SIMULATION EXTRACTION ##############################################

supbase.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))
mihurbase.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))
stclairbase.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))
eriebase.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))

setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/baseline")

for (sim in 1:1000) {
  
  print(sim)
  
  setwd(paste("~/INSERT_WORKING_DIRECTORY_HERE/output/baseline/", sim, sep = ""))
  
  supdata <- data.frame(lapply(read.csv("supforecast.csv", header = F)[2:71,], as.numeric))
  mihurdata <- data.frame(lapply(read.csv("mihurforecast.csv", header = F)[2:71,], as.numeric))
  stclairdata <- data.frame(lapply(read.csv("stclairforecast.csv", header = F)[2:71,], as.numeric))
  eriedata <- data.frame(lapply(read.csv("erieforecast.csv", header = F)[2:71,], as.numeric))
  
  sup.temp <- data.frame(pivot_longer(supdata, colnames(supdata)[2:13], values_to = "water_level"))[,3]
  mihur.temp <- data.frame(pivot_longer(mihurdata, colnames(mihurdata)[2:13], values_to = "water_level"))[,3]
  stclair.temp <- data.frame(pivot_longer(stclairdata, colnames(stclairdata)[2:13], values_to = "water_level"))[,3]
  erie.temp <- data.frame(pivot_longer(eriedata, colnames(eriedata)[2:13], values_to = "water_level"))[,3]
  
  for (i in 1:length(sup.temp)) {
    supbase.Monthly[sim,i] = sup.temp[i]
    mihurbase.Monthly[sim,i] = mihur.temp[i]
    stclairbase.Monthly[sim,i] = stclair.temp[i]
    eriebase.Monthly[sim,i] = erie.temp[i]
  }
}

supbase.linearMonthly <- c(as.matrix(supbase.Monthly))
mihurbase.linearMonthly <- c(as.matrix(mihurbase.Monthly))
stclairbase.linearMonthly <- c(as.matrix(stclairbase.Monthly))
eriebase.linearMonthly <- c(as.matrix(eriebase.Monthly))

supbase.simframe <- matrix(supbase.linearMonthly, nrow = 1000, ncol = 840)
mihurbase.simframe <- matrix(mihurbase.linearMonthly, nrow = 1000, ncol = 840)
stclairbase.simframe <- matrix(stclairbase.linearMonthly, nrow = 1000, ncol = 840)
eriebase.simframe <- matrix(eriebase.linearMonthly, nrow = 1000, ncol = 840)

supbasefull <- cbind(matrix(nrow = 1000, ncol = 852), supbase.simframe)
mihurbasefull <- cbind(matrix(nrow = 1000, ncol = 852), mihurbase.simframe)
stclairbasefull <- cbind(matrix(nrow = 1000, ncol = 852), stclairbase.simframe)
eriebasefull <- cbind(matrix(nrow = 1000, ncol = 852), eriebase.simframe)

zsupbasesims <- zoo(t(supbase.simframe))
supbasesdsims <- rollapply(zsupbasesims, width = 12, FUN = sd)
supbasesdsims <- t(as.matrix(supbasesdsims))
supbasemediansd <- c(supbasesdsims)
supbasesdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), supbasesdsims, matrix(NA, nrow = 1000, ncol = 6))

zmihurbasesims <- zoo(t(mihurbase.simframe))
mihurbasesdsims <- rollapply(zmihurbasesims, width = 12, FUN = sd)
mihurbasesdsims <- t(as.matrix(mihurbasesdsims))
mihurbasemediansd <- c(mihurbasesdsims)
mihurbasesdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), mihurbasesdsims, matrix(NA, nrow = 1000, ncol = 6))

zstclairbasesims <- zoo(t(stclairbase.simframe))
stclairbasesdsims <- rollapply(zstclairbasesims, width = 12, FUN = sd)
stclairbasesdsims <- t(as.matrix(stclairbasesdsims))
stclairbasemediansd <- c(stclairbasesdsims)
stclairbasesdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), stclairbasesdsims, matrix(NA, nrow = 1000, ncol = 6))

zeriebasesims <- zoo(t(eriebase.simframe))
eriebasesdsims <- rollapply(zeriebasesims, width = 12, FUN = sd)
eriebasesdsims <- t(as.matrix(eriebasesdsims))
eriebasemediansd <- c(eriebasesdsims)
eriebasesdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), eriebasesdsims, matrix(NA, nrow = 1000, ncol = 6))

supbasemaxsims <- apply(supbase.simframe, 1, max)
mihurbasemaxsims <- apply(mihurbase.simframe, 1, max)
stclairbasemaxsims <- apply(stclairbase.simframe, 1, max)
eriebasemaxsims <- apply(eriebase.simframe, 1, max)

supbaseminsims <- apply(supbase.simframe, 1, min)
mihurbaseminsims <- apply(mihurbase.simframe, 1, min)
stclairbaseminsims <- apply(stclairbase.simframe, 1, min)
eriebaseminsims <- apply(eriebase.simframe, 1, min)

supbaseavgsims <- apply(supbase.simframe, 1, mean)
mihurbaseavgsims <- apply(mihurbase.simframe, 1, mean)
stclairbaseavgsims <- apply(stclairbase.simframe, 1, mean)
eriebaseavgsims <- apply(eriebase.simframe, 1, mean)

## CLIMATE CHANGE SIMULATION EXTRACTION ########################################

sup.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))
mihur.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))
stclair.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))
erie.Monthly <- data.frame(matrix(data = NA, nrow = 1000, ncol = 840))

setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")

for (sim in 1:1000) {
  
  print(sim)
  
  setwd(paste("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE/", sim, sep = ""))
  
  supdata <- data.frame(lapply(read.csv("supforecast.csv", header = F)[2:71,], as.numeric))
  mihurdata <- data.frame(lapply(read.csv("mihurforecast.csv", header = F)[2:71,], as.numeric))
  stclairdata <- data.frame(lapply(read.csv("stclairforecast.csv", header = F)[2:71,], as.numeric))
  eriedata <- data.frame(lapply(read.csv("erieforecast.csv", header = F)[2:71,], as.numeric))
  
  sup.temp <- data.frame(pivot_longer(supdata, colnames(supdata)[2:13], values_to = "water_level"))[,3]
  mihur.temp <- data.frame(pivot_longer(mihurdata, colnames(mihurdata)[2:13], values_to = "water_level"))[,3]
  stclair.temp <- data.frame(pivot_longer(stclairdata, colnames(stclairdata)[2:13], values_to = "water_level"))[,3]
  erie.temp <- data.frame(pivot_longer(eriedata, colnames(eriedata)[2:13], values_to = "water_level"))[,3]
  
  for (i in 1:length(sup.temp)) {
    sup.Monthly[sim,i] = sup.temp[i]
    mihur.Monthly[sim,i] = mihur.temp[i]
    stclair.Monthly[sim,i] = stclair.temp[i]
    erie.Monthly[sim,i] = erie.temp[i]
  }
  
}

sup.linearMonthly <- c(as.matrix(sup.Monthly))
mihur.linearMonthly <- c(as.matrix(mihur.Monthly))
stclair.linearMonthly <- c(as.matrix(stclair.Monthly))
erie.linearMonthly <- c(as.matrix(erie.Monthly))

sup.simframe <- matrix(sup.linearMonthly, nrow = 1000, ncol = 840)
mihur.simframe <- matrix(mihur.linearMonthly, nrow = 1000, ncol = 840)
stclair.simframe <- matrix(stclair.linearMonthly, nrow = 1000, ncol = 840)
erie.simframe <- matrix(erie.linearMonthly, nrow = 1000, ncol = 840)

supfull <- cbind(matrix(nrow = 1000, ncol = 852), sup.simframe)
mihurfull <- cbind(matrix(nrow = 1000, ncol = 852), mihur.simframe)
stclairfull <- cbind(matrix(nrow = 1000, ncol = 852), stclair.simframe)
eriefull <- cbind(matrix(nrow = 1000, ncol = 852), erie.simframe)

zsupsims <- zoo(t(sup.simframe))
supsdsims <- rollapply(zsupsims, width = 12, FUN = sd)
supsdsims <- t(as.matrix(supsdsims))
supmediansd <- c(supsdsims)
supsdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), supsdsims, matrix(NA, nrow = 1000, ncol = 6))

zmihursims <- zoo(t(mihur.simframe))
mihursdsims <- rollapply(zmihursims, width = 12, FUN = sd)
mihursdsims <- t(as.matrix(mihursdsims))
mihurmediansd <- c(mihursdsims)
mihursdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), mihursdsims, matrix(NA, nrow = 1000, ncol = 6))

zstclairsims <- zoo(t(stclair.simframe))
stclairsdsims <- rollapply(zstclairsims, width = 12, FUN = sd)
stclairsdsims <- t(as.matrix(stclairsdsims))
stclairmediansd <- c(stclairsdsims)
stclairsdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), stclairsdsims, matrix(NA, nrow = 1000, ncol = 6))

zeriesims <- zoo(t(erie.simframe))
eriesdsims <- rollapply(zeriesims, width = 12, FUN = sd)
eriesdsims <- t(as.matrix(eriesdsims))
eriemediansd <- c(eriesdsims)
eriesdsims <- cbind(matrix(NA, nrow = 1000, ncol = 857), eriesdsims, matrix(NA, nrow = 1000, ncol = 6))

supmaxsims <- apply(sup.simframe, 1, max)
mihurmaxsims <- apply(mihur.simframe, 1, max)
stclairmaxsims <- apply(stclair.simframe, 1, max)
eriemaxsims <- apply(erie.simframe, 1, max)

supminsims <- apply(sup.simframe, 1, min)
mihurminsims <- apply(mihur.simframe, 1, min)
stclairminsims <- apply(stclair.simframe, 1, min)
erieminsims <- apply(erie.simframe, 1, min)

supavgsims <- apply(sup.simframe, 1, mean)
mihuravgsims <- apply(mihur.simframe, 1, mean)
stclairavgsims <- apply(stclair.simframe, 1, mean)
erieavgsims <- apply(erie.simframe, 1, mean)

timeline <- seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-12-01"), by = "1 month")

## STANDARD DEVIATION PLOT #############################################################
setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")

plot_out = 'varchange.png'
if (file.exists(plot_out)) {
  file.remove(plot_out)
}
png(file = plot_out, width = 16, height = 9, units = "in", res = 300)

plot.new()
par(mfcol = c(2,4))
par(mar=c(0, 0, 0, 0))
par(oma=c(5, 5, 5, 5))
rows_to_plot <- sample.int(1000, 10)

plot(x = timeline, y = supsdreal, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(0, 0.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = supsdsims[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
     format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
      labels = FALSE, tck = -0.01)
axis(2, at = seq(0, 0.4, 0.1), cex.axis = 1.25)
axis(2, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
title(main = list('Superior', cex = 1.5), line = -2, xpd=NA)
title(ylab = list(substitute(paste(bold('Standard Deviation (m)'))), cex = 1.5), line = 2.5, xpd=NA)

plot(density(supmediansd, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', xlim = c(0, 0.4))
polygon(density(supmediansd, na.rm = TRUE), col=alpha("lightgrey", 0.5))
lines(density(supbasemediansd, na.rm = TRUE))
polygon(density(supbasemediansd, na.rm = TRUE), col=alpha("lightpink", 0.5))
abline(v=median(supmediansd, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(supbasemediansd, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(0, 0.4, 0.1), cex.axis = 1.25)
axis(1, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
title(xlab = list(substitute(paste(bold('Standard Deviation (m)'))), cex = 1.5), line = 3, xpd=NA)
box()

plot(x = timeline, y = mihursdreal, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(0, 0.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = mihursdsims[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(0, 0.4, 0.1), labels = FALSE)
axis(2, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
title(main = list('Michigan-Huron', cex = 1.5), line = -2, xpd=NA)

plot(density(mihurmediansd, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', xlim = c(0, 0.4))
polygon(density(mihurmediansd, na.rm = TRUE), col=alpha("lightgrey",0.5))
lines(density(mihurbasemediansd, na.rm = TRUE))
polygon(density(mihurbasemediansd, na.rm= TRUE), col=alpha("lightpink",0.5))
abline(v=median(mihurmediansd, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(mihurbasemediansd, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(0, 0.4, 0.1), labels = FALSE)
axis(1, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
box()


plot(x = timeline, y = stclairsdreal, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(0, 0.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = stclairsdsims[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(0, 0.4, 0.1), labels = FALSE)
axis(2, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
title(main = list('St. Clair', cex = 1.5), line = -2, xpd=NA)

plot(density(stclairmediansd, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', xlim = c(0, 0.4))
polygon(density(stclairmediansd, na.rm = TRUE), col=alpha("lightgrey", 0.5))
lines(density(stclairbasemediansd, na.rm = TRUE))
polygon(density(stclairbasemediansd, na.rm = TRUE), col=alpha("lightpink", 0.5))
abline(v=median(stclairmediansd, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(stclairbasemediansd, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(0, 0.4, 0.1), cex.axis = 1.25)
axis(1, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
box()


plot(x = timeline, y = eriesdreal, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(0, 0.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = eriesdsims[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(0, 0.4, 0.1), labels = FALSE)
axis(2, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
title(main = list('Erie', cex = 1.5), line = -2, xpd=NA)

plot(density(eriemediansd, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', xlim = c(0, 0.4))
polygon(density(eriemediansd, na.rm = TRUE), col=alpha("lightgrey",0.5))
lines(density(eriebasemediansd, na.rm = TRUE))
polygon(density(eriebasemediansd, na.rm = TRUE), col=alpha("lightpink",0.5))
abline(v=median(eriemediansd, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(eriebasemediansd, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(0, 0.4, 0.1), labels = FALSE)
axis(1, at = seq(0, 0.4, 0.02), labels = FALSE, tck = -0.01)
box()

dev.off()

## MAXIMUM PLOT #############################################################
setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")

plot_out = 'maxchange.png'
if (file.exists(plot_out)) {
  file.remove(plot_out)
}
png(file = plot_out, width = 16, height = 9, units = "in", res = 300)

plot.new()
par(mfcol = c(2,4))
par(mar=c(0, 0, 0, 0))
par(oma=c(5, 5, 5, 5))
rows_to_plot <- sample.int(1000, 10)

plot(x = timeline, y = old.supfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(182.4, 184.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = supfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(183, 184, 1), cex.axis = 1.25)
axis(2, at = seq(182.4, 184.4, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Superior', cex = 1.5), line = -2, xpd=NA)
title(ylab = list(substitute(paste(bold('Water Level (m)'))), cex = 1.5), line = 2.5, xpd=NA)

plot(density(supbasemaxsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(183.4, 184.4), ylim = c(0, 5))
polygon(density(supbasemaxsims, na.rm = TRUE), col=alpha("lightpink", 0.5))
lines(density(supmaxsims, na.rm = TRUE))
polygon(density(supmaxsims, na.rm = TRUE), col=alpha("lightgrey", 0.5))
abline(v=median(supmaxsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(supbasemaxsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(183, 184, 1), cex.axis = 1.25)
axis(1, at = seq(183.4, 184.4, 0.2), labels = FALSE, tck = -0.01)
title(xlab = list(substitute(paste(bold('Maximum Water Level (m)'))), cex = 1.5), line = 3, xpd=NA)
box()

plot(x = timeline, y = old.mihurfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(174.8, 178.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = mihurfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(175, 178, 1), cex.axis = 1.25)
axis(2, at = seq(174.8, 178.4, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Michigan-Huron', cex = 1.5), line = -2, xpd=NA)

plot(density(mihurbasemaxsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(176.6, 178.4), ylim = c(0, 5))
polygon(density(mihurbasemaxsims, na.rm = TRUE), col=alpha("lightpink",0.5))
lines(density(mihurmaxsims, na.rm = TRUE))
polygon(density(mihurmaxsims, na.rm= TRUE), col=alpha("lightgrey",0.5))
abline(v=median(mihurmaxsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(mihurbasemaxsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(177, 178, 1), cex.axis = 1.25)
axis(1, at = seq(176.6, 178.4, 0.2), labels = FALSE, tck = -0.01)
box()

plot(x = timeline, y = old.stclairfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(173.8, 176.6))
for (row in rows_to_plot) {
  lines(x = timeline, y = stclairfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(174, 176, 1), cex.axis = 1.25)
axis(2, at = seq(173.8, 176.6, 0.2), labels = FALSE, tck = -0.01)
title(main = list('St. Clair', cex = 1.5), line = -2, xpd=NA)

plot(density(stclairbasemaxsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(175.2, 176.6), ylim = c(0, 5))
polygon(density(stclairbasemaxsims, na.rm = TRUE), col=alpha("lightpink", 0.5))
lines(density(stclairmaxsims, na.rm = TRUE))
polygon(density(stclairmaxsims, na.rm = TRUE), col=alpha("lightgrey", 0.5))
abline(v=median(stclairmaxsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(stclairbasemaxsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(175, 176, 1), cex.axis = 1.25)
axis(1, at = seq(175.2, 176.6, 0.2), labels = FALSE, tck = -0.01)
box()

plot(x = timeline, y = old.eriefull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(173, 175.8))
for (row in rows_to_plot) {
  lines(x = timeline, y = eriefull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(173, 175, 1), cex.axis = 1.25)
axis(2, at = seq(173, 175.8, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Erie', cex = 1.5), line = -2, xpd=NA)

plot(density(eriebasemaxsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(174.4, 175.8), ylim = c(0, 5))
polygon(density(eriebasemaxsims, na.rm = TRUE), col=alpha("lightpink",0.5))
lines(density(eriemaxsims, na.rm = TRUE))
polygon(density(eriemaxsims, na.rm = TRUE), col=alpha("lightgrey",0.5))
abline(v=median(eriemaxsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(eriebasemaxsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(173, 175, 1), cex.axis = 1.25)
axis(1, at = seq(174.4, 175.8, 0.2), labels = FALSE, tck = -0.01)
box()

dev.off()

## AVERAGE PLOT #############################################################
setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")

plot_out = 'avgchange.png'
if (file.exists(plot_out)) {
  file.remove(plot_out)
}
png(file = plot_out, width = 16, height = 9, units = "in", res = 300)

plot.new()
par(mfcol = c(2,4))
par(mar=c(0, 0, 0, 0))
par(oma=c(5, 5, 5, 5))
rows_to_plot <- sample.int(1000, 10)

plot(x = timeline, y = old.supfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(182.4, 184.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = supfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(183, 184, 1), cex.axis = 1.25)
axis(2, at = seq(182.4, 184.4, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Superior', cex = 1.5), line = -2, xpd=NA)
title(ylab = list(substitute(paste(bold('Water Level (m)'))), cex = 1.5), line = 2.5, xpd=NA)

plot(density(supbaseavgsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(182.4, 184.4), ylim = c(0, 5.5))
polygon(density(supbaseavgsims, na.rm = TRUE), col=alpha("lightpink", 0.5))
lines(density(supavgsims, na.rm = TRUE))
polygon(density(supavgsims, na.rm = TRUE), col=alpha("lightgrey", 0.5))
abline(v=median(supavgsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(supbaseavgsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(183, 184, 1), cex.axis = 1.25)
axis(1, at = seq(182.4, 184.4, 0.2), labels = FALSE, tck = -0.01)
title(xlab = list(substitute(paste(bold('Average Water Level (m)'))), cex = 1.5), line = 3, xpd=NA)
box()

plot(x = timeline, y = old.mihurfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(174.8, 178.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = mihurfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(175, 178, 1), cex.axis = 1.25)
axis(2, at = seq(174.8, 178.4, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Michigan-Huron', cex = 1.5), line = -2, xpd=NA)

plot(density(mihurbaseavgsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(174.8, 178.4), ylim = c(0, 5))
polygon(density(mihurbaseavgsims, na.rm = TRUE), col=alpha("lightpink",0.5))
lines(density(mihuravgsims, na.rm = TRUE))
polygon(density(mihuravgsims, na.rm= TRUE), col=alpha("lightgrey",0.5))
abline(v=median(mihuravgsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(mihurbaseavgsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(175, 178, 1), cex.axis = 1.25)
axis(1, at = seq(174.8, 178.4, 0.2), labels = FALSE, tck = -0.01)
box()

plot(x = timeline, y = old.stclairfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(173.8, 176.6))
for (row in rows_to_plot) {
  lines(x = timeline, y = stclairfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(174, 176, 1), cex.axis = 1.25)
axis(2, at = seq(173.8, 176.6, 0.2), labels = FALSE, tck = -0.01)
title(main = list('St. Clair', cex = 1.5), line = -2, xpd=NA)

plot(density(stclairbaseavgsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(173.8, 176.6), ylim = c(0, 5))
polygon(density(stclairbaseavgsims, na.rm = TRUE), col=alpha("lightpink", 0.5))
lines(density(stclairavgsims, na.rm = TRUE))
polygon(density(stclairavgsims, na.rm = TRUE), col=alpha("lightgrey", 0.5))
abline(v=median(stclairavgsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(stclairbaseavgsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(174, 176, 1), cex.axis = 1.25)
axis(1, at = seq(173.8, 176.6, 0.2), labels = FALSE, tck = -0.01)
box()

plot(x = timeline, y = old.eriefull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(173, 175.8))
for (row in rows_to_plot) {
  lines(x = timeline, y = eriefull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(173, 175, 1), cex.axis = 1.25)
axis(2, at = seq(173, 175.8, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Erie', cex = 1.5), line = -2, xpd=NA)

plot(density(eriebaseavgsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(173, 175.8), ylim = c(0, 5))
polygon(density(eriebaseavgsims, na.rm = TRUE), col=alpha("lightpink",0.5))
lines(density(erieavgsims, na.rm = TRUE))
polygon(density(erieavgsims, na.rm = TRUE), col=alpha("lightgrey",0.5))
abline(v=median(erieavgsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(eriebaseavgsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(173, 175, 1), cex.axis = 1.25)
axis(1, at = seq(173, 175.8, 0.2), labels = FALSE, tck = -0.01)
box() 

dev.off()

## MINIMUM PLOT #############################################################
setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/INSERT_SIM_NAME_HERE")

plot_out = 'minchange.png'
if (file.exists(plot_out)) {
  file.remove(plot_out)
}
png(file = plot_out, width = 16, height = 9, units = "in", res = 300)

plot.new()
par(mfcol = c(2,4))
par(mar=c(0, 0, 0, 0))
par(oma=c(5, 5, 5, 5))
rows_to_plot <- sample.int(1000, 10)

plot(x = timeline, y = old.supfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(182.4, 184.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = supfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(183, 184, 1), cex.axis = 1.25)
axis(2, at = seq(182.4, 184.4, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Superior', cex = 1.5), line = -2, xpd=NA)
title(ylab = list(substitute(paste(bold('Water Level (m)'))), cex = 1.5), line = 2.5, xpd=NA)

plot(density(supbaseminsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(182.4, 183.4), ylim = c(0, 5.5))
polygon(density(supbaseminsims, na.rm = TRUE), col=alpha("lightpink", 0.5))
lines(density(supminsims, na.rm = TRUE))
polygon(density(supminsims, na.rm = TRUE), col=alpha("lightgrey", 0.5))
abline(v=median(supminsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(supbaseminsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(182, 183, 1), cex.axis = 1.25)
axis(1, at = seq(182.4, 183.4, 0.2), labels = FALSE, tck = -0.01)
title(xlab = list(substitute(paste(bold('Minimum Water Level (m)'))), cex = 1.5), line = 3, xpd=NA)
box()

plot(x = timeline, y = old.mihurfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(174.8, 178.4))
for (row in rows_to_plot) {
  lines(x = timeline, y = mihurfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(175, 178, 1), cex.axis = 1.25)
axis(2, at = seq(174.8, 178.4, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Michigan-Huron', cex = 1.5), line = -2, xpd=NA)

plot(density(mihurbaseminsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(174.8, 176.6), ylim = c(0, 5))
polygon(density(mihurbaseminsims, na.rm = TRUE), col=alpha("lightpink",0.5))
lines(density(mihurminsims, na.rm = TRUE))
polygon(density(mihurminsims, na.rm= TRUE), col=alpha("lightgrey",0.5))
abline(v=median(mihurminsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(mihurbaseminsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(174, 177, 1), cex.axis = 1.25)
axis(1, at = seq(174.8, 176.6, 0.2), labels = FALSE, tck = -0.01)
box()

plot(x = timeline, y = old.stclairfull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(173.8, 176.6))
for (row in rows_to_plot) {
  lines(x = timeline, y = stclairfull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          format = "%Y", cex.axis = 1.25)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(174, 176, 1), cex.axis = 1.25)
axis(2, at = seq(173.8, 176.6, 0.2), labels = FALSE, tck = -0.01)
title(main = list('St. Clair', cex = 1.5), line = -2, xpd=NA)

plot(density(stclairbaseminsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(173.8, 175.2), ylim = c(0, 5))
polygon(density(stclairbaseminsims, na.rm = TRUE), col=alpha("lightpink", 0.5))
lines(density(stclairminsims, na.rm = TRUE))
polygon(density(stclairminsims, na.rm = TRUE), col=alpha("lightgrey", 0.5))
abline(v=median(stclairminsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(stclairbaseminsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(173, 176, 1), cex.axis = 1.25)
axis(1, at = seq(173.8, 175.2, 0.2), labels = FALSE, tck = -0.01)
box()

plot(x = timeline, y = old.eriefull, type = 'l', col = "red", xaxt = 'n', yaxt = 'n', ylim = c(173, 175.8))
for (row in rows_to_plot) {
  lines(x = timeline, y = eriefull[row,])
}
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "40 years"), 
          labels = FALSE)
axis.Date(3, at = seq.Date(from = as.Date("1950-01-01"), to = as.Date("2090-01-01"), by = "10 years"), 
          labels = FALSE, tck = -0.01)
axis(2, at = seq(173, 175, 1), cex.axis = 1.25)
axis(2, at = seq(173, 175.8, 0.2), labels = FALSE, tck = -0.01)
title(main = list('Erie', cex = 1.5), line = -2, xpd=NA)

plot(density(eriebaseminsims, na.rm = TRUE), main = " ", xaxt = 'n', yaxt = 'n', 
     xlim = c(173, 174.4), ylim = c(0, 5))
polygon(density(eriebaseminsims, na.rm = TRUE), col=alpha("lightpink",0.5))
lines(density(erieminsims, na.rm = TRUE))
polygon(density(erieminsims, na.rm = TRUE), col=alpha("lightgrey",0.5))
abline(v=median(erieminsims, na.rm = TRUE), col = "black", lwd = 2)
abline(v=median(eriebaseminsims, na.rm = TRUE), col = "red", lwd= 2)
axis(1, at = seq(173, 175, 1), cex.axis = 1.25)
axis(1, at = seq(173, 174.4, 0.2), labels = FALSE, tck = -0.01)
box() 

dev.off()
