library("tidyverse")

old.data <- read.csv("~/INSERT_WORKING_DIRECTORY_HERE/output/GLHYD_data_metric.csv")[997:1248,3:7]

supdata <- read.csv("supforecast.csv", header = F)
mihurdata <- read.csv("mihurforecast.csv", header = F)
stclairdata <- read.csv("stclairforecast.csv", header = F)
eriedata <- read.csv("erieforecast.csv", header = F)

old.sup.linearMonthly <- as.numeric(old.data[,1])
old.mihur.linearMonthly <- as.numeric(old.data[,2])
old.stclair.linearMonthly <- as.numeric(old.data[,3])
old.erie.linearMonthly <- as.numeric(old.data[,4])

sup.linearMonthly <- data.frame(pivot_longer(supdata, colnames(supdata)[2:13], values_to = "water_level"))[,3]
mihur.linearMonthly <- data.frame(pivot_longer(mihurdata, colnames(mihurdata)[2:13], values_to = "water_level"))[,3]
stclair.linearMonthly <- data.frame(pivot_longer(stclairdata, colnames(stclairdata)[2:13], values_to = "water_level"))[,3]
erie.linearMonthly <- data.frame(pivot_longer(eriedata, colnames(eriedata)[2:13], values_to = "water_level"))[,3]

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

sup.dataAverage <- as.numeric(sum(c(old.sup.linearMonthly, sup.linearMonthly))/(length(old.sup.linearMonthly) + length(sup.linearMonthly)))
mihur.dataAverage <- as.numeric(sum(c(old.mihur.linearMonthly, mihur.linearMonthly))/(length(old.mihur.linearMonthly) + length(mihur.linearMonthly)))
stclair.dataAverage <- as.numeric(sum(c(old.stclair.linearMonthly, stclair.linearMonthly))/(length(old.stclair.linearMonthly) + length(stclair.linearMonthly)))
erie.dataAverage <- as.numeric(sum(c(old.erie.linearMonthly, erie.linearMonthly))/(length(old.erie.linearMonthly) + length(erie.linearMonthly)))

old.mihur.linearMonthly <- old.mihur.linearMonthly+4
old.mihur.linearAnnual <- old.mihur.linearAnnual+4
mihur.linearMonthly <- mihur.linearMonthly+4
mihur.linearAnnual <- mihur.linearAnnual+4
mihur.dataAverage <- mihur.dataAverage+4

old.stclair.linearMonthly <- old.stclair.linearMonthly+2
old.stclair.linearAnnual <- old.stclair.linearAnnual+2
stclair.linearMonthly <- stclair.linearMonthly+2
stclair.linearAnnual <- stclair.linearAnnual+2
stclair.dataAverage <- stclair.dataAverage+2

months <- seq.Date(from = as.Date("2000-01-01"), to = as.Date("2070-12-01"), by = "1 month")
oldmonths <- length(old.sup.linearMonthly)
newmonths <- length(sup.linearMonthly)
oldyears <- length(old.sup.linearAnnual)
newyears <- length(sup.linearAnnual)

png("waterLevels.png", width = 7, height = 7, units = "in", res = 300)
par(mar = c(5, 5, 1, 1))
par(xpd=FALSE)

plot.new()

plot(months, c(old.sup.linearMonthly, rep(NA, newmonths)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), sup.linearMonthly), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="skyblue", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.sup.linearAnnual, rep(NA, newmonths)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), sup.linearAnnual), axes=FALSE, 
     col="darkblue", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
segments(x0 = as.Date("2000-01-01"), x1 = as.Date("2070-12-01"), y0 = sup.dataAverage, col = "red", lwd = 0.85)
par(new=TRUE)

plot(months, c(old.mihur.linearMonthly, rep(NA, newmonths)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), mihur.linearMonthly), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="skyblue", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.mihur.linearAnnual, rep(NA, newmonths)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), mihur.linearAnnual), axes=FALSE, 
     col="darkblue", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
segments(x0 = as.Date("2000-01-01"), x1 = as.Date("2070-12-01"), y0 = mihur.dataAverage, col = "red", lwd = 0.85)
par(new=TRUE)


plot(months, c(old.stclair.linearMonthly, rep(NA, newmonths)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), stclair.linearMonthly), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="skyblue", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.stclair.linearAnnual, rep(NA, newmonths)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), stclair.linearAnnual), axes=FALSE, 
     col="darkblue", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
segments(x0 = as.Date("2000-01-01"), x1 = as.Date("2070-12-01"), y0 = stclair.dataAverage, col = "red", lwd = 0.85)
par(new=TRUE)


plot(months, c(old.erie.linearMonthly, rep(NA, newmonths)), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="lightgrey", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), erie.linearMonthly), type="p", xlab="Year", ylab="Water Surface Elevation (meters)", 
     col="skyblue", cex=0.3, pch=20, lwd=1, axes=FALSE, ylim = c(173.5, 184.5))
par(new=TRUE)
plot(months, c(old.erie.linearAnnual, rep(NA, newmonths)), axes=FALSE, 
     col="black", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
plot(months, c(rep(NA, oldmonths), erie.linearAnnual), axes=FALSE, 
     col="darkblue", cex=0.3, pch=20, ylim = c(173.5, 184.5), xlab="", ylab="")
par(new=TRUE)
segments(x0 = as.Date("2000-01-01"), x1 = as.Date("2070-12-01"), y0 = erie.dataAverage, col = "red", lwd = 0.85)
par(new=TRUE)

axis.Date(1, at = seq.Date(as.Date("2000-06-01"), as.Date("2070-06-01"), by = "10 years"),
          format = "%Y", cex.axis = 0.75)
axis.Date(1, at = seq.Date(as.Date("2000-06-01"), as.Date("2070-06-01"), by = "2 years"), 
          labels = FALSE, tck = -0.01)

axis.Date(3, at = seq.Date(as.Date("2000-06-01"), as.Date("2070-06-01"), by = "10 years"), 
          labels = FALSE, cex.axis = 0.75)
axis.Date(3, at = seq.Date(as.Date("2000-06-01"), as.Date("2070-06-01"), by = "2 years"), 
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