library("plotly")

setwd("~/INSERT_WORKING_DIRECTORY_HERE/output/LongTermWLs")
eri <- data.matrix(read.csv("Erie_LongTerm_PEChange_WL.csv", row.names = 1))
mih <- data.matrix(read.csv("MichiganHuron_LongTerm_PEChange_WL.csv", row.names = 1))
stc <- data.matrix(read.csv("StClair_LongTerm_PEChange_WL.csv", row.names = 1))
sup <- data.matrix(read.csv("Superior_LongTerm_PEChange_WL.csv", row.names = 1))

##### ERIE #####

wls <- c(eri)
precip <- rep(seq(0.86, 1.14, 0.02), times = 15)
evap <- rep(seq(0.86, 1.14, 0.02), each = 15)
model <- lm(wls ~ precip + evap)
model$fitted.values
fitted <- matrix(data = model$fitted.values, nrow = 15, ncol = 15)

fig <- plot_ly()
fig <- fig %>% add_markers(x = ~evap, 
                           y = ~precip, 
                           z = ~wls, 
                           marker = list(color = ~wls, colorscale = "Viridis"))
fig <- fig %>% add_surface(x = ~seq(0.86, 1.14, 0.02), 
                           y = ~seq(0.86, 1.14, 0.02), 
                           z = ~fitted,
                           opacity = 0.5)

axx <- list(title = "% Change in Evaporation")
axy <- list(title = "% Change in Precipitation")
axz <- list(title = "Long Term Avg. WL")

fig <- fig %>% layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))

fig

##### MICHIGAN-HURON #####

wls <- c(mih)
precip <- rep(seq(0.86, 1.14, 0.02), times = 15)
evap <- rep(seq(0.86, 1.14, 0.02), each = 15)
model <- lm(wls ~ precip + evap)
model$fitted.values
fitted <- matrix(data = model$fitted.values, nrow = 15, ncol = 15)

fig <- plot_ly()
fig <- fig %>% add_markers(x = ~evap, 
                           y = ~precip, 
                           z = ~wls, 
                           marker = list(color = ~wls, colorscale = "Viridis"))
fig <- fig %>% add_surface(x = ~seq(0.86, 1.14, 0.02), 
                           y = ~seq(0.86, 1.14, 0.02), 
                           z = ~fitted,
                           opacity = 0.5)

axx <- list(title = "% Change in Evaporation")
axy <- list(title = "% Change in Precipitation")
axz <- list(title = "Long Term Avg. WL")

fig <- fig %>% layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))

fig

##### ST. CLAIR #####

wls <- c(stc)
precip <- rep(seq(0.86, 1.14, 0.02), times = 15)
evap <- rep(seq(0.86, 1.14, 0.02), each = 15)
model <- lm(wls ~ precip + evap)
model$fitted.values
fitted <- matrix(data = model$fitted.values, nrow = 15, ncol = 15)

fig <- plot_ly()
fig <- fig %>% add_markers(x = ~evap, 
                           y = ~precip, 
                           z = ~wls, 
                           marker = list(color = ~wls, colorscale = "Viridis"))
fig <- fig %>% add_surface(x = ~seq(0.86, 1.14, 0.02), 
                           y = ~seq(0.86, 1.14, 0.02), 
                           z = ~fitted,
                           opacity = 0.5)

axx <- list(title = "% Change in Evaporation")
axy <- list(title = "% Change in Precipitation")
axz <- list(title = "Long Term Avg. WL")

fig <- fig %>% layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))

fig

##### SUPERIOR #####

wls <- c(sup)
precip <- rep(seq(0.86, 1.14, 0.02), times = 15)
evap <- rep(seq(0.86, 1.14, 0.02), each = 15)
model <- lm(wls ~ precip + evap)
model$fitted.values
fitted <- matrix(data = model$fitted.values, nrow = 15, ncol = 15)

fig <- plot_ly()
fig <- fig %>% add_markers(x = ~evap, 
                           y = ~precip, 
                           z = ~wls, 
                           marker = list(color = ~wls, colorscale = "Viridis"))
fig <- fig %>% add_surface(x = ~seq(0.86, 1.14, 0.02), 
                           y = ~seq(0.86, 1.14, 0.02), 
                           z = ~fitted,
                           opacity = 0.5)

axx <- list(title = "% Change in Evaporation")
axy <- list(title = "% Change in Precipitation")
axz <- list(title = "Long Term Avg. WL")

fig <- fig %>% layout(scene = list(xaxis = axx, yaxis = axy, zaxis = axz))

fig
