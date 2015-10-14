## Install packages for the tutorial
#install.packages(c("demography","StMoMo","rgl","fanplot", "ggplo2", "gridExtra", "reshape2"))

## Install own package from Triple A repository on GitHub
# library(devtools)
# devtools::install_github("TARF/insureR")

## Load required libraries
library(demography)
library(StMoMo)
library(rgl)
library(fanplot)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Own package
library(insureR)

## Source functions
#source("Case_study/temp_functions.R")

# skDemo<-hmd.mx("NLD", username=username, password=password)
load("Case_study/nlDemo.RData")
years <- 1950:2012
ages<- skDemo$age
Dxt <- skDemo$rate[[3]] * skDemo$pop[[3]]
E0xt <- skDemo$pop[[3]] + 0.5 * Dxt
Ecxt <- skDemo$pop[[3]]
Dxt <- Dxt[,as.character(1950:2012)]
E0xt <- E0xt[,as.character(1950:2012)]
qxt <- log(Dxt/E0xt)

## get data
forecastTime <- 120
ages.fit <- 45:90

## Total polulation 3d plot
persp3d(ages[0:90], years, qxt[0:90,], col="green", shade=TRUE,xlab="Ages (0-90)", ylab="Years",zlab="Mortality rate (log)")

## Male 3d plot
Dxtm <- skDemo$rate$male * skDemo$pop$male
E0xtm <- skDemo$pop$male + 0.5 * Dxtm
Dxtm  <- Dxtm[,as.character(1950:2012)]
E0xtm <- E0xtm[,as.character(1950:2012)]
qxtm <- log(Dxtm/E0xtm)

persp3d(ages[0:90], years, qxtm[0:90,], col="skyblue", shade=TRUE,xlab="Ages (0-90)",
        ylab="Years",zlab="Mortality rate (log)")

## Female 3d plot
Dxtf <- skDemo$rate$female * skDemo$pop$female
E0xtf <- skDemo$pop$female + 0.5 * Dxtf
Dxtf  <- Dxtf[,as.character(1950:2012)]
E0xtf <- E0xtf[,as.character(1950:2012)]
qxtf <- log(Dxtf/E0xtf)

persp3d(ages[0:90], years, qxtf[0:90,], col="pink", shade=TRUE,xlab="Ages (0-90)",
        ylab="Years",zlab="Mortality rate (log)")

## Setup weights for fitting
weights <- genWeightMat(ages.fit, years, 3)
  
## Modelling
modelFit <- array(data = NA, c(8, 2))
colnames(modelFit) <- c("AIC", "BIC")
modelNames <- c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT")
rownames(modelFit) <- modelNames 

## LC model under a Binomial setting - M1
LC <- lc(link = "logit")
LCfit <- fit(LC, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

## plot Lee Carter model fit
plot(LCfit, nCol = 3)

## get residual fit
LCres <- residuals(LCfit)
plot(LCres)
plot(LCres, type = "colourmap", reslim = c(-3.5, 3.5))

modelFit[1, 1] <- AIC(LCfit)
modelFit[1, 2] <- BIC(LCfit)

# Forecast
LCfor <- forecast(LCfit, h = forecastTime)
LCqxt <- cbind(LCfor$fitted, LCfor$rates)

## APC Currie (2006) - M3
APC <- apc("logit")
APCfit <- fit(APC, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights, start.ax = LCfit$ax,
              start.bx = LCfit$bx, start.kt = LCfit$kx)

modelFit[2, 1] <- AIC(APCfit)
modelFit[2, 2] <- BIC(APCfit)

APCres <- residuals(APCfit)
plot(APCres)
plot(APCres, type = "colourmap", reslim = c(-3.5, 3.5))

APCfor <- forecast(APCfit, h = forecastTime, gc.order = c(1, 1, 0))
APCqxt <- cbind(APCfor$fitted, APCfor$rates)

## RH
RH <- rh("logit", cohortAgeFun = "1")
RHfit <- fit(RH, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights, start.ax = LCfit$ax,
              start.bx = LCfit$bx, start.kt = LCfit$kx)

modelFit[3, 1] <- AIC(RHfit)
modelFit[3, 2] <- BIC(RHfit)

RHres <- residuals(RHfit)
plot(RHres)
plot(RHres, type = "colourmap", reslim = c(-3.5, 3.5))

RHfor <- forecast(RHfit, h = forecastTime, gc.order = c(1, 1, 0))
RHqxt <- cbind(RHfor$fitted, RHfor$rates)

## CBD model Cairns (2009) under a Binomial distribution of deaths Haberman and Renshaw (2011) - M5
CBD <- cbd()
CBDfit <- fit(CBD, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[4, 1] <- AIC(CBDfit)
modelFit[4, 2] <- BIC(CBDfit)

CBDres <- residuals(CBDfit)
plot(CBDres)
plot(CBDres, type = "colourmap", reslim = c(-3.5, 3.5))

CBDfor <- forecast(CBDfit, h = forecastTime)
CBDqxt <- cbind(CBDfor$fitted, CBDfor$rates)

## M6
M6 <- m6()
M6fit <- fit(M6, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[5, 1] <- AIC(M6fit)
modelFit[5, 2] <- BIC(M6fit)

M6res <- residuals(M6fit)
plot(M6res)
plot(M6res, type = "colourmap", reslim = c(-3.5, 3.5))

M6for <- forecast(M6fit, h = forecastTime, gc.order = c(2, 0, 0))
M6qxt <- cbind(M6for$fitted, M6for$rates)

## M7 under Binomial setting
M7 <- m7(link = "logit")
M7fit <- fit(M7, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[6, 1] <- AIC(M7fit)
modelFit[6, 2] <- BIC(M7fit)

M7res <- residuals(M7fit)
plot(M7res)
plot(M7res, type = "colourmap", reslim = c(-3.5, 3.5))

M7for <- forecast(M7fit, h = forecastTime, gc.order = c(2, 0, 0))
M7qxt <- cbind(M7for$fitted, M7for$rates)

## M8
M8 <- m8(link = "logit", xc = 65)
M8fit <- fit(M8, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)
modelFit[7, 1] <- AIC(M8fit)
modelFit[7, 2] <- BIC(M8fit)

M8res <- residuals(M8fit)
plot(M8res)
plot(M8res, type = "colourmap", reslim = c(-3.5, 3.5))

M8for <- forecast(M8fit, h = forecastTime, gc.order = c(2, 0, 0))
M8qxt <- cbind(M8for$fitted, M8for$rates)

## PLAT
f2 <- function(x, ages) mean(ages) - x
constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
  nYears <- dim(wxt)[2]
  x <- ages
  t <- 1:nYears
  c <- (1 - tail(ages, 1)):(nYears - ages[1])
  xbar <- mean(x)
  #\sum g(c)=0, \sum cg(c)=0, \sum c^2g(c)=0
  phiReg <- lm(gc ~ 1 + c + I(c^2), na.action = na.omit)
  
  phi <- coef(phiReg)
  gc <- gc - phi[1] - phi[2] * c - phi[3] * c^2
  kt[2, ] <- kt[2, ] + 2 * phi[3] * t
  kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t^2 - 2 * xbar * t)
  ax <- ax + phi[1] - phi[2] * x + phi[3] * x^2
  #\sum kt[i, ] = 0
  ci <- rowMeans(kt, na.rm = TRUE)
  ax <- ax + ci[1] + ci[2] * (xbar - x)
  kt[1, ] <- kt[1, ] - ci[1]
  kt[2, ] <- kt[2, ] - ci[2]
  list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
}
PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE,
               periodAgeFun = c("1", f2), cohortAgeFun = "1",
               constFun = constPlat)
PLATfit <- fit(PLAT, Dxt = Dxt, Ext= E0xt, ages = ages, years = years, ages.fit = ages.fit, wxt = weights)

modelFit[8, 1] <- AIC(PLATfit)
modelFit[8, 2] <- BIC(PLATfit)

PLATres <- residuals(PLATfit)
plot(PLATres)
plot(PLATres, type = "colourmap", reslim = c(-3.5, 3.5))

PLATfor <- forecast(PLATfit, h = forecastTime, gc.order = c(2, 0, 0))
PLATqxt <- cbind(PLATfor$fitted, PLATfor$rates)

## Collect fitted models for simulation purposes
modelsFitted <- list(LC = LCfit, APC = APCfit,  RH = RHfit, 
               CBD = CBDfit,  M6 = M6fit, M7 = M7fit, M8 = M8fit, PLAT = PLATfit)

## Model fit criteria AIC and BIC
modelFit

## Ploting - inspection of forecasted mortality rates for 65 year old
## Plots
years_chart <- c(years, (years[length(years)]+1):(years[length(years)]+forecastTime))
plot_df <- data.frame(years = years_chart, LC = LCqxt["65",], CBD = CBDqxt["65",], APC = APCqxt["65",], RH = RHqxt["65",],
                      M6 = M6qxt["65",], M7 = M7qxt["65",], M8 = M8qxt["65",], PLAT = PLATqxt["65",])
plot_df2 <- data.frame(years=years, Observed = (Dxt/E0xt)["65", ])
ncols <- ncol(plot_df)
plot_df <- melt(plot_df, id="years")
plot_q <- ggplot(data=plot_df, aes(x=years, y=value, colour=variable)) +
  geom_line(na.rm=TRUE, alpha=1, size=1) + 
  geom_point(data=plot_df2, aes(x=years, y=Observed), group=1, pch=21, colour="black", size=3, fill="red") +
  theme_minimal() +
  xlab("Years") + ylab("Mortality rates at age 65") + ggtitle("Mortality projection") +
  theme(panel.grid.major.x=element_blank(), axis.line = element_line(colour = "black", size =0.5), axis.title.x=element_text(family="sans",size=rel(1)),
        axis.title.y=element_text(family="sans",size=rel(1), angle=90), plot.title = element_text(family="sans", size=20, face="bold", vjust=2))+
  coord_cartesian(ylim = c(0, 0.023)) 

plot_q

## Extrapolate data to omega age = 120
LCextrapolate <- kannistoExtrapolation(LCqxt, ages.fit, years_chart)
LCqxtExtr <- LCextrapolate$qxt

APCextrapolate <- kannistoExtrapolation(APCqxt, ages.fit, years_chart)
APCqxtExtr <- APCextrapolate$qxt

RHextrapolate <- kannistoExtrapolation(RHqxt, ages.fit, years_chart)
RHqxtExtr <- RHextrapolate$qxt

CBDextrapolate <- kannistoExtrapolation(CBDqxt, ages.fit, years_chart)
CBDqxtExtr <- CBDextrapolate$qxt

M6extrapolate <- kannistoExtrapolation(M6qxt, ages.fit, years_chart)
M6qxtExtr <- M6extrapolate$qxt

M7extrapolate <- kannistoExtrapolation(M7qxt, ages.fit, years_chart)
M7qxtExtr <- M7extrapolate$qxt

M8extrapolate <- kannistoExtrapolation(M8qxt, ages.fit, years_chart)
M8qxtExtr <- M8extrapolate$qxt

extrapolate <- kannistoExtrapolation(PLATqxt, ages.fit, years_chart)
PLATqxtExtr <- extrapolate$qxt

## Read in AG table
ag <- read.csv("Case_study/ag-total.csv", row.names = 1)
colnames(ag) <- as.character(2014:(2014+ncol(ag)-1))
ag <- as.matrix(ag[46:121,], rownames.force=T)

models <- list(LCqxtExtr = LCqxtExtr, APCqxtExtr = APCqxtExtr, RHqxtExtr = RHqxtExtr, 
               CBDqxtExtr = CBDqxtExtr,  M6qxtExtr = M6qxtExtr, M7qxtExtr = M7qxtExtr,
               M8qxtExtr = M8qxtExtr, PLATqxtExtr = PLATqxtExtr, AG = ag)

## Annuity projection
## Assumption annuity due deffered at 65

## Read in portfolio data
portfolio <- read.csv('Case_study/portfolio.csv')

## Read in experience data
experience.factors <- read.csv('Case_study/experience-factors.csv')

ages.fit <- 45:120
valyear <- 2015
pension <- 1000

# calculate ages of insured
portfolio$age <- valyear - portfolio$YoB
experience.factors$total <- (experience.factors$Male + experience.factors$Female)/2

expF <- experience.factors$total[ages.fit]

BEL <- array(NA, c(9,1))
rownames(BEL) <- c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT", "AG")
colnames(BEL) <- "BEL"

for (m in 1:length(models)){
  output <- list()
  output2 <- list()
  for (i in 1:nrow(portfolio)){
    output[[i]] <- PVcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, 
                              valyear = valyear, ir = 0.02, type = "annuity")*pension
    output2[[i]] <- PVcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, 
                               valyear = valyear, ir = 0.02, type = "premium")*portfolio$Premium[i]
  }
  BEL[m, 1] <- do.call(sum, output)-do.call(sum, output2)
}

## Show BEL per model
df = data.frame(models = factor(c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT", "AG"), levels=c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT", "AG")),
                BEL = BEL)

ggplot(data=df, aes(x=models, y=BEL, fill=models)) +
  geom_bar(stat="identity") + theme_minimal() +
  xlab("Models") + ylab("Amount in EUR") + ggtitle("BEL for a portfolio") +
  theme(panel.grid.major.x=element_blank(), axis.line = element_line(colour = "black", size =0.5), legend.position="none", axis.title.x=element_text(family="sans",size=rel(1)),
        axis.title.y=element_text(family="sans",size=rel(1), angle=90), plot.title = element_text(family="sans", size=20, face="bold", vjust=2))

## SIMULATIONS for SCR calculation purposes
set.seed(1234)

# number of simulations
nsim <- 200
models2run <- length(modelsFitted)
ages.fit <- 45:90

modelSim <- list()
selectBEL <- list()

#modelSim <- sapply(modelsFitted, simulate, nsim = nsim, h = forecastTime)

for (m in 1:(models2run)){
  modelSim[[m]] <- simulate(modelsFitted[[m]], nsim = nsim, h = forecastTime)
  collectBEL <- array(NA, c(200,1))
  for (s in 1:nsim){
    prem_s <- 0
    ben_s <- 0
    qx <- cbind(modelSim[[m]]$fitted[, , s], modelSim[[m]]$rates[, , s])
    extrapolate <- kannistoExtrapolation(qx, ages.fit, years_chart)
    for (i in 1:nrow(portfolio)){
      ben_s <- ben_s + PVcashflow(extrapolate$qxt*expF, ageStart = portfolio$age[i], omegaAge = 120, 
                                  pensionAge = 65, valyear = valyear, ir = 0.02, type = "annuity")*pension
      prem_s <- prem_s + PVcashflow(extrapolate$qxt*expF, ageStart = portfolio$age[i], omegaAge = 120, 
                                    pensionAge = 65, valyear = valyear, ir = 0.02, type = "premium")*portfolio$Premium[i]
    }
    collectBEL[s, 1] <- ben_s - prem_s
  }
  selectBEL[[m]] <- quantile(collectBEL, probs = 0.995, type = 1)
}

SCR <- as.numeric(selectBEL) - BEL[1:models2run]
names(SCR) <- "SCR"

## Plot simulations for LC model
qxt <- Dxt/E0xt
plot(LCfit$years, qxt["65", ], xlim = c(1950, 2132), ylim = range(PLATqxt["65",]),
     xlab = "Years", ylab = "Mortality rates", main = "Mortality rates at age 65",
     pch = 20, log = "y", type = "l")
matlines(modelSim[[1]]$years, modelSim[[1]]$rates["65", , 1:20], type = "l", lty = 1, col = 1:20)

## Plot model uncertainity
probs <- c(0.5, 2.5, 10, 25, 75, 90, 97.5, 99.5)
plot(LCfit$years, qxt["65", ], xlim = c(1950, 2132), ylim = c(0.0025, 0.04),
     xlab = "Years", ylab = "Mortality rates at age 65", main = "Uncertainity associated with a model forecast",
     pch = 20, log = "y")
fan(t(modelSim[[1]]$rates["65", , ]), start = 2012, probs = probs, n.fan = 4,
    fan.col = colorRampPalette(c("yellow", "darkgreen")), ln = NULL)

## Show BEL and SCR per model
df = data.frame(models = factor(c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT"), levels=c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT")), BEL = BEL[1:8], SCR = SCR[1:8])
plot_df <- melt(df, id="models")

ggplot(data=plot_df, aes(x=models, y=value, fill=variable)) +
  geom_bar(stat="identity") + theme_minimal() +
  xlab("Models") + ylab("Amount in EUR") + ggtitle("BEL and SCR for a portfolio") +
  theme(panel.grid.major.x=element_blank(), axis.line = element_line(colour = "black", size =0.2), axis.title.x=element_text(family="sans",size=rel(1)),
        axis.title.y=element_text(family="sans",size=rel(1), angle=90), plot.title = element_text(family="sans", size=20, face="bold", vjust=2))