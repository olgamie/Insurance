## Exercise 1

ages.fit <- 45:120
valyear <- 2015
pension <- 1000

# calculate ages of insured
portfolio$age <- ## adjust the portfolio data

BEL_ex1 <- array(NA, c(9,1))
rownames(BEL_ex1) <- c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT", "AG")
colnames(BEL_ex1) <- "BEL"

for (m in 1:length(models)){
  output <- list()
  output2 <- list()
  for (i in 1:nrow(portfolio)){
    output[[i]] <- PVcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, 
                              valyear = valyear, ir = 0.02, type = "annuity")*pension
    output2[[i]] <- PVcashflow(models[[m]]*expF, ageStart = portfolio$age[i], omegaAge = 120, pensionAge = 65, 
                               valyear = valyear, ir = 0.02, type = "premium")*portfolio$Premium[i]
  }
  BEL_ex1[m, 1] <- do.call(sum, output)-do.call(sum, output2)
}

## Show BEL per model
df = data.frame(models = factor(c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT", "AG"), levels=c("LC", "APC", "RH", "CBD", "M6", "M7", "M8", "PLAT", "AG")),
                BEL = BEL_ex1)

ggplot(data=df, aes(x=models, y=BEL, fill=models)) +
  geom_bar(stat="identity") + theme_minimal() +
  xlab("Models") + ylab("Amount in EUR") + ggtitle("BEL for a portfolio") +
  theme(panel.grid.major.x=element_blank(), axis.line = element_line(colour = "black", size =0.5), legend.position="none", axis.title.x=element_text(family="sans",size=rel(1)),
        axis.title.y=element_text(family="sans",size=rel(1), angle=90), plot.title = element_text(family="sans", size=20, face="bold", vjust=2))

## Compare BEL results from the case study with the results from exercise 1
## Create a data.frame or matrix containing the orginal and exercies results
