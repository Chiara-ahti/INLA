rm(list=ls())

###Data from the tot count
setwd("")
hiv_agesex <- read.csv("hiv_agesex2011.csv")

hiv_msm <- hiv_agesex[,c("SexTransRegion","count")]
hiv_aggMSM <- aggregate(count ~ SexTransRegion, data=hiv_msm, FUN=sum)
#head(hiv_aggMSM)


#setwd("/Thu INLA/")
Amsterdam_MSM <- read.table("ahti_sex_22_gebieden.txt",header=TRUE, sep="\t", dec = ",")
Amsterdam_pc4 <- read.table("22_gebieden_pc4.txt",header=TRUE, sep="\t")

##Put data frame together to merge the gebied and the pc4
Amsterdam <- merge(Amsterdam_MSM, Amsterdam_pc4, by.x = c("Gebied"), by.y = c("Gebied"), all=T)
Amsterdam<- Amsterdam[-c(1), ]
library(reshape2)
Amsterdam1<- melt(Amsterdam, id.vars = c("Gebied", "PC4"), measure.vars = c("LHBT", "Hetero"))

pop <- read.csv("pop_age&sex_pc4.csv")
pop_agg <- aggregate(pop ~ pc4, data=pop, FUN=sum)

Amsterdam_pop <- merge(Amsterdam1, pop_agg, by.x = c("PC4"), by.y = c("pc4"), all=T)
Amsterdam_pop$homoPop <- Amsterdam_pop$pop * (Amsterdam_pop$value / 100) #number of LGBT per gebied

##calculate amount tot LGBT per tot population
Ams_hiv_h<- aggregate(homoPop ~ variable, data = Amsterdam_pop, FUN=sum)
Ams_hiv_h<- subset(Ams_hiv_h, Ams_hiv_h$variable == "LHBT")
Ams_hiv_h$SexTransRegion<- rep(c(11))
Ams_hiv_rate<- merge(Ams_hiv_h, hiv_aggMSM, by.x = c("SexTransRegion"), by.y = c("SexTransRegion"), all=T)
Ams_hiv_rate$rate<- Ams_hiv_rate$count/Ams_hiv_rate$homoPop #rate of hiv gay withing the gay pop


##merge tot rate LGBT with population per PC4
data_LGBT <- merge(Amsterdam_pop, Ams_hiv_rate, by.x=c("variable"),by.y=c("variable"),all=TRUE)
data_LGBT$E_j <- data_LGBT$homoPop.x*data_LGBT$rate
data_LGBT <- subset(data_LGBT, data_LGBT$variable == "LHBT")

data_LGBT_agg <- aggregate(E_j ~ PC4, data=data_LGBT, FUN=sum)

#Data on total HIV diagnoses since 2011
hiv_pc4 <- read.csv("/hiv_pc4_year.csv")
hiv_pc4_since2011 <- hiv_pc4[,2:8]
head(hiv_pc4_since2011)
hiv_pc4_since2011$y_i <- rowSums(hiv_pc4_since2011[,2:7],na.rm=TRUE)

#merge expected case table with HIV diagnose table
data_hiv <- merge(data_LGBT_agg, hiv_pc4_since2011,  by.x=c("PC4"),by.y=c("pc4"),all=TRUE, na.rm=T)
data_hiv <- data_hiv[,c("PC4","E_j","y_i")]

completeFun <- function(data_hiv, desiredCols) {
  completeVec <- complete.cases(data_hiv[, desiredCols])
  return(data_hiv[completeVec, ])
}

data_hiv<- completeFun(data_hiv, "y_i")
data_hiv[is.na(data_hiv)] <- 0

names(data_hiv) <- c("PC4","E","y")
write.csv(data_hiv, file="data_hiv_LGBT.csv")





