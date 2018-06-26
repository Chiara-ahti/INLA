rm(list=ls())

#################### Aggregate HIV diagnoses per age & sex combination and Calculate expected cases based on popolation composition
setwd("")
hiv_agesex <- read.csv("hiv_agesex2011.csv")

hiv_agesex <- hiv_agesex[,c("sex","age_link","count")]
hiv_agg <- aggregate(count ~ sex+age_link, data=hiv_agesex, FUN=sum)
head(hiv_agg)
hiv_agg

pop <- read.csv("pop_age&sex_pc4.csv")
pop_agg <- aggregate(pop ~ sex+age_group, data=pop, FUN=sum)

#Merge the two tables to calculate the age and sex standardzed reference rate r_j
pop_hiv <- merge(hiv_agg, pop_agg, by.x= c("age_link","sex"), by.y=c("age_group","sex"),all=TRUE)
pop_hiv$rate <- pop_hiv$count/pop_hiv$pop

pop_hiv <- pop_hiv[,c("age_link","sex","rate")]

#Merge the above table with the population table to calculate the expected cases per PC4
data <- merge(pop, pop_hiv, by.x=c("age_group","sex"),by.y=c("age_link","sex"),all=TRUE)
data$E_j <- data$pop*data$rate
data_agg <- aggregate(E_j ~ pc4, data=data, FUN=sum)

#Data on total HIV diagnoses since 2011
hiv_pc4 <- read.csv("/hiv_pc4_year.csv")
hiv_pc4_since2011 <- hiv_pc4[,2:8]
head(hiv_pc4_since2011)
hiv_pc4_since2011$y_i <- rowSums(hiv_pc4_since2011[,2:7],na.rm=TRUE)

#merge expected case table with HIV diagnose table
data_hiv <- merge(data_agg, hiv_pc4_since2011[,c("pc4","y_i")])
names(data_hiv) <- c("PC4","E","y")
write.csv(data_hiv, file="data_hiv.csv")
