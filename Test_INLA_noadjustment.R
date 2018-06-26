rm(list=ls())
## Install INLA, maptools, spdep package
setwd("Thu INLA/")

install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable", dep=TRUE)
source("http://www.math.ntnu.no/inla/givemeINLA.R")
inla.upgrade(testing=TRUE)
install.packages("maptools")
install.packages("spdep")
install.packages("lattice")

library("INLA")
library("maptools")
library("spdep")
library("lattice")

## Upload Amsterdam shape file
amsmap <- readShapePoly("amsmap_hiv.shp")

## Transform the shape file into an adjacency matrix and to make it compatible with the R-INLA format
temp <- poly2nb(amsmap)
nb2INLA("AMS.graph",temp)
AMS.adj <- paste(getwd(),"/AMS.graph",sep="")

## Import the AMS.adj graph in R-INLA and obtain the adjacency matrix
H <- inla.read.graph(filename="AMS.graph")
image(inla.graph2matrix(H),xlab="",ylab="")

## Specify formula for the model
## ID represents the identifiers for the boroughs and through the 'graph' option
## we include the name of the object containing the neighborhood structure
formula <- y ~ 1 + f(ID, model="bym",graph=AMS.adj, scale.model=TRUE,
                     hyper=list(prec.unstruct=list(prior="loggamma",param=c(1,0.001)), prec.spatial=list(prior="loggamma",param=c(1,0.001))))
#@In the above formular, we specify prior distribution of u and v in 'hyper' argument.

data.hiv <- read.csv("data_hiv_Popcount.csv")
data.hiv <- data.hiv[,-1]    ##remove the row index colume
Nareas <- length(data.hiv[,1])

# The order of the areas needs to be the same between the data and the spatial polygon object obtained importing the shapefile, so we re-order the data
data.boroughs <- attr(amsmap, "data")
data.boroughs <- data.boroughs[data.boroughs$PC4 %in% data.hiv$PC4,]
order <- match(data.boroughs$PC4,data.hiv$PC4)
data.hiv <- data.hiv[order,]
data.hiv$ID <- seq(1,Nareas)

# Include in the amsmap dataframe the ID variable  
attr(amsmap, "data") <- merge(data.boroughs,data.hiv,by=c("PC4"))

mod.hiv <- inla(formula,family="poisson",
                data=data.hiv, E=data.hiv$Pop,
                control.compute=list(dic=TRUE))

round(mod.hiv$summary.fixed,3) 
round(head(mod.hiv$summary.random$ID),3) #partial output
round(mod.hiv$summary.random$ID,3) 

exp.b0.mean <- inla.emarginal(exp,mod.hiv$marginals.fixed[[1]])
exp.b0.mean
exp.b0.95CI <- inla.qmarginal(c(0.025,0.975), inla.tmarginal(exp,mod.hiv$marginals.fixed[[1]]))
exp.b0.95CI


# *** Code for Figure 6.6 left - map of the posterior mean for the borough-specific relative risks of acquiring HIV, compared to the whole of Amsterdam
csi <- mod.hiv$marginals.random$ID[1:Nareas]
zeta <- lapply(csi,function(x) inla.emarginal(exp,x))
# Define the cutoff for zeta
summary(unlist(zeta)) ##get an overview of zeta's values to decide the cutoff 
zeta.cutoff <- c(0.4, 0.8, 1.2, 1.6, 2.0, 3.4)
# Transform zeta in categorical variable
cat.zeta <- cut(unlist(zeta),breaks=zeta.cutoff,include.lowest=TRUE)
# Create a dataframe with all the information needed for the map
maps.cat.zeta <- data.frame(ID=data.hiv$ID, cat.zeta=cat.zeta)
# Add the categorized zeta to the spatial polygon
data.boroughs <- attr(amsmap, "data")
attr(amsmap, "data") <- merge(data.boroughs, maps.cat.zeta, by="ID")

trellis.par.set(axis.line=list(col=NA))
spplot(obj=amsmap, zcol= "cat.zeta", col.regions=gray(seq(0.9,0.1,length=5)), asp=1)
# ***

# *** Code for Figure 6.6 right - The uncertainty associated with the posterior means
a <- 0
prob.csi <- lapply(csi, function(x) {1 - inla.pmarginal(a, x)})
prob.csi.cutoff <- c(0,0.2,0.8,1)
cat.prob.csi <- cut(unlist(prob.csi),breaks=prob.csi.cutoff, include.lowest=TRUE)
#Create a dataframe with all the information needed for the map
maps.cat.prob.csi <- data.frame(ID=data.hiv$ID, cat.prob.csi=cat.prob.csi)

## Optional: Create zeta data to put in Tableau
#data.boroughs <- attr(amsmap, "data")
#maps.zeta <- data.frame(ID=data.hiv$ID, zeta=unlist(zeta))
#attr(amsmap, "data") <- merge(data.boroughs, maps.zeta, by="ID")
#write.csv(attr(amsmap, "data"), file="modelresult_zeta_nosexage.csv")

#Add the categorized zetea to the spatial polygon
data.boroughs <- attr(amsmap, "data")
attr(amsmap, "data") <- merge(data.boroughs, maps.cat.prob.csi, by="ID")

#Map zeta
spplot(obj=amsmap, zcol= "cat.prob.csi", col.regions=gray(seq(0.9,0.1,length=3)))
# ***



## Finally, it could be interesting to evaluate the proportion of variance explained by the 
#@ structured spatial component. 
## To obtain this in INLA, a simulation-based approach is used. For each area we extract a
#@ large enough number of values (e.g., 100 000) from the corresponding marginal posterior 
#@ distribution of v_i and save the simulated values in a matrix with rows equal to the 
#@ number of areas and 100 000 columns. Then we calculate the emperical variance (over column)

mat.marg <- matrix(NA, nrow=Nareas, ncol=100000)
m <- mod.hiv$marginals.random$ID
for (i in 1:Nareas){
  # Remember that the first Nareas values of the random effects
  # are u+v, while u values are stored in the Nareas+1 to 2*Nareas elements.
  u <- m[[Nareas+i]]
  mat.marg[i,] <- inla.rmarginal(100000, u)
}
var.u <- apply(mat.marg, 2, var)  ##the emperical variance (over column) of u (the spatially structured component)

##expected value of the variance for v (the unstructured component) is calculated as
var.v <- inla.rmarginal(100000,inla.tmarginal(function(x) 1/x,
                                              mod.hiv$marginals.hyper$"Precision for ID (iid component)"))

##Finally, we compute the spatial fractional variance as
perc.var.u <- mean(var.u/(var.u+var.v))  
perc.var.u   #0.6396386

## In the case of HIV in Amsterdam, the proportion of spatial variance is about 0.16,
# @suggesting that only a small part of the variability is explained by the spatial structure.

marg.hyper <- inla.hyperpar.sample(100000,mod.hiv)
perc.var.u1 <- mean(marg.hyper[,1] / (marg.hyper[,1]+marg.hyper[,2]))
perc.var.u1




