### this code inputs a presence-absence matrix of OTUs by host individual and
### outputs species accumulation curves
### species == virus OTUs
### sites == host individuals

setwd("Sjodin_etal_XXXX_JOURNAL/Code")

library(vegan) # needed to get true sampling curve
library(iNEXT) # needed for extrapolation

# load in dataset - dataset is the site (host) by species (OTU) matrix output
## from the code OTU_Classification.R PLUS all of the bats/hosts/sites with
## no species/OTUs
raw.data <- read.csv("../Data/OTU_by_bat_pres_abs_mat.csv", header=TRUE)
# take a look at it:
raw.data[1:6,1:6]

# Get rid of uninformative OTU names
raw.data <- raw.data[,-1]

# rename rows to OTU numeration
rownames(raw.data) <- c("OTU 1", "OTU 10", "OTU 11", "OTU 12", "OTU 13", "OTU 14", "OTU 15",
                        "OTU 16", "OTU 17", "OTU 18", "OTU 19", "OTU 2", "OTU 20", "OTU 21",
                        "OTU 22", "OTU 23", "OTU 24", "OTU 25", "OTU 26", "OTU 27", "OTU 28",
                        "OTU 29", "OTU 3", "OTU 30", "OTU 31", "OTU 32", "OTU 33", "OTU 34",
                        "OTU 35", "OTU 36", "OTU 37", "OTU 38", "OTU 39", "OTU 4", "OTU 40",
                        "OTU 41", "OTU 42", "OTU 43", "OTU 5", "OTU 6", "OTU 7", "OTU 8", "OTU 9")

# make data a list to read into iNext
raw.data <- list(as.matrix(raw.data))

########---------------- for all host species in both caves ----------------########

# calculate chao richness estimate
chaorich <- ChaoRichness(raw.data[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.all <-  iNEXT(x = raw.data, q=0, datatype = "incidence_raw",
                    se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.all <- specaccum(t(raw.data[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.all <- iNEXT.all$iNextEst[[1]]$t[which(iNEXT.all$iNextEst[[1]]$t < 1087)]
interpy.all <- iNEXT.all$iNextEst[[1]]$qD[which(iNEXT.all$iNextEst[[1]]$t < 1087)]
extrapx.all <- iNEXT.all$iNextEst[[1]]$t[which(iNEXT.all$iNextEst[[1]]$t > 1086)]
extrapy.all <- iNEXT.all$iNextEst[[1]]$qD[which(iNEXT.all$iNextEst[[1]]$t > 1086)]

# create plot
plot(iNEXT.all$iNextEst[[1]]$t, iNEXT.all$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,60))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.all$iNextEst[[1]]$t, rev(iNEXT.all$iNextEst[[1]]$t)), 
        y = c(iNEXT.all$iNextEst[[1]]$qD.LCL, rev(iNEXT.all$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.all, interpy.all, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.all, extrapy.all, 
      col = "deeppink4", lwd = 3, lty = 2)
# add empirical collection curve
lines(x = as.vector(spp.accum.all[[3]]), y = as.vector(spp.accum.all[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data[[1]]), y = chaorich$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext("All Bats", side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich$Estimator
#chaorich$Observed/chaorich$Estimator
#add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 55.089\n 78% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- All host species in Culebrones Cave ----------------########

# subset data to remove Artibeus jamaicensis and Eptesicus fuscus
raw.data.Cul <- list(raw.data[[1]][, -grep("Artibeus", colnames(raw.data[[1]]))])
raw.data.Cul <- list(raw.data.Cul[[1]][, -grep("Eptesicus", colnames(raw.data.Cul[[1]]))])
# calculate chao richness estimate
chaorich.cul <- ChaoRichness(raw.data.Cul[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Culebrones <-  iNEXT(x = raw.data.Cul, q=0, datatype = "incidence_raw",
                           se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.Cul <- specaccum(t(raw.data.Cul[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.cul <- iNEXT.Culebrones$iNextEst[[1]]$t[which(iNEXT.Culebrones$iNextEst[[1]]$t < (ncol(raw.data.Cul[[1]])+1))]
interpy.cul <- iNEXT.Culebrones$iNextEst[[1]]$qD[which(iNEXT.Culebrones$iNextEst[[1]]$t < (ncol(raw.data.Cul[[1]])+1))]
extrapx.cul <- iNEXT.Culebrones$iNextEst[[1]]$t[which(iNEXT.Culebrones$iNextEst[[1]]$t > ncol(raw.data.Cul[[1]]))]
extrapy.cul <- iNEXT.Culebrones$iNextEst[[1]]$qD[which(iNEXT.Culebrones$iNextEst[[1]]$t > ncol(raw.data.Cul[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Culebrones$iNextEst[[1]]$t, iNEXT.Culebrones$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,55))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Culebrones$iNextEst[[1]]$t, rev(iNEXT.Culebrones$iNextEst[[1]]$t)), 
        y = c(iNEXT.Culebrones$iNextEst[[1]]$qD.LCL, rev(iNEXT.Culebrones$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.cul, interpy.cul, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.cul, extrapy.cul, 
      col = "deeppink4", lwd = 2, lty = 3)
# add empirical collection curve
lines(x = as.vector(spp.accum.Cul[[3]]), y = as.vector(spp.accum.Cul[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Cul[[1]]), y = chaorich.cul$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext("Culebrones Cave", side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.cul$Estimator
#chaorich.cul$Observed/chaorich.cul$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 47.325\n 82% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- All host species in Larva Cave ----------------########
### Only one bat species (A. jamaicensis) tested positive in Larva Cave

# subset data to include only jamaicensis and Eptesicus fuscus
raw.data.Lar <- list(raw.data[[1]][, -grep("Pteronotus", colnames(raw.data[[1]]))])
raw.data.Lar <- list(raw.data.Lar[[1]][, -grep("Mormoops", colnames(raw.data.Lar[[1]]))])
raw.data.Lar <- list(raw.data.Lar[[1]][, -grep("Monophyllus", colnames(raw.data.Lar[[1]]))])
raw.data.Lar <- list(raw.data.Lar[[1]][, -grep("Erophylla", colnames(raw.data.Lar[[1]]))])
raw.data.Lar <- list(raw.data.Lar[[1]][, -grep("Brachyphylla", colnames(raw.data.Lar[[1]]))])
# calculate chao richness estimate
chaorich.lar <- ChaoRichness(raw.data.Lar[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Larva <-  iNEXT(x = raw.data.Lar, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.Lar <- specaccum(t(raw.data.Lar[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.lar <- iNEXT.Larva$iNextEst[[1]]$t[which(iNEXT.Larva$iNextEst[[1]]$t < (ncol(raw.data.Lar[[1]])+1))]
interpy.lar <- iNEXT.Larva$iNextEst[[1]]$qD[which(iNEXT.Larva$iNextEst[[1]]$t < (ncol(raw.data.Lar[[1]])+1))]
extrapx.lar <- iNEXT.Larva$iNextEst[[1]]$t[which(iNEXT.Larva$iNextEst[[1]]$t > ncol(raw.data.Lar[[1]]))]
extrapy.lar <- iNEXT.Larva$iNextEst[[1]]$qD[which(iNEXT.Larva$iNextEst[[1]]$t > ncol(raw.data.Lar[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Larva$iNextEst[[1]]$t, iNEXT.Larva$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,15))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Larva$iNextEst[[1]]$t, rev(iNEXT.Larva$iNextEst[[1]]$t)), 
        y = c(iNEXT.Larva$iNextEst[[1]]$qD.LCL, rev(iNEXT.Larva$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.lar, interpy.lar, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.lar, extrapy.lar, 
      col = "deeppink4", lwd = 2, lty = 3)
# add empirical collection curve
lines(x = as.vector(spp.accum.Lar[[3]]), y = as.vector(spp.accum.Lar[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Lar[[1]]), y = chaorich.lar$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext("Larva Cave", side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.lar$Estimator
#chaorich.lar$Observed/chaorich.lar$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 11.434\n 61% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Pteronotus portoricensis ----------------########
# only 3 OTUs discovered in this species

# subset data
raw.data.Ptepor <- list(raw.data[[1]][, grep("Pteronotus_portoricensis", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.ptepor <- ChaoRichness(raw.data.Ptepor[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Ptepor <- iNEXT(x = raw.data.Ptepor, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.Pp <- specaccum(t(raw.data.Ptepor[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.ptepor <- iNEXT.Ptepor$iNextEst[[1]]$t[which(iNEXT.Ptepor$iNextEst[[1]]$t < (ncol(raw.data.Ptepor[[1]])+1))]
interpy.ptepor <- iNEXT.Ptepor$iNextEst[[1]]$qD[which(iNEXT.Ptepor$iNextEst[[1]]$t < (ncol(raw.data.Ptepor[[1]])+1))]
extrapx.ptepor <- iNEXT.Ptepor$iNextEst[[1]]$t[which(iNEXT.Ptepor$iNextEst[[1]]$t > ncol(raw.data.Ptepor[[1]]))]
extrapy.ptepor <- iNEXT.Ptepor$iNextEst[[1]]$qD[which(iNEXT.Ptepor$iNextEst[[1]]$t > ncol(raw.data.Ptepor[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Ptepor$iNextEst[[1]]$t, iNEXT.Ptepor$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,6))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Ptepor$iNextEst[[1]]$t, rev(iNEXT.Ptepor$iNextEst[[1]]$t)), 
        y = c(iNEXT.Ptepor$iNextEst[[1]]$qD.LCL, rev(iNEXT.Ptepor$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.ptepor, interpy.ptepor, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.ptepor, extrapy.ptepor, 
      col = "deeppink4", lwd = 2, lty = 3)
# add empirical collection curve
lines(x = as.vector(spp.accum.Pp[[3]]), y = as.vector(spp.accum.Pp[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Ptepor[[1]]), y = chaorich.ptepor$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("P. portoricensis")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.ptepor$Estimator
#chaorich.ptepor$Observed/chaorich.ptepor$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 3.242\n 93% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Pteronotus quadridens ----------------########

# subset data
raw.data.Ptequa <- list(raw.data[[1]][, grep("Pteronotus_quadridens", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.ptequa <- ChaoRichness(raw.data.Ptequa[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Ptequa <- iNEXT(x = raw.data.Ptequa, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.Pq <- specaccum(t(raw.data.Ptequa[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.ptequa <- iNEXT.Ptequa$iNextEst[[1]]$t[which(iNEXT.Ptequa$iNextEst[[1]]$t < (ncol(raw.data.Ptequa[[1]])+1))]
interpy.ptequa <- iNEXT.Ptequa$iNextEst[[1]]$qD[which(iNEXT.Ptequa$iNextEst[[1]]$t < (ncol(raw.data.Ptequa[[1]])+1))]
extrapx.ptequa <- iNEXT.Ptequa$iNextEst[[1]]$t[which(iNEXT.Ptequa$iNextEst[[1]]$t > ncol(raw.data.Ptequa[[1]]))]
extrapy.ptequa <- iNEXT.Ptequa$iNextEst[[1]]$qD[which(iNEXT.Ptequa$iNextEst[[1]]$t > ncol(raw.data.Ptequa[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Ptequa$iNextEst[[1]]$t, iNEXT.Ptequa$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,30))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Ptequa$iNextEst[[1]]$t, rev(iNEXT.Ptequa$iNextEst[[1]]$t)), 
        y = c(iNEXT.Ptequa$iNextEst[[1]]$qD.LCL, rev(iNEXT.Ptequa$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.ptequa, interpy.ptequa, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.ptequa, extrapy.ptequa, 
      col = "deeppink4", lwd = 2, lty = 3)
# add empirical collection curve
lines(x = as.vector(spp.accum.Pq[[3]]), y = as.vector(spp.accum.Pq[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Ptequa[[1]]), y = chaorich.ptequa$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("P. quadridens")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.ptequa$Estimator
#chaorich.ptequa$Observed/chaorich.ptequa$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 27.207\n 55% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Monophyllus redmani ----------------########

# subset data
raw.data.Monred <- list(raw.data[[1]][, grep("Monophyllus", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.monred <- ChaoRichness(raw.data.Monred, conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Monred <- iNEXT(x = raw.data.Monred, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.Mr <- specaccum(t(raw.data.Monred[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.monred <- iNEXT.Monred$iNextEst[[1]]$t[which(iNEXT.Monred$iNextEst[[1]]$t < (ncol(raw.data.Monred[[1]])+1))]
interpy.monred <- iNEXT.Monred$iNextEst[[1]]$qD[which(iNEXT.Monred$iNextEst[[1]]$t < (ncol(raw.data.Monred[[1]])+1))]
extrapx.monred <- iNEXT.Monred$iNextEst[[1]]$t[which(iNEXT.Monred$iNextEst[[1]]$t > ncol(raw.data.Monred[[1]]))]
extrapy.monred <- iNEXT.Monred$iNextEst[[1]]$qD[which(iNEXT.Monred$iNextEst[[1]]$t > ncol(raw.data.Monred[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Monred$iNextEst[[1]]$t, iNEXT.Monred$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,25))#, cex.lab = 1.5)
# calculate confidence boundaries
polygon(x = c(iNEXT.Monred$iNextEst[[1]]$t, rev(iNEXT.Monred$iNextEst[[1]]$t)), 
        y = c(iNEXT.Monred$iNextEst[[1]]$qD.LCL, rev(iNEXT.Monred$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation line
lines(interpx.monred, interpy.monred, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.monred, extrapy.monred, 
      col = "deeppink4", lwd = 2, lty = 3)
# add empirical collection curve
lines(x = as.vector(spp.accum.Mr[[3]]), y = as.vector(spp.accum.Mr[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Monred[[1]]), y = chaorich.monred$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("Mon. redmani")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.monred$Estimator
#chaorich.monred$Observed/chaorich.monred$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 16.151\n 74% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Mormoops blainvillii ----------------########

# subset data
raw.data.Morbla <- list(raw.data[[1]][, grep("Mormoops", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.morbla <- ChaoRichness(raw.data.Morbla[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Morbla <- iNEXT(x = raw.data.Morbla, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# create empirical accumulation curve
spp.accum.Mb <- specaccum(t(raw.data.Morbla[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.morbla <- iNEXT.Morbla$iNextEst[[1]]$t[which(iNEXT.Morbla$iNextEst[[1]]$t < (ncol(raw.data.Morbla[[1]])+1))]
interpy.morbla <- iNEXT.Morbla$iNextEst[[1]]$qD[which(iNEXT.Morbla$iNextEst[[1]]$t < (ncol(raw.data.Morbla[[1]])+1))]
extrapx.morbla <- iNEXT.Morbla$iNextEst[[1]]$t[which(iNEXT.Morbla$iNextEst[[1]]$t > ncol(raw.data.Morbla[[1]]))]
extrapy.morbla <- iNEXT.Morbla$iNextEst[[1]]$qD[which(iNEXT.Morbla$iNextEst[[1]]$t > ncol(raw.data.Morbla[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Morbla$iNextEst[[1]]$t, iNEXT.Morbla$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,20))#, cex.lab = 1.5)
# calculate confidence boundaries
polygon(x = c(iNEXT.Morbla$iNextEst[[1]]$t, rev(iNEXT.Morbla$iNextEst[[1]]$t)), 
        y = c(iNEXT.Morbla$iNextEst[[1]]$qD.LCL, rev(iNEXT.Morbla$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.morbla, interpy.morbla, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.morbla, extrapy.morbla, 
      col = "deeppink4", lwd = 3, lty = 2)
# add empirical collection curve
lines(x = as.vector(spp.accum.Mb[[3]]), y = as.vector(spp.accum.Mb[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Morbla[[1]]), y = chaorich.morbla$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("Mor. blainvillii")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.morbla$Estimator
#chaorich.morbla$Observed/chaorich.morbla$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 13.485\n 67% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Erophylla bombifrons ----------------########

# subset data
raw.data.Erobom <- list(raw.data[[1]][, grep("Erophylla", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.Erobom <- ChaoRichness(raw.data.Erobom[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Erobom <- iNEXT(x = raw.data.Erobom, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# calculate empirical accumulation curve
spp.accum.Es <- specaccum(t(raw.data.Erobom[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.Erobom <- iNEXT.Erobom$iNextEst[[1]]$t[which(iNEXT.Erobom$iNextEst[[1]]$t < (ncol(raw.data.Erobom[[1]])+1))]
interpy.Erobom <- iNEXT.Erobom$iNextEst[[1]]$qD[which(iNEXT.Erobom$iNextEst[[1]]$t < (ncol(raw.data.Erobom[[1]])+1))]
extrapx.Erobom <- iNEXT.Erobom$iNextEst[[1]]$t[which(iNEXT.Erobom$iNextEst[[1]]$t > ncol(raw.data.Erobom[[1]]))]
extrapy.Erobom <- iNEXT.Erobom$iNextEst[[1]]$qD[which(iNEXT.Erobom$iNextEst[[1]]$t > ncol(raw.data.Erobom[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Erobom$iNextEst[[1]]$t, iNEXT.Erobom$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,20))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Erobom$iNextEst[[1]]$t, rev(iNEXT.Erobom$iNextEst[[1]]$t)), 
        y = c(iNEXT.Erobom$iNextEst[[1]]$qD.LCL, rev(iNEXT.Erobom$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.Erobom, interpy.Erobom, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.Erobom, extrapy.Erobom, 
      col = "deeppink4", lwd = 3, lty = 2)
# add empirical collection curve
lines(x = as.vector(spp.accum.Es[[3]]), y = as.vector(spp.accum.Es[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Erobom[[1]]), y = chaorich.Erobom$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("E. bombifrons")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.Erobom$Estimator
#chaorich.Erobom$Observed/chaorich.Erobom$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 14.892\n 47% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Brachyphylla cavernarum ----------------########

# subset data
raw.data.Bracav <- list(raw.data[[1]][, grep("Brachyphylla", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.bracav <- ChaoRichness(raw.data.Bracav[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Bracav <- iNEXT(x = raw.data.Bracav, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# calculate empirical accumulation curve
spp.accum.Bc <- specaccum(t(raw.data.Bracav[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.bracav <- iNEXT.Bracav$iNextEst[[1]]$t[which(iNEXT.Bracav$iNextEst[[1]]$t < (ncol(raw.data.Bracav[[1]])+1))]
interpy.bracav <- iNEXT.Bracav$iNextEst[[1]]$qD[which(iNEXT.Bracav$iNextEst[[1]]$t < (ncol(raw.data.Bracav[[1]])+1))]
extrapx.bracav <- iNEXT.Bracav$iNextEst[[1]]$t[which(iNEXT.Bracav$iNextEst[[1]]$t > ncol(raw.data.Bracav[[1]]))]
extrapy.bracav <- iNEXT.Bracav$iNextEst[[1]]$qD[which(iNEXT.Bracav$iNextEst[[1]]$t > ncol(raw.data.Bracav[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Bracav$iNextEst[[1]]$t, iNEXT.Bracav$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,10))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Bracav$iNextEst[[1]]$t, rev(iNEXT.Bracav$iNextEst[[1]]$t)), 
        y = c(iNEXT.Bracav$iNextEst[[1]]$qD.LCL, rev(iNEXT.Bracav$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.bracav, interpy.bracav, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.bracav, extrapy.bracav, 
      col = "deeppink4", lwd = 3, lty = 2)
# add empirical collection curve
lines(x = as.vector(spp.accum.Bc[[3]]), y = as.vector(spp.accum.Bc[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Bracav[[1]]), y = chaorich.bracav$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("B.cavernarum")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.bracav$Estimator
#chaorich.bracav$Observed/chaorich.bracav$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 7.492\n 93% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)

########---------------- Artibeus jamaicensis ----------------########

# subset data
raw.data.Artjam <- list(raw.data[[1]][, grep("Artibeus", colnames(raw.data[[1]]))])
# calculate chao richness estimate
chaorich.artjam <- ChaoRichness(raw.data.Artjam[[1]], conf = 0.95, datatype="incidence_raw")
# generate species accumulation curve data
iNEXT.Artjam <- iNEXT(x = raw.data.Artjam, q=0, datatype = "incidence_raw",
                      se = TRUE, conf = 0.95, nboot=1000)
# calculate empirical accumulation curve
spp.accum.Aj <- specaccum(t(raw.data.Artjam[[1]]), method="collector")

## plot the curve

# subset data to interpolation and extrapolation parts
interpx.artjam <- iNEXT.Artjam$iNextEst[[1]]$t[which(iNEXT.Artjam$iNextEst[[1]]$t < (ncol(raw.data.Artjam[[1]])+1))]
interpy.artjam <- iNEXT.Artjam$iNextEst[[1]]$qD[which(iNEXT.Artjam$iNextEst[[1]]$t < (ncol(raw.data.Artjam[[1]])+1))]
extrapx.artjam <- iNEXT.Artjam$iNextEst[[1]]$t[which(iNEXT.Artjam$iNextEst[[1]]$t > ncol(raw.data.Artjam[[1]]))]
extrapy.artjam <- iNEXT.Artjam$iNextEst[[1]]$qD[which(iNEXT.Artjam$iNextEst[[1]]$t > ncol(raw.data.Artjam[[1]]))]

# for setting plot margins, if necessary: bottom, left, top, right
# change margins if you want to add ylab = "OTU Richness"
#par(mar = c(5,5,4,2) + 0.1) ## default is c(5,4,4,2) + 0.1

# create plot
plot(iNEXT.Artjam$iNextEst[[1]]$t, iNEXT.Artjam$iNextEst[[1]]$qD, main = " ", cex.axis = 1.75, type = "l",
     xlab = "", ylab = "", ylim = c(0,15))#, cex.lab = 1.5)
# add confidence boundaries
polygon(x = c(iNEXT.Artjam$iNextEst[[1]]$t, rev(iNEXT.Artjam$iNextEst[[1]]$t)), 
        y = c(iNEXT.Artjam$iNextEst[[1]]$qD.LCL, rev(iNEXT.Artjam$iNextEst[[1]]$qD.UCL)),
        col = "pink", border = NA)
# add interpolation curve
lines(interpx.artjam, interpy.artjam, 
      col = "deeppink4", lwd = 3)
# add extrapolation curve
lines(extrapx.artjam, extrapy.artjam, 
      col = "deeppink4", lwd = 3, lty = 2)
lines(x = as.vector(spp.accum.Aj[[3]]), y = as.vector(spp.accum.Aj[[4]]), lwd = 3)
# add a point for empirical number of OTUs discovered in number of bats samples
points(x = ncol(raw.data.Artjam[[1]]), y = chaorich.artjam$Observed, 
       col = "red4", pch = 17, cex = 2.5)
# add a description to the plot
mtext(bquote(italic("A. jamaicensis")), side = 3, line = -2, cex = 1.75, adj = 0.05)
# find chao2 estimate and percent completeness of sampling:
#chaorich.artjam$Estimator
#chaorich.artjam$Observed/chaorich.artjam$Estimator
# add chao2 estimate and percent completeness of sampling to graph
mtext("Chao2 = 11.421\n 61% Complete", side = 1, line = -2, cex = 1.5, adj = 0.95)










