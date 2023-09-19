### this code inputs a presence-absence matrix of OTUs by host individual and
### outputs rank frequency distributions (Figure 4)


library(dplyr)

# load in dataset - dataset is the site (host) by species (OTU) matrix output
## from the code OTU_Classification.R PLUS all of the bats/hosts/sites with
## no species/OTUs
raw.data <- read.csv("../Data/OTU_by_bat_pres_abs_mat.csv", header=TRUE)
# take a look at it:
raw.data[1:6,1:6]

# Get rid of uninformative names
raw.data <- raw.data[,-1]

# rename rows to OTU numeration
rownames(raw.data) <- c("OTU 1", "OTU 10", "OTU 11", "OTU 12", "OTU 13", "OTU 14", "OTU 15",
                        "OTU 16", "OTU 17", "OTU 18", "OTU 19", "OTU 2", "OTU 20", "OTU 21",
                        "OTU 22", "OTU 23", "OTU 24", "OTU 25", "OTU 26", "OTU 27", "OTU 28",
                        "OTU 29", "OTU 3", "OTU 30", "OTU 31", "OTU 32", "OTU 33", "OTU 34",
                        "OTU 35", "OTU 36", "OTU 37", "OTU 38", "OTU 39", "OTU 4", "OTU 40",
                        "OTU 41", "OTU 42", "OTU 43", "OTU 5", "OTU 6", "OTU 7", "OTU 8", "OTU 9")

# sum and sort OTU frequencies
all.viral.freq <- rowSums(raw.data)
all.viral.freq <- sort(all.viral.freq, decreasing=TRUE)

# set plot margins to have room for increased font size
par(mar=c(7,4,4,2))  #bottom, left, top, right

### all bats
barplot(all.viral.freq, ylim = c(0,70), las=2, cex.axis = 2, cex.names = 1.5,
        main = "All Bats")
abline(h = (1/nrow(raw.data))*sum(all.viral.freq), col = "red", lwd = 3)

# determine how many OTUs are considered common, rare
length(which(all.viral.freq > (1/nrow(raw.data))*sum(all.viral.freq)))
length(all.viral.freq) - length(which(all.viral.freq > (1/nrow(raw.data))*sum(all.viral.freq)))

### Culebrones
Culebrones <- raw.data[, -grep("Eptesicus", colnames(raw.data))]
Culebrones <- Culebrones[, -grep("Artibeus", colnames(Culebrones))]
ncol(Culebrones)

cul.viral.freq <- rowSums(Culebrones)
cul.viral.freq <- sort(cul.viral.freq, decreasing=TRUE)
# remove zeros
cul.viral.freq <- cul.viral.freq[which(cul.viral.freq != 0)]

barplot(cul.viral.freq, ylim=c(0,70), 
        las=2, cex.axis = 2, cex.names = 1.25)
abline(h = (1/(length(cul.viral.freq)))*sum(cul.viral.freq), col = "red", lwd = 3)

length(cul.viral.freq)
length(which(cul.viral.freq > (1/length(cul.viral.freq))*sum(cul.viral.freq)))
length(cul.viral.freq) - length(which(cul.viral.freq > (1/length(cul.viral.freq))*sum(cul.viral.freq)))

### Larvas - will be the same as A. jamaicensis
##### only A. jamaicensis included in publication figure
Larva <- raw.data[, -grep("Pteronotus", colnames(raw.data))]
Larva <- Larva[, -grep("Mormoops", colnames(Larva))]
Larva <- Larva[, -grep("Monophyllus", colnames(Larva))]
Larva <- Larva[, -grep("Brachyphylla", colnames(Larva))]
Larva <- Larva[, -grep("Erophylla", colnames(Larva))]


lar.viral.freq <- rowSums(Larva)
lar.viral.freq <- sort(lar.viral.freq, decreasing=TRUE)
# remove zeros
lar.viral.freq <- lar.viral.freq[which(lar.viral.freq != 0)]

barplot(lar.viral.freq, ylim=c(0,17), las = 2,
        cex.axis = 2, cex.names = 1.25)
abline(h = (1/(length(lar.viral.freq)))*sum(lar.viral.freq), col = "red", lwd = 3)

length(lar.viral.freq)
length(which(lar.viral.freq > (1/length(lar.viral.freq))*sum(lar.viral.freq)))
length(lar.viral.freq) - length(which(lar.viral.freq > (1/length(lar.viral.freq))*sum(lar.viral.freq)))

######## Pte qua ##########

Ptequa <- raw.data[, grep("quadridens", colnames(raw.data))]

Ptequa.freq <- rowSums(Ptequa)
Ptequa.freq <- sort(Ptequa.freq, decreasing=TRUE)
# remove zeros
Ptequa.freq <- Ptequa.freq[which(Ptequa.freq != 0)]

barplot(Ptequa.freq, ylim=c(0,30), las = 2, main = "Pte qua", 
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Ptequa.freq)))*sum(Ptequa.freq), col = "red", lwd = 3)

ncol(Ptequa)
length(Ptequa.freq)
length(which(Ptequa.freq > (1/length(Ptequa.freq))*sum(Ptequa.freq)))
length(Ptequa.freq) - length(which(Ptequa.freq > (1/length(Ptequa.freq))*sum(Ptequa.freq)))

######## Pte por ########## - NOT USEFUL, ONLY 3 OTUs
## not included in the publication

Ptepor <- raw.data[, grep("portoricensis", colnames(raw.data))]

Ptepor.freq <- rowSums(Ptepor)
Ptepor.freq <- sort(Ptepor.freq, decreasing=TRUE)
# remove zeros
Ptepor.freq <- Ptepor.freq[which(Ptepor.freq != 0)]

barplot(Ptepor.freq, ylim=c(0,3), las = 2, main = "Pte por", 
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Ptepor.freq)))*sum(Ptepor.freq), col = "red", lwd = 3)

ncol(Ptepor)
length(Ptepor.freq)
length(which(Ptepor.freq > (1/length(Ptepor.freq))*sum(Ptepor.freq)))
length(Ptepor.freq) - length(which(Ptepor.freq > (1/length(Ptepor.freq))*sum(Ptepor.freq)))

######## Mor bla ##########

Morbla <- raw.data[, grep("Mormoops", colnames(raw.data))]

Morbla.freq <- rowSums(Morbla)
Morbla.freq <- sort(Morbla.freq, decreasing=TRUE)

# remove zeros
Morbla.freq <- Morbla.freq[which(Morbla.freq != 0)]

barplot(Morbla.freq, ylim=c(0,60), las = 2, main = "Morbla",
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Morbla.freq)))*sum(Morbla.freq), col = "red", lwd = 3)

ncol(Morbla)
length(Morbla.freq)
length(which(Morbla.freq > (1/length(Morbla.freq))*sum(Morbla.freq)))
length(Morbla.freq) - length(which(Morbla.freq > (1/length(Morbla.freq))*sum(Morbla.freq)))

######## Mon red ##########

Monred <- raw.data[, grep("Monophyllus", colnames(raw.data))]

Monred.freq <- rowSums(Monred)
Monred.freq <- sort(Monred.freq, decreasing=TRUE)

#remove zeros
Monred.freq <- Monred.freq[which(Monred.freq != 0)]

barplot(Monred.freq, ylim=c(0,70), las = 2, main = "Monophyllus", 
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Monred.freq)))*sum(Monred.freq), col = "red", lwd = 3)

ncol(Monred)
length(Monred.freq)
length(which(Monred.freq > (1/length(Monred.freq))*sum(Monred.freq)))
length(Monred.freq) - length(which(Monred.freq > (1/length(Monred.freq))*sum(Monred.freq)))


######## Ero bom ##########

Erobom <- raw.data[, grep("Erophylla", colnames(raw.data))]

Erobom.freq <- rowSums(Erobom)
Erobom.freq <- sort(Erobom.freq, decreasing=TRUE)

# remove zeros
Erobom.freq <- Erobom.freq[which(Erobom.freq != 0)]

barplot(Erobom.freq, ylim=c(0,8), las = 2, main = "Erophylla", 
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Erobom.freq)))*sum(Erobom.freq), col = "red", lwd = 3)

ncol(Erobom)
length(Erobom.freq)
length(which(Erobom.freq > (1/length(Erobom.freq))*sum(Erobom.freq)))
length(Erobom.freq) - length(which(Erobom.freq > (1/length(Erobom.freq))*sum(Erobom.freq)))

######## Bra cav ##########

Bracav <- raw.data[, grep("Brachyphylla", colnames(raw.data))]

Bracav.freq <- rowSums(Bracav)
Bracav.freq <- sort(Bracav.freq, decreasing=TRUE)

# remove zeros
Bracav.freq <- Bracav.freq[which(Bracav.freq != 0)]

barplot(Bracav.freq, ylim=c(0,15), las = 2, main = "Brachyphylla",
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Bracav.freq)))*sum(Bracav.freq), col = "red", lwd = 3)

ncol(Bracav)
length(Bracav.freq)
length(which(Bracav.freq > (1/length(Bracav.freq))*sum(Bracav.freq)))
length(Bracav.freq) - length(which(Bracav.freq > (1/length(Bracav.freq))*sum(Bracav.freq)))

######## Art jam ##########

Artjam <- raw.data[, grep("Artibeus", colnames(raw.data))]

Artjam.freq <- rowSums(Artjam)
Artjam.freq <- sort(Artjam.freq, decreasing=TRUE)

# remove zeros
Artjam.freq <- Artjam.freq[which(Artjam.freq != 0)]

barplot(Artjam.freq, ylim=c(0,17), las = 2, main = "Artibeus",
        cex.names = 2, cex.axis = 2)
abline(h = (1/(length(Artjam.freq)))*sum(Artjam.freq), col = "red", lwd = 3)

ncol(Artjam)
length(Artjam.freq)
length(which(Artjam.freq > (1/length(Artjam.freq))*sum(Artjam.freq)))
length(Artjam.freq) - length(which(Artjam.freq > (1/length(Artjam.freq))*sum(Artjam.freq)))


### chi-square simulation to compare RFDs

# Morbla and Ptequa
MbPq <- t(rbind.fill(as.data.frame(t(unname(Morbla.freq))), 
           as.data.frame(t(unname(Ptequa.freq)))))
MbPq[is.na(MbPq)] <- 0
# take a look
MbPq
chisq.test(x = MbPq, simulate.p.value = T, B = 10000)

# Morbla and Monred
MbMr <- t(rbind.fill(as.data.frame(t(unname(Morbla.freq))), 
                     as.data.frame(t(unname(Monred.freq)))))
MbMr[is.na(MbMr)] <- 0
# take a look
MbMr
chisq.test(x = MbMr, simulate.p.value = T, B = 10000)

# Morbla and Erobom
MbEb <- t(rbind.fill(as.data.frame(t(unname(Morbla.freq))), 
                     as.data.frame(t(unname(Erobom.freq)))))
MbEb[is.na(MbEb)] <- 0
# take a look
MbEb
chisq.test(x = MbEb, simulate.p.value = T, B = 10000)

# Morbla and Bracav
MbBc <- t(rbind.fill(as.data.frame(t(unname(Morbla.freq))), 
                     as.data.frame(t(unname(Bracav.freq)))))
MbBc[is.na(MbBc)] <- 0
# take a look
MbBc
chisq.test(x = MbBc, simulate.p.value = T, B = 10000)

# Morbla and Artjam
MbAj <- t(rbind.fill(as.data.frame(t(unname(Morbla.freq))), 
                     as.data.frame(t(unname(Artjam.freq)))))
MbAj[is.na(MbAj)] <- 0
# take a look
MbAj
chisq.test(x = MbAj, simulate.p.value = T, B = 10000)

# Monred and Ptequa
MrPq <- t(rbind.fill(as.data.frame(t(unname(Monred.freq))), 
                     as.data.frame(t(unname(Ptequa.freq)))))
MrPq[is.na(MrPq)] <- 0
# take a look
MrPq
chisq.test(x = MrPq, simulate.p.value = T, B = 10000)

# Erobom and Ptequa
EbPq <- t(rbind.fill(as.data.frame(t(unname(Erobom.freq))), 
                     as.data.frame(t(unname(Ptequa.freq)))))
EbPq[is.na(EbPq)] <- 0
# take a look
EbPq
chisq.test(x = EbPq, simulate.p.value = T, B = 10000)

# Bracav and Ptequa
BcPq <- t(rbind.fill(as.data.frame(t(unname(Bracav.freq))), 
                     as.data.frame(t(unname(Ptequa.freq)))))
BcPq[is.na(BcPq)] <- 0
# take a look
BcPq
chisq.test(x = BcPq, simulate.p.value = T, B = 10000)

# Artjam and Ptequa
AjPq <- t(rbind.fill(as.data.frame(t(unname(Artjam.freq))), 
                     as.data.frame(t(unname(Ptequa.freq)))))
AjPq[is.na(AjPq)] <- 0
# take a look
AjPq
chisq.test(x = AjPq, simulate.p.value = T, B = 10000)

# Monred and Erobom
MrEb <- t(rbind.fill(as.data.frame(t(unname(Monred.freq))), 
                     as.data.frame(t(unname(Erobom.freq)))))
MrEb[is.na(MrEb)] <- 0
# take a look
MrEb
chisq.test(x = MrEb, simulate.p.value = T, B = 10000)

# Monred and Bracav
MrBc <- t(rbind.fill(as.data.frame(t(unname(Monred.freq))), 
                     as.data.frame(t(unname(Bracav.freq)))))
MrBc[is.na(MrBc)] <- 0
# take a look
MrBc
chisq.test(x = MrBc, simulate.p.value = T, B = 10000)

# Monred and Artjam
MrAj <- t(rbind.fill(as.data.frame(t(unname(Monred.freq))), 
                     as.data.frame(t(unname(Artjam.freq)))))
MrAj[is.na(MrAj)] <- 0
# take a look
MrAj
chisq.test(x = MrAj, simulate.p.value = T, B = 10000)

# Bracav and Erobom
BcEb <- t(rbind.fill(as.data.frame(t(unname(Bracav.freq))), 
                     as.data.frame(t(unname(Erobom.freq)))))
BcEb[is.na(BcEb)] <- 0
# take a look
BcEb
chisq.test(x = BcEb, simulate.p.value = T, B = 10000)

# Artjam and Erobom
AjEb <- t(rbind.fill(as.data.frame(t(unname(Artjam.freq))), 
                     as.data.frame(t(unname(Erobom.freq)))))
AjEb[is.na(AjEb)] <- 0
# take a look
AjEb
chisq.test(x = AjEb, simulate.p.value = T, B = 10000)

# Artjam and Bracav
AjBc <- t(rbind.fill(as.data.frame(t(unname(Artjam.freq))), 
                     as.data.frame(t(unname(Bracav.freq)))))
AjBc[is.na(AjBc)] <- 0
# take a look
AjBc
chisq.test(x = AjBc, simulate.p.value = T, B = 10000)


## culebrones vs larvas

BetweenCaves <- t(rbind.fill(as.data.frame(t(unname(cul.viral.freq))), 
                     as.data.frame(t(unname(lar.viral.freq)))))
BetweenCaves[is.na(BetweenCaves)] <- 0
# take a look
BetweenCaves
chisq.test(x = BetweenCaves, simulate.p.value = T, B = 10000)

#### to calculate experimentwise error rate

# k = 15 (all species combos, except Pte por)
# alpha = 0.05

# alpha' = 1 - (1-alpha)^1/k 
## sokal and rolf

1 - ((1-0.05)^(1/15)) # = 0.0034








