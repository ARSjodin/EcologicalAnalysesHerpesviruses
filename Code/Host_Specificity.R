# this code inputs a presence-absence matrix of OTUs by host individual and
# outputs 1) host specificy indices for each OTU infecting >1 host (Table 3) and
# 2) abundance-occupancy relationships of OTUs (Figure 5)

library(picante)

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


# create an empty matrix that is OTU by host species
hosts <- c("Pteronotus_quadridens", "Pteronotus_portoricensis", 
           "Mormoops_blainvillii", "Artibeus_jamaicensis", 
           "Brachyphylla_cavernarum", "Monophyllus_redmani", 
           "Erophylla_sezekorni")
freq.mat <- matrix(nrow = nrow(raw.data), ncol=length(hosts))
rownames(freq.mat) <- rownames(raw.data)
colnames(freq.mat) <- hosts
freq.mat[is.na(freq.mat)] <- 0

# for each row, subset the data to only contain the positives
# for each host name, see if the name is in any of the subsetted column names
# if it is, add up how many columns have that host name, and put that value
# in the freq.mat column that pertains to the correct row and host name

for(i in 1:nrow(freq.mat)) {
  data.subset <- as.data.frame(raw.data[i,which(raw.data[i,] == 1)])
  colnames(data.subset) <- colnames(raw.data[which(raw.data[i,] == 1)])
  for(j in 1:length(hosts)){
    x <- grep(hosts[j], colnames(data.subset), value=FALSE)
    freq.mat[i,j] <- length(x)
  }
}

# make a phylogenetic distance matrix of hosts
## same species = 0
## same genus = 1
## same family = 2
## different family = 3

hpd.mat <- matrix(nrow = length(hosts), ncol = length(hosts))
rownames(hpd.mat) <- hosts
colnames(hpd.mat) <- hosts

hpd.mat[,1] <- c(0,1,2,3,3,3,3)
hpd.mat[,2] <- c(1,0,2,3,3,3,3)
hpd.mat[,3] <- c(2,2,0,3,3,3,3)  
hpd.mat[,4] <- c(3,3,3,0,2,2,2)
hpd.mat[,5] <- c(3,3,3,2,0,2,2)             
hpd.mat[,6] <- c(3,3,3,2,2,0,2)
hpd.mat[,7] <- c(3,3,3,2,2,2,0)

# calculate Std of one of the OTUs to see if it makes sense

# OTU 12 infects Pte qua, Mor bla, Art jam

2*((sum(hpd.mat[1,3], hpd.mat[1,4], hpd.mat[3,4]))/(3*2))
# answer=2.667 .... this makes sense

## make a prevalence matrix for each OTU by host
prev.mat <- matrix(nrow=nrow(raw.data), ncol=length(hosts))
colnames(prev.mat) <- hosts
rownames(prev.mat) <- row.names(raw.data)
num.sampled <- matrix(nrow = 1, ncol=length(hosts))
colnames(num.sampled) <- hosts
# hard code in the number of individuals sampled from each species
num.sampled[1,] <- c(288,31,299,57,64,262,74)

# for each row, divide the number in freq.mat by the number in num.sampled
# and place that number in prev.mat
for(i in 1:nrow(freq.mat)) {
  prev.mat[i,] <- freq.mat[i,]/num.sampled
}

## figure out which OTUs infect +1 host
### can't calculate Std* if only one host is infected

# for each row, how many columns contain non-zeros?

host.infected <- matrix(ncol=1, nrow=nrow(freq.mat))
rownames(host.infected) <- rownames(freq.mat)

for(i in 1:nrow(host.infected)){
  host.infected[i,] <- length(which(freq.mat[i,] != 0))
}

# find the placement of all those OTUs infecting more than one host
## this also names the positions w/ OTU name (i.e. multi contains OTU name and position)
multi <- which(host.infected[,1] > 1)

### now i need to calculate Std* for each OTU in multi
std.star <- matrix(ncol=1, nrow=length(multi))
rownames(std.star) <- names(multi)

# for each OTU in std.star
# for each host species pair 
# multiply the distance between pair (from hpd.mat) by the prevalence in each pair (from prev.mat)
### w * p1 * p2 
# sum the above products for each host species pair

# pte qua and pte por
std.star[1,] <- (hpd.mat[1,2]*prev.mat[1,1]*prev.mat[1,2])/
  (prev.mat[1,1]*prev.mat[1,2])
# pte qua, mor bla, and art jam
std.star[2,] <- (sum((hpd.mat[1,3]*prev.mat[4,1]*prev.mat[4,3]), 
                     (hpd.mat[1,4]*prev.mat[4,1]*prev.mat[4,4]),
                     (hpd.mat[3,4]*prev.mat[4,3]*prev.mat[4,4])))/
  (sum((prev.mat[4,1]*prev.mat[4,3]), 
       (prev.mat[4,1]*prev.mat[4,4]),
       (prev.mat[4,3]*prev.mat[4,4])))
# pte qua and mor bla
std.star[3,] <- (hpd.mat[1,3]*prev.mat[5,1]*prev.mat[5,3])/
  (prev.mat[5,1]*prev.mat[5,3])
# bra cav and mon red
std.star[4,] <- (hpd.mat[5,6]*prev.mat[7,5]*prev.mat[7,6])/
  (prev.mat[7,5]*prev.mat[7,6])
# art jam and mon red
std.star[5,] <- (hpd.mat[4,6]*prev.mat[9,4]*prev.mat[9,6])/
  (prev.mat[9,4]*prev.mat[9,6])
# pte qua, mor bla, and mon red
std.star[6,] <- (sum((hpd.mat[1,3]*prev.mat[12,1]*prev.mat[12,3]), 
                     (hpd.mat[1,6]*prev.mat[12,1]*prev.mat[12,6]),
                     (hpd.mat[3,6]*prev.mat[12,3]*prev.mat[12,6])))/
  (sum((prev.mat[12,1]*prev.mat[12,3]), 
       (prev.mat[12,1]*prev.mat[12,6]),
       (prev.mat[12,3]*prev.mat[12,6])))
# pte qua and mor bla
std.star[7,] <- (hpd.mat[1,3]*prev.mat[18,1]*prev.mat[18,3])/
  (prev.mat[18,1]*prev.mat[18,3])
# pte qua and mor bla
std.star[8,] <- (hpd.mat[1,3]*prev.mat[19,1]*prev.mat[19,3])/
  (prev.mat[19,1]*prev.mat[19,3]) 
# pte qua, mor bla, and mon red
std.star[9,] <- (sum((hpd.mat[1,3]*prev.mat[21,1]*prev.mat[21,3]), 
                     (hpd.mat[1,6]*prev.mat[21,1]*prev.mat[21,6]),
                     (hpd.mat[3,6]*prev.mat[21,3]*prev.mat[21,6])))/
  (sum((prev.mat[21,1]*prev.mat[21,3]), 
       (prev.mat[21,1]*prev.mat[21,6]),
       (prev.mat[21,3]*prev.mat[21,6])))
# pte qua and bra cav
std.star[10,] <- (hpd.mat[1,5]*prev.mat[24,1]*prev.mat[24,5])/
  (prev.mat[24,1]*prev.mat[24,5]) 
# pte qua and ero sez
std.star[11,] <- (hpd.mat[1,7]*prev.mat[28,1]*prev.mat[28,7])/
  (prev.mat[28,1]*prev.mat[28,7]) 
# pte qua, mor bla, and mon red
std.star[12,] <- (sum((hpd.mat[1,3]*prev.mat[34,1]*prev.mat[34,3]), 
                      (hpd.mat[1,6]*prev.mat[34,1]*prev.mat[34,6]),
                      (hpd.mat[3,6]*prev.mat[34,3]*prev.mat[34,6])))/
  (sum((prev.mat[34,1]*prev.mat[34,3]), 
       (prev.mat[34,1]*prev.mat[34,6]),
       (prev.mat[34,3]*prev.mat[34,6])))
# art jam and mon red
std.star[13,] <- (hpd.mat[4,6]*prev.mat[39,4]*prev.mat[39,6])/
  (prev.mat[39,4]*prev.mat[39,6])

std.star <- as.data.frame(std.star)
colnames(std.star)[1] <- "Std.star"
std.star$Subfam <- c("Beta", "Beta", "Beta", "Beta", "Beta",
                     "Beta", "Beta", "Beta", "Gamma", "Gamma",
                     "Gamma", "Beta", "Beta")
multi <- as.matrix(multi)
std.star$NumHost <- host.infected[multi,]
std.star$Pte_qua <- prev.mat[multi,1]
std.star$Pte_por <- prev.mat[multi,2]
std.star$Mor_bla <- prev.mat[multi,3]
std.star$Art_jam <- prev.mat[multi,4]
std.star$Bra_cav <- prev.mat[multi,5]
std.star$Mon_red <- prev.mat[multi,6]
std.star$Ero_sez <- prev.mat[multi,7]

std.star
write.csv(std.star, "../Data/Output/Std.Star.csv")

# compare specificity between species

# is OTU 1 prevalence significantly greater in Pte qua?
f1 <- fisher.test(x = matrix(c(freq.mat[1,1], freq.mat[1,2], 
                               num.sampled[1,1], num.sampled[1,2]), 
                             byrow = TRUE, 2, 2), alternative = "greater")
f1$p.value #no

# is OTU 2 prevalence significantly greater in Pte qua? (row 12)
f2 <- fisher.test(x = matrix(c(freq.mat[12,1], sum(freq.mat[12,3], freq.mat[12,6]),
                               num.sampled[1,1], sum(num.sampled[1,3], num.sampled[1,6])),
                             byrow = TRUE, 2, 2), alternative = "greater")

f2$p.value #yes

# is OTU 4 prevalence significantly greater in Pte qua (row 34)
row.names(freq.mat)[34]
f4 <- fisher.test(x = matrix(c(freq.mat[34,1], sum(freq.mat[34,3], freq.mat[34,6]),
                               num.sampled[1,1], sum(num.sampled[1,3], num.sampled[1,6])),
                             byrow = TRUE, 2, 2), alternative = "greater")
f4$p.value #yes

# is OTU 7 prevalence significantly greater in Mon red (row 39) vs. art jam
row.names(freq.mat)[41]
f7 <- fisher.test(x = matrix(c(freq.mat[41,6], freq.mat[41,4], 
                               num.sampled[1,6], num.sampled[1,4]), 
                             byrow = TRUE, 2,2), alternative = "greater")
f7$p.value #no

# is OTU 12 prevalence significantly greater in Mor bla (row 4) vs. pte qua, art jam
row.names(freq.mat)[4]
f12 <- fisher.test(x = matrix(c(freq.mat[4,3], sum(freq.mat[4,1], freq.mat[4,4]),
                                num.sampled[1,3], sum(num.sampled[1,1], num.sampled[1,4])),
                              byrow = TRUE, 2, 2), alternative = "greater")
f12$p.value #yes

# is OTU 13 prevalence significantly greater in Mor bla (row 5) vs. pte qua
f13 <- fisher.test(x = matrix(c(freq.mat[5,3], freq.mat[5,1], 
                                num.sampled[1,3], num.sampled[1,1]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f13$p.value #yes

# is OTU 15 prevalence significantly greater in Mon red (row 7) vs. Bra cav
f15 <- fisher.test(x = matrix(c(freq.mat[7,6], freq.mat[7,5], 
                                num.sampled[1,6], num.sampled[1,5]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f15$p.value #no

# is OTU 17 prevalence significantly greater in Art jam (row 9) vs. Mon red
f17 <- fisher.test(x = matrix(c(freq.mat[9,4], freq.mat[9,6], 
                                num.sampled[1,4], num.sampled[1,6]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f17$p.value #yes

# is OTU 25 prevalence significantly greater in Mor bla (row 18) vs. pte qua
f25 <- fisher.test(x = matrix(c(freq.mat[18,3], freq.mat[18,1], 
                                num.sampled[1,3], num.sampled[1,1]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f25$p.value #yes

# is OTU 26 prevalence significantly greater in Mor bla (row 19) vs. pte qua
f26 <- fisher.test(x = matrix(c(freq.mat[19,3], freq.mat[19,1], 
                                num.sampled[1,3], num.sampled[1,1]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f26$p.value #no

# is OTU 28 prevalence significantly greater in Mon red (row 21) vs. pte qua, mor bla
f28 <- fisher.test(x = matrix(c(freq.mat[21,6], sum(freq.mat[21,1], freq.mat[21,3]),
                                num.sampled[1,6], sum(num.sampled[1,1], num.sampled[1,3])),
                              byrow = TRUE, 2, 2), alternative = "greater")
f28$p.value #yes

# is OTU 30 prevalence significantly greater in Bra cav (row 24) vs. pte qua
row.names(freq.mat)[24]
f30 <- fisher.test(x = matrix(c(freq.mat[24,5], freq.mat[24,1], 
                                num.sampled[1,5], num.sampled[1,1]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f30$p.value #yes

# is OTU 34 prevalence significantly greater in Ero sez (row 28) vs. pte qua
f34 <- fisher.test(x = matrix(c(freq.mat[28,7], freq.mat[28,1], 
                                num.sampled[1,7], num.sampled[1,1]), 
                              byrow = TRUE, 2,2), alternative = "greater")
f34$p.value #yes

#### to look at abundance-occupancy relationships
# add a column to prev.mat that says how many hosts each otu infects
prev.mat2 <- as.data.frame(prev.mat)
avg.prev <- rowMeans(prev.mat)
num.hosts <- rep(NA, 43)
for(i in 1:nrow(prev.mat2)){
  num.hosts[i] <- as.numeric(length(which(prev.mat2[i,] > 0)))
  prev.mat2$num_hosts <- num.hosts
}
for(i in 1:nrow(prev.mat2)){
  data.subset <- prev.mat[i,]
  prev.mat2$AvgPrev[i] <- mean(data.subset[which(data.subset > 0)])
}
prev.mat2$AvgPrevAllInds <- avg.prev

# Figure 5a
plot(x = prev.mat2$num_hosts, y = prev.mat2$AvgPrevAllInds, ylab = "Average prevalence \n(all hosts)", 
     xlab = "Number of infected host species", pch = 19, 
     cex.lab = 1.5, cex.axis = 1.75, cex = 1.5)
cor.test(prev.mat2$num_hosts, prev.mat2$AvgPrevAllInds)
# Figure 5b
plot(x = prev.mat2$num_hosts, y = prev.mat2$AvgPrev, ylab = "Average prevalence \n(infected hosts)", 
     xlab = "Number of infected host species", pch = 19, 
     cex.lab = 1.5, cex.axis = 1.75, cex = 1.5)
cor.test(prev.mat2$num_hosts, prev.mat2$AvgPrev)
