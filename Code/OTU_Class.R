### this code inputs a sequence-by-sequence percent identity matrix
### and outputs:
#### 1) Host (i.e., site) by OTU (i.e., species) presence/absence (1/0) matrix
#### 2) Affinity propagation clusters: columns == clusters w/ row 1 == exemplar sequence

### contains functions written by Cale Basaraba


setwd("Sjodin_etal_XXXX_JOURNAL/Code")

library(tidyverse)

# input data
oral_swabs <- read.csv("../Data/Emp_HV_PID_Mat.csv")
oral_swabs$X <- as.character(oral_swabs$X)
rownames(oral_swabs) <- oral_swabs$X
oral_swabs <- select(oral_swabs, -X)
oral_swabs <- as.matrix(oral_swabs)

### the following set of code was written by Cale Basaraba 
#### Cale can be found at https://github.com/calebasaraba

##### Cale's code - creates function to sort PID matrix based on given PID#####

create_groups <- function(idMat, cutoff){
  loweroriginal <- idMat
  diag(loweroriginal)<-0
  over<-loweroriginal>=cutoff
  potential<-which(over, arr.ind = TRUE)
  columnvals<-c(potential[,2])
  columnvals<-unique(columnvals)
  taxalist<-list()
  for(j in columnvals){
    sametaxa<-c(potential[which(potential[,2]==j),1])
    finaltaxa<-c(j)
    for(i in sametaxa){
      finaltaxa[length(finaltaxa)+1]<-i
      nextlevel<-c(potential[which(potential[,2]==i),1])
      sametaxa<-c(nextlevel,sametaxa)
      sametaxa<-unique(sametaxa)
    }
    if(j==columnvals[1]){
      taxalist[[length(taxalist)+1]]<-finaltaxa
    }
    else if(!is.element(finaltaxa,taxalist[[length(taxalist)]])[1]){
      taxalist[[length(taxalist)+1]]<-finaltaxa}
  }
  
  alltaxarows<-c()
  for (i in taxalist){
    alltaxarows<-c(alltaxarows, i)
  }
  lonesomerows<-c()
  for (i in (1:length(colnames(idMat)))){
    if (i %in% alltaxarows){
    }else{
      lonesomerows <- c(lonesomerows, i)
    }
  }
  for (i in lonesomerows){
    taxalist[[length(taxalist)+1]]<-c(i)
  }
  #this is a reduction function that makes sure there is no overlap between groups
  reduceT <- function(x, input){
    b<-c()
    for(i in 1:length(input)){
      a<-intersect(x, input[[i]])
      if(length(a)>0){
        if(identical(x, a)){
          b<-c(b,x)
          b<-unique(b)
        }else{
          b<- c(b,x, input[[i]])
          b<-unique(b)
        }
      }
    }
    return(b)
  }
  
  #this is first round of reduction, just in case there is overlap
  scyther<-taxalist
  arcanine<-list()
  for(i in 1:length(scyther)){
    g <- reduceT(scyther[[i]], scyther)
    arcanine[[length(arcanine)+1]]<-g
  }
  
  # this reduces taxalists again
  blastoise<-list()
  for(i in 1:length(arcanine)){
    g <- reduceT(arcanine[[i]], arcanine)
    blastoise[[length(blastoise)+1]]<-g
  }
  blastoise<-lapply(blastoise, sort)
  blastoise<-unique(blastoise)
  taxalist<-blastoise
  # this creates the taxanames you see together
  taxanames<-c()
  for(i in 1:length(taxalist)){
    name<-paste0(cutoff,"group",i)
    taxanames[length(taxanames)+1]<-name
  }
  # this creates the output dataframe
  taxanames
  olaf<-matrix(c(0), nrow=length(taxanames),ncol=length(colnames(idMat)))
  for(i in 1:length(taxalist)){
    x<-taxalist[[i]]
    for(j in 1:length(x)){
      y<-x[j]
      olaf[i,y]<-1
    }
  }
  olaf<-data.frame(olaf)
  rownames(olaf)<-taxanames
  colnames(olaf)<-colnames(loweroriginal)
  final_groups<-olaf[,order(colnames(olaf))]
  final_groups
  results <- list(taxalist, final_groups)
}

######## end Cale's code ##########

### cutoff percentage determined using histogram

os90 <- create_groups(oral_swabs, 90)
os90.mat <- os90[[2]]

os74 <- create_groups(oral_swabs, 74)
os74.mat <- os74[[2]]

## condense names to be just unique ID number - code written by Cale Basaraba (https://github.com/calebasaraba)
os90.mat$group_name <- row.names(os90.mat)
gathered90 <- os90.mat %>%
  gather(key = sample, value = presence, starts_with("ARS"))
mutated90 <- gathered90 %>%
  mutate(sample_id = str_match(sample, "PPR\\.([0-9]{4})\\.([0-9]{2})")[,2]) %>%
  mutate(clone_id = str_match(sample, "PPR\\.([0-9]{4})\\.([0-9]{2})")[,3])
grouped90 <- mutated90 %>%
  group_by(sample_id, group_name) %>%
  summarize(presence = ifelse(sum(presence) > 0,1,0))
spreaded90 <- grouped90 %>%
  spread(key = group_name, value = presence)
final_osmat90 <- as.matrix(spreaded90[,-1])
row.names(final_osmat90) <- spreaded90$sample_id
final_osmat90 <- t(final_osmat90)
# take a look at it
final_osmat90[c(1:6), c(1:6)]
write.csv(final_osmat90, "../Data/Output/Emp_Herp_Site_x_Spp_mat_OTU90.csv")

os74.mat$group_name <- row.names(os74.mat)
gathered74 <- os74.mat %>%
  gather(key = sample, value = presence, starts_with("ARS"))
mutated74 <- gathered74 %>%
  mutate(sample_id = str_match(sample, "PPR\\.([0-9]{4})\\.([0-9]{2})")[,2]) %>%
  mutate(clone_id = str_match(sample, "PPR\\.([0-9]{4})\\.([0-9]{2})")[,3])
grouped74 <- mutated74 %>%
  group_by(sample_id, group_name) %>%
  summarize(presence = ifelse(sum(presence) > 0,1,0))
spreaded74 <- grouped74 %>%
  spread(key = group_name, value = presence)
final_osmat74 <- as.matrix(spreaded74[,-1])
row.names(final_osmat74) <- spreaded74$sample_id
final_osmat74 <- t(final_osmat74)
# take a look at it
final_osmat74[c(1:6), c(1:6)]
write.csv(final_osmat74, "../Data/Output/Emp_Herp_Site_x_Spp_mat_OTU74.csv")

########### To test out the established Herpesvirus species

# load in percent identity (PID) matrix for known herpesvirus species
### trimmed to same 175 basepair (bp) region as used 
### in the manuscript (VanDeVanter et al. 2006)
All175Pol <- read.csv("../Data/Known_HV_PID_Mat_175bp.csv")
All175Pol$X <- as.character(All175Pol$X)
rownames(All175Pol) <- All175Pol$X
All175Pol <- select(All175Pol, -X)
All175Pol <- as.matrix(All175Pol)

# load in percent identity (PID) matrix for known herpesvirus species
### using full polymerase gene sequence
AllFullPol <- read.csv("../Data/Known_HV_PID_Mat_FullPol.csv")
AllFullPol$X <- as.character(AllFullPol$X)
rownames(AllFullPol) <- AllFullPol$X
AllFullPol <- select(AllFullPol, -X)
AllFullPol <- as.matrix(AllFullPol)

# break each PID matrix into groups using same cutoff percent
### as determined with empirical sequences
All175Pol.groups90 <- create_groups(All175Pol, 90)
All175Pol.mat90 <- All175Pol.groups90[[2]]
write.csv(All175Pol.mat90, "../Data/Output/175Pol_groups90.csv")

All175Pol.groups74 <- create_groups(All175Pol, 74)
All175Pol.mat74 <- All175Pol.groups74[[2]]
write.csv(All175Pol.mat74, "../Data/Output/175Pol_groups74.csv")

AllFullPol.groups90 <- create_groups(AllFullPol, 90)
AllFullPol.mat90 <- AllFullPol.groups90[[2]]
write.csv(AllFullPol.mat90, "../Data/Output/FullPol_groups90.csv")

AllFullPol.groups74 <- create_groups(AllFullPol, 74)
AllFullPol.mat74 <- AllFullPol.groups74[[2]]
write.csv(AllFullPol.mat74, "../Data/Output/FullPol_groups74.csv")


############### affinity propagation ###############

library(apcluster)

herpes <- oral_swabs[,-c(3132,3133,3134)]
herpes[is.na(herpes)] <- 100 # gotta change diags to 100, otherwise the NAs will screw it up
herpes <- as.matrix(herpes)


osAP <- apcluster(herpes)
summary(osAP) # 45 clusters

# names(osAP@exemplars[i]) # this pulls out the exemplar of the cluster
# names(osAP@clusters[[i]]) # this will pull out all of the sequences that belong to the cluster
# str(osAP@clusters) # will show me each cluster with the number of sequences in the cluster

exemplars <- names(osAP@exemplars)
clusters <- osAP@clusters
names(clusters) <- names(osAP@exemplars)
str(clusters)

for(i in 1:length(clusters)){
  if(i==1){maxlength <- length(clusters[[i]])}
  if(i > 1){current_length <- length(clusters[[i]])
  if(current_length > maxlength){maxlength <- current_length}}
  
}

for(i in 1:length(clusters)){
  col1 <- names(clusters[[i]])
  length(col1) <- maxlength
  if(i==1){
    df <- col1
  }
  if(i > 1){
    df <- cbind(df, col1)
  }
  if(i == length(clusters)){colnames(df) <- names(clusters)}
}

head(df)
write.csv(df, file = "../Data/Output/Emp_OTUs_AffProp_45.csv")

# make it a loop to do 25 times
nboots <- 25
for(j in 1:nboots){
  osAP <- apcluster(herpes)
  exemplars <- names(osAP@exemplars)
  clusters <- osAP@clusters
  names(clusters) <- names(osAP@exemplars)
  for(i in 1:length(clusters)){
    if(i==1){maxlength <- length(clusters[[i]])}
    if(i > 1){current_length <- length(clusters[[i]])
    if(current_length > maxlength){maxlength <- current_length}}
    
  }
  for(i in 1:length(clusters)){
    col1 <- names(clusters[[i]])
    length(col1) <- maxlength
    if(i==1){
      df <- col1
    }
    if(i > 1){
      df <- cbind(df, col1)
    }
    if(i == length(clusters)){colnames(df) <- names(clusters)}
  }
  
  write.csv(df, 
            file = paste0("../Data/Output/AP_Emp_Bootstrap/OTUs_AffProp_", 
                          j, ".csv"))
}

# load them back in and see if they're the same
ap1 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_1.csv", header = T)
ap2 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_2.csv", header = T)
ap3 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_3.csv", header = T)
ap4 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_4.csv", header = T)
ap5 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_5.csv", header = T)
ap6 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_6.csv", header = T)
ap7 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_7.csv", header = T)
ap8 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_8.csv", header = T)
ap9 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_9.csv", header = T)
ap10 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_10.csv", header = T)
ap11 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_11.csv", header = T)
ap12 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_12.csv", header = T)
ap13 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_13.csv", header = T)
ap14 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_14.csv", header = T)
ap15 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_15.csv", header = T)
ap16 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_16.csv", header = T)
ap17 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_17.csv", header = T)
ap18 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_18.csv", header = T)
ap19 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_19.csv", header = T)
ap20 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_20.csv", header = T)
ap21 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_21.csv", header = T)
ap22 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_22.csv", header = T)
ap23 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_23.csv", header = T)
ap24 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_24.csv", header = T)
ap25 <- read.csv("../Data/Output/AP_Emp_Bootstrap//OTUs_AffProp_25.csv", header = T)

idmat <- matrix(nrow = 25, ncol = 25, data = NA)
idmat <- as.data.frame(idmat)
for(i in 1:nrow(idmat)){
  rownames(idmat)[i] <- paste0("ap", i)
  colnames(idmat)[i] <- paste0("ap", i)
}

for(i in 1:dim(idmat)[1]){
  for(j in 1:dim(idmat)[2]){
    mat1 <- paste0("ap", i)
    mat2 <- paste0("ap", j)
    idmat[i,j] <- identical(mat1, mat2)
  }
}

idmat <- as.matrix(idmat)
diag(idmat) <- NA
library(gdata)
tri <- upperTriangle(idmat, diag = F)
table(tri) # none are identical

exempmat <- matrix(nrow = 25, ncol = 25, data = NA)
exempmat <- as.data.frame(exempmat)
for(i in 1:nrow(exempmat)){
  rownames(exempmat)[i] <- paste0("ap", i)
  colnames(exempmat)[i] <- paste0("ap", i)
}

for(i in 1:dim(exempmat)[1]){
  for(j in 1:dim(exempmat)[2]){
    mat1 <- paste0("ap", i)
    mat2 <- paste0("ap", j)
    exempmat[i,j] <- identical(colnames(mat1), colnames(mat2))
  }
}

exempmat <- as.matrix(exempmat)
tri <- upperTriangle(exempmat, diag = F)
table(tri) # all are identical


## test default parameters of apcluster

osAP.conservative <- apcluster(herpes, q=0)
summary(osAP.conservative)
# more conservative setting gives 44 OTUs

osAP.loose <- apcluster(herpes, q=0.8, convits = 500, maxits=2000)
summary(osAP.loose)
# gives 49 OTUs w/ q=0.8
# gives 87 OTUs w/ q=0.9

