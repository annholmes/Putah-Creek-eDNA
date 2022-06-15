#convert GapeDNA download to RDB format

#set working directory
setwd("/Users/aholmes/Desktop/github/Putah-Creek-eDNA")

#import the csv file downloaded from GAPeDNA
FW_all <-read.csv("Freshwater_Miya_12S_World_with_sequences.csv", header = TRUE)

#remove species without sequences
FW_sequences<- FW_all[-which(is.na(FW_all$sequence)), ]

#retain only taxonomy (Genus + species) and sequences columns
FW <- subset(FW_sequences, select = c(Species,sequence))

#rename columns for format consistent with other RDB files
names(FW)[names(FW) == "Species"] <- "Binomial"
names(FW)[names(FW) == "sequence"] <- "Barcode_trimmed"

library(tidyverse)
#convert binomial with non-numberic separater into separate columns for Genus and species
FW_tax <- separate(data = FW, col = Binomial, into = c("genus", "species"), sep = "_")
head(FW_tax)

library(taxize)
#add family to genus
genus <- as.vector(unique(FW_tax$genus))
num <- as.vector(1:length(genus))
fam <- vector("character", length(genus))

#find family names using taxise seach of NCBI
for(i in num){
  t <- tax_name(query = genus[i], get = "family", db = "ncbi")
  fam[i] <- print(t$family)
}

#create df with family and genus information
fam <- unlist(fam)
temp <- data.frame(
  fam,
  genus
)

#merge df with fish df by genus
FW_fam <- full_join(FW_tax, temp, by = "genus")

#copy to family variable
FW_fam$family <- FW_fam$fam
#remove redundant column
FW_fam <- subset(FW_fam, select = -c(fam))

#add order to family
fam1 <- as.vector(unique(FW_fam$family))
num <- as.vector(1:length(fam1))
ord <- vector("character", length(fam1))

#find family names using taxise seach of NCBI
for(i in num){
  t <- tax_name(query = fam1[i], get = "order", db = "ncbi")
  ord[i] <- print(t$order)
}

#create df with family and genus information
ord <- unlist(ord)
temp1 <- data.frame(
  ord,
  fam1
)

#rename column
names(temp1)[names(temp1) == "fam1"] <- "family"
names(temp1)[names(temp1) == "ord"] <- "order"

#merge df with fish df by genus
FW_RDB <- full_join(FW_fam, temp1, by = "family")
#somehow I didn't have to do the next 2 steps?
#copy to family variable
#FW_fam$family <- FW_fam$fam
#remove redundant column
#FW_fam <- subset(FW_fam, select = -c(fam))

FW_RDB$Class <- 'Actinopterygii'
FW_RDB$Phylum <- 'Chordata'
FW_RDB$Kingdom <- '<Animalia'

#rename columns
names(FW_RDB)[names(FW_RDB) == "order"] <- "Order"
names(FW_RDB)[names(FW_RDB) == "family"] <- "Family"
names(FW_RDB)[names(FW_RDB) == "genus"] <- "Genus"
names(FW_RDB)[names(FW_RDB) == "species"] <- "Species"

#reorder columns
FW_final <- FW_RDB[,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Barcode_trimmed")]
head(FW_final)
write.csv(FW_final,"12S_FW_GAPeDNA_21Apr22.csv", row.names = FALSE)
