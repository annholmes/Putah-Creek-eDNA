#Putah Creek eDNA Project
#Spatiotemporal stability of fish communities in a regulated stream: Insights from environmental DNA (Holmes et al.)

#bioinformatics pipeline for eDNA metabarcoding data
#taxonomic assignment using custom RDB
#some basic plots
#some R script is modified from DADA2 tutorial (Callahan et al. 2016)
#python script for cutadapt (Martin 2011) is called here but needs to be installed separately

library(here)
library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(writexl)
library(cowplot)
library(vegan)
library(viridis)

path <- "~/Raw_data"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

#provide primer sequences for cutadapt
FWD <- "GTCGGTAAAACTCGTGCCAGC" # MiFish-U-F forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # MiFish-U-R reverse primer sequence 27 bp

allOrients <- function(primer) {
  #creates all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

#calculate number of reads containing forward and reverse primer sequences (only exact matches)
#check the first set of paired end fastq.gz files 
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

for(a in 1:length(fnFs)) {
  b <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[a]]), 
             FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[a]]), 
             REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[a]]), 
             REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[a]]))
  print(b)
}

#trim primers by sequence using cutadapt
cutadapt <- Sys.getenv("CUTADAPT_PATH") #path is in .Renviron file
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, #default allowed error rate of 0.1
                             "-m", 1, # -m 1 discards reads having a length of zero bp after primers removed
                             "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

#count primers in the first cutdapt-ed sample to check cutadapt worked as intended
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#dada2
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")

#extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

#generate quality profile plots as a quality check
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#preliminary filterAndTrim with no length trimming retains off target amplification (mostly >195 bp, but some <150 bp)
out_prelim <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, 
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out_prelim

#final filtering parameter is set to max length 195 and min length 150 
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, 
                     minLen = 150, maxLen = 195, 
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)

#check number of reads after filtering
out

x <- filtFs
y <- filtRs

#calculate error rates
set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#visualize estimated error rate
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")

#inspect the returned dada-class object of the first sample
#denoising results, i.e., the number of ASVs inferred from the sequence data
dadaFs[[1]]
dadaRs[[1]]

#merge forward and reverse reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

#inspect the merged data.frame of the first sample
head(mergers[[1]])

#create the ASV table and remove chimeric sequences
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim.PC <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim.PC)
table(nchar(getSequences(seqtab.nochim.PC)))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim.PC))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

#save DADA2 output
saveRDS(seqtab.nochim.PC, "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/seqtab.nochim.PC.rds")
#can import here later if needed
seqtab.nochim.PC <- readRDS("~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/seqtab.nochim.PC.rds")

#optimizing taxonomic assignment 
#first testing bootstrap value 50 (default), then test 80
#min boot 50 using custom regional reference sequence database (modified from Nagarajan et al. 2023)
set.seed(100)
taxa_minboot50 <- assignTaxonomy(seqtab.nochim.PC, 
                                 "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/12S_MiFish_RDB_PutahCK_eDNA.txt", 
                                 minBoot = 50, #default value
                                 multithread=TRUE,
                                 outputBootstraps=TRUE,
                                 tryRC=TRUE)
#inspect the min boot 50 assignment results
taxa_minboot50.print <- taxa_minboot50  #removing sequence row names for display only
rownames(taxa_minboot50.print) <- NULL
head(taxa_minboot50.print)
#save min boot 50 taxonomy results to Excel
df_minboot50 <- as.data.frame(taxa_minboot50)
df_minboot50$ASV <- rownames(df_minboot50)
#relative abundance
df_minboot50$Count <- colSums(seqtab.nochim.PC)
df_minboot50$Relative_abundance <- df_minboot50$Count / sum(df_minboot50$Count)
write_xlsx(df_minboot50, "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/PC_taxonomy_minboot50_26Oct24.xlsx")
#assign taxonomy at min boot 80
set.seed(100)
taxa_minboot80 <- assignTaxonomy(seqtab.nochim.PC, 
                       "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/12S_MiFish_RDB_PutahCK_eDNA.txt", 
                       minBoot = 80,
                       multithread=TRUE,
                       outputBootstraps=TRUE,
                       tryRC=TRUE)
#inspect the min boot 80 assignment results
taxa_minboot80.print <- taxa_minboot80  #removing sequence row names for display only
rownames(taxa_minboot80.print) <- NULL
head(taxa_minboot80.print)
#save min boot 80 taxonomy results to Excel
df_minboot80 <- as.data.frame(taxa_minboot80)
df_minboot80$ASV <- rownames(df_minboot80)
#relative abundance
df_minboot80$Count <- colSums(seqtab.nochim.PC)
df_minboot80$Relative_abundance <- df_minboot80$Count / sum(df_minboot80$Count)
write_xlsx(df_minboot80, "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/PC_taxonomy_minboot80.xlsx")

#final taxonomic assignment
#assign taxonomy based results above (minboot50 generally provides accurate assignment for this dataset with this RDB)
#can't use previous output as taxa file because including the bootstrap values in the output changes the file structure
set.seed(100)
taxa <- assignTaxonomy(seqtab.nochim.PC, 
                       "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/12S_MiFish_RDB_PutahCK_eDNA.txt", 
                       minBoot = 50, #default value
                       multithread=TRUE,
                       tryRC=TRUE)
# Inspect the assignment result
taxa.print <- taxa  
rownames(taxa.print) <- NULL #remove sequence row names for display only
head(taxa.print)

#save dada2 taxonomic assignment results to Excel
df <- as.data.frame(taxa)
df$ASV <- rownames(df)
# relative abundance
df$Count <- colSums(seqtab.nochim.PC)
df$Relative_abundance <- df$Count / sum(df$Count)
write_xlsx(df, "~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/PC_taxonomy_all_taxa.xlsx")

#import decontaminated ASV and taxa tables
seqtab.nochim.decontam <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/Manuscript/final/seqtab.nochim.decontam.csv", row.names = 1)
seqtab.nochim.decontam <- data.matrix(seqtab.nochim.decontam)
taxa.decontam <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/Manuscript/final/taxa.decontam.csv", row.names = 1)
taxa.decontam <- as.matrix(taxa.decontam)

#phyloseq
#read sample metadata
metadata <- read.csv("~/Desktop/github/Putah-Creek-eDNA/Manuscript/final/PC_metadata_all.csv")
#remove negative control samples in metadata so it will match decontaminated OTU table
metadata_field <- metadata[-c(60:69), ]
rownames(metadata_field) <- metadata_field$sample_id

#create phyloseq object
ps.decontam <- phyloseq(otu_table(seqtab.nochim.decontam, taxa_are_rows=FALSE), 
               sample_data(metadata_field), 
               tax_table(taxa.decontam))

#save phyloseq object
saveRDS(ps.decontam, "~/Supporting_Data_and_Resources/ps.decontam.rds")

#read in RDS
ps.decontam <- readRDS(here("Supporting_Data_and_Resources", "ps.decontam.rds"))

#check read sums for all samples
sample_sums(ps.decontam)

#merge ASVs that have the same taxonomy rank (Species)
ps.decontam.glom <- tax_glom(ps.decontam, taxrank = "Species") 

#save as an excel file
PC_species_table=as(otu_table(ps.decontam.glom), "matrix")
if(taxa_are_rows(ps.decontam.glom)){PC_speces_table<-t(PC_species_table)}
PC_species_table_df<-as.data.frame(PC_species_table)
write_xlsx(PC_species_table_df, "~/Desktop/github/Putah-Creek-eDNA/Manuscript/final/ps_species_table.xlsx")
