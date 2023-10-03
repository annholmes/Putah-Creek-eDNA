#Putah Creek Run all
#with negative controls included

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(psadd) #phyloseq add ons, requires package remotes to install

library(xlsx)
library(writexl)

path <- "~/Desktop/github/Putah-Creek-eDNA/PC_all"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# cutadapt
FWD <- "GTCGGTAAAACTCGTGCCAGC" # forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # reverse primer sequence

allOrients <- function(primer) {
  # Create all orientations of the input sequence
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

# Calculate number of reads containing forward and reverse primer sequences (Only exact matches are found.)
# Only one set of paired end fastq.gz files will be checked (first sample in this case)
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

cutadapt <- "/Users/aholmes/miniconda3/envs/cutadaptenv/bin/cutadapt"
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
  system2(cutadapt, args = c(R1.flags, R2.flags, # we do not change the default allowed error rate of 0.1
                             "-m", 1, # -m 1 discards reads having a length of zero bp after cutadapting to avoid problem
                             "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}

# Count the presence of primers in the first cutdapt-ed sample to check if cutadapt worked as intended
#see loop in MiSebastes R script
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# dada2
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")

# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

# Generate quality profile plots for our reads
# plotQualityProfile(cutFs[1:2])
# plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))
# Set filtering parameter max length 195 min length 150
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, 
                     maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, 
                     minLen = 150, maxLen = 195, 
                     rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

x <- filtFs
y <- filtRs

# 0 read -> delete any files that don't exist eg file a bile b
#indices <- c(file a, file b)
#filtFs <- filtFs[-indices]
#filtRs <- filtRs[-indices]

set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates:
# plotErrors(errF, nominalQ = TRUE)
# plotErrors(errR, nominalQ = TRUE)

#if there is an error in the dada command below about files not existing because reads got filtered out, run this:
filtFs <- filtFs[file.exists(filtFs)]
filtRs <- filtRs[file.exists(filtRs)]

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")
# Inspecting the returned dada-class object of the first sample:
dadaFs[[1]]
dadaRs[[1]]

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
# Inspect the merger data.frame of the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

# delete 13, 42 from out & sample.names
#z <- out
#out <- out[-indices,]
#w <- sample.names
#sample.names <- sample.names[-indices]

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track
# 0 read: 14, 37

# create RDB if needed
# RDB <- read.xlsx("12S_PC_RDB.xlsx", 1)
# write.table(RDB, file = "12S_PC_RDB.txt", sep = ";",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)
# go to RDB_final-formatting.Rmd

set.seed(100)
taxa <- assignTaxonomy(seqtab.nochim, 
                       "~/Desktop/github/Putah-Creek-eDNA/RunD/MiFishU_SFE_RDB_21Sept23.txt", 
                       multithread=TRUE)
# Inspect the assignment result
taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel
df <- as.data.frame(taxa)
df$ASV <- rownames(df)
# relative abundance
df$Count <- colSums(seqtab.nochim)
df$Relative_abundance <- df$Count / sum(df$Count)
write_xlsx(df, "~/Desktop/github/Putah-Creek-eDNA/PC_all/PC_23Sept23_taxonomy.xlsx")

# manually fill in unidentified ASVs if needed
# master <- read.xlsx("PC_RunAB_tax.xlsx", 1)
# master <- as.matrix(master)
# rownames(master) <- rownames(taxa)

# phylogenetic tree 
library(DECIPHER)
library(phangorn)
theme_set(theme_bw())

# Extract sequences from DADA2 output
sequences <- getSequences(seqtab.nochim)
names(sequences) <- sequences

# Build a neighbor-joining tree
# then fit a maximum likelihood tree
# using the neighbor-joining tree as a starting point

# Run sequence alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
#Change sequence alignment output into a phyDat structure
phang_align <- phyDat(as(alignment, 'matrix'), type='DNA')
# Create distance matrix
dm <- dist.ml(phang_align)
# Perform neighbor joining
treeNJ <- NJ(dm)  # note, tip order != sequence order
# Internal maximum likelihood
fit = pml(treeNJ, data=phang_align)
# negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2)

# phyloseq
metadata <- read.csv("~/Desktop/github/Putah-Creek-eDNA/PC_all/PutahCk_all_metadata_23Sept23.csv")
rownames(metadata) <- metadata$sample_id
rownames(seqtab.nochim) <- metadata$sample_id
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), tax_table(taxa), phy_tree(fitGTR$tree))

# subset taxa to fish only
ps.fish <- subset_taxa(ps, Class=="Actinopterygii")
#subset to field samples only
ps.field <- subset_samples(ps.fish, type=="Sample")
#remove some samples from sites that were oversampled relative to others (>3 replicates per week)
ps.samples <- subset_samples(ps.field, sample_name != "No")
#remove ASVs that don't have at least 100 reads in the entire
library(MicEco)
ps.filt.100 <- ps_prune(ps.samples, min.reads = 100)
# merge ASVs that have the same taxonomy rank (Species)
ps.glom <- tax_glom(ps.filt.100, taxrank = "Species")

#decontaminate ASV table (aka OTU table)
#remove species that don't have >500 reads in 1 rep or >10 reads in at least 2 reps in same sample
#all of the species removed are present in samples from other projects on the same run
#detections presumed to be from tag jumping or contamination 
ps.decontam <- subset_taxa(ps.glom, Genus!="Sebastes" & Genus!="Alosa" & Genus!="Lavinia" & Genus!="Pognocihthys" )
#Sebastes (rockfish) is a marine species
#Alosa and Pogonithys are more estuarine and not recorded in Putah Ck long-term data from Jacinto et al 2023
#Lavinia (Sac hitch) recorded in e-fishing data from 2020 (n=15) and closely related CA roach recorded in Jacinto et al
#However, Lavinia detections do not pass decontamination criteria (detected in 4 reps from different samples: 64, 194, 17, 51 @ SW1D, MACD, SW2D, SW2D)


# Species Relative abundance
# calculate relative abundance within each sample
ps.rel1 <- transform_sample_counts(ps.decontam, function(x) x / sum(x))

custom_cbPalette <- c("#c5f0ee", "#c30010", "#ffffff", "#601d2c", 
                      "#006ddb", "#22cf22", "#db6d00", "#ffdf4d",
                      "#004949", "#009999", "#fff9ba", "#ff66d1",
                      "#666600", "#b0c0c9", "#caca00", "#876090", 
                      "#00008B", "#cfa382", "#f8766d", "#c7f6b6", 
                      "#8f4e00", "#444444", "#ffe6ee", "#920000", 
                      "#b66dff", "#490092", "#fbbf77", "#E73009",
                      "#96E709", "#30D5C8")
plot_bar(ps.final.glom, x="sample_name", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~sampling_year+sampling_event+sampling_season, space = "free", scales = "free")
sample_sums(ps.final.glom)

plot_bar(ps.rel1, x="sample_name", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~sampling_year+sampling_event+sampling_season, space = "free", scales = "free")

OTU1=as(otu_table(ps.final.glom), "matrix")
if(taxa_are_rows(ps.final.glom)){OTU1<-t(OTU1)}
OTUdf<-as.data.frame(OTU1)
write_xlsx(OTUdf, "~/Desktop/github/Putah-Creek-eDNA/PC_all/PC_OTU_table_23Sept23.xlsx")

OTU2=as(otu_table(ps), "matrix")
if(taxa_are_rows(ps)){OTU2<-t(OTU2)}
OTU2df<-as.data.frame(OTU2)
write_xlsx(OTU2df, "~/Desktop/github/Putah-Creek-eDNA/PC_all/PC_ps_OTU_table_23Sept23.xlsx")

sample_sums(ps)
sample_sums(ps.fish)
sample_sums(ps.filt.100)
#look at only field negative controls (FNCs, n=4)
ps.nc <- subset_samples(ps.filt, type=="NC")
plot_bar(ps.nc, x="sample_id", fill="Species") #visualize abundance of sequences in FNCs
#this is fine, <50 sequences in all
#loook at only field samples
ps.samples <- subset_samples(ps.filt, type=="Sample")
plot_bar(ps.samples, x="sample_id", fill="Species")
# prune to 100 reads
ps.100 <- ps_prune(ps.samples, min.reads = 100)



#troubleshooting SP20
#loook at only field samples
ps.filt.2020 <- subset_samples(ps.filt, sampling_year=="2020")
ps.filt.SP20 <- subset_samples(ps.filt.2020, sampling_season=="Spring")
plot_bar(ps.filt.SP20, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette2) +
  facet_grid(~Fix, space = "free", scales = "free")



custom_cbPalette2 <- c("#c5f0ee", "#c30010", "#ffffff", "#601d2c", 
                      "#006ddb", "#22cf22", "#db6d00", "#ffdf4d",
                      "#004949", "#009999", "#fff9ba", "#ff66d1",
                      "#666600", "#b0c0c9", "#caca00", "#876090", 
                      "#00008B", "#cfa382", "#f8766d", "#c7f6b6", 
                      "#8f4e00", "#444444", "#ffe6ee", "#920000", 
                      "#b66dff", "#490092", "#30D5C8", "#E73009",
                      "#96E709", "#3C469E", "#9E3C4D", "#c5f0ee", "#c30010", "#ffffff", "#601d2c", 
                      "#006ddb", "#22cf22")



ps.nc.100 <- ps_prune(ps.nc, min.reads = 100)
plot_bar(ps.nc.100, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~sampling_year+sampling_season+site_code, space = "free", scales = "free")

ps.unk <- subset_samples(ps.filt, type=="Unknown")
ps.unk.100 <- ps_prune(ps.unk, min.reads = 100)
plot_bar(ps.unk.100, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~sampling_year+sampling_season+site_code, space = "free", scales = "free")
sample_sums(ps.unk.100)

plot_bar(ps.100.glom, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~river_mile+sampling_year+site_code+Season, space = "free", scales = "free")


sample_sums(ps.samples)
sample_sums(ps.100.glom)


######
# if I need to restart R w/o losing environment: .rs.restartR()
#skip this
ntaxa(ps.fnc) #141 ASVs
ps.fnc.species <-tax_glom(ps.fnc, taxrank = "Species")
ntaxa(ps.fnc.species) #34 fish species in dataset, but some are zero in FNC subset
ps.fnc.species.detects <- subset_samples_no_zero(ps.fnc.species)
ntaxa(ps.fnc.species.detects) #23 fish species in FNC data
tax_fnc_species <- tax_table(ps.fnc.species.detects)
df_tax_fnc_species <- as.data.frame(tax_fnc_species)
fnc_species <- otu_table(ps.fnc.species)
df_fnc_species <- as.data.frame(fnc_species)
colnames(df_fnc_species) <- df_tax_fnc_species$Species
#write a file to summarize FNC species detections for decontamination
write_xlsx(df_fnc_species, "~/Desktop/github/Putah-Creek-eDNA/RunsAB/species_fnc.xlsx")

#skip this
#paired samples only; min 50,000 reads, max 200,000 reads
ps.paired <- subset_samples(ps.filt, 
                            paired %in% c("FYK", "BEL", "RST") & type %in% c("sample", "maybetest") )
ps.paired.detects <- subset_samples_no_zero(ps.paired)
plot_bar(ps.paired.detects, x="sample_id", fill="Species") #visualize abundance of sequences in paired samples
ntaxa(ps.paired.detects) #87 
ps.paired.species <-tax_glom(ps.paired.detects, taxrank = "Species")
ntaxa(ps.paired.species) #28
tax_paired_species <- tax_table(ps.paired.species)
df_tax_paired_species <- as.data.frame(tax_paired_species)
paired_species <- otu_table(ps.paired.species)
df_paired_species <- as.data.frame(paired_species)
colnames(df_paired_species) <- df_tax_paired_species$Species
#write a file to summarize detection in samples paired with FYK, BEL and RST
write_xlsx(df_paired_species, "~/Desktop/github/Putah-Creek-eDNA/RunsAB/species_paired.xlsx")

#write a file to summarize detection of taxa >100 reads in all samples
tax_100 <- tax_table(ps.100)
df_tax_100 <- as.data.frame(tax_100)
otu_100 <- otu_table(ps.100)
df_otu_100 <- as.data.frame(otu_100)
colnames(df_otu_100) <- df_tax_100$Species
write_xlsx(df_otu_100, "~/Desktop/github/Putah-Creek-eDNA/RunsAB/otu_100.xlsx")

# convert categorical data (site_details) from character variable to factor variable
sample_data(ps.rel1)$Season <- as.factor(sample_data(ps.rel1)$Season)
sample_data(ps.rel1)$km_from_MD <- as.factor(sample_data(ps.rel1)$km_from_MD)


# because some sites include more than one sample, merge samples by sites
#ps.rel2 <- merge_samples(ps.rel1, "site_date")

# repair the merged values associated with each site after merge
#sample_data(ps.rel2)$site_details <- levels(sample_data(ps.rel1)$site_details)[get_variable(ps.rel2, "site_details")]
#sample_data(ps.rel2)$sampling_date <- levels(sample_data(ps.rel1)$sampling_date)[get_variable(ps.rel2, "sampling_date")]

# transform to percentages of total available.
#ps.rel2 <- transform_sample_counts(ps.rel1, function(x) x / sum(x) * 100)




#22Sept23

#Krustal Wallis
ps.species <-taxa_level(ps.100.glom, "Species")
install.packages("devtools")

library(devtools)

install_github("umerijaz/microbiomeSeq") 
library(microbiomeSeq)


#write a file to summarize rel abund of species >100 reads in all samples
tax_rel1 <- tax_table(ps.rel1)
df_tax_rel1 <- as.data.frame(tax_rel1)
otu_rel1 <- otu_table(ps.rel1)
df_otu_rel1 <- as.data.frame(otu_rel1)
colnames(df_otu_rel1) <- df_tax_rel1$Species
write_xlsx(df_otu_rel1, "~/Desktop/github/Putah-Creek-eDNA/RunsAB/otu_rel1.xlsx")

plot_bar(ps.rel2, x="site_details", fill="Species") +
  facet_wrap(~river_mile~factor(sampling_date, levels=c('12/19/2019', '5/21/2020', '5/25/2020', '5/28/2020', '10/19/2020',
                                                        '11/17/2020', '4/12/2021', '5/27/2021')), ncol=1, scales="free") +
  labs(title = "PC RunAB (filtered) v2", x = "Site Details", y = "Relative Abundance (%)") + coord_flip() +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(angle=0))

library(randomcoloR)
n <- 41
palette <- distinctColorPalette(n)

plot_bar(ps.rel2, x="site_details", fill="Species") +
  facet_grid(~river_mile~factor(sampling_date, levels=c('12/19/2019', '5/21/2020', '5/25/2020', '5/28/2020', '10/19/2020',
                                                        '11/17/2020', '4/12/2021', '5/27/2021')), scales="free") +
  labs(title = "PC RunAB (filtered) v3", x = "Site Details", y = "Relative Abundance (%)") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(angle=90)) + 
  scale_fill_manual(values=palette)

#Relative abundance added 13Mar23
# transform to percentages of total available.
ps.rel <- transform_sample_counts(ps.filt, function(x) x / sum(x) * 100)

custom_palette_41 <- c("#FFD580", "#c30010", "#ffffff", "#601d2c", 
                      "#006ddb", "#22cf22", "#db6d00", "#ffdf4d",
                      "#004949", "#009999", "#fff9ba", "#ff66d1",
                      "#666600", "#b0c0c9", "#caca00", "#876090", 
                      "#00008B", "#cfa382", "#f8766d", "#c7f6b6", 
                      "#8f4e00", "#444444", "#ffe6ee", "#920000", 
                      "#b66dff", "#490092","coral1", "cornflowerblue",
                      "cadetblue","bisque3","brown2","azure2",
                      "deeppink4","gold1","springgreen4","yellow2",
                      "palegreen4","snow2","khaki","darkolivegreen",
                      "darkorange")

plot_bar(ps.rel, x="sample_id", fill="Species") +
  facet_wrap(~site_code, nrow=1, scales="free") +
  labs(x = "Samples", y = "Relative sequnece abundance")  +
  labs(title = "Putah Creek eDNA") +
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5)) +
  scale_fill_manual(values=custom_palette_41)

plot_bar(ps.rel, x="sample_id", fill="Species") +
  facet_grid(river_mile ~ sampling_year, scales="free", space = "free") +
  labs(x = "Samples", y = "Relative sequnece abundance")  +
  labs(title = "Putah Creek eDNA") +
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5)) +
  scale_fill_manual(values=custom_palette_41)

#summary stats
library(microbiome)
library(knitr)
summarize_phyloseq(ps.filt)

#subset e-fishing sites
# ps.filt is what you get after removing Homo sapiens
efish <- subset_samples(ps.filt, paired=="Electrofishing")

# Relative abundance
# merges ASVs that have the same taxonomy rank (Species)
efish.rel1 <- tax_glom(efish, taxrank = "Species")

# calculate relative abundance within each sample
efish.rel1 <- transform_sample_counts(efish.rel1, function(x) x / sum(x))

# convert categorical data (site_details) from character variable to factor variable
sample_data(efish.rel1)$site_details <- as.factor(sample_data(efish.rel1)$site_details)
sample_data(efish.rel1)$sampling_date <- as.factor(sample_data(efish.rel1)$sampling_date)

sample_data(efish.rel1)$site_date <- paste0(sample_data(efish.rel1)$site_details, sample_data(efish.rel1)$sampling_date)

### because some sites include more than one sample, merge samples by sites
#efish.rel2 <- merge_samples(efish.rel1, "site_date")

### repair the merged values associated with each site after merge
#sample_data(ps.rel2)$site_details <- levels(sample_data(ps.rel1)$site_details)[get_variable(ps.rel2, "site_details")]
#sample_data(ps.rel2)$sampling_date <- levels(sample_data(ps.rel1)$sampling_date)[get_variable(ps.rel2, "sampling_date")]

### transform to percentages of total available.
#ps.rel2 <- transform_sample_counts(ps.rel2, function(x) x / sum(x) * 100)

#plot_bar(ps.rel2, x="site_details", fill="Species") +
  #facet_grid(~river_mile~sampling_date, scales="free") +
  facet_wrap(~river_mile, ncol=1, scales="free") +
  labs(title = "PC RunAB (filtered) v1", x = "Site Details", y = "Relative Abundance") + coord_flip() +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(angle=0))

#plot_bar(ps.rel2, x="site_details", fill="Species") + facet_wrap(~river_mile~factor(sampling_date, levels=c('12/19/2019', '5/21/2020', '5/25/2020', '5/28/2020', '10/19/2020',
                                                        '11/17/2020', '4/12/2021', '5/27/2021')), ncol=1, scales="free") +
  labs(title = "PC RunAB (filtered) v2", x = "Site Details", y = "Relative Abundance (%)") + coord_flip() +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(angle=0))

library(randomcoloR)
n <- 41
palette <- distinctColorPalette(n)

plot_bar(ps.rel2, x="site_details", fill="Species") +
  facet_grid(~river_mile~factor(sampling_date, levels=c('12/19/2019', '5/21/2020', '5/25/2020', '5/28/2020', '10/19/2020',
                                                        '11/17/2020', '4/12/2021', '5/27/2021')), scales="free") +
  labs(title = "PC RunAB (filtered) v3", x = "Site Details", y = "Relative Abundance (%)") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(angle=90)) + 
  scale_fill_manual(values=palette)

#Relative abundance added 13Mar23
# transform to percentages of total available.
ps.rel <- transform_sample_counts(ps.filt, function(x) x / sum(x) * 100)

custom_palette_41 <- c("#FFD580", "#c30010", "#ffffff", "#601d2c", 
                       "#006ddb", "#22cf22", "#db6d00", "#ffdf4d",
                       "#004949", "#009999", "#fff9ba", "#ff66d1",
                       "#666600", "#b0c0c9", "#caca00", "#876090", 
                       "#00008B", "#cfa382", "#f8766d", "#c7f6b6", 
                       "#8f4e00", "#444444", "#ffe6ee", "#920000", 
                       "#b66dff", "#490092","coral1", "cornflowerblue",
                       "cadetblue","bisque3","brown2","azure2",
                       "deeppink4","gold1","springgreen4","yellow2",
                       "palegreen4","snow2","khaki","darkolivegreen",
                       "darkorange")

plot_bar(efish.rel1, x="sample_id", fill="Species") +
  facet_wrap(~site_code, nrow=1, scales="free") +
  labs(x = "Samples", y = "Relative sequence abundance")  +
  labs(title = "Putah Creek eDNA") +
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5)) +
  scale_fill_manual(values=custom_palette_41)

summarize_phyloseq(efish.rel1) #error

plot_tree(efish.rel1,
          color = "site_code",
          justify = "yes please", 
          ladderize = "left", label.tips = "Species")

plot_tree(efish.rel1,
          color = "site_code",
          justify = "yes please", 
          ladderize = "left", label.tips = "Order")

plot_bar(ps.rel, x="sample_id", fill="Species") +
  facet_grid(river_mile ~ sampling_year, scales="free", space = "free") +
  labs(x = "Samples", y = "Relative sequence abundance")  +
  labs(title = "Putah Creek eDNA") +
  theme(plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5)) +
  scale_fill_manual(values=custom_palette_41)
