#Putah Creek eDNA study (Holmes et al. 2025)

#metabarcoding bioinformatics pipeline, taxonomic assignment, and bar plots
#R script modified from DADA2 R package tutorial (Callahan et al. 2016)

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(writexl)
library(cowplot)
library(vegan)
library(viridis)
#cutadapt python script is called from R and needs to be installed separately

path <- "~/Desktop/github/Putah-Creek-eDNA/PC_all"
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
saveRDS(ps.decontam, "~/Desktop/github/Putah-Creek-eDNA/Manuscript/final/ps.decontam.rds")
ps.decontam <- readRDS("~/Desktop/github/Putah-Creek-eDNA/PC_all/final_bioinformatics_files/ps.PC.all.rds")

#check read sums for all samples
sample_sums(ps.decontam)

#merge ASVs that have the same taxonomy rank (Species)
ps.decontam.glom <- tax_glom(ps.decontam, taxrank = "Species") 

#save as an excel file
PC_species_table=as(otu_table(ps.decontam.glom), "matrix")
if(taxa_are_rows(ps.decontam.glom)){PC_speces_table<-t(PC_species_table)}
PC_species_table_df<-as.data.frame(PC_species_table)
write_xlsx(PC_species_table_df, "~/Desktop/github/Putah-Creek-eDNA/Manuscript/final/ps_species_table.xlsx")

#data visualization
#create custom color palette
custom_cbPalette <- c("#601d2c", #A. catus
                      "#fbbf77", #A. melas/nebulosus
                      "#ffdf4d", #C. auratus
                      "#009999", #C. occidentalis
                      "#006ddb", #C. asper
                      "#fff9ba", #C. carpio
                      "#b0c0c9", #D. petenense
                      "#db6d00", #G. affinis
                      "#c5f0ee", #G. aculeatus
                      "#22cf22", #H traskii
                      "#c7f6b6", #L. cyanellus
                      "#666600", #L. gulosus
                      "#b66dff", #L. macrochirus
                      "#920000", #L. microlophus
                      "#cfa382", #M. beryllina
                      "#876090", #M. dolomieu
                      "#f8766d", #M. salmonoides
                      "#caca00", #N. crysoleucas
                      "#00008B", #O. mykiss
                      "#444444", #O. nerka
                      "#004949", #O. tshawytscha
                      "#30D5C8", #O. microlepidotus
                      "#ffe6ee", #P. macrolepida
                      "#ff66d1", #P. promelas
                      "#8f4e00", #P. nigromaculatus
                      "#006400") #P. grandis

#bar plots of sequence counts by sample (not relative abunance)
plot_bar(ps.decontam.glom, x="sample_name", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~sampling_year+sampling_event+sampling_season, space = "free", scales = "free")

## Species Relative abundance
# calculate relative abundance within each sample
ps.rel1 <- transform_sample_counts(ps.decontam.glom, function(x) x / sum(x))

#bar plots of sequences by sample, relative abunance
plot_bar(ps.rel1, x="sample_name", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_wrap(~sampling_year+sampling_event+sampling_season, scales = "free", ncol=7)+
  ylab("Relative sequence abundance")+
  xlab("Sample site and replicate by collecting event")

#for visualization, remove D. petenense and O. nerka, which only appear in 1 collecting event (#1) at low read counts
#these species are most likely eDNA transported from Lake Berryessa
ps.rel1_filter <- subset_taxa(ps.rel1, Species!= "Dorosoma petenense" & Species!= "Oncorhynchus nerka")  

#redo the palette removing 2 species
custom_cbPalette_filter <- c("#601d2c", #A. catus
                      "#fbbf77", #A. melas/nebulosus
                      "#ffd700", #C. auratus
                      "#009999", #C. occidentalis
                      "#006ddb", #C. asper
                      "#fff9ba", #C. carpio
                      "#db6d00", #G. affinis
                      "#c5f0ee", #G. aculeatus
                      "#22cf22", #H traskii
                      "#c7f6b6", #L. cyanellus
                      "#666600", #L. gulosus
                      "#b66dff", #L. macrochirus
                      "#920000", #L. microlophus
                      "#cfa382", #M. beryllina
                      "#876090", #M. dolomieu
                      "#f8766d", #M. salmonoides
                      "#caca00", #N. crysoleucas
                      "#00008B", #O. mykiss
                      "#004949", #O. tshawytscha
                      "#30D5C8", #O. microlepidotus
                      "#ffe6ee", #P. macrolepida
                      "#ff66d1", #P. promelas
                      "#8f4e00", #P. nigromaculatus
                      "#006400") #P. grandis

PC_rel_abund_bar_plots <- plot_bar(ps.rel1_filter, x="sample_name", fill="Species") +
  scale_fill_manual(values=custom_cbPalette_filter) +
  facet_wrap(~sampling_year+sampling_event+sampling_season, scales = "free", ncol=7)+
  ylab("Relative sequence abundance")+
  xlab("Sample site and replicate by collecting event")+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  guides(fill = guide_legend(ncol = 1))

#save barplots for supplementary figure
save_plot("Putah_Ck_eDNA_rel_abund_bar_plots.jpg", 
          PC_rel_abund_bar_plots, 
          dpi = 300,
          base_height = 8, #height in inches
          base_width = 11, #width in inches
          bg = "white")

#using phyloseq object from bioinformatics script
#for visualization, remove D. petenense and O. nerka, which only appear in 1 collecting event (#1) at low read counts
#these species are most likely eDNA transported from Lake Berryessa
ps.decontam.glom.filter <- subset_taxa(ps.decontam.glom, Species!= "Dorosoma petenense" & Species!= "Oncorhynchus nerka")  

#species accumulation curves
#make a matrix for the species (not really OTUs but that's what the code uses)
otu_matrix <-as.matrix(otu_table(ps.decontam.glom.filter))
df_otu <-as.data.frame(otu_matrix)

#using vegan
S <- specnumber(df_otu) # observed number of species in each of 59 field samples

#select 58 colors from viridis
colors <- viridis_pal(option = "D")(59)

#create the species accumulation curves
rarecurve(df_otu, 
            col = colors, 
            step=20, 
            lwd=3, 
            xlab= "Number of sequences", 
            ylab="Number of species detected by sample", 
            main="Species accumulaton by sequencing depth",
            label=FALSE)

#capture the base R plot
vegan_plot <- recordPlot()
#convert the recorded plot to a grob
vegan_grob <- as_grob(vegan_plot)

#save species accumulation curve plot
#Supplementary Figure S2
save_plot("Putah_Ck_eDNA_SAC.jpg", 
          vegan_grob, 
          dpi = 300,
          base_height = 5, #height in inches
          base_width = 7, #width in inches
          bg = "white")

#non-metric multidimensional scaling (NMDS) plot with stress value
m <- betadiver(df_otu)
plot(m)
## The indices
betadiver(help=TRUE)
## The basic Whittaker index
d <- betadiver(df_otu, "w")
## This should be equal to Sorensen index (binary Bray-Curtis in vegan)
range(d - vegdist(df_otu, binary=TRUE))

ordination_palette <- c("#009999", "#fbbf77", "#6f2da8", "#006ddb", "#22cf22",
                        "#db6d00", "#b0c0c9", "#ffdf4d", "#30D5C8", "#ff66d1", 
                        "#cfa382", "#601d2c", 
                        "#00008B", "#f8766d", "#c7f6b6", 
                        "#8f4e00", "#920000", "#ffe6ee", "#444444", "#b66dff", 
                        "#006400", "#c5f0ee")
sample_data(ps.decontam.glom.filter)$river_km <- as.factor(sample_data(ps.decontam.glom.filter)$river_km)

samples.ord <- ordinate(ps.decontam.glom.filter, "NMDS", "bray")

#see this link to produce a stress value
#https://ourcodingclub.github.io/tutorials/ordination/
#also here https://rpubs.com/CPEL/NMDS

dist <- vegdist(df_otu,  method = "bray")

#define a function NMDS.scree()
#performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

#use the function to choose the optimal nr of dimensions
NMDS.scree(dist)

#the final result depends on the initial random placement of the points so set a seed to make the results reproducible
set.seed(2)

#perform the final analysis and check the result
NMDS1 <- metaMDS(dist, k = 2, trymax = 100, trace = F)

#results
NMDS1
#stress value = 0.15 (stress <0.20 is ok; stress <0.10 is better)

stressplot(NMDS1)
#Results: Non-metric fit, R2 = 0.976, Linear fit, R2 = 0.883

legend_values <- c(10.6, 11.1, 16.9, 19.2, 20.8, 26.7, 34.1, 34.3, 35.1, 36.5,51.0)
formatted_values <- format(legend_values, nsmall = 1)

#plot with stress value


PC_NMDS <-plot_ordination(ps.decontam.glom.filter, samples.ord, 
                          type = "sample_id", 
                          color = "river_km", 
                          shape = "season") + 
  scale_shape_manual(values = c(17,15),
                     name="Season") +
  scale_color_manual(values = ordination_palette,
                     labels = format(as.numeric(levels(ps.decontam.glom.filter@sam_data$river_km)), nsmall = 1),#values so 51 is 51.0 in legend
                     name="Distance from confluence (km)") +
  geom_jitter(size = 8) +
  theme_bw() +
  annotate("text", 
           x = .8, y = .85, 
           size = 5,
           label = "Stress = 0.15")

save_plot("Putah_Ck_eDNA_NMDS.jpg", 
          PC_NMDS, 
          dpi = 300,
          base_height = 5.5, #height in inches
          base_width = 8, #width in inches
          bg = "white")

#PERMANOVA using adonis

sample_data(ps.decontam.glom.filter)$river_km <- as.numeric(sample_data(ps.decontam.glom.filter)$river_km)

ps.decontam_bray <-phyloseq::distance(ps.decontam.glom.filter, method="bray")
df_ps.decontam <- data.frame(sample_data(ps.decontam.glom.filter))

adonis2(ps.decontam_bray~river_km, data = df_ps.decontam)
#significant p=0.001

adonis2(ps.decontam_bray~season, data = df_ps.decontam)
#significant at p<.05 (p=0.025)

