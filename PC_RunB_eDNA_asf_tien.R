library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library("phyloseq")

# included: 49, 51, 52, 53, 55-66
# negative controls is sample 60 (field negative control)
path <- "~/Downloads/eDNA/Putah-Creek-eDNA/PC_RunB_fastq"
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# cutadapt
FWD <- "GTCGGTAAAACTCGTGCCAGC" # forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # reverse primer sequence
#Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
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
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

cutadapt <- "/Library/Frameworks/Python.framework/Versions/3.9/bin/cutadapt"
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
                             "-m 1", # -m 1 discards reads having a length of zero bp after cutadapting to avoid problem
                             "-n 2", # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs[i], fnRs[i])) # input files
}
# Count the presence of primers in the first cutdapt-ed sample to check if cutadapt worked as intended
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

filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))
# Set filtering parameter
# out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
#                      truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, maxLen = 215, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")

# Inspecting the returned dada-class object of the first sample:
dadaFs
dadaRs

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 10, maxMismatch = 1, verbose = TRUE)

# Inspect the merger data.frame of the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
# How many sequence variants were inferred?
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/eDNA/Putah-Creek-eDNA/12S_RDB.txt", multithread=TRUE)

taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel
# df <- as.data.frame(taxa)
# library("writexl")
# write_xlsx(df, "~/Downloads/eDNA/Putah-Creek-eDNA/PC_RunB_taxonomy.xlsx")

# phyloseq
metadata <- read.csv("PC_RunB_map.csv")
rownames(metadata) <- metadata$sample_id
rownames(seqtab.nochim) <- metadata$sample_id

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), 
               tax_table(taxa))

# filter out low abundance OTUs
get_taxa_unique(ps, 'Species', errorIfNULL=FALSE) # view species present before filtering
minTotRelAbun <- 1e-3 #set filtering threshold to 0.1%
x <- taxa_sums(ps) 
keepTaxa <- taxa_names(ps)[which((x / sum(x)) > minTotRelAbun)]
ps_filtered <- prune_taxa(keepTaxa, ps)
get_taxa_unique(ps_filtered, 'Species', errorIfNULL=FALSE) # view species present after filtering

ps_sub <- subset_samples(ps_filtered, site_code != "FDNC")

# Figures
dna <- Biostrings::DNAStringSet(taxa_names(ps_sub))
names(dna) <- taxa_names(ps_sub)
ps_sub <- merge_phyloseq(ps_sub, dna)
taxa_names(ps_sub) <- paste0("ASV", seq(ntaxa(ps_sub)))

# merges ASVs that have the same taxonomy rank (Species)
ps_species = tax_glom(ps_sub, taxrank = "Species")

# calculate relative abundance within each sample
rb = transform_sample_counts(ps_species, function(x) x / sum(x))

# convert categorical data (site_details) from character variable to factor variable
sample_data(rb)$site_details <- as.factor(sample_data(rb)$site_details)

# because some sites include more than one sample, merge samples by sites
ps_sites <- merge_samples(rb, "site_details")

# repair the merged values associated with each site after merge
sample_data(ps_sites)$site_details <- levels(sample_data(rb)$site_details)

# transform to percentages of total available.
ps_sites <- transform_sample_counts(ps_sites, function(x) x / sum(x) * 100)


# Custom palate based off some options here 
palette <- c("chartreuse4", "firebrick3", "dodgerblue1", "purple3",
             "ivory", "chocolate4", "gold1", "darkorange2",
             "darkcyan", "orchid3", "mediumvioletred", "mediumblue",
             "dimgray")

# Plot for asf
plot_bar(ps_sites, fill="Species") +
  facet_wrap(~river_mile, ncol=1, scales="free_y") + coord_flip() +
  scale_fill_manual(values=palette) +
  labs(x = "Site Details", y = "Relative Abundance (%)")

# Extra figures

custom_cbPalette <- c("#FFD580", "#cfa382", "#ffffff", "#601d2c", 
                      "#b66dff", "#490092", "#006ddb", "#22cf22",
                      "#444444", "#8f4e00", "#db6d00", "#ffdf4d",
                      "#004949", "#009999", "#fff9ba", "#ff66d1",
                      "#666600", "#68789c", "#876090", "#b0c0c9", 
                      "#00008B", "#caca00", "#f8766d", "#c7f6b6",
                      "#ffe6ee", "#920000") 

plot_bar(ps, x="site_code", fill="Species") +
  coord_flip() +
  facet_wrap(~river_mile, ncol=1, strip.position = "right", scales="free") + 
  scale_fill_manual(values=custom_cbPalette)

plot_bar(ps, x="sample_id", fill="Species") +
  facet_wrap(~river_mile, nrow=1, scales="free_x") + 
  scale_fill_manual(values=custom_cbPalette) +
  labs(x = "Sample ID", y = "Relative Abundance (%)")
 
plot_bar(ps, x="sample_id", fill="Species") +
  facet_grid(site_code~river_mile, scales="free", space="free") + 
  scale_fill_manual(values=custom_cbPalette) +
  labs(x = "Sample ID", y = "Relative Abundance (%)")

plot_bar(ps, "sampling_date", fill="Species") +
  facet_grid(site_code~river_mile, scales="free", space="free_x") + 
  scale_fill_manual(values=custom_cbPalette)

plot_bar(ps, x="sampling_date", fill="Species") + 
  facet_wrap(~site_code, scales="free_x") + 
  scale_fill_manual(values=custom_cbPalette)

ps_sub1 <- subset_samples(ps, site_code=="MBLU")
ps_sub2 <- subset_samples(ps, site_code=="RRAN")
ps_sub3 <- subset_samples(ps, site_code=="ST2D")
ps_sub4 <- subset_samples(ps, site_code=="ST2U")
ps_sub5 <- subset_samples(ps, site_code=="SCCA")
ps_sub5 <- subset_samples(ps, site_code=="FDNC")

plot_bar(ps_sub1, x="sampling_date", fill="Species") + 
  scale_fill_manual(values=custom_cbPalette)
plot_bar(ps_sub2, x="sampling_date", fill="Species") + 
  scale_fill_manual(values=custom_cbPalette)
plot_bar(ps_sub3, x="sampling_date", fill="Species") + 
  scale_fill_manual(values=custom_cbPalette)
plot_bar(ps_sub4, x="sampling_date", fill="Species") + 
  scale_fill_manual(values=custom_cbPalette)
plot_bar(ps_sub5, x="sampling_date", fill="Species") + 
  scale_fill_manual(values=custom_cbPalette)

# Shiny app
install.packages("shiny")
shiny::runGitHub("shiny-phyloseq","joey711")
