#QC, pair, merge reads; assign taxonomy; preliminary figures
#this is only Putah Creek eDNA samples from Run B (June 2021)

library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)

#dada2 part of code from Tien's PC.R file
#set path for my computer
#also change cutadapt on line 64
path <- "/Users/aholmes/desktop/github/Putah-Creek-eDNA/PC_22Apr22"
#organize files from Run B
#Putah Ck PC_47 to PC_70, except that not all are included
#Samples deleted before run due to low concentration: 47, 48
#Negative controls not sequenced due to too many: 67, 68, 69, 70
#Putah Ck sample deleted but Mock Comm sample had the same barcodes: 54
#included: 49, 51, 52, 53, 55-66
#negative controls is sample 60 (field negative control)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Remove primers with cutadapt #
FWD <- "GTCGGTAAAACTCGTGCCAGC" # forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # reverse primer sequence
#Verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients

# Calculate number of reads containing forward and reverse primer sequences (Only exact matches are found.)
# Only one set of paired end fastq.gz files will be checked (first sample in this case)
# This is sufficient, assuming all the files were created using the same library preparation.
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))
# Output interpretation
# FWD primer is found in the forward reads in its forward orientation
# FWD primer is found in the reverse reads in its reverse-complement orientation 
# REV primer is found in the reverse reads in its forward orientation 
# REV primer is found in the forward reads in its reverse-complement orientation

# Use cutadapt for primer removal: #

# This requires prior installation of cutadapt on your machine via python, anaconda, etc.
# You can do this via python's pip function:
# Install python on your machine
# Run "pip install cutadapt" in command line

# Define the path to cutadapt.exe file:
cutadapt <- "/Users/aholmes/miniconda3/envs/cutadaptenv/bin/cutadapt"
# To see if this worked and cutadapt has indeed been found, check installed version of cutadapt:
system2(cutadapt, args = "--version") # Run shell commands from R. See if R recognizes cutadapt and shows its version

# Create output filenames for the cupadapt-ed files
# Define the parameters for the cutadapt command
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

# Run cutadapt
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
# R1.flags <- paste("-g", paste0("\"", FWD, "...", REV.RC, ";optional\""))
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
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

## The primer-devoid sequence files are now ready to be analyzed. ##

# Similar to the earlier steps of reading in FASTQ files, read in the names of the cutadapt-ed FASTQ files.
# Get the matched lists of forward and reverse fastq files.

# Forward and reverse fastq filenames have the format (make sure file name patterns are correct):
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Check if forward and reverse files match:
if(length(cutFs) == length(cutRs)) print("Forward and reverse files match. Go forth and explore")
if(length(cutFs) != length(cutRs)) print("Forward and reverse files do no match. Better go back and have a check")
# Extract sample names, assuming filenames have format: SAMPLENAME.XXX.fastq.gz
sample.names <- sapply(strsplit(basename(cutFs), "\\."), `[`, 1)
sample.names

# Generate quality profile plots for our reads
# In case we have more than 20 fastq files, the following command will randomly choose 20 files to be plotted
if(length(cutFs) <= 20) {
  fwd_qual_plots <- plotQualityProfile(cutFs) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
} else {
  rand_samples <- sample(size = 20, 1:length(cutFs)) #grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(cutFs[rand_samples]) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
  rev_qual_plots <- plotQualityProfile(cutRs[rand_samples]) +
    scale_x_continuous(breaks = seq(0, 300, 20)) +
    scale_y_continuous(breaks = seq(0, 40, 5)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = 30)
}
fwd_qual_plots
rev_qual_plots

# Filter and trim #
# Assigning the directory for the filtered reads to be stored in
filtFs <- file.path(path, "filtered", basename(cutFs))
filtRs <- file.path(path, "filtered", basename(cutRs))
# Set filtering parameter
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

# Error model generation #
# dada2 learns the specific error-signature of our dataset
# This is why files from separate sequencing runs have to be processed separately till later on
set.seed(100) # set seed to ensure that randomized steps are replicable
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates:
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# Apply the dada2's core sequence-variant inference algorithm: #
# Set pool = TRUE to allow information to be shared across samples
# This makes it easier to resolve rare variants which occur just once or twice in one sample but a few more times across samples
# This will increase computational time (most likely problematic once large data sets are analyzed)
# An alternative is pseudo-pooling, an intermediate solution with medium computational time
dadaFs <- dada(filtFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(filtRs, err = errR, multithread = TRUE, pool = "pseudo")

# Inspecting the returned dada-class object of the first sample:
dadaFs
dadaRs

# Merge the forward and reverse reads together to obtain the full denoised sequence #
# Adjust the minimum overlap (default = 12) and maximum mismatch allowed (e.g = 1) if necessary
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 10, maxMismatch = 1, verbose = TRUE)

# Inspect the merger data.frame of the first sample
head(mergers[[1]])

# Construct an amplicon sequence variant table (ASV) table #
seqtab <- makeSequenceTable(mergers)
# How many sequence variants were inferred?
dim(seqtab)

# Remove chimeras #
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)
# Frequency of chimeric sequences
sum(seqtab.nochim)/sum(seqtab)

# Create a table to track read numbers throughout the pipeline: #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

# Assign taxonomy #
taxa <- assignTaxonomy(seqtab.nochim, "/Users/aholmes/desktop/github/Putah-Creek-eDNA/PC_22Apr22/12S_RDB_7Feb22.txt", multithread=TRUE)

# Inspect the assignment result:
taxa.print <- taxa  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

df <- as.data.frame(taxa)
library("writexl")
write_xlsx(df, "/Users/aholmes/desktop/github/Putah-Creek-eDNA/PC_22Apr22/putahcreek_taxonomy.xlsx")

####################move on to phyloseq
library("phyloseq")

setwd("/Users/aholmes/desktop/github/Putah-Creek-eDNA/PC_22Apr22")

#import sample metadata
PC <- read.csv("PC_RunB_map_22Apr22.csv")

samples.out <- rownames(seqtab.nochim)
rownames(PC) <- samples.out

#construct a phyloseq object from the dada2 outputs
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(PC), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="river_mile", measures=c("Shannon", "Simpson"), color="sampling_date")
plot_richness(ps, x="site_code", measures=c("Shannon", "Simpson"), color="sampling_date")

#custom palate based off some options here https://lospec.com/palette-list/vibrant-banana-8
custom_cbPalette <- c("#FFD580", "#cfa382", "#ffffff", "#601d2c", 
               "#b66dff", "#490092", "#006ddb", "#22cf22",
               "#444444", "#8f4e00", "#db6d00", "#ffdf4d",
               "#004949", "#009999", "#fff9ba", "#ff66d1",
               "#666600", "#68789c", "#876090", "#b0c0c9", 
               "#00008B", "#caca00", "#f8766d", "#c7f6b6",
               "#ffe6ee", "#920000")

#relative abundance by species
plot_bar(ps, 
         x="sample_id", fill="Species") + 
  facet_wrap(~site_code, scales="free_x") + 
  scale_fill_manual(values=custom_cbPalette)

#and some other options for relative abundance 
plot_bar(ps, 
         x="sample_id", fill="Genus") + 
        facet_wrap(~site_code, scales="free_x") + 
        scale_fill_manual(values=custom_cbPalette)

top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)
plot_bar(ps.top50, x="river_mile", fill="Species") + 
  facet_wrap(~site_code, scales="free_x") + 
  scale_fill_manual(values=custom_cbPalette)

plot_bar(ps, x="sampling_date", fill="Order") + 
        facet_wrap(~site_code, scales="free_x") +
        scale_fill_manual(values=custom_cbPalette)
