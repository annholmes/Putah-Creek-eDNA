library(dada2)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)
library(psadd)

library(xlsx)
library(writexl)

path <- "~/Desktop/Putah_Creek/PC_RunAB/PC_RunAB_fastq"
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# cutadapt
FWD <- "GTCGGTAAAACTCGTGCCAGC" # forward primer sequence
REV <- "CATAGTGGGGTATCTAATCCCAGTTTG" # reverse primer sequence

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna),
               Reverse = reverse(dna), RevComp = reverseComplement(dna))
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

cutadapt <- "/Users/tienly/.local/bin/cutadapt"
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
# Set filtering parameter
out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
# Check how many reads remain after filtering
out

x <- filtFs
y <- filtRs
# 0 read -> delete 13, 42 from the file path because the files do not exist
indices <- c(13, 42)
filtFs <- filtFs[-indices]
filtRs <- filtRs[-indices]

set.seed(100)
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
# Visualize the estimated error rates:
# plotErrors(errF, nominalQ = TRUE)
# plotErrors(errR, nominalQ = TRUE)

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
z <- out
out <- out[-indices,]
w <- sample.names
sample.names <- sample.names[-indices]

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
taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/Putah_Creek/PC_RunAB/12S_PC_RDB.txt", multithread=TRUE)
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
write_xlsx(df, "~/Desktop/Putah_Creek/PC_RunAB/PC_RunAB_tax1.xlsx")

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
# delete 13, 42 from the map
metadata <- read.csv("PC_RunAB_map_6Oct22.csv")
rownames(metadata) <- metadata$sample_id
rownames(seqtab.nochim) <- metadata$sample_id

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), tax_table(taxa), phy_tree(fitGTR$tree))
# Filter samples by number of reads
# 10, 17, 32, , 46, 60
ps.filt <- prune_samples(sample_sums(ps) >= 5000, ps)
# Filter ASVs by length
ps.filt <- prune_taxa(nchar(taxa_names(ps)) >= 100 & nchar(taxa_names(ps)) <= 215, ps)
# Filter ASVs by number of reads
ps.filt <- prune_taxa(taxa_sums(ps.filt) >= 100, ps.filt)
# Remove negative control
ps.filt <- subset_samples_no_zero(ps.filt, river_mile != "NA")
# Remove Homo sapiens
ps.filt <- subset_taxa(ps.filt, Species != "Homo sapiens")

# Relative abundance
# merges ASVs that have the same taxonomy rank (Species)
ps.rel1 <- tax_glom(ps.filt, taxrank = "Species")

# calculate relative abundance within each sample
ps.rel1 <- transform_sample_counts(ps.rel1, function(x) x / sum(x))

# convert categorical data (site_details) from character variable to factor variable
sample_data(ps.rel1)$site_details <- as.factor(sample_data(ps.rel1)$site_details)
sample_data(ps.rel1)$sampling_date <- as.factor(sample_data(ps.rel1)$sampling_date)

sample_data(ps.rel1)$site_date <- paste0(sample_data(ps.rel1)$site_details, sample_data(ps.rel1)$sampling_date)

# because some sites include more than one sample, merge samples by sites
ps.rel2 <- merge_samples(ps.rel1, "site_date")

# repair the merged values associated with each site after merge
sample_data(ps.rel2)$site_details <- levels(sample_data(ps.rel1)$site_details)[get_variable(ps.rel2, "site_details")]
sample_data(ps.rel2)$sampling_date <- levels(sample_data(ps.rel1)$sampling_date)[get_variable(ps.rel2, "sampling_date")]

# transform to percentages of total available.
ps.rel2 <- transform_sample_counts(ps.rel2, function(x) x / sum(x) * 100)

plot_bar(ps.rel2, x="site_details", fill="Species") +
  #facet_grid(~river_mile~sampling_date, scales="free") +
  facet_wrap(~river_mile, ncol=1, scales="free") +
  labs(title = "PC RunAB (filtered) v1", x = "Site Details", y = "Relative Abundance") + coord_flip() +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5),
        axis.text.x = element_text(angle=0))

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

