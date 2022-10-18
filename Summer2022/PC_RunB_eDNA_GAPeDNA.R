library(dada2)

library(ShortRead)
library(Biostrings)
library(ggplot2)
library(phyloseq)

library(ape)

library(xlsx)
library(writexl)

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

set.seed(100)
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/eDNA/Putah-Creek-eDNA/12S_RDB_without-GAPeDNA.txt", multithread=TRUE)

taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel
df <- as.data.frame(taxa)
df$ASV <- rownames(taxa)
write_xlsx(df, "~/Downloads/eDNA/Putah-Creek-eDNA/results.xlsx")

# relative abundance
rownames(df) <- taxa_names(ps)
df$Count <- taxa_sums(ps)
df$Relative_abundance <- df$Count / sum(df$Count)
write_xlsx(df, "~/Downloads/eDNA/Putah-Creek-eDNA/rel_abun.xlsx")

# manually fill in unidentified ASVs
master <-read.xlsx("results.xlsx", 1)
master <- as.matrix(master)
rownames(master) <- rownames(taxa)

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
metadata <- read.csv("PC_RunB_map.csv")
rownames(metadata) <- metadata$sample_id
rownames(seqtab.nochim) <- metadata$sample_id

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(metadata), tax_table(master),
               phy_tree(fitGTR$tree))

# use short names for ASVs (e.g. ASV21) rather than the full DNA sequence 
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

plot_tree(ps, color = "sample_id", shape = "site_code", label.tips = "Species")

# filter out low abundance OTUs
# across all samples
get_taxa_unique(ps, 'Species', errorIfNULL=FALSE) # view species present before filtering
minTotRelAbun <- 1e-3 #set filtering threshold to 0.1%
x <- taxa_sums(ps) 
keepTaxa <- taxa_names(ps)[which((x / sum(x)) > minTotRelAbun)]
ps.filtered.all.samples <- prune_taxa(keepTaxa, ps)
get_taxa_unique(ps.filtered.all.samples, 'Species', errorIfNULL=FALSE) # view species present after filtering

# per sample
ps.filtered.per.sample = filter_taxa(transform_sample_counts(ps, function(x) x / sum(x)),
                   function(x) sum(x) > 1e-3, TRUE)
get_taxa_unique(ps.filtered.per.sample, 'Species', errorIfNULL=FALSE)

# root tree

# no filter
ps.rooted <- ps
phy_tree(ps.rooted) <- root(phy_tree(ps.rooted), taxa_names(ps.rooted)[51], resolve.root = TRUE)
is.rooted(phy_tree(ps.rooted))

ps.rooted <- subset_samples(ps.rooted, site_code != "FDNC")

plot_tree(ps.rooted, color="site_details", shape = "site_details", label.tips = "Species",
          ladderize = TRUE, plot.margin=0.1) +
  labs(title = "Rooted Phylogenetic Tree (No Filtering)", color="Side Details", fill="Side Details", shape = "Side Details") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.6))

# across all samples
ps.rooted.filtered.all <- ps.filtered.all.samples
phy_tree(ps.rooted.filtered.all) <- root(phy_tree(ps.rooted.filtered.all), taxa_names(ps.rooted.filtered.all)[17], resolve.root = TRUE)

ps.rooted.filtered.all <- subset_samples(ps.rooted.filtered.all, site_code != "FDNC")

plot_tree(ps.rooted.filtered.all, color="site_details", shape = "site_details", label.tips = "Species",
          ladderize = TRUE, plot.margin=0.1) +
  labs(title = "Rooted Phylogenetic Tree (Filter Across All Samples By Relative Abundance)", color="Side Details", fill="Side Details", shape = "Side Details") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5))

# per sample
ps.rooted.filtered.per <- ps.filtered.per.sample
phy_tree(ps.rooted.filtered.per) <- root(phy_tree(ps.rooted.filtered.per), taxa_names(ps.rooted.filtered.per)[48], resolve.root = TRUE)

ps.rooted.filtered.per <- subset_samples(ps.rooted.filtered.per, site_code != "FDNC")

plot_tree(ps.rooted.filtered.per, color="site_details", shape = "site_details", label.tips = "Species",
          ladderize = TRUE, plot.margin=0.1) +
  labs(title = "Rooted Phylogenetic Tree (Filter Per Sample By Relative Abundance)", color="Side Details", fill="Side Details", shape = "Side Details") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.5))

# abundance
custom_cbPalette <- c("#FFD580", "#cfa382", "#ffffff", "#601d2c", 
                      "#b66dff", "#490092", "#006ddb", "#22cf22",
                      "#444444", "#8f4e00", "#db6d00", "#ffdf4d",
                      "#004949", "#009999", "#fff9ba", "#ff66d1",
                      "#666600", "#68789c", "#876090", "#b0c0c9", 
                      "#00008B", "#caca00", "#f8766d", "#c7f6b6",
                      "#ffe6ee", "#920000") 

plot_bar(ps, x="sample_id", fill="Species") +
  facet_wrap(~site_code, nrow=1, scales="free") +
  labs(x = "Sample ID", y = "Abundance")  +
  labs(title = "No Filtering") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.7))

plot_bar(ps.filtered.all.samples, x="sample_id", fill="Species") +
  facet_wrap(~site_code, nrow=1, scales="free") +
  labs(x = "Sample ID", y = "Abundance")  +
  labs(title = "Filter Across All Samples By Relative Abundance") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.7))

plot_bar(ps.filtered.per.sample, x="sample_id", fill="Species") +
  facet_wrap(~site_code, nrow=1, scales="free") +
  labs(x = "Sample ID", y = "Relative Abundance")  +
  labs(title = "Filter Per Sample By Relative Abundance") +
  theme(plot.title = element_text(color="blue", size=18, face="bold", hjust = 0.7))
