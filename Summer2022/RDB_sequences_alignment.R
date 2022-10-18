library(msa)
library(microseq)
library(writexl)

mySeqs <- readAAStringSet("alignment.fasta")

myAlignment <- msa(mySeqs, method = "ClustalW")
print(myAlignment, show="complete")
writeXStringSet(unmasked(myAlignment), file="msa.fasta")

myAlignment <- readFasta("msa.fasta")
print(str_length(myAlignment$Sequence))
myAlignment.trimmed <- msaTrim(myAlignment, gap.end = 0, gap.mid = 0)
print(str_length(myAlignment.trimmed$Sequence))
write_xlsx(myAlignment.trimmed, "~/Downloads/eDNA/Putah-Creek-eDNA/msa_trimmed.xlsx")

temp <- read.xlsx("PC_RunB_GAPeDNA-comparisons.xlsx", 1)
write.table(temp, file = "temp.txt", sep = ";",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

set.seed(100)
taxa <- assignTaxonomy(seqtab.nochim, "~/Downloads/eDNA/Putah-Creek-eDNA/temp.txt", multithread=TRUE)

taxa.print <- taxa  # Removing sequence row names for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# save dada2 results to Excel
df <- as.data.frame(taxa)
df$ASV <- rownames(taxa)
write_xlsx(df, "~/Downloads/eDNA/Putah-Creek-eDNA/results.xlsx")

library(scales)
library(ggtree)
ggtree(ps.rooted.filtered.per, ladderize = FALSE) +
  geom_text2(aes(subset=!isTip, label=label), hjust=-.2, size=4) +
  geom_tiplab(aes(label=Species), hjust=-.3) +
  geom_point(aes(x=x+hjust, color=site_code, size=Abundance),na.rm=TRUE) +
  scale_size_continuous(trans=log_trans(5)) +
  theme(legend.position="right")

ggtree(ps.rooted.filtered.per, ladderize = FALSE) + xlim(0, 1) +
  geom_tiplab(aes(label=Species), size=3)


library(metacoder)
library(phyloseq)

a <- read.csv("PC_RunB_map.csv")
rownames(a) <- a$sample_id
a$sample_id <- NULL
b <- seqtab.nochim

c <- phyloseq(otu_table(b, taxa_are_rows=FALSE), 
               sample_data(a), tax_table(master),
               phy_tree(fitGTR$tree))

phy_tree(c) <- root(phy_tree(c), taxa_names(c)[51], resolve.root = TRUE)

x <- parse_phyloseq(c)

# Plot data
set.seed(100)
heat_tree(x,
          node_label = taxon_names,
          node_size = n_obs,
          node_color = n_obs,
          node_color_axis_label = "OTU count",
          layout = "davidson-harel", initial_layout = "reingold-tilford")

