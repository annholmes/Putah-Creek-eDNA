setwd("/Users/cristina/Desktop/GVL/putah_creek")

# load phyloseq object
ps <- readRDS("ps.rds")

# Plots
library(MicEco)
ps_venn(
  ps,
  'site_code',
  fraction = 0,
  weight = FALSE, # if true, the overlaps are weighted by abundance
  relative = TRUE,
  plot = TRUE,
)

library(ggVennDiagram)
library(ggplot2)
sample_variables(ps)
rank_names(ps)

eDNA <- c('ST2D', 'ST2U')
step_weir_taxa <- subset_samples(ps, site_code %in% eDNA) # subset taxa at step weir
step_weir_taxa <- prune_taxa(taxa_sums(step_weir_taxa) > 0, step_weir_taxa) 
step_weir_vec <- get_taxa_unique(step_weir_taxa, 'species', errorIfNULL=FALSE) 
step_weir_vec <- step_weir_vec[!is.na(step_weir_vec)]

mblu_taxa <- subset_samples(ps, site_code == 'MBLU') # subset taxa at mace blvd
mblu_taxa <- prune_taxa(taxa_sums(mblu_taxa) > 0, mblu_taxa) 
mblu_vec <- get_taxa_unique(mblu_taxa, 'species', errorIfNULL=FALSE)
mblu_vec <- mblu_vec[!is.na(mblu_vec)]

rran_taxa <- subset_samples(ps, site_code == 'RRAN') # subset taxa at russell ranch
rran_taxa <- prune_taxa(taxa_sums(rran_taxa) > 0, rran_taxa) 
rran_vec <- get_taxa_unique(rran_taxa, 'species', errorIfNULL=FALSE)
rran_vec <- rran_vec[!is.na(rran_vec)]

scca_taxa <- subset_samples(ps, site_code == 'SCCA') # subset taxa at stebbins cold canyon
scca_taxa <- prune_taxa(taxa_sums(scca_taxa) > 0, scca_taxa) 
scca_vec <- get_taxa_unique(scca_taxa, 'species', errorIfNULL=FALSE)
scca_vec <- scca_vec[!is.na(scca_vec)]

jpeg("site_venn.jpeg", width = 1500, height = 1100, res = 480)
site_list <- list(step_weir_vec, mblu_vec, rran_vec, scca_vec)
ggVennDiagram(site_list, label_alpha = 0, label = "count", label_size = 4, set_size = 2.5, edge_size = 0.5, category.names = c("Step\nWeir 2", "Mace Blvd", "Russell Ranch", "Stebbins \nCold \nCanyon")) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  theme(legend.position = "none")
dev.off()

# venn diagram comparing eDNA and screw trap 
# read in screw trap data
library(readxl)
library(janitor)
library(tidyr)
screw_trap <- read_excel("screw_trap.xlsx")

# create species vector for eDNA
eDNA_taxa <- subset_samples(ps, site_code %in% eDNA) # subset taxa at step weir
eDNA_taxa <- prune_taxa(taxa_sums(eDNA_taxa) > 0, eDNA_taxa) 
eDNA_vec <- get_taxa_unique(eDNA_taxa, 'species', errorIfNULL=FALSE) 
eDNA_vec <- eDNA_vec[!is.na(eDNA_vec)]
# create species vector for screw trap
screw_vec <- screw_trap %>%
  drop_na(species) %>% 
  pull(species)
screw_vec <- as.character(screw_vec)
class(screw_vec)

# Make list of vectors comparing screw trap and eDNA
method_list <- list(eDNA_vec, screw_vec)

pdf("method_venn.pdf")
ggVennDiagram(method_list, label_alpha = 0, label = "count", label_size = 8, set_size = 7, edge_size = 0.5, category.names = c("eDNA", "Screw Trap")) +
  scale_color_manual(values = c("black", "black", "black", "black")) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# alpha diversity
jpeg("alpha.jpeg", res = 300, height = 1000, width = 1000)
plot_richness(ps, x="river_mile", measures=c("Shannon", "Simpson"), color="site_code")+
  theme_linedraw()+
  labs(x = "River Mile",
       color = "Site") +
  geom_point(size=0.8) 
dev.off()
