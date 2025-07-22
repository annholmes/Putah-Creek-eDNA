#data visualization

library(here)
library(phyloseq)
library(writexl)
library(cowplot)
library(vegan)
library(viridis)
library(gridGraphics)
library(tidyverse) #includes ggplot2 and dpylr

#read in file
ps.decontam <- readRDS(here("Supporting_Data_and_Resources", "ps.decontam.rds")) #decontaminated fish reads

#read sums by sample
sample_reads <-sample_sums(ps.decontam)
sample_reads_df <- data.frame(SampleID = names(sample_reads),
                              ReadCount = as.numeric(sample_reads))
write.csv(sample_reads_df, "decontaminated_fish_sample_reads.csv", row.names = FALSE)
mean(sample_reads)
median(sample_reads)
min(sample_reads)
max(sample_reads)
sum(sample_sums(ps.decontam))

#compare with
raw_ASVs <- readRDS(here("Supporting_Data_and_Resources", "seqtab.nochim.PC.rds"))
sum(raw_ASVs)

#merge ASVs that have the same taxonomy rank (Species)
ps.decontam.glom <- tax_glom(ps.decontam, taxrank = "Species") 

#Species Accumulation Curves (SAC)
#make a matrix for the species (not really OTUs but that's what the code uses)
otu_matrix <-as.matrix(otu_table(ps.decontam.glom))
df_otu <-as.data.frame(otu_matrix)

#using vegan
S <- specnumber(df_otu) # observed number of species in each of 59 field samples

#select 58 colors from viridis
colors <- viridis_pal(option = "D")(59)

#create SAC plot
rarecurve(df_otu, 
          col = colors, 
          step=20, 
          lwd=3, 
          xlab= "Number of sequences", 
          ylab="Cumulative species richness", 
          label=FALSE,
          yaxt = "n", #suppress default y-axis labels
          ylim = c(0, 20)) #y limit to 20 since curves go to 19

#customize y-axis ticks
axis(2, at = c(5,10,15,20), las = 1) 

#capture the base R plot
vegan_plot <- recordPlot()
#convert the recorded plot to a grob
vegan_grob <- as_grob(vegan_plot)

#save species accumulation curve plot
#Supplementary Figure S2
save_plot("Putah_Ck_eDNA_SAC.jpg", 
          vegan_grob, 
          dpi = 300,
          base_height = 6, #height in inches
          base_width = 6, #width in inches
          bg = "white")

#bar plots
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
plot_bar(ps.decontam.glom, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_grid(~sampling_event_no, space = "free", scales = "free")

#calculate relative abundance within each sample
ps.rel1 <- transform_sample_counts(ps.decontam.glom, function(x) x / sum(x))

#bar plots of sequences by sample, relative abunance
plot_bar(ps.rel1, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_wrap(~sampling_event_no, scales = "free", ncol=7)+
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

#using phyloseq object from bioinformatics script
#for visualization, remove D. petenense and O. nerka, which only appear in 1 collecting event (#1) at low read counts
#these species are most likely eDNA transported from Lake Berryessa
ps.decontam.glom.filter <- subset_taxa(ps.decontam.glom, Species!= "Dorosoma petenense" & Species!= "Oncorhynchus nerka")  

#bar plots of sequences by sample, relative abunance
plot_bar(ps.rel1_filter, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette) +
  facet_wrap(~sampling_event_no, scales = "free", ncol=7)+
  ylab("Relative sequence abundance")+
  xlab("Sample site and replicate by collecting event")

PC_rel_abund_bar_plots <- plot_bar(ps.rel1_filter, x="sample_id", fill="Species") +
  scale_fill_manual(values=custom_cbPalette_filter) +
  facet_wrap(~sampling_event_no, scales = "free", ncol=7)+
  ylab("Relative read abundance")+
  xlab("eDNA replicates grouped by sampling event")+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  guides(fill = guide_legend(ncol = 1))

PC_rel_abund_bar_plots

#save barplots for supplementary figure
save_plot("Putah_Ck_eDNA_rel_abund_bar_plots.jpg", 
          PC_rel_abund_bar_plots, 
          dpi = 300,
          base_height = 10, #height in inches
          base_width = 8, #width in inches
          bg = "white")

#get species richness, shannon diversity and evenness metrics based on merged replicate data

asv <- read_csv(here("Supporting_Data_and_Resources", "Putah_Creek_ASV_table_glom_metrics.csv"))

meta <- read_csv(here("Supporting_Data_and_Resources", "metadata_sampling_event.csv")) %>%
  mutate(sampling_event_no = as.integer(sampling_event_no))

status <- read_csv(here("Supporting_Data_and_Resources", "metadata_status.csv")) %>%
  mutate(species = as.character(species), 
         status = as.character(status))

glimpse(asv)
glimpse(meta)
glimpse(status)

#separate sample_id and species columns
asv_long <- asv %>%
  pivot_longer(-sample_id, names_to = "species", values_to = "reads")

#join with sampling event metadata
asv_long <- asv_long %>%
  left_join(meta, by = "sample_id")

#sum technical replicates by sampling_event_no and species
asv_merged <- asv_long %>%
  group_by(sampling_event_no, species) %>%
  summarise(reads = sum(reads), .groups = "drop")

#pivot back to wide format
asv_wide <- asv_merged %>%
  pivot_wider(names_from = species, values_from = reads, values_fill = 0)

#matrix of species read counts
asv_matrix <- asv_wide %>%
  column_to_rownames("sampling_event_no") %>%
  as.matrix()

#species richness (number of species with > 0 reads)
richness <- specnumber(asv_matrix)

#Shannon diversity
shannon <- diversity(asv_matrix, index = "shannon")

#Pielou’s evenness = Shannon / log(species richness)
evenness <- shannon / log(richness)

#combine into a single data frame
metrics <- data.frame(
  sampling_event_no = as.numeric(rownames(asv_matrix)),
  richness = richness,
  shannon = shannon,
  evenness = evenness)

#reshape species status to make it easier to join
asv_long_status <- asv_merged %>%
  left_join(status, by = "species")

#calculate richness by group
status_richness <- asv_long_status %>%
  filter(reads > 0) %>%
  group_by(sampling_event_no, status) %>%
  summarise(richness = n(), .groups = "drop") %>%
  pivot_wider(names_from = status, values_from = richness, values_fill = 0)

#join with previous metrics
metrics_all <- left_join(metrics, status_richness, by = "sampling_event_no")
metrics_all <- metrics_all[, c(1,2,5,6,3,4)] 
write_xlsx(metrics_all, "Putah_Creek_eDNA_metrics.xlsx")


#non-metric multidimensional scaling (NMDS) plot with stress value

#uses eDNA Index values, not read counts
df_edna_index <- read.csv(here("Supporting_Data_and_Resources", "PutahCk_eDNA_Index.csv"), 
                          row.names = 1)
df_edna_index <- t(df_edna_index)

#metadata 
metadata <- read.csv(here("Supporting_Data_and_Resources", "Putah_Ck_eDNA_metadata_field.csv"))

#check sample ID matching
stopifnot(all(rownames(df_edna_index) %in% metadata$sample_id))

#reorder metadata rows to match eDNA index samples exactly
metadata <- metadata %>% filter(sample_id %in% rownames(df_edna_index)) %>%
  arrange(match(sample_id, rownames(df_edna_index)))

#calculate Bray-Curtis distance matrix
dist_bray <- vegdist(df_edna_index, method = "bray")

#NMDS
set.seed(100)
nmds_result <- metaMDS(dist_bray, k = 2, trymax = 100, autotransform = FALSE, trace = FALSE)

#check stress
cat("NMDS stress:", nmds_result$stress, "\n")

#prepare NMDS scores & metadata for plotting
nmds_scores <- as.data.frame(scores(nmds_result))
nmds_scores$sample_id <- rownames(nmds_scores)

#join with metadata
plot_data <- left_join(nmds_scores, metadata, by = "sample_id")

#plot NMDS using colors for river_km (and shapes for season)
ordination_palette <- c("#009999", "#fbbf77", "#6f2da8", "#006ddb", "#22cf22",
                        "#db6d00", "#b0c0c9", "#ffdf4d", "#30D5C8", "#ff66d1", 
                        "#cfa382", "#601d2c", "#00008B", "#f8766d", "#c7f6b6", 
                        "#8f4e00", "#920000", "#ffe6ee", "#444444", "#b66dff", 
                        "#006400", "#c5f0ee")


NMDS <- ggplot(plot_data, aes(x = NMDS1, y = NMDS2, color = as.factor(river_km), shape = as.factor(season))) +
  geom_point(size = 6, alpha = 0.8) +   # bigger points here
  scale_color_manual(
    values = ordination_palette, 
    name = "Distance from\nconfluence (km)",
    labels = function(x) format(as.numeric(x), nsmall = 1)  #force a decimal place
  ) +
  scale_shape_manual(values = c(17, 15), name = "Season") +
  theme_bw() +
  theme(
    legend.position = "right",
    text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  annotate(
    "text",
    x = max(plot_data$NMDS1, na.rm = TRUE),
    y = max(plot_data$NMDS2, na.rm = TRUE),
    label = paste0("Stress: ", round(nmds_result$stress, 3)),
    hjust = 1, vjust = 1,
    size = 6
  )

save_plot(filename = "PutahCk_eDNA_NMDS.png", 
          plot = NMDS, 
          base_height = 5, 
          base_width = 8, 
          dpi = 300)

#effects of river_km and season with PERMANOVA
#convert river_km and season to numeric and factor
metadata$river_km <- as.numeric(metadata$river_km)
metadata$season <- as.factor(metadata$season)

#test terms and interaction
adonis2(dist_bray ~ river_km + season + river_km:season, data = metadata, permutations = 999, by = "terms")
