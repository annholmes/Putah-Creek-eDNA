#R version 4.4.1
library(here)
library(cowplot) 
library(dplyr)
library(devtools)
library(tidyverse)
library(writexl) 
library(magick)
library(ggtext)

#installing gradientForest can be a challenge, see https://github.com/z0on/RDA-forest?tab=readme-ov-file
library(gradientForest) #gradientForest version 0.1.37
library(extendedForest) #extendedForest version 1.6.2

#eDNA Index
PC_species_for_eDNA_index <- read.csv("PC_species_table_transp.csv")

rownames(PC_species_for_eDNA_index) <- PC_species_for_eDNA_index$X
PC_species_for_eDNA_index <- PC_species_for_eDNA_index[,-1]

#modified from Kelley et al. 2019 https://doi.org/10.1038/s41598-019-48546-x
#function calculates eDNA index from ASV table (which has been agglomerated to species)
eDNAINDEX <- function(x) { #where x is a dataframe with taxa/OTUs/etc in rows, and samples in columns
  rowMax <- function(x){apply(x, MARGIN = 1, FUN = max)}
  temp <- sweep(x, MARGIN = 2, STATS = colSums(x), FUN = "/") #create proportion for each taxon/OTU within each sample
  sweep(temp, MARGIN = 1, STATS = rowMax(temp), FUN = "/")
}

PC_eDNA_Index_t <- eDNAINDEX(PC_species_for_eDNA_index)
PC_eDNA_Index <- t(PC_eDNA_Index_t)

#save as an excel file to view
PC_eDNA_Index1 <- cbind(" "=rownames(PC_eDNA_Index), PC_eDNA_Index)
write_xlsx(PC_eDNA_Index1, "/Users/aholmes/Desktop/github/Putah-Creek-eDNA/Analysis/PC_eDNA_Index.xlsx")

#save as a csv file for supplementary file
write.csv(PC_eDNA_Index, "/Users/aholmes/Desktop/github/Putah-Creek-eDNA/Analysis/PC_eDNA_Index.csv")


#GF model code modified from Monuki et al. 2021 https://doi.org/10.1371/journal.pone.0253104
#eDNA Index as ASV table
ind_val <- as.data.frame(PC_eDNA_Index)

#make numeric
ind_val <- sapply(ind_val, as.numeric)

#re-apply row and col names
row.names(ind_val) <- row.names(PC_eDNA_Index)
colnames(ind_val) <- colnames(PC_eDNA_Index)

#rename columns
species_names <- colnames(PC_eDNA_Index)
easy_name <- vector(length=length(species_names))

for (i in 1:length(easy_name)) {
  easy_name[i] <- paste0("species", i)
}

#store species name pairs
species_list <- as.list(easy_name,species_names)
names(species_names) <- easy_name

#input and format metadata
metadata_GF <- read.csv("/Users/annholmes/Desktop/github/Putah-Creek-eDNA/Analysis/PC_metadata_GF.csv")
row_namers <- metadata_GF$sample_id
row.names(metadata_GF) <- row_namers
metadata_GF <- metadata_GF[,-1] 

data_GF <-merge(metadata_GF,ind_val,by="row.names")
g <- length(data_GF )
data_GF  <- data_GF [,2:g]
data_GF  <- as.data.frame(data_GF )
row.names(data_GF ) <- row_namers

ind_val <- as.data.frame(ind_val)

#Gradient Forest Model
ind_val <- ind_val %>%
  rename_all(~ make.names(.))
#remove 2 Lake Berryessa species where detections are probably from eDNA transport (n=2 species)
#remove species with 5 or fewer data points (detected in 5 or fewer samples) (n=3 species)
ind_val_21spp <- ind_val[,-c(20:23,26)] 

metadata_GF <- metadata_GF %>%
  rename_all(~ make.names(.))

data_GF  <- data_GF  %>%
  rename_all(~ make.names(.))

#identify predictor columns and response columns
preds <-colnames(metadata_GF) #predictors are km from confluence and julian date
resps <-colnames(ind_val_21spp) #responses are eDNA index values for each species by sample

#forest of 10000 regression trees for each of 21 species
eDNA_gradFor<-gradientForest(data_GF,
                             preds, 
                             resps,
                             ntree = 10000, 
                             transform = NULL, 
                             compact = F,
                             corr.threshold = 0.5, 
                             trace=T)

#plot overall importance of variables
importance(eDNA_gradFor)

#color palette
cols <- c("red3", "cornflowerblue", "gold")

tiff(file="Imp_GF.tiff",
     width = 1000, height = 1000, units = "px", pointsize = 25, res = 75)
plot(eDNA_gradFor,plot.type="Overall.Importance", col = cols, cex.axis = 0.8)
dev.off()

tiff(file="Perf_GF.tiff",
     width = 1000, height = 1000, units = "px", pointsize = 25, res = 75)
plot(eDNA_gradFor, plot.type = "Performance", show.names = F, horizontal = F,
     cex.axis = 1, cex.labels = 0.7) + abline(h=0.5, lty = 2)
dev.off()

#sort species by perfomance
sort(eDNA_gradFor$result, decreasing=T) -> sort_species
sort_species
#there are 13 species before steep drop off of R2 at 0.5
#top 13 species: 0.8626645, 0.8523772, 0.7674909, 0.7510680, 0.7071549, 0.6980784, 0.6820493, 0.6343605, 0.6220526, 0.5933452, 0.5494337, 0.5383996, 0.5933452
#remaining species: 0.4113253, 0.3533982, 0.2337443, 0.1943564, 0.1137715

#top 13 species 
#species 1 Lepomis macrochirus
anova(lm(sqrt(Lepomis.macrochirus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 2 Gasterosteus aculeatus
anova(lm(sqrt(Gasterosteus.aculeatus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 3 Cyprinus carpio
anova(lm(sqrt(Cyprinus.carpio) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 4 Lepomis microlophus
anova(lm(sqrt(Lepomis.microlophus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 5 Oncorhynchus mykiss
anova(lm(sqrt(Oncorhynchus.mykiss) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 6 Ptychocheilus grandis
anova(lm(sqrt(Ptychocheilus.grandis) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 7 Orthodon microlepidotus 
anova(lm(sqrt(Orthodon.microlepidotus ) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 8 Carassius.auratus      
anova(lm(sqrt(Carassius.auratus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 9 Cottus asper
anova(lm(sqrt(Cottus.asper) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 10 Micropterus salmoides
anova(lm(sqrt(Micropterus.salmoides) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 11 Lepomis.cyanellus
anova(lm(sqrt(Lepomis.cyanellus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 12 Catostomus.occidentalis
anova(lm(sqrt(Catostomus.occidentalis) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 13 Notemigonus.crysoleucas
anova(lm(sqrt(Notemigonus.crysoleucas) ~ river_km + julian_date + river_km:julian_date, data = data_GF))

#remaining species 
#species 14 Ameiurus.melas.nebulosus
anova(lm(sqrt(Ameiurus.melas.nebulosus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 15 Menidia.beryllina.complex
anova(lm(sqrt(Menidia.beryllina.complex) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 16 Gambusia affinis
anova(lm(sqrt(Gambusia.affinis) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 17 Hysterocarpus.traskii
anova(lm(sqrt(Hysterocarpus.traskii) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#species 18 Oncorhynchus.tshawytscha
anova(lm(sqrt(Oncorhynchus.tshawytscha) ~ river_km + julian_date + river_km:julian_date, data = data_GF))

#species not included in model
#white catfish, non-native
anova(lm(sqrt(Ameiurus.catus) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#smallmouth bass, non-native
anova(lm(sqrt(Micropterus.dolomieu) ~ river_km + julian_date + river_km:julian_date, data = data_GF))
#fathead minnow, non-native
anova(lm(sqrt(Pimephales.promelas) ~ river_km + julian_date + river_km:julian_date, data = data_GF))

#GF performance plot needs labeling to show species corresponding to numbers on x-axis
#read in plot for gradient forest performance (tiff file) that needs the legend
gf_img <- image_read(here("Manuscript_Outputs", "Putah_Ck_eDNA_Perf_GF_model.tiff"))

#convert to ggplot object
gf_plot <- ggdraw() + draw_image(gf_img)

#custom legend text with species that are significant in model (1–13) in bold
legend_text <- paste0(
  "<b>Significant species:</b><br>",
  "<b>1.</b> Lepomis macrochirus<br>",
  "<b>2.</b> Gasterosteus aculeatus<br>",
  "<b>3.</b> Cyprinus carpio<br>",
  "<b>4.</b> Lepomis microlophus<br>",
  "<b>5.</b> Oncorhynchus mykiss<br>",
  "<b>6.</b> Ptychocheilus grandis<br>",
  "<b>7.</b> Orthodon microlepidotus<br>",
  "<b>8.</b> Carassius auratus<br>",
  "<b>9.</b> Cottus asper<br>",
  "<b>10.</b> Micropterus salmoides<br>",
  "<b>11.</b> Lepomis cyanellus<br>",
  "<b>12.</b> Catostomus occidentalis<br>",
  "<b>13.</b> Notemigonus crysoleucas<br><br>",
  "<b>Not significant:</b><br>",
  "<b>14.</b> Ameiurus melas/nebulosus<br>",
  "<b>15.</b> Menidia beryllina complex<br>",
  "<b>16.</b> Gambusia affinis<br>",
  "<b>17.</b> Hysterocarpus traskii<br>",
  "<b>18.</b> Oncorhynchus tshawytscha"
)

#legend plot
legend_plot <- ggplot() +
  annotate("richtext", x = 0, y = 0.5, hjust = 0, vjust = 0.5, label = legend_text,
           fill = NA, label.color = NA, size = 4.5, lineheight = 1.2) +
  theme_void() +
  xlim(0, 1) +
  ylim(0, 1)

#combine the GF perf plot and legend 
GF_perf_legend_plot <- plot_grid(gf_plot, legend_plot, ncol = 2, rel_widths = c(2.5, 1), align = "v")

GF_perf_legend_plot

#save combined plot
save_plot("Manuscript_Outputs/Putah_Ck_GF_Perf_legend.png",
          GF_perf_legend_plot,
          base_width = 10,
          base_height = 6,
          dpi = 300,
          bg = "white")
