##Rank Abundance Curves for Putah Creek data

#code modified from Aviolo et al. 2019 Fig 6 see https://github.com/mavolio/RACs_paper
#see also Jacinto et al. 2023 for (colors match Fig 5, except rainbow trout is purple instead of blue)

library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(devtools)
library(codyn)

theme_set(theme_bw(12))

#df with columns for species, relative abundance, site comparison (codyn is written for a time variable but site works too), and replicate (unique within site) 
RAC_data <-read.csv("~/Desktop/github/Putah-Creek-eDNA/PC_all/Putah_Ck_spp_rel_abund_tidy.csv")
#remove cases when Rel_abund column is NA (i.e. the corredsponding species in the same row wasn't detected in that sample)
RAC_data <-RAC_dat %>% drop_na(Rel_abund)

### Making Rank Abundance Curves for 27 species
##RAC plots by sampling event

#color by native/non-native, similar to Jacinto et al. 2023 Fig5a
ggplot(data=RAC_data, aes(x=Rank, y=Rel_abund, group = interaction(Year, Season, Replicate)))+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  geom_line(color="darkgray", size=1)+
  geom_point(aes(color=Status), size=3)+
  facet_wrap(vars(RM), nrow = 2, ncol = 6)+
  scale_color_manual(values = c("blue4", "gold1"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

#all samples
#facet by site (river mile)
#color by 2 Native and 2 Non-native, following Jacinto et al. 2023 Fig5b & 5c
#prickly sculpin (C. asper)=green; rainbow trout (O. mykiss)=purple; 
#largemouth bass (M. salmonoides)=orange; common carp (C. carpio)=red
#other native and non-native species are gray #7D7D7D
ggplot(data=RAC_data, aes(x=Rank, y=Rel_abund, group = interaction(Year, Season, Replicate)))+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  geom_line(color="darkgray", size=1)+
  geom_point(aes(color=Fig_color_Jacinto), size=2)+
  facet_wrap(vars(RM), nrow = 2, ncol = 6)+
  scale_color_manual(values = c("springgreen4", "red2", "darkorange", "#7D7D7D", "#7D7D7D", "orchid4"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

#facet by site (river mile)
#color by top 2 native and top 2 non-native species detected using eDNA (following Aviolo et al. Fig 6)
#Sacramento pikeminnow (P. grandis)=beige; Sacramento sucker (C. occidentalis)=pink; 
#largemouth bass (M. salmonoides)=orange, as above; bluegill (L. macrochirus)=light blue
#other native and non-native species are gray #7D7D7D
ggplot(data=RAC_data, aes(x=Rank, y=Rel_abund, group = interaction(Year, Season, Replicate)))+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans='log10')+
  geom_line(color="darkgray", size=1)+
  geom_point(aes(color=Fig_color_top4), size=2)+
  facet_wrap(vars(RM), nrow = 2, ncol = 6)+
  scale_color_manual(values = c("deeppink", "skyblue", "darkorange", "#7D7D7D", "#7D7D7D", "chartreuse"))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())


#comparisions in ecological measures and RACs

abund_ch <- abundance_change(df = RAC_data, 
                     time.var = "sample_id", 
                     species.var = "Species",
                     abundance.var = "Rel_abund")

RAC_ch <- RAC_change(df = RAC_data, 
                     time.var = "sample_id", 
                     species.var = "Species",
                     abundance.var = "Rel_abund",
                     reference.time = "PC_007")

curve_ch <- curve_change(df = RAC_data, 
                         time.var = "sample_id", 
                         species.var = "Species",
                         abundance.var = "Rel_abund")

MRS <- rank_shift(df = RAC_data, 
                  time.var = "RM", 
                  species.var = "Species",
                  abundance.var = "Rel_abund",
                  replicate.var = "Rep_unique")


##get RAC changes
rac <- RAC_change(testyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID")

##doing curve change
cc <- curve_change(testyears, time.var = "RecYear", species.var = "GenusSpecies", abundance.var = "Abundance", replicate.var = "PlotID")

