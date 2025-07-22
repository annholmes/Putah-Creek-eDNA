#allluvial plots for gear comp

library(ggplot2)
library(ggalluvial)
library(dplyr)
library(cowplot)
setwd("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp")

fyke <- read.csv("~/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_fyke.csv")
RST_April <- read.csv("~/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_RST_April.csv")
RST_May <- read.csv("~/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_RST_May.csv")

MACU_FA20$Names_by_family <- factor(MACU_FA20$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))

#fyke alluvial plot
fyke_plot <-ggplot(fyke,
       aes(x = Method, 
           stratum = Names, 
           alluvium = Names,
           y = Relative_abundance,
           fill = Names, 
           label = Common_name_label)) + 
  scale_fill_manual(values = c("bisque",#bigscale logperch BLP non-native"
                               "deeppink4", #bluegill BGL non-native
                               "aliceblue", #Chinook salmon CHN  native
                               "salmon", #common carp CRP non-native
                               "goldenrod", #green sunfish GSL non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "purple", #Pacific lamprey PLM native
                               "blue", #prickly sculpin PSC native
                               "darkcyan", #rainbow trout RBT native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "lightblue4", #Sacramento sucker SKR native
                               "deeppink", #smallmouth bass SMB non-native
                               "aquamarine3", #threespine stickleback TSS native
                               "cyan"))+ #tule perch TUP native
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16)) +
  labs(y = "Relative abundance") 
  
fyke_plot

#RST alluvial plot - April

RST_April_plot <- ggplot(RST_April,
       aes(x = Method, 
           stratum = Names, 
           alluvium = Names,
           y = Relative_abundance,
           fill = Names, 
           label = Common_name_label)) + 
  scale_fill_manual(values = c("bisque",#bigscale logperch BLP non-native"
                               "deeppink4", #bluegill BGL non-native
                               "aliceblue", #Chinook salmon CHN  native
                               "salmon", #common carp CRP non-native
                               "goldenrod", #green sunfish GSL non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "purple", #Pacific lamprey PLM native
                               "blue", #prickly sculpin PSC native
                               "darkcyan", #rainbow trout RBT native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "lightblue4", #Sacramento sucker SKR native
                               "deeppink", #smallmouth bass SMB non-native
                               "aquamarine3", #threespine stickleback TSS native
                               "cyan"))+ #tule perch TUP native
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width = 2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16)) +
  labs(y = "Relative abundance") 

RST_April_plot

#RST alluvial plot - May

RST_May_plot <- ggplot(RST_May,
                       aes(x = Method, 
                           stratum = Names, 
                           alluvium = Names,
                           y = Relative_abundance,
                           fill = Names, 
                           label = Common_name_label)) + 
  scale_fill_manual(values = c("bisque",#bigscale logperch BLP non-native"
                               "deeppink4", #bluegill BGL non-native
                               "aliceblue", #Chinook salmon CHN  native
                               "salmon", #common carp CRP non-native
                               "goldenrod", #green sunfish GSL non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "purple", #Pacific lamprey PLM native
                               "blue", #prickly sculpin PSC native
                               "darkcyan", #rainbow trout RBT native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "lightblue4", #Sacramento sucker SKR native
                               "deeppink", #smallmouth bass SMB non-native
                               "aquamarine3", #threespine stickleback TSS native
                               "cyan"))+ #tule perch TUP native
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

RST_May_plot

#to extract legend
legend <- ggplot(RST_May,
       aes(x = Method, 
           stratum = Names, 
           alluvium = Names,
           y = Relative_abundance,
           fill = Names, 
           label = Common_name_label)) + 
  scale_fill_manual(values = c("bisque",#bigscale logperch BLP non-native"
                               "deeppink4", #bluegill BGL non-native
                               "cadetblue1", #Chinook salmon CHN  native
                               "salmon", #common carp CRP non-native
                               "darkgoldenrod4", #green sunfish GSL non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "purple", #Pacific lamprey PLM native
                               "blue", #prickly sculpin PSC native
                               "darkcyan", #rainbow trout RBT native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "lightblue4", #Sacramento sucker SKR native
                               "orange", #smallmouth bass SMB non-native
                               "aquamarine3", #threespine stickleback TSS native
                               "cyan"))+ #tule perch TUP native
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.key = element_rect(color="black", fill = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        text=element_text(size=16)) +
  guides(override.aes = list(linetype = c(0,1,1,1)))

#extract legend from last plot
grobs <- ggplotGrob(legend)$grobs
legend_plot <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

# build grid without legends
gear_comp_traps <-plot_grid(fkye_plot, RST_April_plot, RST_May_plot, 
                            labels = c('A', 'B', 'C'), 
                            label_size = 18,
                            ncol = 1,
                            scale = 0.95)
gear_comp_traps
# add legend
gear_comp_traps_w_legend <- plot_grid(gear_comp_traps, legend_plot, 
                                      ncol = 2,
                                      rel_widths = c(1, .6))
gear_comp_traps_w_legend

save_plot("~/Desktop/github/Putah-Creek-eDNA/gear_comp/gear_comp_traps_w_legend.jpg", 
          gear_comp_traps_w_legend, 
          dpi = 300,
          base_height = 12, #height in inches, equal to 30.48cm 
          base_width = 13.25, 
          bg = "white")

#electrofishing survey comps

#MACU 2020
MACU_FA20 <- read.csv("~/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_MACUFA20.csv")

#to order the species names for the legend
MACU_FA20$Names_by_family <- factor(MACU_FA20$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                 "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                 "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                 "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                 "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                 "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                 "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                 "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                 "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                 "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                 "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                 "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                 "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                 "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                 "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                 "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                 "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                 "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                 "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                 "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                 "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                 "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                 "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                 "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                 "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                 "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
MACU_FA20_plot <-ggplot(MACU_FA20, aes(x = Method, 
                     stratum = Names_by_family, 
                     alluvium = Names_by_family,
                     y = Relative_abundance,
                     fill = Names_by_family, 
                     label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

MACU_FA20_plot
  
#MACU 2021
MACU_FA21 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_MACUFA21.csv")

#to order the species names for the legend
MACU_FA21$Names_by_family <- factor(MACU_FA21$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
MACU_FA21_plot <-ggplot(MACU_FA21, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+ 
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

MACU_FA21_plot

#OLDR20
OLDR_FA20 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_OLDRFA20.csv")

#to order the species names for the legend
OLDR_FA20$Names_by_family <- factor(OLDR_FA20$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
OLDR_FA20_plot <-ggplot(OLDR_FA20, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+ 
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

OLDR_FA20_plot

#FRPG21
FRPG_FA21 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_FRPGFA21.csv")

#to order the species names for the legend
FRPG_FA21$Names_by_family <- factor(FRPG_FA21$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
FRPG_FA21_plot <-ggplot(FRPG_FA21, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "gray", #undetermined sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+ 
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

FRPG_FA21_plot

#PED20
PEDR_FA20 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_PEDRFA20.csv")

#to order the species names for the legend
PEDR_FA20$Names_by_family <- factor(PEDR_FA20$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
PEDR_FA20_plot <-ggplot(PEDR_FA20, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "gray", #undetermined sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+ 
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

PEDR_FA20_plot

#PED21
PEDR_FA21 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_PEDRFA21.csv")

#to order the species names for the legend
PEDR_FA21$Names_by_family <- factor(PEDR_FA21$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
PEDR_FA21_plot <-ggplot(PEDR_FA21, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "gray", #undetermined sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+ 
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

PEDR_FA21_plot

#I505
I505_FA21 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_I505FA21.csv")

#to order the species names for the legend
I505_FA21$Names_by_family <- factor(I505_FA21$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
I505_FA21_plot <-ggplot(I505_FA21, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "gray", #undetermined sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+  
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

I505_FA21_plot

#WPCP
WPCP_FA21 <- read.csv("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/gear_comp/Putah_alluvial_plot_efish_WPCPFA21.csv")

#to order the species names for the legend
WPCP_FA21$Names_by_family <- factor(WPCP_FA21$Names_by_family, levels=c("Atherinopsidae: Menidia beryllina complex (silverside) [Introduced]",
                                                                        "Catostomidae: Catostomus occidentalis (Sacramento sucker) [Native]",
                                                                        "Centrarchidae: Lepomis cyanellus (green sunfish) [Introduced]",
                                                                        "Centrarchidae: L. gulosus (warmouth) [Introduced]",
                                                                        "Centrarchidae: L. macrochirus (bluegill) [Introduced]",
                                                                        "Centrarchidae: L. microlophus (redear sunfish) [Introduced]",
                                                                        "Centrarchidae: Lepomis sp. (undetermined sunfish) [Introduced]",
                                                                        "Centrarchidae: Micropterus dolomieu (smallmouth bass) [Introduced]",
                                                                        "Centrarchidae: M. salmoides (largemouth bass)  [Introduced]",
                                                                        "Centrarchidae: Pomoxis nigromaculatus (black crappie) [Introduced]",
                                                                        "Cottidae: Cottus asper (prickly sculpin) [Native]",
                                                                        "Cyprinidae: Carassius auratus (goldfish) [Introduced]",
                                                                        "Cyprinidae: Cyprinus carpio (common carp) [Introduced]",
                                                                        "Cyprinidae: Lavinia exilicauda (hitch) [Native]",
                                                                        "Cyprinidae: Notemigonus crysoleucas (golden shiner) [Introduced]",
                                                                        "Cyprinidae: Orthodon microlepidotus (Sacramento blackfish) [Native]",
                                                                        "Cyprinidae: Pimephales promelas (fathead minnow) [Introduced]",
                                                                        "Cyprinidae: Ptychocheilus grandis (Sacramento pikeminnow) [Native]",
                                                                        "Embiotocidae: Hysterocarpus traskii (tule perch) [Native]",
                                                                        "Gasterosteidae: Gasterosteus aculeatus (threespine stickleback) [Native]",
                                                                        "Ictaluridae: Ameiurus catus (white catfish) [Introduced]",
                                                                        "Ictaluridae: A. melas/nebulosus (bullhead) [Introduced]",
                                                                        "Percidae: Percina macrolepida (bigscale logperch) [Introduced]",
                                                                        "Poeciliidae: Gambusia affinis (western mosquitofish) [Introduced]",
                                                                        "Salmonidae: Oncorhynchus mykiss (rainbow trout) [Native]",
                                                                        "Salmonidae: O. tshawytscha (Chinook salmon) [Native]"))
WPCP_FA21_plot <-ggplot(WPCP_FA21, aes(x = Method, 
                      stratum = Names_by_family, 
                      alluvium = Names_by_family,
                      y = Relative_abundance,
                      fill = Names_by_family, 
                      label = Common_name_label)) + 
  scale_fill_manual(values = c("red",#silverside non-native
                               "chartreuse", #Sacramento sucker SKR native
                               "goldenrod", #green sunfish GSL non-native  
                               "pink", #warmouth non-native
                               "deeppink4", #bluegill BGL non-native
                               "darkorange4", #redear sunfish non-native
                               "gray", #undetermined sunfish non-native
                               "deeppink", #smallmouth bass SMB non-native
                               "darkorange2", #largemouth bass LMB  non-native
                               "black", #black crappie non-native
                               "blue", #prickly sculpin native
                               "gold", # goldfish non-native
                               "salmon", #common carp CRP non-native
                               "springgreen3", #hitch native
                               "yellow", #goldenshiner non-native
                               "darkgreen", #sacramento blackfish native
                               "firebrick1", #fatheadminnow non-native
                               "skyblue", #Sacramento pikeminnow PKM native
                               "cyan", #tule perch native
                               "aquamarine3", #threespine stickleback TSS native
                               "ivory", # white catfish non-native
                               "coral", #black/brown bullhead non-native
                               "bisque",#bigscale logperch BLP non-native"
                               "orchid", #western mosquitofish non-native
                               "darkcyan", #rainbow trout RBT native
                               "aliceblue" #Chinook salmon CHN  native
                               ))+  
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width=2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", 
            size = 4) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        text=element_text(size=16))+
  labs(y = "Relative abundance") 

WPCP_FA21_plot

#make a plot only for the legend
legend <- ggplot(WPCP_FA21, aes(x = Method, 
                                stratum = Names_by_family, 
                                alluvium = Names_by_family,
                                y = Relative_abundance,
                                fill = Names_by_family, 
                                label = Common_name_label)) + 
  scale_fill_manual(values = c(
    "red", "chartreuse", "goldenrod", "pink", "deeppink4", 
    "darkorange4", "gray", "deeppink", "darkorange2", "black", 
    "blue", "gold", "salmon", "springgreen3", "yellow", 
    "darkgreen", "firebrick1", "skyblue", "cyan", "aquamarine3", 
    "ivory", "coral", "bisque", "orchid", "darkcyan", "aliceblue"
  )) + 
  geom_flow(curve_type = "sigmoid") +
  geom_stratum(width = 2/3) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_text(stat = "stratum", size = 4) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.key = element_rect(color = "black", fill = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    text = element_text(size = 16)) +
  guides(
    fill = guide_legend(
      ncol = 1)) #display legend in one column

#extract the legend from last plot
grobs <- ggplotGrob(legend)$grobs
legend_plot <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

#build grid of plot without legends
gear_comp_efish <-plot_grid(MACU_FA20_plot, MACU_FA21_plot, 
                            OLDR_FA20_plot, FRPG_FA21_plot,
                            PEDR_FA20_plot, PEDR_FA21_plot,
                            I505_FA21_plot, WPCP_FA21_plot,
                            labels = c('A', 'B','C', 'D','E', 'F','G', 'H'), 
                            label_size = 18,
                            ncol = 1,
                            scale = 1)
gear_comp_efish

#add the legend to the plots
gear_comp_efish_w_legend <- plot_grid(gear_comp_efish, 
                                      plot_grid(legend_plot, NULL, ncol = 1, rel_heights = c(3.8, 10)), #adjust height of legend
                                      ncol = 2, 
                                      rel_widths = c(1, 0.8))

gear_comp_efish_w_legend

#save the plot
save_plot("~/Desktop/github/Putah-Creek-eDNA/gear_comp/gear_comp_efish_w_legend.jpg", 
          gear_comp_efish_w_legend, 
          dpi = 300,
          base_height = 25, #height in inches, equal to 30.48cm 
          base_width = 15, 
          bg = "white")

#Wilcoxon signed-rank test: if more species in eDNA samples vs electrofishing (one-tailed hypothesis)
#7 paired sites

eDNA <- c(14, 18, 10, 9, 13, 9, 11)
efish <- c(10, 7, 6, 8, 9, 5, 7)

wilcox.test(eDNA, efish, 
            alternative = "greater", 
            paired = TRUE)

#mean for each method
mean(eDNA)
mean(efish)

#mean difference (eDNA - electrofishing)
mean(eDNA - efish)

#mean for trap sample comparison (no statistics due to small sample size)
eDNA_trap <- c(9, 8, 10)
traps <- c(2, 6, 3)

#mean for each method
mean(eDNA_trap)
mean(traps)

#mean difference (eDNA - electrofishing)
mean(eDNA_trap - traps)
