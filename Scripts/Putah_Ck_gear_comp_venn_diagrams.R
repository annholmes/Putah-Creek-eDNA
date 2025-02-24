#Putah Creek eDNA study (Holmes et al. 2025)

#Venn diagrams comparing species detected by eDNA vs. conventional methods
#there are 3 trap sites and 8 electrofishing sites
#note not all eDNA sites were sampled by a conventional method

library(ggplot2)
library(cowplot)

#fyke trap and rotary screw trap
#fyke 
RST_fyke20_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightpink", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Fyke trap"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "Sacramento sucker\nprickly sculpin\ntule perch\ngreen sunfish\nsmallmouth bass\nrainbow trout\nbigscale logperch"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by trap only
  geom_text(aes(x = 1.43, y = 0.58, label = ""), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "bluegill\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

RST_fyke20_venn

#RST_Apr20
RST_Apr21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "wheat", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Rotary screw trap"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "threespine stickleback\ntule perch\nrainbow trout"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by trap only
  geom_text(aes(x = 1.43, y = 0.58, label = "Pacific lamprey"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "Sacramento sucker\nprickly sculpin\nbluegill\nChinook salmon\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

RST_Apr21_venn

#RST_May 
RST_May21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "wheat", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Rotary screw trap"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "Sacramento sucker\nprickly sculpin\ncommon carp\nthreespine stickleback\ntule perch\nlargemouth bass\nbigscale logperch"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by trap only
  geom_text(aes(x = 1.43, y = 0.58, label = ""), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "rainbow trout\nChinook salmon\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

RST_May21_venn

#build grid of venn diagrams for trap sites
gear_comp_traps_venn <-plot_grid(RST_fyke20_venn, RST_Apr21_venn, RST_May21_venn,
                                 labels = c('A', 'B','C'), 
                                 label_size = 18,
                                 ncol = 2,
                                 scale = 1.01)

gear_comp_traps_venn

#save the plot
save_plot("~/Desktop/github/Putah-Creek-eDNA/gear_comp/gear_comp_traps_venn.jpg", 
          gear_comp_traps_venn, 
          dpi = 300,
          base_height = 12.5, #width in inches, equal to 30.48cm, half of efishing plot
          base_width = 20, # same as efishing plot
          bg = "white")

#electrofishing sites and multimethod survey (WFC120 class)
#MACU20
MACU_FA20_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "white catfish\ncommon carp\nwarmouth\nChinook salmon\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = "bigscale logperch"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "black bulllhead\nprickly sculpin\nwestern mosquitofish\ngreen sunfish\nbluegill\nredear sunfish\nsilverside\nlargemouth bass\ngolden shiner"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#MACU_FA21
MACU_FA21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "white catfish\nblack/brown bullhead\ngoldfish\nSacramento sucker\ncommon carp\nwestern mosquitofish\ntule perch\nwarmouth\nbigscale logperch\nblack crappie\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = ""), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "green sunfish\nbluegill\nredear sunfish\nsilverside\nsmallmouth bass\nlargemouth bass\ngolden shiner"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#OLDR_FA20
OLDR_FA20_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "Sacramento sucker\ncommon carp\nwestern mosquitofish\nredear sunfish\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = "bigscale logperch"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "white catfish\nprickly sculpin\ngreen sunfish\nbluegill\nlargemouth bass"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#FRPG_FA21
FRPG_FA21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  #eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "palegreen", alpha = 0.5) +  #Multiple methods circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Multiple methods"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "(black/brown bullhead)\ncommon carp\nredear sunfish\nwestern mosquitofish\n(Chinook salmon)"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by conventional methods only
  geom_text(aes(x = 1.43, y = 0.58, label = "undetermined sunfish\nsilverside\nsmallmouth bass\ngolden shiner\nfathead minnow"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "Sacramento sucker\nprickly sculpin\ntule perch\nbluegill\nredear sunfish\nlargemouth bass\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#PEDR_FA20
#plot circles manually
PEDR_FA20_venn <- ggplot() +
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "bluegill\nredear sunfish\nsmallmouth bass"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = "hitch\nbigscale logperch"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "Sacramento sucker\nprickly sculpin\nwestern mosquitofish\ntule perch\nlargemouth bass\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#PEDR_FA21
PEDR_FA21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  # eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  # Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "common carp\nredear sunfish\nsmallmouth bass\nChinook salmon"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = "hitch"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "Sacramento sucker\nprickly sculpin\nwestern mosquitofish\ntule perch\nbluegill\nsilverside\nlargemouth bass\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#I505_FA21
I505_FA21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  #eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  #Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "nwestern mosquitofish\nthreespine stickleback\nsmallmouth bass\nbluegill\nlargemouth bass\nChinook salmon"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = ""), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "Sacramento sucker\nprickly sculpin\ntule perch\nrainbow trout\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)

#WPCP_FA21  
WPCP_FA21_venn <- ggplot() +
  #plot circles manually
  geom_point(aes(x = 0.75, y = 0.6), size = 175, shape = 21, fill = "lightblue", alpha = 0.5) +  #eDNA circle
  geom_point(aes(x = 1.25, y = 0.6), size = 175, shape = 21, fill = "lightyellow", alpha = 0.5) +  #Electrofishing circle
  #add labels for gear type
  geom_text(aes(x = 0.75, y = 1.1, label = "eDNA"), size = 6, fontface = "bold") +  #left circle
  geom_text(aes(x = 1.25, y = 1.1, label = "Electrofishing"), size = 6, fontface = "bold") +  #right circle
  #add species names for species detected by eDNA only
  geom_text(aes(x = 0.60, y = 0.58, label = "goldfish\nwestern mosquitofish\nlargemouth bass\nChinook salmon"), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by electrofishing only
  geom_text(aes(x = 1.43, y = 0.58, label = ""), 
            size = 6, hjust = 0.5) +
  #add species names for species detected by both methods
  geom_text(aes(x = 1, y = 0.6, label = "Sacramento sucker\nprickly sculpin\nthreespine stickleback\ntule perch\nbluegill\nrainbow trout\nSacramento pikeminnow"), 
            size = 6, hjust = 0.5) +
  theme_minimal() +
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), 
        legend.position = "none", 
        panel.grid = element_blank()) +
  xlim(0.25, 1.75) +
  ylim(0, 1.25)


WPCP_FA21_venn

#build grid of venn diagrams for electrofishing sites
gear_comp_efish_venn <-plot_grid(MACU_FA20_venn, MACU_FA21_venn, 
                            OLDR_FA20_venn, FRPG_FA21_venn,
                            PEDR_FA20_venn, PEDR_FA21_venn,
                            I505_FA21_venn, WPCP_FA21_venn,
                            labels = c('A', 'B','C', 'D','E', 'F','G', 'H'), 
                            label_size = 18,
                            ncol = 2,
                            scale = 1.01)

gear_comp_efish_venn

#save the plot
save_plot("~/Desktop/github/Putah-Creek-eDNA/gear_comp/gear_comp_efish_venn.jpg", 
          gear_comp_efish_venn, 
          dpi = 300,
          base_height = 25, #height in inches, equal to 30.48cm 
          base_width = 20, 
          bg = "white")
