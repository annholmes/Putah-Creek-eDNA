library(tidyverse)
library(ggtree)
library(ggplot2)
library(rphylopic)
library(phytools)
library(magick)
library(reshape2)
library(cowplot)

tree <- read.tree("/Users/aholmes/Desktop/github/Putah-Creek-eDNA/PC_Tree_UPGMA.newick")

#show the basic tree structure
ggtree(tree)

#tree structure with assigned node numbers
ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

#add species names from imported tree file
ggtree(tree) + geom_text2(aes(subset=! isTip, label=node), hjust=-.3) + geom_tiplab()

#import Phylopic images to add to tree
#see tutorial https://rphylopic.palaeoverse.org/articles/b-advanced-ggplot.html

#create df with species
#UUID is the code used to retrieve the phylopic images
fish_pics <- data.frame(species = tree$tip.label, uuid = NA)

#get the phylopic UUIDs for the species names
fish_pics$uuid <- sapply(tree$tip.label,
                         function(x) {
                           tryCatch(get_uuid(x), error = function(e) NA)
                         })

#only pulling phylopic images for bluegill and pumpkinseed 
#add the missing pics manually using the UUID
#some aren't pulling a pic because there isn't one, so add an image as a stand in

#Cyprinidae
fish_pics$uuid[1] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio
fish_pics$uuid[2] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio
fish_pics$uuid[3] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio
fish_pics$uuid[4] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio 
fish_pics$uuid[5] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio
fish_pics$uuid[6] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #C. carpio
fish_pics$uuid[7] <- "de187ba5-0498-4d2e-a688-70245c354c5f" #C. auratus 
fish_pics$uuid[8] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio
fish_pics$uuid[9] <- "25a76141-cca2-4b56-bd9f-9f5fdcd8efa1" #stand in: C. carpio

fish_pics$uuid[10] <- "c657b67e-db70-4f8d-a66f-e6b287b148a6" #stand in: C. commersonii -ok?

fish_pics$uuid[11] <- "1df5b970-3302-4e5b-8df1-a3d9fb5744d3" #O. mykiss
fish_pics$uuid[12] <- "3c098bb8-4158-4777-b567-80e48049409c" #O. nerka
fish_pics$uuid[13] <- "9150be88-6910-4374-aa54-a7f8f3d79fb6" #O. tshawytscha

fish_pics$uuid[14] <- "4e31c12d-1c23-4a64-bcc6-0730ab6d0617" #stand in: Odontesthes sp.

fish_pics$uuid[15] <- "6a920565-e1fc-43f0-952a-83b0243ee97a" #stand in: Ameiurus sp.
fish_pics$uuid[16] <- "7a6448e5-09c4-40c8-8378-599d7f974bfe" #A. melas
fish_pics$uuid[17] <- "4317e447-87bb-42b9-94bc-8c910afdbe4b" #I.punctatus

#Centrarchidae
fish_pics$uuid[18] <- "172f60da-2129-4d3d-b680-686cff376d24" #M. salmoides
fish_pics$uuid[19] <- "ab0e38bc-f514-46fb-8255-17adbc6c8d37" #M. dolomieu
fish_pics$uuid[20] <- "172f60da-2129-4d3d-b680-686cff376d24" #stand in: M. salmoides
fish_pics$uuid[21] <- "ccee99b3-052a-4e2d-9db2-80a3a75484f0" #P. nigromaculatus
fish_pics$uuid[22] <- "ae674be3-6a35-4c96-bdb5-31f2363a8bab" #A. interuptus stand in: A. rupestris
#23 bluegill
fish_pics$uuid[24] <- "986a42ec-51f6-4c0e-9e62-2fcbad2657c3" #stand in: L. macrochirus
fish_pics$uuid[25] <- "986a42ec-51f6-4c0e-9e62-2fcbad2657c3" #stand in: L. macrochirus
fish_pics$uuid[26] <- "986a42ec-51f6-4c0e-9e62-2fcbad2657c3" #stand in: L. macrochirus
#27 pumpkinseed
fish_pics$uuid[28] <- "addc129c-d78f-4250-aba7-eb229f0b57e6" #G. aculeatus
fish_pics$uuid[29] <- "966455e3-87b9-4bb3-86ef-7856e8cc3456" #stand in for bs lp: A. flavimanus
fish_pics$uuid[30] <- "592e3594-8bdf-4bd1-b581-56cbdc3e9873" #G. affinis

#other families
fish_pics$uuid[31] <- "0b9cdf1f-ccbc-4922-8cf6-60f90d07107e" #stand in: Hyperprosopon argenteum - close?
fish_pics$uuid[32] <- "bae30310-1dbd-4057-b673-abc8ebd9f7a9" #stand in: M. chrysops - not quite right
fish_pics$uuid[33] <- "a3f9f381-998b-4459-b505-2946e6cb928a" #C. asper stand in: C. aleuticus
fish_pics$uuid[34] <- "999337e0-ecf8-458a-8da1-f2489cf780b6" #stand in: D. cepedianum

#flip this branch so lamprey (outgroup) are on the bottom
fish_pics$uuid[35] <- "63953094-ac64-42c3-920e-53ee87ab188f" #stand in: P. marinus
fish_pics$uuid[36] <- "966455e3-87b9-4bb3-86ef-7856e8cc3456" #A. flavimanus

#all common names with 2 words had '' around them, so add a label column
fish_pics$label2[1] <- "hitch"
fish_pics$label2[2] <- "California roach"
fish_pics$label2[3] <- "Sacramento pikeminnow"
fish_pics$label2[4] <- "Sacramento blackfish"
fish_pics$label2[5] <- "golden shiner" 
fish_pics$label2[6] <- "common carp"
fish_pics$label2[7] <- "goldfish"
fish_pics$label2[8] <- "fathead minnow"
fish_pics$label2[9] <- "red shiner"
fish_pics$label2[10] <- "Sacramento sucker" 
fish_pics$label2[11] <- "rainbow trout" 
fish_pics$label2[12] <- "kokanee salmon"
fish_pics$label2[13] <- "Chinook salmon" 
fish_pics$label2[14] <- "inland silverside" 
fish_pics$label2[15] <- "white catfish"
fish_pics$label2[16] <- "black/brown bullhead" 
fish_pics$label2[17] <- "channel catfish"
fish_pics$label2[18] <- "largemouth bass"
fish_pics$label2[19] <- "smallmouth bass" 
fish_pics$label2[20] <- "spotted bass"
fish_pics$label2[21] <- "black crappie"
fish_pics$label2[22] <- "Sacramento perch" 
fish_pics$label2[23] <- "bluegill"
fish_pics$label2[24] <- "green sunfish"
fish_pics$label2[25] <- "warmouth" 
fish_pics$label2[26] <- "redear sunfish"
fish_pics$label2[27] <- "pumkinseed"
fish_pics$label2[28] <- "three-spined stickleback"
fish_pics$label2[29] <- "bigscale logperch" 
fish_pics$label2[30] <- "western mosquitofish"
fish_pics$label2[31] <- "tule perch" 
fish_pics$label2[32] <- "striped bass" 
fish_pics$label2[33] <- "prickly sculpin"
fish_pics$label2[34] <- "threadfin shad" 
#flip this branch so lamprey (outgroup) are on the bottom
fish_pics$label2[35] <- "Pacific lamprey"
fish_pics$label2[36] <- "yellowfin goby"

#retrieve images
fish_pics$svg <- sapply(fish_pics$uuid,
                        function(x) {
                          tryCatch(get_phylopic(x), error = function(e) NA)
                        })

#add native and non-native status
fish_pics <- fish_pics %>%
  mutate(status=c("native", "native", "native", "native", "non-native",
                  "non-native", "non-native", "non-native", "non-native", "native",
                  "native", "non-native", "native", "non-native", "non-native",
                  "non-native", "non-native", "non-native", "non-native", "non-native",
                  "non-native", "native", "non-native", "non-native", "non-native",
                  "non-native", "non-native", "native", "non-native", "non-native",
                  "native", "non-native", "native", "non-native", "native", "non-native"))

#create phylogenetic tree with fish images
#use geom_tiplab to adjust placement of phylopics
#native species are blue
#non-native species are orange
phylo_tree <- 
  ggtree(tree) %<+% 
  fish_pics + 
  geom_tiplab(aes(label=label2), align=TRUE, linesize=.2, offset = 0.1)+
  geom_tiplab(aes(image=uuid, color=status), geom = "phylopic", 
              size=c(0.05, 0.05, 0.05, 0.045, 0.05, #1 tule perch #2 largemouth bass #3 spotted bass #4 pumpkinseed #5 rainbow trout
                     0.05, 0.05, 0.05, 0.05, 0.05, #1 CA roach #2 hitch #3 Sac bf #4 pikeminnow #5 golden shiner
                     0.05, 0.05, 0.05, 0.07, 0.065, #1 common carp #2 fathead minnow #3 red shiner #4 kokanee #5channel cat
                     0.05, 0.05, 0.05, 0.06, 0.06, #1 inland ss #2 western mf #Pac lamprey #4 black bh #5 wh cf
                     0.07, 0.05, 0.05, 0.045, 0.045, #1 chinook #2 yf goby #3 bigscale logperch  #4 bluegill #5 green sf 
                     0.045, 0.045, 0.05, 0.05, 0.05, #1 warmouth #2 redear  #3 threadfin #4 sculpin #5 sm bass
                     0.05, 0.06, 0.05, 0.06, 0.05, 0.05), #1 stickleback #2 Sac perch #3 str bass #4 Sac sucker #4 #5 blk crappie #6 goldfish
              offset = 0.04)+
  scale_color_manual("status", values = c("blue", "orange"))+
  theme(legend.position = "top",
        legend.title=element_text(size=10),
        legend.text=element_text(size=10))+
  xlim(0,0.5)

phylo_tree

#add heatmap of presence absence of each species at each site
fish_by_site <- read.csv("all_species_occ_by_site_fam.csv")
fish_by_site$Species_order <- as.factor(fish_by_site$Species_order)

#bright tiles are present
#pale tiles are absent
tiles<-
  ggplot(melt(fish_by_site), aes(Species_order, variable, fill = Status, alpha = value)) + 
  scale_fill_manual(values = c("blue", "orange"))+
  geom_tile(colour = "white") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank())+
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  coord_flip()

tiles

tree_with_occurence <- phylo_tree + tiles

save_plot("fig2_tree_with_occurence.jpg", 
          tree_with_occurence, 
          base_width = 8,
          base_height = 8,
          dpi=300)