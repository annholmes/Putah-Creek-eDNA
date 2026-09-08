library(tidyverse) # v2.0.0; for ggplot and dplyr
library(ape) # v5.8.1
library(DECIPHER) # v3.8.1
library(ggtree) # v4.2.0
library(reshape2)
library(grImport2) # to read local svgs
library(rsvg) # to normalize svg to Cairo format
library(cowplot)

# load data
# reference sequences
seqs <- readDNAStringSet("Supporting_Data_and_Resources/12S_MiFish_PutahCK_species_list.txt", format = "fasta")
# path to locally stored svg files
svg_dir <- "fish_pics/svg"
# import df showing local files with species silhouettes
fish_pics <- read.csv("Supporting_Data_and_Resources/fish_silhouettes.csv", stringsAsFactors = FALSE) %>%
  select(-image_source)

# set up dataframes
# extract binomials from reference sequence file
species_names <- sapply(strsplit(names(seqs), ";"), function(x) tail(x, 1))
# add the binomials to the fish images files
fish_pics$binomial <- trimws(sub("\\(.*", "", fish_pics$species))
# check that files align
setdiff(fish_pics$binomial, species_names)   # fish_pics with no FASTA match
setdiff(species_names, fish_pics$binomial)   # FASTA files with no fish_pics match
# finish adding species names
new_names <- fish_pics$common_name[match(species_names, fish_pics$binomial)]
names(seqs) <- new_names
sum(is.na(new_names))

# align reference sequences
aligned <- AlignSeqs(seqs, anchor = NA)
# inspect reference sequences
BrowseSeqs(aligned)

# check genetic distances by generating a distance matrix for MiFish barcodes
aln_dna <- as.DNAbin(aligned)
dist_matrix <- dist.dna(aln_dna, model = "raw", pairwise.deletion = TRUE)
dist_mat <- as.matrix(dist_matrix)
diag(dist_mat) <- NA
low_dist_pairs <- which(dist_mat < 0.03, arr.ind = TRUE) # set >97% similarity as warranting further consideration
low_dist_df <- data.frame(
  species1 = rownames(dist_mat)[low_dist_pairs[, 1]],
  species2 = colnames(dist_mat)[low_dist_pairs[, 2]],
  distance = dist_mat[low_dist_pairs])
low_dist_df_derep <- low_dist_df %>%
  mutate(pair_key = paste(pmin(species1, species2), pmax(species1, species2))) %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  select(-pair_key) %>%
  mutate(distance = round(distance, 4),
         pct_similarity = round((1 - distance) * 100, 2))

# view and save the output
low_dist_df_derep
write.csv(low_dist_df_derep, "Supporting_Data_and_Resources/Table_species_pairs_dist_mat.csv", row.names = FALSE)

# create a phylogenetic tree from the distance matrix
tree <- as.phylo(hclust(dist_matrix, method = "average"))
tree <- root(tree, outgroup = "Pacific lamprey", resolve.root = TRUE)
tree <- ladderize(tree, right = FALSE)
ggtree(tree)

# tree structure with assigned node numbers
ggtree(tree) + geom_text(aes(label=node), hjust=-.3)

# add species names from imported tree file
ggtree(tree) + geom_text2(aes(subset=! isTip, label=node), hjust=-.3) + geom_tiplab()

# add fish silhouettes
# check that common names match tree tip labels
setdiff(tree$tip.label, fish_pics$common_name)

# create helper function to load local svg files into an object that ggtree and rphylopic can use
load_local_phylopic <- function(file) {
  path <- file.path(svg_dir, file)
  tryCatch({
    tmp <- tempfile(fileext = ".svg")
    rsvg::rsvg_svg(path, tmp)   #normalize to Cairo SVG
    grImport2::readPicture(tmp)
  }, error = function(e) {
    warning(paste("Could not load silhouette for:", path))
    NA
  })
}

# add each silhouette
fish_pics$img <- lapply(fish_pics$svg, load_local_phylopic)

# check if any fish pics failed to load, should return 0
sum(is.na(fish_pics$img))

# build tree with silhouettes
fish_pics_join <- fish_pics %>% 
  select(common_name, everything()) %>%
  mutate(svg_path = file.path(svg_dir, svg))  #build file paths

# create a column to facilitate minor adjustments in fish image sizes
fish_pics_join$img_size <- 0.05 # default size

# fine tune the size of each fish species; easier to track if unchanged species are also listed
# roach unchanged
# hitch unchanged
fish_pics_join$img_size[fish_pics_join$common_name == "Sacramento pikeminnow"] <- 0.07
# Sac blackfish unchanged
# golden shiner unchanged
fish_pics_join$img_size[fish_pics_join$common_name == "goldfish"] <- 0.045
# common carp unchanged
# fathead minnow unchanged
# red shiner unchanged
fish_pics_join$img_size[fish_pics_join$common_name == "Sacramento sucker"] <- 0.07
fish_pics_join$img_size[fish_pics_join$common_name == "black bullhead/brown bullhead"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "white catfish"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "channel catfish"] <- 0.075
fish_pics_join$img_size[fish_pics_join$common_name == "inland silverside"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "kokanee"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "rainbow trout"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "Chinook salmon"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "green sunfish"] <- 0.045
fish_pics_join$img_size[fish_pics_join$common_name == "bluegill"] <- 0.04
fish_pics_join$img_size[fish_pics_join$common_name == "warmouth"] <- 0.045
fish_pics_join$img_size[fish_pics_join$common_name == "pumpkinseed"] <- 0.045
fish_pics_join$img_size[fish_pics_join$common_name == "redear sunfish"] <- 0.045
# smb unchanged
fish_pics_join$img_size[fish_pics_join$common_name == "largemouth bass"] <- 0.07
# spotted bass unchanged
fish_pics_join$img_size[fish_pics_join$common_name == "Sacramento perch"] <- 0.045
fish_pics_join$img_size[fish_pics_join$common_name == "black crappie"] <- 0.045
fish_pics_join$img_size[fish_pics_join$common_name == "bigscale logperch"] <- 0.06
# stickleback unchanged
# mosquitofish unchaged
fish_pics_join$img_size[fish_pics_join$common_name == "stripped bass"] <- 0.075
# tule perch unchanged
# prickly sculpin unchanged
# threadfin shad unchanged
fish_pics_join$img_size[fish_pics_join$common_name == "yellowfin goby"] <- 0.06
fish_pics_join$img_size[fish_pics_join$common_name == "Pacific lamprey"] <- 0.06

# create the tree
phylo_tree <- 
  ggtree(tree) %<+% 
  fish_pics_join + 
  geom_tiplab(align = TRUE, linetype = "22", linewidth = 0.3, offset = 0.09) +
  geom_tiplab(aes(image = svg_path, color = status, size = img_size), geom = "image",
              offset = 0.04, show.legend = FALSE) +
  scale_size_identity() +
  geom_point(data = data.frame(status = c("native", "non-native")),
             aes(x = -Inf, y = -Inf, color = status),
             alpha = 0, size = 0, show.legend = TRUE) +
  scale_color_manual("Status", values = c("native" = "blue", "non-native" = "orange")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 16))) +
  theme(legend.position = c(0.02, 1.05),
        legend.justification = c(0, 1),
        legend.direction = "horizontal",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.background = element_rect(fill = scales::alpha("white", 0.6)),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 0, 0)) +
  xlim(0, 0.5)

phylo_tree

# save the tree to check formatting
save_plot("fig2_test.jpg",
          phylo_tree,
          base_width = 6,
          base_height = 8,
          dpi = 300)

#create the other half of the figure showing sites where species were detected by eDNA
# create df for +/- at each site
sites <- c("MACD", "MACU", "OLDR", "APOP", "PEDR", "RRAN", "SW2D", "SW1D", "I505", "WPCP", "SCCA")
species <- c("hitch", "California roach", "Sacramento pikeminnow", "Sacramento blackfish", "golden shiner","goldfish", 
             "common carp", "fathead minnow", "red shiner", "Sacramento sucker",
             "black/brown bullhead", "white catfish", "channel catfish", "inland silverside", "kokanee",
             "rainbow trout", "Chinook salmon", "green sunfish", "bluegill", "warmouth",
             "pumpkinseed", "redear sunfish", "smallmouth bass", "largemouth bass", "spotted bass",
             "Sacramento perch", "black crappie", "bigscale logperch", "three-spined stickleback", "western mosquitofish",
             "striped bass", "tule perch", "prickly sculpin", "threadfin shad", "yellowfin goby", "Pacific lamprey")

# start all species with 0 (not detected)
species_by_site <- data.frame(Species_order = species)
for (s in sites) {
  species_by_site[[s]] <- 0
}

# for species detected at all sites, set entire row to 1
all_sites_detected <- c("prickly sculpin", "largemouth bass", "Sacramento sucker", "Sacramento pikeminnow")
species_by_site[species_by_site$Species_order %in% all_sites_detected, sites] <- 1

# other species by sites
species_by_site[species_by_site$Species_order == "Sacramento blackfish", c("MACD", "OLDR")] <- 1
species_by_site[species_by_site$Species_order %in% c("golden shiner", "goldfish", "white catfish", "inland silverside", "black crappie"), c("MACD", "MACU", "OLDR")] <- 1
species_by_site[species_by_site$Species_order == "goldfish", "WPCP"] <- 1
species_by_site[species_by_site$Species_order == "inland silverside", "APOP"] <- 1
species_by_site[species_by_site$Species_order == "common carp", setdiff(sites, c("RRAN", "SW1D", "WPCP"))] <- 1
species_by_site[species_by_site$Species_order == "fathead minnow", "OLDR"] <- 1
species_by_site[species_by_site$Species_order == "black/brown bullhead", c("MACU", "OLDR")] <- 1
species_by_site[species_by_site$Species_order == "rainbow trout", setdiff(sites, c("MACD", "MACU", "OLDR", "APOP", "PEDR")) ] <- 1
species_by_site[species_by_site$Species_order == "Chinook salmon", setdiff(sites, c("OLDR", "APOP")) ] <- 1
species_by_site[species_by_site$Species_order == "green sunfish", setdiff(sites, c("PEDR", "SW2D", "I505", "SCCA")) ] <- 1
species_by_site[species_by_site$Species_order == "bluegill", setdiff(sites, c("RRAN", "SCCA")) ] <- 1
species_by_site[species_by_site$Species_order == "warmouth", "MACU"] <- 1
species_by_site[species_by_site$Species_order == "redear sunfish", setdiff(sites, c("SW1D", "I505", "WPCP", "SCCA")) ] <- 1
species_by_site[species_by_site$Species_order == "smallmouth bass", setdiff(sites, c("SW2D", "SW1D", "WPCP")) ] <- 1
species_by_site[species_by_site$Species_order == "bigscale logperch", c("MACD", "MACU", "RRAN", "SW1D", "SW2D")] <- 1
species_by_site[species_by_site$Species_order == "three-spined stickleback", setdiff(sites, c("PEDR", "OLDR"))] <- 1
species_by_site[species_by_site$Species_order %in% c("tule perch", "western mosquitofish"), setdiff(sites, "SCCA")] <- 1
species_by_site[species_by_site$Species_order %in% c("kokanee", "threadfin shad"), "SCCA"] <- 1

sum(rowSums(species_by_site[sites]) > 0) # checks that 26 taxa were detected at at least one site

# check if tree and species by site list are in the same order
tip_order <- get_taxa_name(phylo_tree)
current_levels <- levels(species_by_site$Species_order)
identical(tip_order, current_levels)
# change to factor and recheck
species_by_site$Species_order <- factor(species_by_site$Species_order, levels = tip_order) # make factor
identical(levels(species_by_site$Species_order), tip_order)

# add heat map
# eDNA detections by site (all years and seasons pooled), colored by native/non-native status
species_by_site$Species_order <- as.factor(species_by_site$Species_order)

# add native/non-native status from fish_pics_join, matched by species name
status_lookup <- fish_pics_join %>% select(common_name, status)

tiles_long <- melt(species_by_site, id.vars = "Species_order") %>%
  left_join(status_lookup, by = c("Species_order" = "common_name")) %>%
  mutate(Species_order = factor(Species_order, levels = rev(tip_order)))

# detected native = blue, detected non-native = orange, not detected = grey
tiles_long <- tiles_long %>%
  mutate(fill_color = case_when(
    value == 1 & status == "native"     ~ "blue",
    value == 1 & status == "non-native" ~ "orange",
    TRUE                                ~ "grey95"
  ))

# format the tiles indicating presence of species where detected by eDNA
tiles <-
  ggplot(tiles_long, aes(Species_order, variable, fill = fill_color)) +
  scale_fill_identity() +
  geom_tile(colour = "white") +
  coord_equal(expand = 0) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, color = "black"),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  coord_flip()

which(is.na(species_by_site$Species_order))
tiles

# combine the tree and the tiles
tree_with_occurrence <- phylo_tree + tiles

# view and save final figure
tree_with_occurrence
save_plot("Manuscript_Outputs/Figures/Fig2_PutahCk_eDNA_tree_with_tiles.jpg", 
          tree_with_occurrence, 
          base_width = 7.5,
          base_height = 8,
          dpi = 300)
