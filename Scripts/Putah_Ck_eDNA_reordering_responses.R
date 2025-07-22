##Rank Abundance Curves for Putah Creek data

#some code is modified from Aviolo et al. 2019 Fig 6 see https://github.com/mavolio/RACs_paper
#see also Jacinto et al. 2023 for color matches to that study, except rainbow trout is purple instead of blue)

library(here)
library(tidyverse)
library(vegan)
library(gridExtra)
library(grid)
library(gtable)
library(devtools)
library(codyn)
library(scales)
library(stringr)
library(viridis)

#read in the ASV table glom to species with sample IDs as row names
asv_data <- read.csv(here("Supporting_Data_and_Resources", "Putah_Creek_ASV_table_glom.csv"), 
                     row.names = 1,
                     check.names = FALSE)

#transpose so species are rows and samples are columns
asv_t <- t(asv_data)

#convert to relative abundance (each column = one sample)
rel_abund <- sweep(asv_t, 2, colSums(asv_t), "/")

#tidy to long format
RAC_data <- rel_abund %>%
  as.data.frame() %>%
  rownames_to_column(var = "species") %>%  # move row names into a "species" column
  pivot_longer(
    cols = -species,
    names_to = "sample_id",
    values_to = "rel_abund"
  ) %>%
  filter(!is.na(rel_abund) & rel_abund > 0)

#join metadata to get sampling_event_no
RAC_data <- RAC_data %>%
  left_join(metadata, by = "sample_id")

#aggregate replicates by sampling_event_no and species, averaging rel_abund
RAC_data_agg <- RAC_data %>%
  group_by(sampling_event_no, species) %>%
  summarise(rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
  group_by(sampling_event_no) %>%
  arrange(desc(rel_abund)) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

RAC_data_agg <- RAC_data_agg %>%
  left_join(
    metadata %>%
      select(sampling_event_no, river_km) %>%
      distinct(),  # ensures one row per sampling_event_no
    by = "sampling_event_no"
  )

#join species metadata for status (native/non-native) and colors
RAC_data_agg <- RAC_data_agg %>%
  left_join(species_metadata, by = "species")

#convert river_km to numeric with decimal formatting
RAC_data_agg <- RAC_data_agg %>%
  mutate(river_km = sprintf("%.1f", as.numeric(river_km)))  #for label with "51.0" not "51"

#RAC plots 
#plots showing native (blue) and non-native (gold) species
RAC_status <- ggplot(RAC_data_agg, aes(x = Rank, y = rel_abund, group = sampling_event_no)) +
  geom_line(color = "gray80", size = 0.8) +
  geom_point(aes(color = status), size = 3, shape = 17) +
  scale_y_log10(
    name = "Relative Abundance (log scale)",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  scale_x_log10(
    name = "Species Rank (log scale)",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(~river_km) +  
  scale_color_manual(name = "Status", #capitalizes "status"
                       values = c("native" = "blue4", "non-native" = "gold1")) +
  theme_bw() +
  theme(
    legend.position = c(0.88, 0.15),
    legend.background = element_blank(),
    legend.title = element_text(size = 14),  # Title size
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3)
  )

#save plot
save_plot("PutahCK_eDNA_RAC_status.png",
          RAC_status,
          base_width = 8,
          base_height = 5,
          dpi = 300,
          bg = "white")

#plots following color scheme for top 4 species from electrofishing data analyzed in Jacinto et al. 2022
RAC_top4_efish <- ggplot(RAC_data_agg, aes(x = Rank, y = rel_abund, group = sampling_event_no)) +
  geom_line(color = "gray80", size = 0.8) +
  geom_point(aes(color = species), size = 3, shape = 17) +
  scale_y_log10(
    name = "Relative Abundance (log scale)",
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_x_log10(
    name = "Species Rank (log scale)",
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  facet_wrap(~river_km) +
  scale_color_manual(
    name = "Species",
    values = c(
      "Cottus asper" = "springgreen4",
      "Cyprinus carpio" = "red2",
      "Micropterus salmoides" = "darkorange",
      "Oncorhynchus mykiss" = "orchid4",
      "other" = "#7D7D7D"
    ),
    breaks = c("Cottus asper", "Cyprinus carpio", "Micropterus salmoides", "Oncorhynchus mykiss")
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.885, 0.10),  # adjust as needed
    legend.background = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

#save plot
save_plot("PutahCK_eDNA_RAC_top4_efish_highlighted.png",
          RAC_top4_efish,
          base_width = 8,
          base_height = 5,
          dpi = 300,
          bg = "white")


#plots highlighting the top 4 species from the eDNA data
RAC_top4_eDNA <- ggplot(RAC_data_agg, aes(x = Rank, y = rel_abund, group = sampling_event_no)) +
  geom_line(color = "gray80", size = 0.8) +
  geom_point(aes(color = species), size = 3, shape = 17) +
  scale_y_log10(
    name = "Relative Abundance (log scale)",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_x_log10(
    name = "Species Rank (log scale)",
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  facet_wrap(~river_km) +
  scale_color_manual(
    name = "Species",
    values = c(
      "Catostomus occidentalis"  = "#009999",
      "Ptychocheilus grandis"     = "#006400",
      "Micropterus salmoides"    = "darkorange",
      "Lepomis macrochirus"      = "#b66dff",
      "other"                    = "#7D7D7D"
    ),
    breaks = c(
      "Catostomus occidentalis",
      "Ptychocheilus grandis",
      "Micropterus salmoides",
      "Lepomis macrochirus"
    )
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.885, 0.10),
    legend.background = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

#save plot
save_plot("PutahCK_eDNA_RAC_top4_eDNA_highlighted.png",
          RAC_top4_eDNA,
          base_width = 8,
          base_height = 5,
          dpi = 300,
          bg = "white")


#exploratory MRS values
#there are 9 possible comparisons across sequential years or seasons; put these into a data frame
comparison_table <- tribble(
  ~event_1, ~event_2, ~comparison_type,         ~category,
  3,        7,        "Fall 2019 → Spring 2020", "Season to Season",
  2,        6,        "Fall 2019 → Spring 2020", "Season to Season",
  1,        5,        "Fall 2019 → Spring 2020", "Season to Season",
  5,        13,       "Spring 2020 → Spring 2021", "Year to Year",
  8,        16,       "Fall 2020 → Fall 2021",   "Year to Year",
  9,        17,       "Fall 2020 → Fall 2021",   "Year to Year",
  10,       18,       "Fall 2020 → Fall 2021",   "Year to Year",
  8,        14,       "Fall 2020 → Spring 2021", "Season to Season",
  14,       16,       "Spring 2021 → Fall 2021", "Season to Season"
)

#function to calculate MRS for each comparison
calculate_mrs <- function(ev1, ev2, df) {
  sub_df <- df %>%
    filter(sampling_event_no %in% c(ev1, ev2)) %>%
    group_by(sampling_event_no, species) %>%
    summarise(rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
    mutate(sampling_event_no = as.numeric(sampling_event_no))
  
  if (length(unique(sub_df$sampling_event_no)) == 2) {
    rank_shift(
      df = sub_df,
      time.var = "sampling_event_no",
      species.var = "species",
      abundance.var = "rel_abund",
      replicate.var = NA_character_
    )
  } else {
    NULL
  }
}

#get river_km for each sampling event
event_rkm <- RAC_data %>%
  distinct(sampling_event_no, river_km)

#join river_km info to comparison table and filter for within-site comparisons
comparison_table_rkm <- comparison_table %>%
  left_join(event_rkm, by = c("event_1" = "sampling_event_no")) %>%
  rename(rkm_1 = river_km) %>%
  left_join(event_rkm, by = c("event_2" = "sampling_event_no")) %>%
  rename(rkm_2 = river_km) %>%
  filter(rkm_1 == rkm_2) %>%
  mutate(river_km = rkm_1) %>%
  select(event_1, event_2, comparison_type, category, river_km)

#MRS calculation
MRS_results <- comparison_table_rkm %>%
  mutate(result = map2(event_1, event_2, ~calculate_mrs(.x, .y, RAC_data))) %>%
  filter(!map_lgl(result, is.null)) %>%
  unnest(result) %>%
  select(event_1, event_2, comparison_type, category, MRS, river_km)

#create a new column so the plot displays a mean value instead of just one of the season-to-season values for 11.1
MRS_results <- MRS_results %>%
  group_by(river_km, category) %>%
  mutate(MRS_avg_for_plot = if (n() > 1) mean(MRS) else MRS) %>%
  ungroup()

#update the river_km labels so all have 1 decimal place
MRS_results <- MRS_results %>%
  mutate(river_km_label = case_when(
    river_km == 51 ~ "51.0", #so all axis values have 1 decimal place
    TRUE ~ as.character(river_km)
  ),
  river_km_label = factor(river_km_label, levels = c("10.6","11.1", "16.9", "20.8", "26.7", "34.4", "51.0")))

#wrap the x axis category labels
MRS_results <- MRS_results %>%
  mutate(category_wrapped = str_wrap(category, width = 10))

#heatmap plot
MRS_heatmap <- ggplot(MRS_results, aes(x = river_km_label, y = category_wrapped, fill = MRS_avg_for_plot)) +
  geom_tile(color = "black", width = 0.85, height = 0.85) +
  geom_text(aes(label = round(MRS_avg_for_plot, 1)), color = "white", fontface = "bold", size = 4.5) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "MRS", limits = c(0, NA)) +
  theme_minimal() +
  labs(x = "Distance from confluence (km)", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey80")
  )

#check out the plot (will combine with turnover plot later)
save_plot("MRS_heatmap.png",          
          MRS_heatmap,                          
          base_width = 7,             
          base_height = 3.5,             
          dpi = 300,                  
          bg = "white")

# Function to calculate species turnover for each comparison
calculate_turnover <- function(ev1, ev2, df) {
  sub_df <- df %>%
    filter(sampling_event_no %in% c(ev1, ev2)) %>%
    group_by(sampling_event_no, species) %>%
    summarise(rel_abund = mean(rel_abund, na.rm = TRUE), .groups = "drop") %>%
    mutate(present = ifelse(rel_abund > 0, 1, 0)) %>%
    select(sampling_event_no, species, present) %>%
    pivot_wider(names_from = sampling_event_no, values_from = present, values_fill = 0)
  
  tp <- names(sub_df)[-1]
  sp1 <- sub_df[[tp[1]]]
  sp2 <- sub_df[[tp[2]]]
  
  total_turnover <- sum(sp1 != sp2) / length(sp1)
  
  tibble(turnover = total_turnover)
}

#use the same comparison_table  and filtered for within-site comparisons as with MRS
event_rkm <- RAC_data %>%
  distinct(sampling_event_no, river_km)

comparison_table_rkm <- comparison_table %>%
  left_join(event_rkm, by = c("event_1" = "sampling_event_no")) %>%
  rename(rkm_1 = river_km) %>%
  left_join(event_rkm, by = c("event_2" = "sampling_event_no")) %>%
  rename(rkm_2 = river_km) %>%
  filter(rkm_1 == rkm_2) %>%
  mutate(river_km = rkm_1) %>%
  select(event_1, event_2, comparison_type, category, river_km)

#calculate species turnover
turnover_results <- comparison_table_rkm %>%
  mutate(result = map2(event_1, event_2, ~calculate_turnover(.x, .y, RAC_data))) %>%
  filter(!map_lgl(result, is.null)) %>%
  unnest(result) %>%
  select(event_1, event_2, comparison_type, category, turnover, river_km)

#again get a mean value for 11.1 season to season which has 2 values
turnover_results <- turnover_results %>%
  group_by(river_km, category) %>%
  mutate(turnover_avg_for_plot = if (n() > 1) mean(turnover) else turnover) %>%
  ungroup()

#update labels and wrap
turnover_results <- turnover_results %>%
  mutate(river_km_label = case_when(
    river_km == 51 ~ "51.0",
    TRUE ~ as.character(river_km)
  ),
  river_km_label = factor(river_km_label, levels = c("10.6","11.1", "16.9", "20.8", "26.7", "34.4", "51.0")),
  category_wrapped = str_wrap(category, width = 10))

#make the heatmap plot for turnover
turnover_heatmap <- ggplot(turnover_results, aes(x = river_km_label, y = category_wrapped, fill = turnover_avg_for_plot)) +
  geom_tile(color = "black", width = 0.85, height = 0.85) +
  geom_text(aes(label = round(turnover_avg_for_plot, 2)), color = "white", fontface = "bold", size = 4.5) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "Turnover", limits = c(0, 1)) +
  theme_minimal() +
  labs(x = "Distance from confluence (km)", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_line(color = "grey80")
  )

turnover_heatmap

#can save this plot to check it out, but will use combined plot below
save_plot("turnover_heatmap.png",          
          turnover_heatmap,                          
          base_width = 7,             
          base_height = 3.5,             
          dpi = 300,                  
          bg = "white")

#combine the 2 heatmap plots
combined_heatmap <- plot_grid(
  MRS_heatmap,
  turnover_heatmap,
  ncol = 1,        # Stack vertically
  labels = c("A", "B"),  # Optional: panel labels
  align = "v",     # Align vertically
  rel_heights = c(1, 1)  # Equal height for both
)

#save the combined plot
save_plot("combined_heatmap.png",
          combined_heatmap,
          base_width = 7,
          base_height = 7,
          dpi = 300,
          bg = "white")
