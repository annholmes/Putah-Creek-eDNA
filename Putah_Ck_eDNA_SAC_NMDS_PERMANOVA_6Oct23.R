##Figure 2 Species accumulation curves

library(phyloseq)
library(vegan)
library(ggplot2)

#read in file if needed
ps.decontam <- readRDS("~/Desktop/github/Putah-Creek-eDNA/PC_all/ps.decontam.rds")

otu_matrix <-as.matrix(otu_table(ps.decontam))
df_otu <-as.data.frame(otu_matrix)

S <- specnumber(df_otu) # observed number of species in each of 58 field samples

library(viridis)
#select 58 colors from viridis
colors <- viridis_pal(option = "D")(58)

fish_spp_accum<- rarecurve(df_otu, 
          col = colors, 
          step=20, 
          lwd=3, 
          xlab= "Number of sequences", 
          ylab="Number of species detected", 
          main="Species accumulaton by sequencing depth",
          label=FALSE)
#this phyloseq object has not been fully decontaminated- the manual stuff in excel isn't reflected here
#max # of spp should be 18 but looks to be 21 or 22 in the SAC; ok for now but need to redo


#cowplot, ggsave and png() don't work for the object 
#just saved as pdf for now 6 Oct 23

##Figure 4 NMDS with stress value

m <- betadiver(df_otu)
plot(m)
## The indices
betadiver(help=TRUE)
## The basic Whittaker index
d <- betadiver(df_otu, "w")
## This should be equal to Sorensen index (binary Bray-Curtis in vegan)
range(d - vegdist(df_otu, binary=TRUE))

ordination_palette <- c("#009999", "#fbbf77", "#6f2da8", "#006ddb", "#22cf22",
                      "#db6d00", "#b0c0c9", "#ffdf4d", "#30D5C8", "#ff66d1", 
                      "#cfa382", "#601d2c", 
                      "#00008B", "#f8766d", "#c7f6b6", 
                      "#8f4e00", "#920000", "#ffe6ee", "#444444", "#b66dff", 
                      "#006400", "#c5f0ee")
sample_data(ps.decontam)$river_mile <- as.factor(sample_data(ps.decontam)$river_mile)

samples.ord <- ordinate(ps.decontam, "NMDS", "bray")

#produce a stress value
#https://ourcodingclub.github.io/tutorials/ordination/
#also here https://rpubs.com/CPEL/NMDS

dist <- vegdist(df_otu,  method = "bray")

# define a function NMDS.scree()
# performs a NMDS for 1-10 dimensions and plots the nr of dimensions vs the stress
NMDS.scree <- function(x) { #where x is the name of the data frame variable
  plot(rep(1, 10), replicate(10, metaMDS(x, autotransform = F, k = 1)$stress), xlim = c(1, 10),ylim = c(0, 0.30), xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")
  for (i in 1:10) {
    points(rep(i + 1,10),replicate(10, metaMDS(x, autotransform = F, k = i + 1)$stress))
  }
}

# use the function to choose the optimal nr of dimensions
NMDS.scree(dist)

# Because the final result depends on the initial 
# random placement of the points 
# we`ll set a seed to make the results reproducible
set.seed(2)

# Here, we perform the final analysis and check the result
NMDS1 <- metaMDS(dist, k = 2, trymax = 100, trace = F)
# Do you know what the trymax = 100 and trace = F means?
# Let's check the results
NMDS1
#stress value = 0.15 (stress <0.20 is ok; stress <0.10 is better)

stressplot(NMDS1)
#Non-metric fit, R2 = 0.977
#Linear fit, R2 = 0.89

#plot with stress value

PC_NMDS <- plot_ordination(ps.decontam, samples.ord, 
                           type = "sample_id", 
                           color = "river_mile", 
                           shape = "season") + 
  scale_shape_manual(values = c(17,15)) +
  scale_color_manual(values = ordination_palette) +
  geom_jitter(size = 8) +
  theme_bw() +
  theme(legend.title = NULL,
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)) +
  annotate("text", 
           x = .8, y = .85, 
           size = 7,
           label = "Stress = 0.15")

PC_NMDS

library(cowplot)
save_plot("Putah_Ck_NMDS_6Oct23.png", 
          PC_NMDS, 
          dpi = 300,
          base_height = 5, #height in inches
          base_width = 8, #width in inches
          bg = "white")

#PERMANOVA using adonis

sample_data(ps.decontam)$river_mile <- as.numeric(sample_data(ps.decontam)$river_mile)

ps.decontam_bray <-phyloseq::distance(ps.decontam, method="bray")
df_ps.decontam <- data.frame(sample_data(ps.decontam))

adonis2(ps.decontam_bray~river_mile, data = df_ps.decontam)
adonis2(ps.decontam_bray~site, data = df_ps.decontam)
adonis2(ps.decontam_bray~Region, data = df_ps.decontam)
#all are significant p<.001

adonis2(ps.decontam_bray~season, data = df_ps.decontam)
#significant at p<.05 (p=0.2)

