# PCOA---- 
setwd("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div")

load("core_rel.RData")
load("../phylogeny/host_tree_tri.RData")

# Metrics calcul ----
core_rel <- transform_sample_counts(gut_core, function(x) x / sum(x) )
samp_data <- data.frame(sample_data(core_rel))
save(core_rel, file =  "core_rel.RData")

# __ ** Bray-Curtis ----
# PCOA coda & Permanova
library(vegan)
otu.bc <- vegdist(core_rel@otu_table, method = "bray")
pcoa.sub.bc <- pcoa(otu.bc)
pcoa_coord.bc <- pcoa.sub.bc$vectors[,1:3]
# Contruction of the table for graphic 
library(stringr)
hull.bc <- cbind(pcoa_coord.bc, samp_data)

# What is the percentage of the explicative variance? 
paste("Axis 1 :",percent(pcoa.sub.bc$values$Relative_eig[1])) # 10 %
paste("Axis 2 :",percent(pcoa.sub.bc$values$Relative_eig[2])) # 6 %
paste("Axis 3 :",percent(pcoa.sub.bc$values$Relative_eig[3])) # 5 %
paste("Axis 4 :",percent(pcoa.sub.bc$values$Relative_eig[4])) # 4 %

save(otu.bc, hull.bc, file = "beta.bc.RData")

# __ ** Weighted Unifrac ----
# PCOA coda & Permanova
otu.wu <- UniFrac(core_rel, weighted = T, normalized=F, parallel = F, fast=T)
pcoa.sub.wu <- pcoa(otu.wu)
pcoa_coord.wu <- pcoa.sub.wu$vectors[,1:3]
hull.wu <- cbind(pcoa_coord.wu, samp_data)

# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.wu$values$Relative_eig[1])) # 26 %
paste("Axis 2 :",percent(pcoa.sub.wu$values$Relative_eig[2])) # 19 %
paste("Axis 3 :",percent(pcoa.sub.wu$values$Relative_eig[3])) # 12 %
paste("Axis 4 :",percent(pcoa.sub.wu$values$Relative_eig[4])) # 5 %
save(otu.wu, hull.wu, file = "beta.wu.RData")

# __ ** Unweighted Unifrac ----
# PCOA coda & Permanova
otu.uu <- UniFrac(core_rel, weighted = F, normalized=F, parallel = F, fast=T)
pcoa.sub.uu <- pcoa(otu.uu)
pcoa_coord.uu <- pcoa.sub.uu$vectors[,1:3]
hull.uu <- cbind(pcoa_coord.uu, samp_data)

# What is the percentage of the explicative variance? 
library(scales)
paste("Axis 1 :",percent(pcoa.sub.uu$values$Relative_eig[1])) # 14 %
paste("Axis 2 :",percent(pcoa.sub.uu$values$Relative_eig[2])) # 8 %
paste("Axis 3 :",percent(pcoa.sub.uu$values$Relative_eig[3])) # 7 %
paste("Axis 4 :",percent(pcoa.sub.uu$values$Relative_eig[4])) # 5 %

save(otu.uu, hull.uu, file = "beta.uu.RData")

## PCOA plots -----
# __ ** Bray-Curtis
hull_fam.bc <- subset(hull.bc , hull.bc$family %in% names(which(table(hull.bc$family) >= 3)))
write.table(hull_fam.bc, file = "hull_fam.bc.txt", quote = T, sep = "\t")
hull_fam.bc <- as.data.frame(read_delim("hull_fam.bc.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
hull_fam.bc$family <- factor(hull_fam.bc$family, levels = c("Pomacentridae", "Carangidae","Aulostomidae","Mullidae","Labridae",
                                                            "Scaridae","Acanthuridae","Zanclidae","Siganidae","Chaetodontidae",
                                                            "Haemulidae","Lutjanidae", "Serranidae", "Lethrinidae","Diodontidae",
                                                            "Monacanthidae","Ostraciidae","Holocentridae","Synodontidae"))
hull_fam.bc$diet2 <- str_replace(hull_fam.bc$diet2,"Planktonivorous", "Omnivorous")
hull_fam.bc$diet1 <- str_replace(as.factor(hull_fam.bc$diet1),"PK", "OM")
hull_fam.bc$diet1 <- factor(hull_fam.bc$diet1, levels = c("H","HD","OM", "FC","MI","SI"))
hull_fam.bc$diet2 <- factor(hull_fam.bc$diet2, levels = c("Herbivorous","Detritivorous","Omnivorous", "Piscivorous",
                                                          "Mobile invertebrate","Sessile invertebrate"))

save(hull_fam.bc, file = "hull_fam.bc.RData")

# __ ** W-Uni
hull_fam.wu <- subset(hull.wu , hull.wu$family %in% names(which(table(hull.wu$family) >= 3)))
write.table(hull_fam.wu, file = "hull_fam.wu.txt", quote = T, sep = "\t")
hull_fam.wu <- as.data.frame(read_delim("hull_fam.wu.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
hull_fam.wu$family <- factor(hull_fam.wu$family, levels = c("Pomacentridae", "Carangidae","Aulostomidae","Mullidae","Labridae",
                                                            "Scaridae","Acanthuridae","Zanclidae","Siganidae","Chaetodontidae",
                                                            "Haemulidae","Lutjanidae", "Serranidae", "Lethrinidae","Diodontidae",
                                                            "Monacanthidae","Ostraciidae","Holocentridae","Synodontidae"))
hull_fam.wu$diet2 <- str_replace(hull_fam.wu$diet2,"Planktonivorous", "Omnivorous")
hull_fam.wu$diet1 <- str_replace(as.factor(hull_fam.wu$diet1),"PK", "OM")
hull_fam.wu$diet1 <- factor(hull_fam.wu$diet1, levels = c("H","HD","OM", "FC","MI","SI"))
hull_fam.wu$diet2 <- factor(hull_fam.wu$diet2, levels = c("Herbivorous","Detritivorous","Omnivorous", "Piscivorous",
                                                          "Mobile invertebrate","Sessile invertebrate"))

save(hull_fam.wu, file = "hull_fam.wu.RData")

# __ ** U-Uni
hull_fam.uu <- subset(hull.uu , hull.uu$family %in% names(which(table(hull.uu$family) >= 3)))
write.table(hull_fam.uu, file = "hull_fam.uu.txt", quote = T, sep = "\t")
hull_fam.uu <- as.data.frame(read_delim("hull_fam.uu.txt", "\t", escape_double = FALSE, trim_ws = TRUE))
hull_fam.uu$family <- factor(hull_fam.uu$family, levels = c("Pomacentridae", "Carangidae","Aulostomidae","Mullidae","Labridae",
                                                            "Scaridae","Acanthuridae","Zanclidae","Siganidae","Chaetodontidae",
                                                            "Haemulidae","Lutjanidae", "Serranidae", "Lethrinidae","Diodontidae",
                                                            "Monacanthidae","Ostraciidae","Holocentridae","Synodontidae"))
hull_fam.uu$diet2 <- str_replace(hull_fam.uu$diet2,"Planktonivorous", "Omnivorous")
hull_fam.uu$diet1 <- str_replace(as.factor(hull_fam.uu$diet1),"PK", "OM")
hull_fam.uu$diet1 <- factor(hull_fam.uu$diet1, levels = c("H","HD","OM", "FC","MI","SI"))
hull_fam.uu$diet2 <- factor(hull_fam.uu$diet2, levels = c("Herbivorous","Detritivorous","Omnivorous", "Piscivorous",
                                                          "Mobile invertebrate","Sessile invertebrate"))

save(hull_fam.uu, file = "hull_fam.uu.RData")

n <- length(levels(hull_fam.wu$family))
colfunc <-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

# ... with viridis palette ----
pcoa.fam.bc <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_viridis_d() +
  geom_point(data = hull_fam.bc, aes(x=Axis.1, y=Axis.2, color = family), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.bc$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.bc$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 16, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))#+
  #facet_wrap(~ diet2)
pcoa.fam.bc

pcoa.fam.wu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_viridis_d() +
  geom_point(data = hull_fam.wu, aes(x=Axis.1, y=Axis.2, color = family), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.wu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.wu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 16, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))#+
  #facet_wrap(~ diet2)
pcoa.fam.wu

pcoa.fam.uu <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_viridis_d() +
  geom_point(data = hull_fam.uu, aes(x=Axis.1, y=Axis.2, color = family), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.uu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.uu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 16, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))#+
  #facet_wrap(~ diet2)
pcoa.fam.uu


pdf(file ="pcoa_fam_viridis.pdf", he = 7 , wi = 10)
pcoa.fam.bc
pcoa.fam.wu
pcoa.fam.uu
dev.off()

# ... with rainbow palette ----
colfunc <-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

pcoa.fam.bc2 <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = colfunc(n)) +
  geom_point(data = hull_fam.bc, aes(x=Axis.1, y=Axis.2, color = family), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.bc$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.bc$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 16, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))#+
  #facet_wrap(~ diet2)
pcoa.fam.bc2

pcoa.fam.wu2 <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = colfunc(n)) +
  geom_point(data = hull_fam.wu, aes(x=Axis.1, y=Axis.2, color = family), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.wu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.wu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 16, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))#+
  #facet_wrap(~ diet2)
pcoa.fam.wu2

pcoa.fam.uu2 <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  scale_color_manual(values = colfunc(n)) +
  geom_point(data = hull_fam.uu, aes(x=Axis.1, y=Axis.2, color = family), alpha = 0.7, size = 4, shape = 16) +
  xlab(paste("PCo1 (", round(pcoa.sub.uu$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub.uu$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=18, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=18, family = "serif"), # remove y-axis labels
        axis.text.x = element_text(size=18, family = "serif"),
        axis.text.y = element_text(size=18, family = "serif"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 16, family = "serif"),
        legend.title = element_text(size = 0,family = "serif"),
        strip.text.x = element_text(size=14, family = "serif"))#+
  #facet_wrap(~ diet2)
pcoa.fam.uu2

pdf(file ="pcoa_fam_rainbow_entire.pdf", he = 7 , wi = 10)
pcoa.fam.bc2
pcoa.fam.wu2
pcoa.fam.uu2
dev.off()
