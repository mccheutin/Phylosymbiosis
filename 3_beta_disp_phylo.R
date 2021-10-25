## Dispersion results between species----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div")
load("beta.bc.RData") ; load("beta.uu.RData") ; load("beta.wu.RData")
load("beta.uu.sp.RData") ; load("beta.bc.sp.RData") ; load("beta.wu.sp.RData")
load("core_rel_tri_sp.RData") ; load("core_rel.RData")

# __ ** Species -----
library(vegan)
samp_data_sp <- core_rel_tri_sp %>% sample_data() %>% as.data.frame()
disp_sp.bc <- betadisper(otu.bc.sp, samp_data_sp$tax1)
Distances.bc <- disp_sp.bc$distances
disp_sp.wu <- betadisper(otu.wu.sp, samp_data_sp$tax1)
Distances.wu <- disp_sp.wu$distances
sig.sp.wu <- permutest(disp_sp.wu)
disp_sp.uu <- betadisper(otu.uu.sp, samp_data_sp$tax1)
Distances.uu <- disp_sp.uu$distances
sig.sp.uu <- permutest(disp_sp.uu)
samp_data_sp <- cbind("ID"= rownames(samp_data_sp), samp_data_sp)
disp.sp <- samp_data_sp %>% select(ID, tax1, family)
Distances <- rbind(cbind(disp.sp , "Distance" = Distances.bc, "Index" = rep("Bray-Curtis")), 
                   cbind(disp.sp , "Distance" = Distances.wu, "Index" = rep("W-Uni")),
                   cbind(disp.sp, "Distance" = Distances.uu, "Index" = rep("U-Uni")))

load(file = "../Physeq_objects/levels.tax.RData")
Distances$tax1 <- factor(Distances$tax1 , levels = levels_sp[levels_sp %in% Distances$tax1])
names_tree <- str_replace_all(levels(Distances$tax1) , c(" " = "_"))
host_tree_sp <- keep.tip(host_tree, names_tree)
write.tree(host_tree_sp, file= "host_tree_sp.tre")

write.table(Distances, file ="Distances_beta_sp.txt", sep ="\t", row.names = F)
Distances <- read.table(file ="Distances_beta_sp.txt", sep ="\t", header = T)

Distances$family <- factor(Distances$family , levels = levels_fam[levels_fam %in% Distances$family])

plot.disp.sp =  ggplot(Distances, aes(x = Distance , y = tax1, color = family)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0, size = 3) + theme_bw() +
  geom_boxplot(alpha=0.1, outlier.colour = NA) +
  scale_color_manual(values = colfunc(length(levels(Distances$family))))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=0, family = "serif", face = "italic"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=12, family = "serif"),
        legend.position="none",
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif"))+
  facet_wrap( ~ Index, nrow=1, ncol=3, scales = "free")

pdf(file ="plot.disp.sp.pdf", he = 7 , wi = 7)
plot.disp.sp
dev.off()

# __ ** Family -----
library(vegan)
load("../physeq_objects/gut_core.RData")
db <- sample_data(gut_core) %>% ad.data.frame()
fam_tri <- rownames(which(table(db$family)>=3, T))
fam_beta <- subset_samples(gut_core , family %in% fam_tri)
fam_beta_rel <- transform_sample_counts(fam_beta, function(x) x / sum(x) )
otu.bc.fam <- vegdist(fam_beta_rel@otu_table, method = "bray")
otu.wu.fam <- UniFrac(fam_beta_rel, weighted = T, normalized=F, parallel = F, fast=T)
otu.uu.fam <- UniFrac(fam_beta_rel, weighted = F, normalized=F, parallel = F, fast=T)
save(otu.bc.fam, otu.uu.fam, otu.wu.fam , file ="otu.matrice.fam.RData")


samp_data_fam <- fam_beta_rel %>% sample_data() %>% as.data.frame()
disp_fam.bc <- betadisper(otu.bc.fam, samp_data_fam$family)
Distances.bc <- disp_fam.bc$distances
disp_fam.wu <- betadisper(otu.wu.fam, samp_data_fam$family)
Distances.wu <- disp_fam.wu$distances
sig.fam.wu <- permutest(disp_fam.wu)
disp_fam.uu <- betadisper(otu.uu.fam, samp_data_fam$family)
Distances.uu <- disp_fam.uu$distances
sig.fam.uu <- permutest(disp_fam.uu)
samp_data_fam <- cbind("ID"= rownames(samp_data_fam), samp_data_fam)
disp.fam <- samp_data_fam %>% select(ID, tax1, family)
Distances <- rbind(cbind(disp.fam , "Distance" = Distances.bc, "Index" = rep("Bray-Curtis")), 
                   cbind(disp.fam , "Distance" = Distances.wu, "Index" = rep("W-Uni")),
                   cbind(disp.fam, "Distance" = Distances.uu, "Index" = rep("U-Uni")))


load(file = "../Physeq_objects/levels.tax.RData")
Distances$family <- factor(Distances$family , levels = levels_fam[levels_fam %in% Distances$family])
n <- length(levels(Distances$family))
colfunc <-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))
write.table(Distances , file ="Distances_betadisp_fam.txt", sep ="\t", row.names = F)

plot.disp.fam =  ggplot(Distances, aes(x = Distance , y = family, color = family)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0, size = 3) + theme_bw() +
  geom_boxplot(alpha=0.1, outlier.colour = NA) +
  scale_color_manual(values = colfunc(n))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif", face = "italic"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=12, family = "serif"),
        legend.position="none",
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif")) +
  facet_wrap( ~ Index, nrow=1, ncol=3, scales = "free")

pdf(file ="plot.disp.fam.pdf", he = 7 , wi = 7)
plot.disp.fam
dev.off()


## Mean ----
Distances_sp <- read.table("Distances_beta_sp.txt", sep = "\t", header = T)
mean_betadisp <- Distances_sp %>% filter(Index == "Bray-Curtis")  %>% summarize(mean = mean(Distance))
mean_betadisp %>% arrange(mean)


# Dispersion intra sp versus interspecies -----
library(vegan)
library(phyloseq)
library(ggplot2)

setwd("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div")
load("core_rel.RData")

## __ BC dissimilariy ----
bc_mat <- as.matrix(otu.bc)
rownames(bc_mat) <- samp_data$tax2
colnames(bc_mat) <- samp_data$tax2
bc_vec <- as.vector(bc_mat)
n <- length(colnames(bc_mat))
sp_1 <- rep(colnames(bc_mat), n)
sp_2 <- character()
for(i in 1:n) {
  vec_out <- rep(rownames(bc_mat)[i],n)
  sp_2 <- c(sp_2, vec_out)
}
sp_2

bc_tab <- as.data.frame(cbind("beta"= bc_vec, "sp_1" = sp_1, "sp_2" = sp_2))
intra_sp <- cbind(bc_tab[which(bc_tab$sp_1 == bc_tab$sp_2,T),] , "sp"= rep("Intra species"))
inter_sp <- cbind(bc_tab[which(bc_tab$sp_1 != bc_tab$sp_2,T),] , "sp"= rep("Inter species"))
bc_tab2 <- rbind(intra_sp, inter_sp)
bc_tab2$beta <- as.numeric(as.character(bc_tab2$beta))
bc_tab3 <- bc_tab2[which(bc_tab2$beta > 0,T), ]
save(bc_tab3, file = "bc_tab3.RData")

## __ Supp W-uni ----
wu_mat <- as.matrix(otu.wu)
rownames(wu_mat) <- samp_data$tax2
colnames(wu_mat) <- samp_data$tax2
wu_vec <- as.vector(wu_mat)
n <- length(colnames(wu_mat))
sp_1 <- rep(colnames(wu_mat), n)
sp_2 <- character()
for(i in 1:n) {
  vec_out <- rep(rownames(wu_mat)[i],n)
  sp_2 <- c(sp_2, vec_out)
}
sp_2

wu_tab <- as.data.frame(cbind("beta"= wu_vec, "sp_1" = sp_1, "sp_2" = sp_2))
intra_sp <- cbind(wu_tab[which(wu_tab$sp_1 == wu_tab$sp_2,T),] , "sp"= rep("Intra species"))
inter_sp <- cbind(wu_tab[which(wu_tab$sp_1 != wu_tab$sp_2,T),] , "sp"= rep("Inter species"))
wu_tab2 <- rbind(intra_sp, inter_sp)
wu_tab2$beta <- as.numeric(as.character(wu_tab2$beta))
wu_tab3 <- wu_tab2[which(wu_tab2$beta > 0,T), ]
save(wu_tab3, file = "wu_tab3.RData")

## __ Supp Un-uni ----
uu_mat <- as.matrix(otu.uu)
rownames(uu_mat) <- samp_data$tax2
colnames(uu_mat) <- samp_data$tax2
uu_vec <- as.vector(uu_mat)
n <- length(colnames(uu_mat))
sp_1 <- rep(colnames(uu_mat), n)
sp_2 <- character()
for(i in 1:n) {
  vec_out <- rep(rownames(uu_mat)[i],n)
  sp_2 <- c(sp_2, vec_out)
}
sp_2

uu_tab <- as.data.frame(cbind("beta"= uu_vec, "sp_1" = sp_1, "sp_2" = sp_2))
intra_sp <- cbind(uu_tab[which(uu_tab$sp_1 == uu_tab$sp_2,T),] , "sp"= rep("Intra species"))
inter_sp <- cbind(uu_tab[which(uu_tab$sp_1 != uu_tab$sp_2,T),] , "sp"= rep("Inter species"))
uu_tab2 <- rbind(intra_sp, inter_sp)
uu_tab2$beta <- as.numeric(as.character(uu_tab2$beta))
uu_tab3 <- uu_tab2[which(uu_tab2$beta > 0,T), ]
save(uu_tab3, file = "uu_tab3.RData")

## __ All metrics----
betadisp_tab <- rbind(cbind(bc_tab3 , "Distance" = rep("Bray-Curtis")),
                      cbind(wu_tab3 , "Distance" = rep("Weighted Unifrac")),
                      cbind(uu_tab3, "Distance" = rep ("Unweighted Unifrac")))
save(betadisp_tab , file ="betadisp_tab.RData")

betadisp_sp =  ggplot(betadisp_tab, aes(x = sp , y = beta)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 2) + theme_bw() +
  geom_boxplot(alpha=0.1, outlier.colour = NA) +
  labs(y = "Bacteriome Dissimilarity")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20,  family = "serif"),
        axis.text.y = element_text(size=20, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 20, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none", 
        strip.text.x = element_text(size=14, family = "serif")) + 
  facet_wrap(~ Distance, nrow = 1)+
  geom_signif(comparisons = list(c("Intra species", "Inter species")), 
              map_signif_level = TRUE, textsize=12, color="black", family = "serif", vjust = 0.5)

pdf(file = "betadisp_intra_inter_sp.pdf", he = 7 , wi = 12)
betadisp_sp
dev.off()

library(pgirmess)




