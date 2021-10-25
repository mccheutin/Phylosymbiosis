####### Alpha-diversity ----
setwd( "/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis")
alpha_div <- "./alpha_div"
dir.create(alpha_div, recursive = T)
setwd(alpha_div)

rff_core <- rarefy_even_depth(gut_core, sample.size = 1808)
rff_core
cov_core_rff <- goods(rff_core@otu_table@.Data)
paste0(mean(cov_core_rff$goods)," +- " ,sd(cov_core_rff$goods))
save(rff_core, file = "rff_core.RData")

# __ Metrics calcul ----
library(picante)
library(entropart)
box1 = estimate_richness(rff_core , measures = c("Shannon", "Chao1"))
box1$Shannon = exp(box1$Shannon)
box2 = pd(samp= rff_core@otu_table , tree = rff_core@phy_tree , include.root =  T)
box3 <- as.data.frame(cbind("Sample" = rownames(box1), box1, box2))
db <- rff_core %>% sample_data() %>% as.data.frame()
db <- cbind("Sample" = rownames(db), db)
library(stringr)
db$Sample <- str_replace_all(db$Sample, "-", ".")
box3 <- left_join(box3, db , by = "Sample") 
save(box3, file = "box.alpha.RData")

mean(box3$SR) ; sd(box3$SR)
mean(box3$PD) ; sd(box3$PD)
mean(box3$Shannon) ; sd(box3$Shannon)
mean(box3$Chao1) ; sd(box3$Chao1)

### __ Alpha between dietecies----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/alpha_div")
dir.create("./dietecies")
setwd("./dietecies")

load("../../phylogeny/core.phylosymb.tri.RData")
rff_core_diet_merged <- rarefy_even_depth(core.phylosymb.tri, sample.size = 7670)
cov_core_rff <- goods(rff_core_diet_merged@otu_table@.Data)
paste0(mean(cov_core_rff$goods)," +- " ,sd(cov_core_rff$goods))
save(rff_core_diet_merged, file = "rff_core_diet_merged.RData")

## ** Calcul alpha ----
load("./rff_core_diet_merged")
library(phyloseq)
library(picante)
box1_diet = estimate_richness(rff_core_diet_merged , measures = c("Shannon", "Chao1"))
box1_diet$Shannon = exp(box1_diet$Shannon)
box2_diet = pd(samp= rff_core_diet_merged@otu_table , tree = rff_core_diet_merged@phy_tree , include.root =  T)

alpha_box_diet <- as.data.frame(rbind(cbind("tax2" = rep(rownames(box1_diet)),
                         "variable" = rep("Observed"),
                         "value" = box2_diet$SR),
                   cbind("tax2" = rep(rownames(box1_diet)) ,
                         "variable" = rep("PD"),
                         "value" = box2_diet$PD),
                   cbind("tax2" = rep(rownames(box1_diet)) ,
                         "variable" = rep("Shannon"),
                         "value" = box1_diet$Shannon),
                   cbind("tax2" = rep(rownames(box1_diet)) ,
                         "variable" = rep("Chao1"),
                         "value" = box1_diet$Chao1)))

alpha_box_diet$value <- as.numeric(alpha_box_diet$value)
save(alpha_box_diet, file = "alpha_box_diet.RData")
write.table(alpha_box_diet, file ="alpha_box_diet.txt", sep = "\t", row.names = F, quote = F)

library(entropart)
library(dplyr)
box_obs <- alpha_box_diet %>% filter(variable == 'Observed')  
box_pd <- alpha_box_diet %>% filter(variable == 'PD')  
box_shannon <- alpha_box_diet %>% filter(variable == 'Shannon')
box_chao <- alpha_box_diet %>% filter(variable == 'Chao1')

data_alpha <- cbind('Observed' = as.numeric(box_obs$value), 
                    'PD' = as.numeric(box_pd$value), 
                    'Shannon' = as.numeric(box_shannon$value),
                    'Chao1' = as.numeric(box_chao$value))
rownames(data_alpha) <- box_obs$tax2  
save(data_alpha, file = "data_alpha.RData")

## ** Moran Lippa ----
library(phylobase)
load("../../phylogeny/host_tree_tri.RData") ; load("./data_alpha.RData") 
data_tree <- host_tree_tri 
write.tree(data_tree, file= "data_tree.tre")
data_alpha

g1 <- as(data_tree, "phylo4")
g2 <- phylo4d(g1, data_alpha, missing.data="warn")
head(g2)

library(phylosignal)
barplot.phylo4d(g2)
I_result <- phyloSignal(g2)
table_I <- as.data.frame(cbind("Alpha measure" = rownames(I_result$stat),"Index" = I_result$stat$I , "p-value" = I_result$pvalue$I))
write.table(table_I,file ="table_I.txt", sep = "\t", quote = F, row.names= F)

observed.cg <- phyloCorrelogram(g2, trait = "Observed")
PD.cg <- phyloCorrelogram(g2, trait = "PD")
shannon.cg <- phyloCorrelogram(g2, trait = "Shannon")
chao.cg <- phyloCorrelogram(g2, trait = "Chao1")
plot(observed.cg)
plot(PD.cg)
plot(shannon.cg)
plot(chao.cg)

local.i.Observed <- lipaMoran(g2, trait = "Observed" ,prox.phylo = "nNodes", as.p4d = TRUE)
points.Observed <- lipaMoran(g2, trait = "Observed",prox.phylo = "nNodes")$p.value
points.col.Observed <- ifelse(points.Observed < 0.05, "red", "black")

obs_plot <- dotplot.phylo4d(local.i.Observed, dot.col = points.col.Observed)

local.i.PD <- lipaMoran(g2, trait = "PD" ,prox.phylo = "nNodes", as.p4d = TRUE)
points.PD <- lipaMoran(g2, trait = "PD",prox.phylo = "nNodes")$p.value
points.col.PD <- ifelse(points.PD < 0.05, "red", "black")
PD_plot <- dotplot.phylo4d(local.i.PD, dot.col = points.col.PD)

local.i.shannon <- lipaMoran(g2, trait = "Shannon" ,prox.phylo = "nNodes", as.p4d = TRUE)
points.shannon <- lipaMoran(g2, trait = "Shannon",prox.phylo = "nNodes")$p.value
points.col.shannon <- ifelse(points.shannon < 0.05, "red", "black")
shannon_plot <- dotplot.phylo4d(local.i.shannon, dot.col = points.col.shannon)

local.i.chao <- lipaMoran(g2, trait = "Chao1" ,prox.phylo = "nNodes", as.p4d = TRUE)
points.chao <- lipaMoran(g2, trait = "Chao1",prox.phylo = "nNodes")$p.value
points.col.chao <- ifelse(points.chao < 0.05, "red", "black")
chao_plot <- dotplot.phylo4d(local.i.chao, dot.col = points.col.chao)

points_i <- as.data.frame(rbind(cbind('tax2' = rownames(points.col.Observed),
                        'Variable' = rep('Observed'), 
                        'I-Moran' = as.numeric(local.i.Observed@data$Observed), 
                        'p-value' = as.numeric(points.Observed[,1])),
                  (cbind('tax2' = rownames(points.col.Observed),
                         'Variable' = rep('PD'), 
                         'I-Moran' = as.numeric(local.i.PD@data$PD), 
                         'p-value' = as.numeric(points.PD[,1]))),
                  (cbind('tax2' = rownames(points.col.Observed),
                         'Variable' = rep('Shannon'), 
                         'I-Moran' = as.numeric(local.i.shannon@data$Shannon), 
                         'p-value' = as.numeric(points.shannon[,1]))),
                  cbind('tax2' = rownames(points.col.Observed),
                        'Variable' = rep('Chao1'), 
                        'I-Moran' = as.numeric(local.i.chao@data$Chao1), 
                        'p-value' = as.numeric(points.chao[,1]))))

points_i$'I-Moran' <- round(as.numeric(points_i$'I-Moran'), 3)

write.table(points_i , file = "points_i.txt", row.names = F, quote = F, sep ="\t")
points_i <- read.table(file = "points_i.txt", sep = '\t', header = T)


pdf("Moran_alpha_index.pdf", he = 7, wi = 7)
# run lines
dev.off()

## ** Plot ----
load(file = "../box_alpha.RData")

diet_to_keep <- levels(factor(box3$tax2))[which(levels(factor(box3$tax2)) %in% levels(factor(points_i$tax2)))]
box3_diet <-  subset(box3, tax2 %in% diet_to_keep)
write.table(box3_diet , file = "box3_diet.txt", sep = "\t", quote = T, row.names= F)
box3_diet <- read.table(file ="box3_diet.txt", sep = '\t', header = T)
levels <- data_tree$tip.label[which(data_tree$tip.label %in% levels(factor(box3_diet$tax2)))]
box3_diet$tax2 <- factor(box3_diet$tax2 , levels = levels)
levels <- str_replace(levels, "_", " ")
box3_diet$tax1 <- factor(box3_diet$tax1 , levels = levels)

box3_diet_obs <- box3_diet$SR
box3_diet_PD <- box3_diet$PD
box3_diet_shannon <- box3_diet$Shannon
box3_diet_chao <- box3_diet$Chao1
box3_diet_value <- as.data.frame(rbind(cbind("variable" = rep("Observed"), "value" = box3_diet_obs),
                       cbind("variable" = rep("PD"), "value" = box3_diet_PD),
                       cbind("variable" = rep("Shannon"), "value" = box3_diet_shannon),
                       cbind("variable" = rep("Chao1"), "value" = box3_diet_chao)))

box_diet_plot <- box3_diet %>% select(Sample, tax1, tax2, family, genus, diet1, diet2, diet3, region, ocean)
box3_diet_plot <- cbind(box_diet_plot, box3_diet_value)
box3_diet_plot$value <- as.numeric(box3_diet_plot$value)
box3_diet_plot$variable <- factor(box3_diet_plot$variable , levels = c("Observed", "PD","Shannon","Chao1"))
save(box3_diet_plot, file = "box3_diet_plot.RData")

colfunc <-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

plot.alpha.diet.rainbow =  ggplot(box3_diet_plot, aes(x = value , y = tax1, color = tax1)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0, size = 3) + theme_bw() +
  geom_boxplot(alpha=0.1, outlier.colour = NA) +
  scale_color_manual(values = colfunc(48))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=0, family = "serif", face = "italic"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, family = "serif"),
        legend.position="none",
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif"))+
  facet_wrap( ~ variable, nrow=1, ncol=4, scales = "free")

plot.alpha.diet.rainbow

names =  ggplot(box3_diet_plot[box3_diet_plot$variable == "Observed",], aes(x = value , y = tax1, color = tax1)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0, size = 3) + theme_bw() +
  geom_boxplot(alpha=0.1, outlier.colour = NA) +
  scale_color_viridis_d()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif", face = "italic"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, family = "serif"),
        legend.position="none",
        legend.text = element_text(family = "serif"),
        legend.title = element_text(family = "serif"))
names

plot_tree(tree_data)


pdf(file = "alpha_diet.pdf" , he = 7 , wi = 12)
plot.alpha.diet.rainbow
names
dev.off()

## ** Dunn test ----
load("./box3_diet_plot.RData")
library(dunn.test)
obs.diet <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Observed"] , box3_diet_plot$tax1[box3_diet_plot$variable=="Observed"], method = "bonferroni")
pd.diet <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "PD"] , box3_diet_plot$tax1[box3_diet_plot$variable=="PD"], method = "bonferroni")
shannon.diet <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Shannon"] , box3_diet_plot$tax1[box3_diet_plot$variable=="Shannon"], method = "bonferroni")
chao1.diet <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Chao1"] , box3_diet_plot$tax1[box3_diet_plot$variable=="Chao1"], method = "bonferroni")
comp.diet <- as.data.frame(cbind("dietecies" = obs.diet$comparisons, "Observed" = round(obs.diet$P.adjusted,3) , "Faith" =  round(pd.diet$P.adjusted,3),"Shannon"= round(shannon.diet$P.adjusted,3), "Chao1" = round(chao1.diet$P.adjusted,3)))
write.table(comp.diet, file = "../comp.diet.txt", sep = "\t", quote = F, row.names = F)

### __ Alpha between families----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/alpha_div")
dir.create("./Family")
setwd("./Family")

load("../box.alpha.RData")
box3_fam <- subset(box3, family %in% names(which(table(box3$family) >= 3)))
write.table(box3_fam , file = "box3_fam.txt", sep = "\t", quote = T, row.names= F)
box3_fam <- read.table(file ="box3_fam.txt", sep = '\t', header = T)

## ** Plot ----
box3_fam_obs <- box3_fam$SR
box3_fam_PD <- box3_fam$PD
box3_fam_shannon <- box3_fam$Shannon
box3_fam_chao <- box3_fam$Chao1
box3_fam_value <- as.data.frame(rbind(cbind("variable" = rep("Observed"), "value" = box3_fam_obs),
                                     cbind("variable" = rep("PD"), "value" = box3_fam_PD),
                                     cbind("variable" = rep("Shannon"), "value" = box3_fam_shannon),
                                     cbind("variable" = rep("Chao1"), "value" = box3_fam_chao)))

box_fam_plot <- box3_fam %>% select(Sample, tax1, tax2, family, genus, diet1, diet2, diet3, region, ocean)
box3_fam_plot <- cbind(box_fam_plot, box3_fam_value)
box3_fam_plot$value <- as.numeric(box3_fam_plot$value)
box3_fam_plot$variable <- factor(box3_fam_plot$variable , levels = c("Observed", "PD","Shannon","Chao1"))
save(box3_fam_plot, file = "box3_fam_plot.RData")

## ** Dunn test ----
load("./box3_fam_plot.RData")
library(dunn.test)
obs.fam <- dunn.test(box3_fam_plot$value[box3_fam_plot$variable == "Observed"] , box3_fam_plot$family[box3_fam_plot$variable=="Observed"], method = "bonferroni")
pd.fam <- dunn.test(box3_fam_plot$value[box3_fam_plot$variable == "PD"] , box3_fam_plot$family[box3_fam_plot$variable=="PD"], method = "bonferroni")
shannon.fam <- dunn.test(box3_fam_plot$value[box3_fam_plot$variable == "Shannon"] , box3_fam_plot$family[box3_fam_plot$variable=="Shannon"], method = "bonferroni")
chao1.fam <- dunn.test(box3_fam_plot$value[box3_fam_plot$variable == "Chao1"] , box3_fam_plot$family[box3_fam_plot$variable=="Chao1"], method = "bonferroni")
comp.fam <- as.data.frame(cbind("Family" = obs.fam$comparisons, "Observed" = round(obs.fam$P.adjusted,3) , "Faith" =  round(pd.fam$P.adjusted,3),"Shannon"= round(shannon.fam$P.adjusted,3), "Chao1" = round(chao1.fam$P.adjusted,3)))
write.table(comp.fam, file = "comp.fam.txt", sep = "\t", quote = F, row.names = F)

# __ Alpha between diet ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/alpha_div")
dir.create("./Diet")
setwd("./Diet")
load("../box.alpha.RData")

box3_diet <- box3
levels(box3$diet2)[2] <- "Strict Herbivores"
box3_diet$diet2 <- box3_diet$diet3
box3_diet <- rbind(box3, box3_diet)
levels(box3_diet$diet2) <- c("Detritivores", "Strict Herbivores", "Mobile invertivores", 
                        "Omnivores","Piscivores","Planktonivores","Sessile invertivores","Carnivores",
                        "Herbivores")

box3_diet$diet2 <- factor(box3_diet$diet2, 
                          levels = c("Herbivores", "Detritivores","Strict Herbivores",
                                    "Omnivores", "Carnivores", "Planktonivores",
                                    "Sessile invertivores", "Mobile invertivores", "Piscivores"))
                           

write.table(box3_diet , file ="box3_diet.txt" , sep = "\t" , quote = T, row.names = F)
box3_diet <- read.table("box3_diet.txt", sep ="\t", header = T)

## ** Plot ----
box3_diet_obs <- box3_diet$SR
box3_diet_PD <- box3_diet$PD
box3_diet_shannon <- box3_diet$Shannon
box3_diet_chao <- box3_diet$Chao1
box3_diet_value <- as.data.frame(rbind(cbind("variable" = rep("Observed"), "value" = box3_diet_obs),
                                      cbind("variable" = rep("PD"), "value" = box3_diet_PD),
                                      cbind("variable" = rep("Shannon"), "value" = box3_diet_shannon),
                                      cbind("variable" = rep("Chao1"), "value" = box3_diet_chao)))

box_diet_plot <- box3_diet %>% select(Sample, tax1, tax2, family, genus, diet1, diet2, diet3, region, ocean)
box3_diet_plot <- cbind(box_diet_plot, box3_diet_value)
box3_diet_plot$value <- as.numeric(box3_diet_plot$value)
box3_diet_plot$variable <- factor(box3_diet_plot$variable , levels = c("Observed", "PD","Shannon","Chao1"))
box3_diet_plot$diet2 <- factor(box3_diet_plot$diet2 , levels = c("Detritivores","Strict Herbivores","Herbivores", 
                                                       "Omnivores",  "Planktonivores",
                                                       "Sessile invertivores", "Mobile invertivores", "Piscivores","Carnivores"))
save(box3_diet_plot, file = "box3_diet_plot.RData")


diet_colors <- c("darkkhaki", "green" , "darkgreen" ,"darkblue", "gray", "yellow", "darkorange","tan4","darkred")
pdf(file = "alpha_diet_all.pdf" , he = 7 , wi = 7)
diet_alpha_div =  ggplot(box3_diet_plot, aes(x = value , y = diet2, color = diet2)) +
  geom_jitter(position = position_jitter(width = .20), alpha = 0.5, size = 3) + theme_bw() +
  scale_colour_manual(values= diet_colors) +
  geom_boxplot(alpha=0.1, outlier.colour = NA) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=12, family = "serif"),
        axis.text = element_text(family = "serif",size = 14),
        axis.text.x = element_text(family = "serif",size = 14, angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14, family = "serif"),
        legend.position="none")+
  facet_grid( ~ variable,  scales = "free")
diet_alpha_div

pdf(file = "alpha_plot_diet.pdf" , he = 5 , wi = 10)
diet_alpha_div
dev.off()

# ** Dunn test ----
library(dunn.test)
obs.diet1 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Observed"] , box3_diet_plot$diet1[box3_diet_plot$variable=="Observed"], method = "bonferroni")
pd.diet1 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "PD"] , box3_diet_plot$diet1[box3_diet_plot$variable=="PD"], method = "bonferroni")
shannon.diet1 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Shannon"] , box3_diet_plot$diet1[box3_diet_plot$variable=="Shannon"], method = "bonferroni")
chao1.diet1 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Chao1"] , box3_diet_plot$diet1[box3_diet_plot$variable=="Chao1"], method = "bonferroni")
comp.diet1 <- as.data.frame(cbind("Diet" = obs.diet1$comparisons, "Observed" = round(obs.diet1$P.adjusted,3) , "Faith" =  round(pd.diet1$P.adjusted,3),"Shannon"= round(shannon.diet1$P.adjusted,3), "Chao1" = round(chao1.diet1$P.adjusted,3)))
write.table(comp.diet1, file = "./comp.diet1.txt", sep = "\t", quote = F, row.names = F)

obs.diet3 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Observed"] , box3_diet_plot$diet3[box3_diet_plot$variable=="Observed"], method = "bonferroni")
pd.diet3 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "PD"] , box3_diet_plot$diet3[box3_diet_plot$variable=="PD"], method = "bonferroni")
shannon.diet3 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Shannon"] , box3_diet_plot$diet3[box3_diet_plot$variable=="Shannon"], method = "bonferroni")
chao1.diet3 <- dunn.test(box3_diet_plot$value[box3_diet_plot$variable == "Chao1"] , box3_diet_plot$diet3[box3_diet_plot$variable=="Chao1"], method = "bonferroni")
comp.diet3 <- as.data.frame(cbind("Diet" = obs.diet3$comparisons, "Observed" = round(obs.diet3$P.adjusted,3) , "Faith" =  round(pd.diet3$P.adjusted,3),"Shannon"= round(shannon.diet3$P.adjusted,3), "Chao1" = round(chao1.diet3$P.adjusted,3)))
write.table(comp.diet3, file = "./comp.diet3.txt", sep = "\t", quote = F, row.names = F)

