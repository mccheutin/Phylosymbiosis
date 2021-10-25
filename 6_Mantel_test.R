## Phylosymbiosis (B-divergence) - Mantel test -----
# Hypothesis: the more families diverged, the beta diversity increased.
# 
# Setup ----
setwd("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/phylogeny")
dir.create("./matrices", recursive = T)
dir.create("./Figures", recursive = T)

library("ape")
library("Biostrings")
library("ggplot2")
library(ggtree)

load("../Physeq_objects/gut_core.RData")
fish_tree <- read.tree(file = "actinopt_12k_treePL_test.tre")

samp_data <- sample_data(gut_core)
tax_names <- levels(factor(samp_data$tax2))
out_tax_names <- tax_names[which(!tax_names %in% fish_tree$tip.label)]
# First step, replace missing species in a tree manually thanks to TimeTree.org to 
# get the missing and replace them in a simplified tree in TreeGraph2.
in_tax_names <- tax_names[which(tax_names %in% fish_tree$tip.label)]
simple_tree <- keep.tip(fish_tree, in_tax_names)
write.tree(my_tree, file= "simple_tree.tre")

#Reimport the modified tree 
host_tree <- read.tree(file = "fish_tree_modif.tre")
mat_cophenetic <- cophenetic.phylo(host_tree)/2
host_sp <- tax_names[which(tax_names %in% host_tree$tip.label)]

# Mantel tests ----
# __ Bacteriome ----
# __ ** On all host ----
core.phylosymb <- subset_samples(gut_core, tax2 %in% host_sp)
core.phylosymb <- merge_samples(core.phylosymb, "tax2", fun = mean) #pour rassembler les abondances dans chque grande catégories d'échantillons
sample_data(core.phylosymb)$tax2 <- rownames(sample_data(core.phylosymb))
samp_data2 <- samp_data[which(samp_data$tax2 %in% sample_data(core.phylosymb)$tax2, T),]
samp_data_sp_col <- unique(samp_data2[, c(1,2,4,6,7,8,10,11,12,21)])
samp_data_sp_rows <- samp_data_sp_col[samp_data_sp_col$tax2 %in% rownames(sample_data(core.phylosymb)),]
samp_data_sp <- samp_data_sp_rows[order(samp_data_sp_rows$tax2),]
rownames(samp_data_sp) <- rownames(sample_data(core.phylosymb))
sample_data(core.phylosymb) <- samp_data_sp
save(core.phylosymb, file = "core.phylosymb.RData")
# Calcul of the matrix with species merged 
core.phylosymb.rel <- transform_sample_counts(core.phylosymb, function(x) x / sum(x) )
otu.jacc.phylo <- as.matrix(vegdist(core.phylosymb.rel@otu_table, method = "jaccard"))
otu.bc.phylo <- as.matrix(vegdist(core.phylosymb.rel@otu_table, method = "bray"))
otu.wu.phylo <- as.matrix(UniFrac(core.phylosymb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo <- as.matrix(UniFrac(core.phylosymb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo,otu.bc.phylo,otu.wu.phylo,otu.uu.phylo , file = "matrices/beta_matrices.all.RData")
# subset both matrices in case some species are only in one
dim(mat_cophenetic) ; dim(as.matrix(otu.bc.phylo))
#mat_cophenetic <-mat_cophenetic[which(rownames(mat_cophenetic) %in% rownames(as.matrix(otu.jacc.phylo))),
                               #which(colnames(mat_cophenetic) %in% colnames(as.matrix(otu.jacc.phylo)))]
#otu.jacc.phylo <- as.matrix(otu.jacc.phylo)[which(rownames(as.matrix(otu.jacc.phylo)) %in% rownames(mat_cophenetic)),
                          # which(colnames(as.matrix(otu.jacc.phylo)) %in% colnames(mat_cophenetic))]

# reorder cophenetic matrix' rows and columns so they are in the ame order as in otu.bc.phylo
mat_cophenetic <- mat_cophenetic[match(rownames(otu.jacc.phylo), rownames(mat_cophenetic)),
                               match(colnames(otu.jacc.phylo), colnames(mat_cophenetic))]

save(mat_cophenetic, file = "matrices/mat_cophenetic.RData")
jacc.all.mantel = mantel(otu.jacc.phylo, mat_cophenetic , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.all.mantel = mantel(otu.bc.phylo, mat_cophenetic , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.all.mantel = mantel(otu.wu.phylo, mat_cophenetic , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.all.mantel = mantel(otu.uu.phylo, mat_cophenetic , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_jacc_all <- paste("r =", round(jacc.all.mantel$statistic,3), ";","p =",round(jacc.all.mantel$signif,3))
phylo_bc_all <- paste("r =", round(bc.all.mantel$statistic,3), ";","p =",round(bc.all.mantel$signif,3))
phylo_wu_all <- paste("r =", round(wu.all.mantel$statistic,3), ";","p =",round(wu.all.mantel$signif,3))
phylo_uu_all <- paste("r =", round(uu.all.mantel$statistic,3), ";","p =",round(uu.all.mantel$signif,3))
phylo_all <- cbind("Jaccard" = phylo_jacc_all, "Bray-Curtis" =phylo_bc_all,"Weighted Unifrac" =phylo_wu_all,"Unweighted Unifrac" =phylo_uu_all)

# __ ** Triplicate species ----
#host_tree <- read.tree(file = "fish_tree_modif.tre") #tree modified after replacment of species thanks to the litteracy
sp_tri <- rownames(which(table(samp_data$tax2)>=3,T))
in_tax_names_tri <- sp_tri[which(sp_tri %in% host_tree$tip.label)]
samp_tri <- samp_data_sp[samp_data_sp$tax2 %in% in_tax_names_tri,]

host_tree_tri <- keep.tip(host_tree, in_tax_names_tri)
save(host_tree_tri, file ="host_tree_tri.RData")
mat_cophenetic_tri<- cophenetic.phylo(host_tree_tri)/2
core.phylosymb.tri <- subset_samples(core.phylosymb, tax2 %in% in_tax_names_tri)
save(core.phylosymb.tri, file = "core.phylosymb.tri.RData")
# Calcul of the matrix with species merged 
core.phylosymb.tri.rel <- transform_sample_counts(core.phylosymb.tri, function(x) x / sum(x))
save(core.phylosymb.tri.rel, file = "core.phylosymb.tri.rel.RData")

otu.jacc.phylo.tri <- as.matrix(vegdist(core.phylosymb.tri.rel@otu_table, method = "jaccard"))
otu.bc.phylo.tri <- as.matrix(vegdist(core.phylosymb.tri.rel@otu_table, method = "bray"))
otu.wu.phylo.tri <- as.matrix(UniFrac(core.phylosymb.tri.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.tri <- as.matrix(UniFrac(core.phylosymb.tri.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.tri,otu.bc.phylo.tri,otu.wu.phylo.tri,otu.uu.phylo.tri , file = "matrices/beta_matrices.tri.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_tri) ; dim(otu.bc.phylo.tri)
mat_cophenetic_tri <- mat_cophenetic_tri[match(rownames(otu.jacc.phylo.tri), rownames(mat_cophenetic_tri)),
                                 match(colnames(otu.jacc.phylo.tri), colnames(mat_cophenetic_tri))]

save(mat_cophenetic_tri, file = "matrices/mat_cophenetic_tri.RData")
jacc.tri.mantel = mantel(otu.jacc.phylo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.tri.mantel = mantel(otu.bc.phylo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.tri.mantel = mantel(otu.wu.phylo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.tri.mantel = mantel(otu.uu.phylo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_jacc_tri <- paste("r =", round(jacc.tri.mantel$statistic,3), ";","p =",round(jacc.tri.mantel$signif,3))
phylo_bc_tri <- paste("r =", round(bc.tri.mantel$statistic,3), ";","p =",round(bc.tri.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.tri.mantel$statistic,3), ";","p =",round(wu.tri.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.tri.mantel$statistic,3), ";","p =",round(uu.tri.mantel$signif,3))
phylo_tri <- cbind("Jaccard" = phylo_jacc_tri, "Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Acanthuridae ----
sp_acanth <- names(table(samp_tri[samp_tri$family == "Acanthuridae",]$tax2))
core.phylosymb.acanth.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_acanth)
host_tree_acanth <- keep.tip(host_tree, sp_acanth)
mat_cophenetic_acanth <- cophenetic.phylo(host_tree_acanth)/2

# Calcul of the maacanthx with species merged 
otu.jacc.phylo.acanth <- as.matrix(vegdist(core.phylosymb.acanth.rel@otu_table, method = "jaccard"))
otu.bc.phylo.acanth <- as.matrix(vegdist(core.phylosymb.acanth.rel@otu_table, method = "bray"))
otu.wu.phylo.acanth <- as.matrix(UniFrac(core.phylosymb.acanth.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.acanth <- as.matrix(UniFrac(core.phylosymb.acanth.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.acanth,otu.bc.phylo.acanth,otu.wu.phylo.acanth,otu.uu.phylo.acanth , file = "matrices/beta_matrices.acanth.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_acanth) ; dim(otu.bc.phylo.acanth)
mat_cophenetic_acanth <- mat_cophenetic_acanth[match(rownames(otu.jacc.phylo.acanth), rownames(mat_cophenetic_acanth)),
                                         match(colnames(otu.jacc.phylo.acanth), colnames(mat_cophenetic_acanth))]
save(mat_cophenetic_acanth, file = "matrices/mat_cophenetic_acanth.RData")
jacc.acanth.mantel = mantel(otu.jacc.phylo.acanth, mat_cophenetic_acanth , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.acanth.mantel = mantel(otu.bc.phylo.acanth, mat_cophenetic_acanth , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.acanth.mantel = mantel(otu.wu.phylo.acanth, mat_cophenetic_acanth , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.acanth.mantel = mantel(otu.uu.phylo.acanth, mat_cophenetic_acanth , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_acanth <- paste("r =", round(jacc.acanth.mantel$statistic,3), ";","p =",round(jacc.acanth.mantel$signif,3))
phylo_bc_acanth <- paste("r =", round(bc.acanth.mantel$statistic,3), ";","p =",round(bc.acanth.mantel$signif, 3))
phylo_wu_acanth <- paste("r =", round(wu.acanth.mantel$statistic,3), ";","p =",round(wu.acanth.mantel$signif, 3))
phylo_uu_acanth <- paste("r =", round(uu.acanth.mantel$statistic,3), ";","p =",round(uu.acanth.mantel$signif, 3))
phylo_acanth <- cbind("Jaccard" = phylo_jacc_acanth, "Bray-Curtis" =phylo_bc_acanth,"Weighted Unifrac" =phylo_wu_acanth,"Unweighted Unifrac" =phylo_uu_acanth)


# __ ** Holocentridae ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_holo <- names(table(samp_tri[samp_tri$family == "Holocentridae",]$tax2))
core.phylosymb.holo.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_holo)
host_tree_holo <- keep.tip(host_tree, sp_holo)
mat_cophenetic_holo <- cophenetic.phylo(host_tree_holo)/2
# Calcul of the maholox with species merged 
otu.jacc.phylo.holo <- as.matrix(vegdist(core.phylosymb.holo.rel@otu_table, method = "jaccard"))
otu.bc.phylo.holo <- as.matrix(vegdist(core.phylosymb.holo.rel@otu_table, method = "bray"))
otu.wu.phylo.holo <- as.matrix(UniFrac(core.phylosymb.holo.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.holo <- as.matrix(UniFrac(core.phylosymb.holo.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.holo,otu.bc.phylo.holo,otu.wu.phylo.holo,otu.uu.phylo.holo , file = "matrices/beta_matrices.holo.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_holo) ; dim(otu.bc.phylo.holo)
mat_cophenetic_holo <- mat_cophenetic_holo[match(rownames(otu.jacc.phylo.holo), rownames(mat_cophenetic_holo)),
                                               match(colnames(otu.jacc.phylo.holo), colnames(mat_cophenetic_holo))]
save(mat_cophenetic_holo, file = "matrices/mat_cophenetic_holo.RData")
jacc.holo.mantel = mantel(otu.jacc.phylo.holo, mat_cophenetic_holo , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.holo.mantel = mantel(otu.bc.phylo.holo, mat_cophenetic_holo , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.holo.mantel = mantel(otu.wu.phylo.holo, mat_cophenetic_holo , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.holo.mantel = mantel(otu.uu.phylo.holo, mat_cophenetic_holo , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_holo <- paste("r =", round(jacc.holo.mantel$statistic,3), ";","p =",round(jacc.holo.mantel$signif, 3))
phylo_bc_holo <- paste("r =", round(bc.holo.mantel$statistic, 3), ";","p =",round(bc.holo.mantel$signif, 3))
phylo_wu_holo <- paste("r =", round(wu.holo.mantel$statistic,3), ";","p =",round(wu.holo.mantel$signif,3))
phylo_uu_holo <- paste("r =", round(uu.holo.mantel$statistic,3), ";","p =",round(uu.holo.mantel$signif,3))
phylo_holo <- cbind("Jaccard" = phylo_jacc_holo, "Bray-Curtis" =phylo_bc_holo,"Weighted Unifrac" =phylo_wu_holo,"Unweighted Unifrac" =phylo_uu_holo)

# __ ** Labridae ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_lab <- names(table(samp_tri[samp_tri$family == "Labridae",]$tax2))
core.phylosymb.lab.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_lab)
host_tree_lab <- keep.tip(host_tree, sp_lab)
mat_cophenetic_lab <- cophenetic.phylo(host_tree_lab)/2
# Calcul of the malabx with species merged 
otu.jacc.phylo.lab <- as.matrix(vegdist(core.phylosymb.lab.rel@otu_table, method = "jaccard"))
otu.bc.phylo.lab <- as.matrix(vegdist(core.phylosymb.lab.rel@otu_table, method = "bray"))
otu.wu.phylo.lab <- as.matrix(UniFrac(core.phylosymb.lab.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.lab <- as.matrix(UniFrac(core.phylosymb.lab.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.lab,otu.bc.phylo.lab,otu.wu.phylo.lab,otu.uu.phylo.lab , file = "matrices/beta_matrices.lab.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_lab) ; dim(otu.bc.phylo.lab)
mat_cophenetic_lab <- mat_cophenetic_lab[match(rownames(otu.jacc.phylo.lab), rownames(mat_cophenetic_lab)),
                                               match(colnames(otu.jacc.phylo.lab), colnames(mat_cophenetic_lab))]
save(mat_cophenetic_lab, file = "matrices/mat_cophenetic_lab.RData")
jacc.lab.mantel = mantel(otu.jacc.phylo.lab, mat_cophenetic_lab , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.lab.mantel = mantel(otu.bc.phylo.lab, mat_cophenetic_lab , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.lab.mantel = mantel(otu.wu.phylo.lab, mat_cophenetic_lab , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.lab.mantel = mantel(otu.uu.phylo.lab, mat_cophenetic_lab , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_lab <- paste("r =", round(jacc.lab.mantel$statistic,3), ";","p =",round(jacc.lab.mantel$signif,3))
phylo_bc_lab <- paste("r =", round(bc.lab.mantel$statistic,3), ";","p =", round(bc.lab.mantel$signif , 3))
phylo_wu_lab <- paste("r =", round(wu.lab.mantel$statistic,3), ";","p =",round(wu.lab.mantel$signif,3))
phylo_uu_lab <- paste("r =", round(uu.lab.mantel$statistic,3), ";","p =",round(uu.lab.mantel$signif,3))
phylo_lab <- cbind("Jaccard" = phylo_jacc_lab, "Bray-Curtis" =phylo_bc_lab,"Weighted Unifrac" =phylo_wu_lab,"Unweighted Unifrac" =phylo_uu_lab)

# __ ** Lethrinidae ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_leth <- names(table(samp_tri[samp_tri$family == "Lethrinidae",]$tax2))
core.phylosymb.leth.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_leth)
host_tree_leth <- keep.tip(host_tree, sp_leth)
mat_cophenetic_leth <- cophenetic.phylo(host_tree_leth)/2
# Calcul of the malethx with species merged 
otu.jacc.phylo.leth <- as.matrix(vegdist(core.phylosymb.leth.rel@otu_table, method = "jaccard"))
otu.bc.phylo.leth <- as.matrix(vegdist(core.phylosymb.leth.rel@otu_table, method = "bray"))
otu.wu.phylo.leth <- as.matrix(UniFrac(core.phylosymb.leth.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.leth <- as.matrix(UniFrac(core.phylosymb.leth.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.leth,otu.bc.phylo.leth,otu.wu.phylo.leth,otu.uu.phylo.leth , file = "matrices/beta_matrices.leth.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_leth) ; dim(otu.bc.phylo.leth)
mat_cophenetic_leth <- mat_cophenetic_leth[match(rownames(otu.jacc.phylo.leth), rownames(mat_cophenetic_leth)),
                                           match(colnames(otu.jacc.phylo.leth), colnames(mat_cophenetic_leth))]
save(mat_cophenetic_leth, file = "matrices/mat_cophenetic_leth.RData")

jacc.leth.mantel = mantel(otu.jacc.phylo.leth, mat_cophenetic_leth , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.leth.mantel = mantel(otu.bc.phylo.leth, mat_cophenetic_leth , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.leth.mantel = mantel(otu.wu.phylo.leth, mat_cophenetic_leth , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.leth.mantel = mantel(otu.uu.phylo.leth, mat_cophenetic_leth , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_leth <- paste("r =", round(jacc.leth.mantel$statistic,3), ";","p =",round(jacc.leth.mantel$signif,3))
phylo_bc_leth <- paste("r =", round(bc.leth.mantel$statistic,3), ";","p =",round(bc.leth.mantel$signif,3))
phylo_wu_leth <- paste("r =", round(wu.leth.mantel$statistic,3), ";","p =", round(wu.leth.mantel$signif,3))
phylo_uu_leth <- paste("r =", round(uu.leth.mantel$statistic,3), ";","p =", round(uu.leth.mantel$signif,3))
phylo_leth <- cbind("Jaccard" = phylo_jacc_leth, "Bray-Curtis" =phylo_bc_leth,"Weighted Unifrac" =phylo_wu_leth,"Unweighted Unifrac" =phylo_uu_leth)

# __ ** Lutjanidae ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_lutj <- names(table(samp_tri[samp_tri$family == "Lutjanidae",]$tax2))
core.phylosymb.lutj.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_lutj)
host_tree_lutj <- keep.tip(host_tree, sp_lutj)
mat_cophenetic_lutj <- cophenetic.phylo(host_tree_lutj)/2
# Calcul of the malutjx with species merged 
otu.jacc.phylo.lutj <- as.matrix(vegdist(core.phylosymb.lutj.rel@otu_table, method = "jaccard"))
otu.bc.phylo.lutj <- as.matrix(vegdist(core.phylosymb.lutj.rel@otu_table, method = "bray"))
otu.wu.phylo.lutj <- as.matrix(UniFrac(core.phylosymb.lutj.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.lutj <- as.matrix(UniFrac(core.phylosymb.lutj.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.lutj,otu.bc.phylo.lutj,otu.wu.phylo.lutj,otu.uu.phylo.lutj , file = "matrices/beta_matrices.lutj.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_lutj) ; dim(otu.bc.phylo.lutj)
mat_cophenetic_lutj <- mat_cophenetic_lutj[match(rownames(otu.jacc.phylo.lutj), rownames(mat_cophenetic_lutj)),
                                           match(colnames(otu.jacc.phylo.lutj), colnames(mat_cophenetic_lutj))]
save(mat_cophenetic_lutj, file = "matrices/mat_cophenetic_lutj.RData")
jacc.lutj.mantel = mantel(otu.jacc.phylo.lutj, mat_cophenetic_lutj , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.lutj.mantel = mantel(otu.bc.phylo.lutj, mat_cophenetic_lutj , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.lutj.mantel = mantel(otu.wu.phylo.lutj, mat_cophenetic_lutj , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.lutj.mantel = mantel(otu.uu.phylo.lutj, mat_cophenetic_lutj , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_lutj <- paste("r =", round(jacc.lutj.mantel$statistic,3), ";","p =", round(jacc.lutj.mantel$signif,3))
phylo_bc_lutj <- paste("r =", round(bc.lutj.mantel$statistic,3), ";","p =", round(bc.lutj.mantel$signif,3))
phylo_wu_lutj <- paste("r =", round(wu.lutj.mantel$statistic,3), ";","p =", round(wu.lutj.mantel$signif,3))
phylo_uu_lutj <- paste("r =", round(uu.lutj.mantel$statistic,3), ";","p =", round(uu.lutj.mantel$signif, 3))
phylo_lutj <- cbind("Jaccard" = phylo_jacc_lutj, "Bray-Curtis" =phylo_bc_lutj,"Weighted Unifrac" =phylo_wu_lutj,"Unweighted Unifrac" =phylo_uu_lutj)

# __ ** Pomacentridae ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_poma <- names(table(samp_tri[samp_tri$family == "Pomacentridae",]$tax2))
core.phylosymb.poma.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_poma)
host_tree_poma <- keep.tip(host_tree, sp_poma)
mat_cophenetic_poma <- cophenetic.phylo(host_tree_poma)/2
# Calcul of the mapomax with species merged 
otu.jacc.phylo.poma <- as.matrix(vegdist(core.phylosymb.poma.rel@otu_table, method = "jaccard"))
otu.bc.phylo.poma <- as.matrix(vegdist(core.phylosymb.poma.rel@otu_table, method = "bray"))
otu.wu.phylo.poma <- as.matrix(UniFrac(core.phylosymb.poma.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.poma <- as.matrix(UniFrac(core.phylosymb.poma.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.poma,otu.bc.phylo.poma,otu.wu.phylo.poma,otu.uu.phylo.poma , file = "matrices/beta_matrices.poma.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_poma) ; dim(otu.bc.phylo.poma)
mat_cophenetic_poma <- mat_cophenetic_poma[match(rownames(otu.jacc.phylo.poma), rownames(mat_cophenetic_poma)),
                                           match(colnames(otu.jacc.phylo.poma), colnames(mat_cophenetic_poma))]
save(mat_cophenetic_poma, file = "matrices/mat_cophenetic_poma.RData")
jacc.poma.mantel = mantel(otu.jacc.phylo.poma, mat_cophenetic_poma , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.poma.mantel = mantel(otu.bc.phylo.poma, mat_cophenetic_poma , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.poma.mantel = mantel(otu.wu.phylo.poma, mat_cophenetic_poma , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.poma.mantel = mantel(otu.uu.phylo.poma, mat_cophenetic_poma , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_poma <- paste("r =", round(jacc.poma.mantel$statistic,3), ";","p =", round(jacc.poma.mantel$signif,3))
phylo_bc_poma <- paste("r =", round(bc.poma.mantel$statistic, 3),";","p =", round(bc.poma.mantel$signif,3))
phylo_wu_poma <- paste("r =", round(wu.poma.mantel$statistic,3), ";","p =", round(wu.poma.mantel$signif,3))
phylo_uu_poma <- paste("r =", round(uu.poma.mantel$statistic,3), ";","p =",round(uu.poma.mantel$signif,3))
phylo_poma <- cbind("Jaccard" = phylo_jacc_poma, "Bray-Curtis" =phylo_bc_poma,"Weighted Unifrac" =phylo_wu_poma,"Unweighted Unifrac" =phylo_uu_poma)

# __ ** Scaridae ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_scar <- names(table(samp_tri[samp_tri$family == "Scaridae",]$tax2))
core.phylosymb.scar.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_scar)
host_tree_scar <- keep.tip(host_tree, sp_scar)
mat_cophenetic_scar <- cophenetic.phylo(host_tree_scar)/2
# Calcul of the mascarx with species merged 
otu.jacc.phylo.scar <- as.matrix(vegdist(core.phylosymb.scar.rel@otu_table, method = "jaccard"))
otu.bc.phylo.scar <- as.matrix(vegdist(core.phylosymb.scar.rel@otu_table, method = "bray"))
otu.wu.phylo.scar <- as.matrix(UniFrac(core.phylosymb.scar.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.scar <- as.matrix(UniFrac(core.phylosymb.scar.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.scar,otu.bc.phylo.scar,otu.wu.phylo.scar,otu.uu.phylo.scar , file = "matrices/beta_matrices.scar.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_scar) ; dim(otu.bc.phylo.scar)
mat_cophenetic_scar <- mat_cophenetic_scar[match(rownames(otu.jacc.phylo.scar), rownames(mat_cophenetic_scar)),
                                           match(colnames(otu.jacc.phylo.scar), colnames(mat_cophenetic_scar))]
save(mat_cophenetic_scar, file = "matrices/mat_cophenetic_scar.RData")
jacc.scar.mantel = mantel(otu.jacc.phylo.scar, mat_cophenetic_scar , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.scar.mantel = mantel(otu.bc.phylo.scar, mat_cophenetic_scar , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.scar.mantel = mantel(otu.wu.phylo.scar, mat_cophenetic_scar , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.scar.mantel = mantel(otu.uu.phylo.scar, mat_cophenetic_scar , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_scar <- paste("r =", round(jacc.scar.mantel$statistic,3), ";","p =",round(jacc.scar.mantel$signif,3))
phylo_bc_scar <- paste("r =", round(bc.scar.mantel$statistic,3), ";","p =", round(bc.scar.mantel$signif,3))
phylo_wu_scar <- paste("r =", round(wu.scar.mantel$statistic,3), ";","p =",round(wu.scar.mantel$signif,3))
phylo_uu_scar <- paste("r =", round(uu.scar.mantel$statistic,3), ";","p =", round(uu.scar.mantel$signif,3))
phylo_scar <- cbind("Jaccard" = phylo_jacc_scar, "Bray-Curtis" =phylo_bc_scar,"Weighted Unifrac" =phylo_wu_scar,"Unweighted Unifrac" =phylo_uu_scar)

# __ ** Herbivores ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
core.phylosymb.herb.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_herb)
host_tree_herb <- keep.tip(host_tree, sp_herb)
mat_cophenetic_herb <- cophenetic.phylo(host_tree_herb)/2
# Calcul of the maherbx with species merged 
otu.jacc.phylo.herb <- as.matrix(vegdist(core.phylosymb.herb.rel@otu_table, method = "jaccard"))
otu.bc.phylo.herb <- as.matrix(vegdist(core.phylosymb.herb.rel@otu_table, method = "bray"))
otu.wu.phylo.herb <- as.matrix(UniFrac(core.phylosymb.herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.herb <- as.matrix(UniFrac(core.phylosymb.herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.herb,otu.bc.phylo.herb,otu.wu.phylo.herb,otu.uu.phylo.herb , file = "matrices/beta_matrices.herb.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_herb) ; dim(otu.bc.phylo.herb)
mat_cophenetic_herb <- mat_cophenetic_herb[match(rownames(otu.jacc.phylo.herb), rownames(mat_cophenetic_herb)),
                                           match(colnames(otu.jacc.phylo.herb), colnames(mat_cophenetic_herb))]
save(mat_cophenetic_herb, file = "matrices/mat_cophenetic_herb.RData")
jacc.herb.mantel = mantel(otu.jacc.phylo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.herb.mantel = mantel(otu.bc.phylo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.herb.mantel = mantel(otu.wu.phylo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.herb.mantel = mantel(otu.uu.phylo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_herb <- paste("r =", round(jacc.herb.mantel$statistic,3), ";","p =",round(jacc.herb.mantel$signif,3))
phylo_bc_herb <- paste("r =", round(bc.herb.mantel$statistic,3), ";","p =", round(bc.herb.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.herb.mantel$statistic,3), ";","p =", round(wu.herb.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.herb.mantel$statistic,3), ";","p =", round(uu.herb.mantel$signif,3))
phylo_herb <- cbind("Jaccard" = phylo_jacc_herb, "Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Detritivores ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_detri <- names(table(samp_tri[samp_tri$diet2 == "Detritivorous",]$tax2))
core.phylosymb.detri.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_detri)
host_tree_detri <- keep.tip(host_tree, sp_detri)
mat_cophenetic_detri <- cophenetic.phylo(host_tree_detri)/2
# Calcul of the madetrix with species merged 
otu.jacc.phylo.detri <- as.matrix(vegdist(core.phylosymb.detri.rel@otu_table, method = "jaccard"))
otu.bc.phylo.detri <- as.matrix(vegdist(core.phylosymb.detri.rel@otu_table, method = "bray"))
otu.wu.phylo.detri <- as.matrix(UniFrac(core.phylosymb.detri.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.detri <- as.matrix(UniFrac(core.phylosymb.detri.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.detri,otu.bc.phylo.detri,otu.wu.phylo.detri,otu.uu.phylo.detri , file = "matrices/beta_matrices.detri.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_detri) ; dim(otu.bc.phylo.detri)
mat_cophenetic_detri <- mat_cophenetic_detri[match(rownames(otu.jacc.phylo.detri), rownames(mat_cophenetic_detri)),
                                           match(colnames(otu.jacc.phylo.detri), colnames(mat_cophenetic_detri))]
save(mat_cophenetic_detri, file = "matrices/mat_cophenetic_detri.RData")
jacc.detri.mantel = mantel(otu.jacc.phylo.detri, mat_cophenetic_detri , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.detri.mantel = mantel(otu.bc.phylo.detri, mat_cophenetic_detri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.detri.mantel = mantel(otu.wu.phylo.detri, mat_cophenetic_detri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.detri.mantel = mantel(otu.uu.phylo.detri, mat_cophenetic_detri , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_detri <- paste("r =", round(jacc.detri.mantel$statistic,3), ";","p =",round(jacc.detri.mantel$signif,3))
phylo_bc_detri <- paste("r =", round(bc.detri.mantel$statistic,3), ";","p =", round(bc.detri.mantel$signif,3))
phylo_wu_detri <- paste("r =", round(wu.detri.mantel$statistic,3), ";","p =", round(wu.detri.mantel$signif,3))
phylo_uu_detri <- paste("r =", round(uu.detri.mantel$statistic,3), ";","p =", round(uu.detri.mantel$signif,3))
phylo_detri <- cbind("Jaccard" = phylo_jacc_detri, "Bray-Curtis" =phylo_bc_detri,"Weighted Unifrac" =phylo_wu_detri,"Unweighted Unifrac" =phylo_uu_detri)

# __ ** Strict herbivores ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_str_herb <- names(table(samp_tri[samp_tri$diet2 == "Herbivorous",]$tax2))
core.phylosymb.str_herb.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_str_herb)
host_tree_str_herb <- keep.tip(host_tree, sp_str_herb)
mat_cophenetic_str_herb <- cophenetic.phylo(host_tree_str_herb)/2
# Calcul of the mastr_herbx with species merged 
otu.jacc.phylo.str_herb <- as.matrix(vegdist(core.phylosymb.str_herb.rel@otu_table, method = "jaccard"))
otu.bc.phylo.str_herb <- as.matrix(vegdist(core.phylosymb.str_herb.rel@otu_table, method = "bray"))
otu.wu.phylo.str_herb <- as.matrix(UniFrac(core.phylosymb.str_herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.str_herb <- as.matrix(UniFrac(core.phylosymb.str_herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.str_herb,otu.bc.phylo.str_herb,otu.wu.phylo.str_herb,otu.uu.phylo.str_herb , file = "matrices/beta_matrices.str_herb.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_str_herb) ; dim(otu.bc.phylo.str_herb)
mat_cophenetic_str_herb <- mat_cophenetic_str_herb[match(rownames(otu.jacc.phylo.str_herb), rownames(mat_cophenetic_str_herb)),
                                             match(colnames(otu.jacc.phylo.str_herb), colnames(mat_cophenetic_str_herb))]
save(mat_cophenetic_str_herb, file = "matrices/mat_cophenetic_str_herb.RData")
jacc.str_herb.mantel = mantel(otu.jacc.phylo.str_herb, mat_cophenetic_str_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.str_herb.mantel = mantel(otu.bc.phylo.str_herb, mat_cophenetic_str_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.str_herb.mantel = mantel(otu.wu.phylo.str_herb, mat_cophenetic_str_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.str_herb.mantel = mantel(otu.uu.phylo.str_herb, mat_cophenetic_str_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_str_herb <- paste("r =", round(jacc.str_herb.mantel$statistic,3), ";","p =",round(jacc.str_herb.mantel$signif,3))
phylo_bc_str_herb <- paste("r =", round(bc.str_herb.mantel$statistic,3), ";","p =", round(bc.str_herb.mantel$signif,3))
phylo_wu_str_herb <- paste("r =", round(wu.str_herb.mantel$statistic,3), ";","p =", round(wu.str_herb.mantel$signif,3))
phylo_uu_str_herb <- paste("r =", round(uu.str_herb.mantel$statistic,3), ";","p =", round(uu.str_herb.mantel$signif,3))
phylo_str_herb <- cbind("Jaccard" = phylo_jacc_str_herb, "Bray-Curtis" =phylo_bc_str_herb,"Weighted Unifrac" =phylo_wu_str_herb,"Unweighted Unifrac" =phylo_uu_str_herb)

# __ ** Carnivores ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
core.phylosymb.carn.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_carn)
host_tree_carn <- keep.tip(host_tree, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(host_tree_carn)/2
# Calcul of the macarnx with species merged 
otu.jacc.phylo.carn <- as.matrix(vegdist(core.phylosymb.carn.rel@otu_table, method = "jaccard"))
otu.bc.phylo.carn <- as.matrix(vegdist(core.phylosymb.carn.rel@otu_table, method = "bray"))
otu.wu.phylo.carn <- as.matrix(UniFrac(core.phylosymb.carn.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.carn <- as.matrix(UniFrac(core.phylosymb.carn.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.carn,otu.bc.phylo.carn,otu.wu.phylo.carn,otu.uu.phylo.carn , file = "matrices/beta_matrices.carn.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_carn) ; dim(otu.bc.phylo.carn)
mat_cophenetic_carn <- mat_cophenetic_carn[match(rownames(otu.jacc.phylo.carn), rownames(mat_cophenetic_carn)),
                                           match(colnames(otu.jacc.phylo.carn), colnames(mat_cophenetic_carn))]
save(mat_cophenetic_carn, file = "matrices/mat_cophenetic_carn.RData")
jacc.carn.mantel = mantel(otu.jacc.phylo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.carn.mantel = mantel(otu.bc.phylo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.carn.mantel = mantel(otu.wu.phylo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.carn.mantel = mantel(otu.uu.phylo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_carn <- paste("r =", round(jacc.carn.mantel$statistic,3), ";","p =",round(jacc.carn.mantel$signif,3))
phylo_bc_carn <- paste("r =", round(bc.carn.mantel$statistic,3), ";","p =", round(bc.carn.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.carn.mantel$statistic,3), ";","p =", round(wu.carn.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.carn.mantel$statistic,3), ";","p =", round(uu.carn.mantel$signif,3))
phylo_carn <- cbind("Jaccard" = phylo_jacc_carn, "Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Invertivore mobile ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_mi <- names(table(samp_tri[samp_tri$diet2 == "Mobile invertebrate",]$tax2))
core.phylosymb.mi.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_mi)
host_tree_mi <- keep.tip(host_tree, sp_mi)
mat_cophenetic_mi <- cophenetic.phylo(host_tree_mi)/2
# Calcul of the mamix with species merged 
otu.jacc.phylo.mi <- as.matrix(vegdist(core.phylosymb.mi.rel@otu_table, method = "jaccard"))
otu.bc.phylo.mi <- as.matrix(vegdist(core.phylosymb.mi.rel@otu_table, method = "bray"))
otu.wu.phylo.mi <- as.matrix(UniFrac(core.phylosymb.mi.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.mi <- as.matrix(UniFrac(core.phylosymb.mi.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.mi,otu.bc.phylo.mi,otu.wu.phylo.mi,otu.uu.phylo.mi , file = "matrices/beta_matrices.mi.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_mi) ; dim(otu.bc.phylo.mi)
mat_cophenetic_mi <- mat_cophenetic_mi[match(rownames(otu.jacc.phylo.mi), rownames(mat_cophenetic_mi)),
                                                   match(colnames(otu.jacc.phylo.mi), colnames(mat_cophenetic_mi))]
save(mat_cophenetic_mi, file = "matrices/mat_cophenetic_mi.RData")
jacc.mi.mantel = mantel(otu.jacc.phylo.mi, mat_cophenetic_mi , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.mi.mantel = mantel(otu.bc.phylo.mi, mat_cophenetic_mi , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.mi.mantel = mantel(otu.wu.phylo.mi, mat_cophenetic_mi , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.mi.mantel = mantel(otu.uu.phylo.mi, mat_cophenetic_mi , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_mi <- paste("r =", round(jacc.mi.mantel$statistic,3), ";","p =",round(jacc.mi.mantel$signif,3))
phylo_bc_mi <- paste("r =", round(bc.mi.mantel$statistic,3), ";","p =", round(bc.mi.mantel$signif,3))
phylo_wu_mi <- paste("r =", round(wu.mi.mantel$statistic,3), ";","p =", round(wu.mi.mantel$signif,3))
phylo_uu_mi <- paste("r =", round(uu.mi.mantel$statistic,3), ";","p =", round(uu.mi.mantel$signif,3))
phylo_mi <- cbind("Jaccard" = phylo_jacc_mi, "Bray-Curtis" =phylo_bc_mi,"Weighted Unifrac" =phylo_wu_mi,"Unweighted Unifrac" =phylo_uu_mi)

# __ ** Invertivore sessile ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_si <- names(table(samp_tri[samp_tri$diet2 == "Sessile invertebrate",]$tax2))
core.phylosymb.si.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_si)
host_tree_si <- keep.tip(host_tree, sp_si)
mat_cophenetic_si <- cophenetic.phylo(host_tree_si)/2
# Calcul of the matrix with species merged 
otu.jacc.phylo.si <- as.matrix(vegdist(core.phylosymb.si.rel@otu_table, method = "jaccard"))
otu.bc.phylo.si <- as.matrix(vegdist(core.phylosymb.si.rel@otu_table, method = "bray"))
otu.wu.phylo.si <- as.matrix(UniFrac(core.phylosymb.si.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.si <- as.matrix(UniFrac(core.phylosymb.si.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.si,otu.bc.phylo.si,otu.wu.phylo.si,otu.uu.phylo.si , file = "matrices/beta_matrices.si.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_si) ; dim(otu.bc.phylo.si)
mat_cophenetic_si <- mat_cophenetic_si[match(rownames(otu.jacc.phylo.si), rownames(mat_cophenetic_si)),
                                       match(colnames(otu.jacc.phylo.si), colnames(mat_cophenetic_si))]
save(mat_cophenetic_si, file = "matrices/mat_cophenetic_si.RData")
jacc.si.mantel = mantel(otu.jacc.phylo.si, mat_cophenetic_si , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.si.mantel = mantel(otu.bc.phylo.si, mat_cophenetic_si , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.si.mantel = mantel(otu.wu.phylo.si, mat_cophenetic_si , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.si.mantel = mantel(otu.uu.phylo.si, mat_cophenetic_si , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_si <- paste("r =", round(jacc.si.mantel$statistic,3), ";","p =",round(jacc.si.mantel$signif,3))
phylo_bc_si <- paste("r =", round(bc.si.mantel$statistic,3), ";","p =", round(bc.si.mantel$signif,3))
phylo_wu_si <- paste("r =", round(wu.si.mantel$statistic,3), ";","p =", round(wu.si.mantel$signif,3))
phylo_uu_si <- paste("r =", round(uu.si.mantel$statistic,3), ";","p =", round(uu.si.mantel$signif,3))
phylo_si <- cbind("Jaccard" = phylo_jacc_si, "Bray-Curtis" =phylo_bc_si,"Weighted Unifrac" =phylo_wu_si,"Unweighted Unifrac" =phylo_uu_si)

# __ ** Piscivorous ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_fc <- names(table(samp_tri[samp_tri$diet2 == "Piscivorous",]$tax2))
core.phylosymb.fc.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_fc)
host_tree_fc <- keep.tip(host_tree, sp_fc)
mat_cophenetic_fc <- cophenetic.phylo(host_tree_fc)/2
# Calcul of the matrix with species merged 
otu.jacc.phylo.fc <- as.matrix(vegdist(core.phylosymb.fc.rel@otu_table, method = "jaccard"))
otu.bc.phylo.fc <- as.matrix(vegdist(core.phylosymb.fc.rel@otu_table, method = "bray"))
otu.wu.phylo.fc <- as.matrix(UniFrac(core.phylosymb.fc.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.fc <- as.matrix(UniFrac(core.phylosymb.fc.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.fc,otu.bc.phylo.fc,otu.wu.phylo.fc,otu.uu.phylo.fc , file = "matrices/beta_matrices.fc.RData")

# Mantel test on divrgence time 
dim(mat_cophenetic_fc) ; dim(otu.bc.phylo.fc)
mat_cophenetic_fc <- mat_cophenetic_fc[match(rownames(otu.jacc.phylo.fc), rownames(mat_cophenetic_fc)),
                                       match(colnames(otu.jacc.phylo.fc), colnames(mat_cophenetic_fc))]
save(mat_cophenetic_fc, file = "matrices/mat_cophenetic_fc.RData")
jacc.fc.mantel = mantel(otu.jacc.phylo.fc, mat_cophenetic_fc , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.fc.mantel = mantel(otu.bc.phylo.fc, mat_cophenetic_fc , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.fc.mantel = mantel(otu.wu.phylo.fc, mat_cophenetic_fc , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.fc.mantel = mantel(otu.uu.phylo.fc, mat_cophenetic_fc , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_fc <- paste("r =", round(jacc.fc.mantel$statistic,3), ";","p =",round(jacc.fc.mantel$signif,3))
phylo_bc_fc <- paste("r =", round(bc.fc.mantel$statistic,3), ";","p =", round(bc.fc.mantel$signif,3))
phylo_wu_fc <- paste("r =", round(wu.fc.mantel$statistic,3), ";","p =", round(wu.fc.mantel$signif,3))
phylo_uu_fc <- paste("r =", round(uu.fc.mantel$statistic,3), ";","p =", round(uu.fc.mantel$signif,3))
phylo_fc <- cbind("Jaccard" = phylo_jacc_fc, "Bray-Curtis" =phylo_bc_fc,"Weighted Unifrac" =phylo_wu_fc,"Unweighted Unifrac" =phylo_uu_fc)


# __ ** Omnivores ----
#host_tree <- read.tree(file = "fish_tree_modif.tre")
#load("core.phylosymb.tri.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
core.phylosymb.omni.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_omni)
host_tree_omni <- keep.tip(host_tree, sp_omni)
mat_cophenetic_omni <- cophenetic.phylo(host_tree_omni)/2
# Calcul of the maomnix with species merged 
otu.jacc.phylo.omni <- as.matrix(vegdist(core.phylosymb.omni.rel@otu_table, method = "jaccard"))
otu.bc.phylo.omni <- as.matrix(vegdist(core.phylosymb.omni.rel@otu_table, method = "bray"))
otu.wu.phylo.omni <- as.matrix(UniFrac(core.phylosymb.omni.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.omni <- as.matrix(UniFrac(core.phylosymb.omni.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.omni,otu.bc.phylo.omni,otu.wu.phylo.omni,otu.uu.phylo.omni , file = "matrices/beta_matrices.omni.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_omni) ; dim(otu.bc.phylo.omni)
mat_cophenetic_omni <- mat_cophenetic_omni[match(rownames(otu.jacc.phylo.omni), rownames(mat_cophenetic_omni)),
                                           match(colnames(otu.jacc.phylo.omni), colnames(mat_cophenetic_omni))]
save(mat_cophenetic_omni, file = "matrices/mat_cophenetic_omni.RData")
jacc.omni.mantel = mantel(otu.jacc.phylo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.omni.mantel = mantel(otu.bc.phylo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.omni.mantel = mantel(otu.wu.phylo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.omni.mantel = mantel(otu.uu.phylo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_omni <- paste("r =", round(jacc.omni.mantel$statistic,3), ";","p =", round(jacc.omni.mantel$signif,3))
phylo_bc_omni <- paste("r =", round(bc.omni.mantel$statistic,3), ";","p =",round(bc.omni.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.omni.mantel$statistic,3), ";","p =",round(wu.omni.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.omni.mantel$statistic,3), ";","p =", round(uu.omni.mantel$signif,3))
phylo_omni <- cbind("Jaccard" = phylo_jacc_omni, "Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

# __ ** Perciformes ----
sp_perci <- names(table(samp_tri[samp_tri$order == "Perciformes",]$tax2))
core.phylosymb.perci.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_perci)
host_tree_perci <- keep.tip(host_tree, sp_perci)
mat_cophenetic_perci <- cophenetic.phylo(host_tree_perci)/2

# Calcul of the matrix with species merged 
otu.jacc.phylo.perci <- as.matrix(vegdist(core.phylosymb.perci.rel@otu_table, method = "jaccard"))
otu.bc.phylo.perci <- as.matrix(vegdist(core.phylosymb.perci.rel@otu_table, method = "bray"))
otu.wu.phylo.perci <- as.matrix(UniFrac(core.phylosymb.perci.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.perci <- as.matrix(UniFrac(core.phylosymb.perci.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.perci,otu.bc.phylo.perci,otu.wu.phylo.perci,otu.uu.phylo.perci , file = "matrices/beta_matrices.perci.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_perci) ; dim(otu.bc.phylo.perci)
mat_cophenetic_perci <- mat_cophenetic_perci[match(rownames(otu.jacc.phylo.perci), rownames(mat_cophenetic_perci)),
                                               match(colnames(otu.jacc.phylo.perci), colnames(mat_cophenetic_perci))]
save(mat_cophenetic_perci, file = "matrices/mat_cophenetic_perci.RData")
jacc.perci.mantel = mantel(otu.jacc.phylo.perci, mat_cophenetic_perci , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.perci.mantel = mantel(otu.bc.phylo.perci, mat_cophenetic_perci , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.perci.mantel = mantel(otu.wu.phylo.perci, mat_cophenetic_perci , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.perci.mantel = mantel(otu.uu.phylo.perci, mat_cophenetic_perci , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_perci <- paste("r =", round(jacc.perci.mantel$statistic,3), ";","p =",round(jacc.perci.mantel$signif,3))
phylo_bc_perci <- paste("r =", round(bc.perci.mantel$statistic,3), ";","p =",round(bc.perci.mantel$signif, 3))
phylo_wu_perci <- paste("r =", round(wu.perci.mantel$statistic,3), ";","p =",round(wu.perci.mantel$signif, 3))
phylo_uu_perci <- paste("r =", round(uu.perci.mantel$statistic,3), ";","p =",round(uu.perci.mantel$signif, 3))
phylo_perci <- cbind("Jaccard" = phylo_jacc_perci, "Bray-Curtis" =phylo_bc_perci,"Weighted Unifrac" =phylo_wu_perci,"Unweighted Unifrac" =phylo_uu_perci)

# __ ** Labridae/Scaridae ----
sp_labsca <- names(table(samp_tri[samp_tri$family %in% c("Labridae", "Scaridae"),]$tax2))
core.phylosymb.labsca.rel <- subset_samples(core.phylosymb.tri.rel, tax2 %in% sp_labsca)
host_tree_labsca <- keep.tip(host_tree, sp_labsca)
mat_cophenetic_labsca <- cophenetic.phylo(host_tree_labsca)/2

# Calcul of the matrix with species merged 
otu.jacc.phylo.labsca <- as.matrix(vegdist(core.phylosymb.labsca.rel@otu_table, method = "jaccard"))
otu.bc.phylo.labsca <- as.matrix(vegdist(core.phylosymb.labsca.rel@otu_table, method = "bray"))
otu.wu.phylo.labsca <- as.matrix(UniFrac(core.phylosymb.labsca.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.phylo.labsca <- as.matrix(UniFrac(core.phylosymb.labsca.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.phylo.labsca,otu.bc.phylo.labsca,otu.wu.phylo.labsca,otu.uu.phylo.labsca , file = "matrices/beta_matrices.labsca.RData")
# Mantel test on divrgence time 
dim(mat_cophenetic_labsca) ; dim(otu.bc.phylo.labsca)
mat_cophenetic_labsca <- mat_cophenetic_labsca[match(rownames(otu.jacc.phylo.labsca), rownames(mat_cophenetic_labsca)),
                                             match(colnames(otu.jacc.phylo.labsca), colnames(mat_cophenetic_labsca))]
save(mat_cophenetic_labsca, file = "matrices/mat_cophenetic_labsca.RData")
jacc.labsca.mantel = mantel(otu.jacc.phylo.labsca, mat_cophenetic_labsca , method = "spearman", permutations = 9999, na.rm = TRUE)
bc.labsca.mantel = mantel(otu.bc.phylo.labsca, mat_cophenetic_labsca , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.labsca.mantel = mantel(otu.wu.phylo.labsca, mat_cophenetic_labsca , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.labsca.mantel = mantel(otu.uu.phylo.labsca, mat_cophenetic_labsca , method = "spearman", permutations = 9999, na.rm = TRUE)
phylo_jacc_labsca <- paste("r =", round(jacc.labsca.mantel$statistic,3), ";","p =",round(jacc.labsca.mantel$signif,3))
phylo_bc_labsca <- paste("r =", round(bc.labsca.mantel$statistic,3), ";","p =",round(bc.labsca.mantel$signif, 3))
phylo_wu_labsca <- paste("r =", round(wu.labsca.mantel$statistic,3), ";","p =",round(wu.labsca.mantel$signif, 3))
phylo_uu_labsca <- paste("r =", round(uu.labsca.mantel$statistic,3), ";","p =",round(uu.labsca.mantel$signif, 3))
phylo_labsca <- cbind("Jaccard" = phylo_jacc_labsca, "Bray-Curtis" =phylo_bc_labsca,"Weighted Unifrac" =phylo_wu_labsca,"Unweighted Unifrac" =phylo_uu_labsca)


## Table : Divergence time ----
Mantel_phylo <-rbind(phylo_all, phylo_tri, phylo_acanth, phylo_holo, 
                     phylo_lab,phylo_leth,phylo_lutj, phylo_poma, phylo_scar,
                     phylo_herb,phylo_carn,phylo_omni)
rownames(Mantel_phylo) <- c("All","Triplicate","Acanthuridae","Holocentridae","Labridae","Lethrinidae", 
                            "Lutjanidae","Pomacentridae", "Scaridae","Herbivores","Carnivores","Omnivores")

write.table(Mantel_phylo, file = "Mantel_phylo_results.txt", sep = "\t", quote = F)

## Graphiques ----
## __ A) For triplicate host ----
load("matrices/beta_matrices.tri.RData")
vec.jacc.phylo = as.vector(as.matrix(otu.jacc.phylo.tri))
vec.bc.phylo = as.vector(as.matrix(otu.bc.phylo.tri))
vec.wu.phylo = as.vector(as.matrix(otu.wu.phylo.tri))
vec.uu.phylo = as.vector(as.matrix(otu.uu.phylo.tri))
vec.div = as.vector(mat_cophenetic_tri)

mat.phylo = as.data.frame(rbind(cbind(vec.jacc.phylo, rep("Jaccard", length(vec.jacc.phylo))),
                                cbind(vec.bc.phylo, rep("Bray-Curtis", length(vec.bc.phylo))),
                                cbind(vec.wu.phylo, rep("W-Unifrac", length(vec.wu.phylo))),
                                cbind(vec.uu.phylo, rep("U-Unifrac", length(vec.uu.phylo)))))
mat = cbind(mat.phylo, rep(vec.div,4))
colnames(mat) <- c("Dissimilarity" , "Index", "Divergence")
mat <- mat[mat$Dissimilarity !=0,]
mat$Dissimilarity <- as.numeric(as.character(mat$Dissimilarity))

Fig.A1 = ggplot(mat[mat$Index == "Bray-Curtis",], aes(y = Dissimilarity, x = Divergence)) + 
  geom_point(size = 4, alpha = 0.2, colour = "black",shape = 16) + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
  labs(x = "Host Divergence Time (My)", y = "Bacteriome Dissimilarity (Bray-Curtis)") + 
  theme( axis.text.x = element_text(colour = "black", size = 16,family = "serif"), 
         axis.text.y = element_text(size = 16, colour = "black",family = "serif"), 
         axis.title= element_text(size = 20, colour = "black", family = "serif"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 3, family = "serif"),
         legend.title = element_text(size = 8, family = "serif"))#,
#strip.text.x = element_text(size=14, family = "serif")) +
#facet_wrap(~ Index)
Fig.A1

pdf(file = "Mantel_divergence_BC_tri.pdf", he = 7 , wi = 7)
Fig.A
dev.off()

Fig.A2 = ggplot(mat[mat$Index != "Bray-Curtis",], aes(y = Dissimilarity, x = Divergence)) + 
  geom_point(size = 4, alpha = 0.2, colour = "black",shape = 16) + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
  labs(x = "Host Divergence Time (My)", y = "Bacteriome Dissimilarity") + 
  theme( axis.text.x = element_text(colour = "black", size = 16,family = "serif"), 
         axis.text.y = element_text(size = 16, colour = "black",family = "serif"), 
         axis.title= element_text(size = 20, colour = "black", family = "serif"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 3, family = "serif"),
         legend.title = element_text(size = 8, family = "serif"),
         strip.text.x = element_text(size=14, family = "serif")) +
  facet_wrap(~ Index , nrow = 1, ncol = 3)
Fig.A2

pdf(file = "Supp_Mantel_divergence_all.pdf", he = 7 , wi = 12)
Fig.A2
dev.off()

## __ B) Diets BC ----
load("matrices/beta_matrices.herb.RData")
load("matrices/beta_matrices.carn.RData")
load("matrices/beta_matrices.omni.RData")

vec.bc.phylo.herb = as.vector(as.matrix(otu.bc.phylo.herb))
vec.bc.phylo.carn = as.vector(as.matrix(otu.bc.phylo.carn))
vec.wu.phylo.omni = as.vector(as.matrix(otu.bc.phylo.omni))
mat.phylo = as.data.frame(rbind(cbind(vec.bc.phylo.herb, rep("Herbivores", length(vec.bc.phylo.herb))),
                                cbind(vec.bc.phylo.carn, rep("Carnivores", length(vec.bc.phylo.carn))),
                                cbind(vec.wu.phylo.omni, rep("Omnivores", length(vec.wu.phylo.omni)))))


vec.div.herb = as.vector(mat_cophenetic_herb)
vec.div.carn = as.vector(mat_cophenetic_carn)
vec.div.omni = as.vector(mat_cophenetic_omni)
vec.div <- c(vec.div.herb, vec.div.carn, vec.div.omni)

mat = cbind(mat.phylo, vec.div)
colnames(mat) <- c("Dissimilarity" , "Diet", "Divergence")
mat <- mat[mat$Dissimilarity !=0,]
mat$Dissimilarity <- as.numeric(as.character(mat$Dissimilarity))

Fig.B1 = ggplot(mat, aes(y = Dissimilarity, x = Divergence)) + 
  geom_point(size = 4, alpha = 0.2, colour = "black",shape = 16) + 
  geom_smooth(method = "lm", colour = "red", alpha = 0.2) + 
  labs(x = "Host Divergence Time (My)", y = "Bacteriome Dissimilarity (Bray-Curtis)") + 
  theme( axis.text.x = element_text(colour = "black", size = 12,family = "serif"), 
         axis.text.y = element_text(size = 12, colour = "black",family = "serif"), 
         axis.title= element_text(size = 0, colour = "black", family = "serif"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"),
         legend.position = "top",
         legend.text = element_text(size = 3, family = "serif"),
         legend.title = element_text(size = 8, family = "serif"),
         strip.text.x = element_text(size=12, family = "serif")) +
  facet_wrap(~ Diet, ncol = 1)
Fig.B1

pdf(file = "Mantel_divergence_BC_diet.pdf", he = 7 , wi = 3)
Fig.B1
dev.off()




### DISTRIBUTION ------

## All files have been proceeded from the "6_topology_congruency.R" script
library("devtools")
install_github("WGS-TB/RFDistributionR")
library("rfdistr")
library("ape")
rfdistr::ntt_polynomial(rtree(5),6)
library("rfdistr")
library("ape")
rfdistr::polynomial(rtree(5),6)


# CTU - Firmicutes ----
setwd("./Firmicutes")

# __ ** Triplicate species ----
load("../matrices/beta_matrices.tri.RData")
load("tree.firmi.RData")
load("matrices/mat_cophenetic_firm")

bc.firm.mantel = mantel(otu.bc.firmi.tri, mat_cophenetic_firm , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.firm.mantel = mantel(otu.wu.firmi.tri, mat_cophenetic_firm , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.firm.mantel = mantel(otu.uu.firmi.tri, mat_cophenetic_firm , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.firm.mantel$statistic,3), ";","p =",round(bc.firm.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.firm.mantel$statistic,3), ";","p =",round(wu.firm.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.firm.mantel$statistic,3), ";","p =",round(uu.firm.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("tree.firmi.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.firmi.herb.rel.RData")

sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
tree_firmi_herb <- keep.tip(tree_firmi, sp_herb)

bc.firm.mantel = mantel(otu.bc.firmi.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.firm.mantel = mantel(otu.wu.firmi.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.firm.mantel = mantel(otu.uu.firmi.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.firm.mantel$statistic,3), ";","p =",round(bc.firm.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.firm.mantel$statistic,3), ";","p =",round(wu.firm.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.firm.mantel$statistic,3), ";","p =",round(uu.firm.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("tree.firmi.RData")
load("matrices/mat_cophenetic_carn.RData")
load("core.phylosymb.firmi.carn.rel.RData")

sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_firmi$tip.label)]
tree_firmi_carn <- keep.tip(tree_firmi, sp_carn)

bc.firm.mantel = mantel(otu.bc.firmi.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.firm.mantel = mantel(otu.wu.firmi.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.firm.mantel = mantel(otu.uu.firmi.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.firm.mantel$statistic,3), ";","p =",round(bc.firm.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.firm.mantel$statistic,3), ";","p =",round(wu.firm.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.firm.mantel$statistic,3), ";","p =",round(uu.firm.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Omnivores ----
load("matrices/mat_cophenetic_omni.RData")
load("core.phylosymb.firmi.omni.rel.RData")

sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_firmi$tip.label)]
tree_firmi_omni <- keep.tip(tree_firmi, sp_omni)

bc.firm.mantel = mantel(otu.bc.firmi.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.firm.mantel = mantel(otu.wu.firmi.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.firm.mantel = mantel(otu.uu.firmi.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.firm.mantel$statistic,3), ";","p =",round(bc.firm.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.firm.mantel$statistic,3), ";","p =",round(wu.firm.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.firm.mantel$statistic,3), ";","p =",round(uu.firm.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_tri)

firmi_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(firmi_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(firmi_phylo, file = "firmi_phylo_results.txt", sep = "\t", quote = F)

# CTU - Clostridiales ----
setwd("./Clostridiales")
load("core.clostri.RData")
load("tree.clostri.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")

bc.clostri.mantel = mantel(otu.bc.clostri.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri.mantel = mantel(otu.wu.clostri.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri.mantel = mantel(otu.uu.clostri.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.clostri.mantel$statistic,3), ";","p =",round(bc.clostri.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.clostri.mantel$statistic,3), ";","p =",round(wu.clostri.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.clostri.mantel$statistic,3), ";","p =",round(uu.clostri.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.clostri.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
tree_clostri_herb <- keep.tip(tree_clostri, sp_herb)

bc.clostri.mantel = mantel(otu.bc.clostri.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri.mantel = mantel(otu.wu.clostri.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri.mantel = mantel(otu.uu.clostri.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.clostri.mantel$statistic,3), ";","p =",round(bc.clostri.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.clostri.mantel$statistic,3), ";","p =",round(wu.clostri.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.clostri.mantel$statistic,3), ";","p =",round(uu.clostri.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.RData")
load("core.clostri.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_clostri$tip.label)]
tree_clostri_carn <- keep.tip(tree_clostri, sp_carn)

bc.clostri.mantel = mantel(otu.bc.clostri.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri.mantel = mantel(otu.wu.clostri.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri.mantel = mantel(otu.uu.clostri.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.clostri.mantel$statistic,3), ";","p =",round(bc.clostri.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.clostri.mantel$statistic,3), ";","p =",round(wu.clostri.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.clostri.mantel$statistic,3), ";","p =",round(uu.clostri.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.RData")
load("core.clostri.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_clostri$tip.label)]
tree_clostri_omni <- keep.tip(tree_clostri, sp_omni)

bc.clostri.mantel = mantel(otu.bc.clostri.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri.mantel = mantel(otu.wu.clostri.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri.mantel = mantel(otu.uu.clostri.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.clostri.mantel$statistic,3), ";","p =",round(bc.clostri.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.clostri.mantel$statistic,3), ";","p =",round(wu.clostri.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.clostri.mantel$statistic,3), ";","p =",round(uu.clostri.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_tri)

clostri_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(clostri_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(clostri_phylo, file = "clostri_phylo_results.txt", sep = "\t", quote = F)

# CTU - Clostridium_sensu_stricto_1 ----
setwd("./Clostridium_sensu_stricto_1")
load("core.clostri1.RData")
load("tree.clostri1.RData")

# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")
load("core.phylosymb.clostri1.tri.rel.RData")

bc.clostri1.mantel = mantel(otu.bc.clostri1.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri1.mantel = mantel(otu.wu.clostri1.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri1.mantel = mantel(otu.uu.clostri1.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.clostri1.mantel$statistic,3), ";","p =",round(bc.clostri1.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.clostri1.mantel$statistic,3), ";","p =",round(wu.clostri1.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.clostri1.mantel$statistic,3), ";","p =",round(uu.clostri1.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.clostri1.herb.rel.RData")

sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_clostri1$tip.label)]
tree_clostri1_herb <- keep.tip(tree_clostri1, sp_herb)

bc.clostri1.mantel = mantel(otu.bc.clostri1.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri1.mantel = mantel(otu.wu.clostri1.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri1.mantel = mantel(otu.uu.clostri1.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.clostri1.mantel$statistic,3), ";","p =",round(bc.clostri1.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.clostri1.mantel$statistic,3), ";","p =",round(wu.clostri1.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.clostri1.mantel$statistic,3), ";","p =",round(uu.clostri1.mantel$signif,3))
phylo_herb <- cbind("Jaccard" = phylo_jacc_herb, "Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.RData")
load("core.clostri1.carn.rel.RData")

sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_clostri1$tip.label)]
tree_clostri1_carn <- keep.tip(tree_clostri1, sp_carn)

bc.clostri1.mantel = mantel(otu.bc.clostri1.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri1.mantel = mantel(otu.wu.clostri1.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri1.mantel = mantel(otu.uu.clostri1.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.clostri1.mantel$statistic,3), ";","p =",round(bc.clostri1.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.clostri1.mantel$statistic,3), ";","p =",round(wu.clostri1.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.clostri1.mantel$statistic,3), ";","p =",round(uu.clostri1.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ **Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.RData")
load("core.clostri1.omni.rel.RData")

sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_clostri1$tip.label)]
tree_clostri1_omni <- keep.tip(tree_clostri1, sp_omni)

bc.clostri1.mantel = mantel(otu.bc.clostri1.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.clostri1.mantel = mantel(otu.wu.clostri1.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.clostri1.mantel = mantel(otu.uu.clostri1.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.clostri1.mantel$statistic,3), ";","p =",round(bc.clostri1.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.clostri1.mantel$statistic,3), ";","p =",round(wu.clostri1.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.clostri1.mantel$statistic,3), ";","p =",round(uu.clostri1.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

clostri1_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(clostri1_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(clostri1_phylo, file = "clostri1_phylo_results.txt", sep = "\t", quote = F)

# CTU - Lachnospiraceae ----
setwd("../lachnospiraceae")
load("core.lachno.RData")
load("tree.lachno.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_lachno.RData")
load("core.phylosymb.lachno.tri.rel.RData")
bc.lachno.mantel = mantel(otu.bc.lachno.tri, mat_cophenetic_lachno , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.lachno.mantel = mantel(otu.wu.lachno.tri, mat_cophenetic_lachno , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.lachno.mantel = mantel(otu.uu.lachno.tri, mat_cophenetic_lachno , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.lachno.mantel$statistic,3), ";","p =",round(bc.lachno.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.lachno.mantel$statistic,3), ";","p =",round(wu.lachno.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.lachno.mantel$statistic,3), ";","p =",round(uu.lachno.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.lachno.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_lachno$tip.label)]

tree_lachno_herb <- keep.tip(tree_lachno, sp_herb)

bc.lachno.mantel = mantel(otu.bc.lachno.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.lachno.mantel = mantel(otu.wu.lachno.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.lachno.mantel = mantel(otu.uu.lachno.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.lachno.mantel$statistic,3), ";","p =",round(bc.lachno.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.lachno.mantel$statistic,3), ";","p =",round(wu.lachno.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.lachno.mantel$statistic,3), ";","p =",round(uu.lachno.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.phylosymb.lachno.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_lachno$tip.label)]
tree_lachno_carn <- keep.tip(tree_lachno, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(tree_lachno_carn)/2

bc.lachno.mantel = mantel(otu.bc.lachno.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.lachno.mantel = mantel(otu.wu.lachno.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.lachno.mantel = mantel(otu.uu.lachno.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.lachno.mantel$statistic,3), ";","p =",round(bc.lachno.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.lachno.mantel$statistic,3), ";","p =",round(wu.lachno.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.lachno.mantel$statistic,3), ";","p =",round(uu.lachno.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.lachno.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_lachno$tip.label)]
tree_lachno_omni <- keep.tip(tree_lachno, sp_omni)

bc.lachno.mantel = mantel(otu.bc.lachno.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.lachno.mantel = mantel(otu.wu.lachno.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.lachno.mantel = mantel(otu.uu.lachno.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.lachno.mantel$statistic,3), ";","p =",round(bc.lachno.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.lachno.mantel$statistic,3), ";","p =",round(wu.lachno.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.lachno.mantel$statistic,3), ";","p =",round(uu.lachno.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

lachno_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(lachno_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(lachno_phylo, file = "lachno_phylo_results.txt", sep = "\t", quote = F)

# CTU - Epulopiscium ----
setwd("../Epulopiscium")
load("core.epulo.RData")
load("tree.epulo.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")

bc.epulo.mantel = mantel(otu.bc.epulo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.epulo.mantel = mantel(otu.wu.epulo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.epulo.mantel = mantel(otu.uu.epulo.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.epulo.mantel$statistic,3), ";","p =",round(bc.epulo.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.epulo.mantel$statistic,3), ";","p =",round(wu.epulo.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.epulo.mantel$statistic,3), ";","p =",round(uu.epulo.mantel$signif,3))
phylo_tri <- cbind("Jaccard" = phylo_jacc_tri, "Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.epulo.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_epulo$tip.label)]
core.epulo.herb.rel <- subset_samples(core.epulo, tax2 %in% sp_herb)
tree_epulo_herb <- keep.tip(tree_epulo, sp_herb)

bc.epulo.mantel = mantel(otu.bc.epulo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.epulo.mantel = mantel(otu.wu.epulo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.epulo.mantel = mantel(otu.uu.epulo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.epulo.mantel$statistic,3), ";","p =",round(bc.epulo.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.epulo.mantel$statistic,3), ";","p =",round(wu.epulo.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.epulo.mantel$statistic,3), ";","p =",round(uu.epulo.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.phylosymb.epulo.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_epulo$tip.label)]
tree_epulo_carn <- keep.tip(tree_epulo, sp_carn)

bc.epulo.mantel = mantel(otu.bc.epulo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.epulo.mantel = mantel(otu.wu.epulo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.epulo.mantel = mantel(otu.uu.epulo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.epulo.mantel$statistic,3), ";","p =",round(bc.epulo.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.epulo.mantel$statistic,3), ";","p =",round(wu.epulo.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.epulo.mantel$statistic,3), ";","p =",round(uu.epulo.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)


# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.phylosymb.epulo.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_epulo$tip.label)]
tree_epulo_omni <- keep.tip(tree_epulo, sp_omni)

bc.epulo.mantel = mantel(otu.bc.epulo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.epulo.mantel = mantel(otu.wu.epulo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.epulo.mantel = mantel(otu.uu.epulo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.epulo.mantel$statistic,3), ";","p =",round(bc.epulo.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.epulo.mantel$statistic,3), ";","p =",round(wu.epulo.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.epulo.mantel$statistic,3), ";","p =",round(uu.epulo.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

epulo_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(epulo_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(epulo_phylo, file = "epulo_phylo_results.txt", sep = "\t", quote = F)

# CTU - Rhodobacterales ----
setwd("../../../Rhodobacterales")
load("core.rhodobacterales.RData")
load("tree.rhodobacterales.RData")

# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")
load("core.phylosymb.rhodobacterales.tri.rel.RData")

bc.rhodobacterales.mantel = mantel(otu.bc.rhodobacterales.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacterales.mantel = mantel(otu.wu.rhodobacterales.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacterales.mantel = mantel(otu.uu.rhodobacterales.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.rhodobacterales.mantel$statistic,3), ";","p =",round(bc.rhodobacterales.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.rhodobacterales.mantel$statistic,3), ";","p =",round(wu.rhodobacterales.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.rhodobacterales.mantel$statistic,3), ";","p =",round(uu.rhodobacterales.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ **Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.rhodobacterales.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_rhodobacterales$tip.label)]
tree_rhodobacterales_herb <- keep.tip(tree_rhodobacterales, sp_herb)

bc.rhodobacterales.mantel = mantel(otu.bc.rhodobacterales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacterales.mantel = mantel(otu.wu.rhodobacterales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacterales.mantel = mantel(otu.uu.rhodobacterales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.rhodobacterales.mantel$statistic,3), ";","p =",round(bc.rhodobacterales.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.rhodobacterales.mantel$statistic,3), ";","p =",round(wu.rhodobacterales.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.rhodobacterales.mantel$statistic,3), ";","p =",round(uu.rhodobacterales.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ **Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.RData")
load("core.phylosymb.rhodobacterales.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_rhodobacterales$tip.label)]
tree_rhodobacterales_carn <- keep.tip(tree_rhodobacterales, sp_carn)

bc.rhodobacterales.mantel = mantel(otu.bc.rhodobacterales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacterales.mantel = mantel(otu.wu.rhodobacterales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacterales.mantel = mantel(otu.uu.rhodobacterales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.rhodobacterales.mantel$statistic,3), ";","p =",round(bc.rhodobacterales.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.rhodobacterales.mantel$statistic,3), ";","p =",round(wu.rhodobacterales.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.rhodobacterales.mantel$statistic,3), ";","p =",round(uu.rhodobacterales.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ **Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.RData")
load("core.phylosymb.rhodobacterales.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_rhodobacterales$tip.label)]
tree_rhodobacterales_omni <- keep.tip(tree_rhodobacterales, sp_omni)

bc.rhodobacterales.mantel = mantel(otu.bc.rhodobacterales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacterales.mantel = mantel(otu.wu.rhodobacterales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacterales.mantel = mantel(otu.uu.rhodobacterales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.rhodobacterales.mantel$statistic,3), ";","p =",round(bc.rhodobacterales.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.rhodobacterales.mantel$statistic,3), ";","p =",round(wu.rhodobacterales.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.rhodobacterales.mantel$statistic,3), ";","p =",round(uu.rhodobacterales.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

rhodobacterales_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(rhodobacterales_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(rhodobacterales_phylo, file = "rhodobacterales_phylo_results.txt", sep = "\t", quote = F)

# CTU - Rhodobacteraceae ----
setwd("./rhodobacteraceae")
load("core.rhodobacteraceae.RData")
load("tree.rhodobacteraceae.RData")

# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")
load("core.phylosymb.rhodobacteraceae.tri.rel.RData")

bc.rhodobacteraceae.mantel = mantel(otu.bc.rhodobacteraceae.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacteraceae.mantel = mantel(otu.wu.rhodobacteraceae.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacteraceae.mantel = mantel(otu.uu.rhodobacteraceae.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.rhodobacteraceae.mantel$statistic,3), ";","p =",round(bc.rhodobacteraceae.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(wu.rhodobacteraceae.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(uu.rhodobacteraceae.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ **Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.rhodobacteraceae.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_rhodobacteraceae$tip.label)]
tree_rhodobacteraceae_herb <- keep.tip(tree_rhodobacteraceae, sp_herb)

bc.rhodobacteraceae.mantel = mantel(otu.bc.rhodobacteraceae.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacteraceae.mantel = mantel(otu.wu.rhodobacteraceae.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacteraceae.mantel = mantel(otu.uu.rhodobacteraceae.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.rhodobacteraceae.mantel$statistic,3), ";","p =",round(bc.rhodobacteraceae.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(wu.rhodobacteraceae.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(uu.rhodobacteraceae.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ **Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.RData")
load("core.phylosymb.rhodobacteraceae.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_rhodobacteraceae$tip.label)]
tree_rhodobacteraceae_carn <- keep.tip(tree_rhodobacteraceae, sp_carn)

bc.rhodobacteraceae.mantel = mantel(otu.bc.rhodobacteraceae.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacteraceae.mantel = mantel(otu.wu.rhodobacteraceae.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacteraceae.mantel = mantel(otu.uu.rhodobacteraceae.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.rhodobacteraceae.mantel$statistic,3), ";","p =",round(bc.rhodobacteraceae.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(wu.rhodobacteraceae.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(uu.rhodobacteraceae.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ **Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.RData")
load("core.phylosymb.rhodobacteraceae.omni.rel.RData")

sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_rhodobacteraceae$tip.label)]
tree_rhodobacteraceae_omni <- keep.tip(tree_rhodobacteraceae, sp_omni)

bc.rhodobacteraceae.mantel = mantel(otu.bc.rhodobacteraceae.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhodobacteraceae.mantel = mantel(otu.wu.rhodobacteraceae.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhodobacteraceae.mantel = mantel(otu.uu.rhodobacteraceae.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.rhodobacteraceae.mantel$statistic,3), ";","p =",round(bc.rhodobacteraceae.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(wu.rhodobacteraceae.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.rhodobacteraceae.mantel$statistic,3), ";","p =",round(uu.rhodobacteraceae.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

rhodobacteraceae_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(rhodobacteraceae_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(rhodobacteraceae_phylo, file = "rhodobacteraceae_phylo_results.txt", sep = "\t", quote = F)

# CTU - Rhizobiales ----
setwd("../Rhizobiales")
load("core.rhizo.RData")
load("tree.rhizo.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")
load("core.phylosymb.rhizo.tri.rel.RData")
bc.rhizo.mantel = mantel(otu.bc.rhizo.tri, mat_cophenetic_rhizo , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhizo.mantel = mantel(otu.wu.rhizo.tri, mat_cophenetic_rhizo , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhizo.mantel = mantel(otu.uu.rhizo.tri, mat_cophenetic_rhizo , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.rhizo.mantel$statistic,3), ";","p =",round(bc.rhizo.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.rhizo.mantel$statistic,3), ";","p =",round(wu.rhizo.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.rhizo.mantel$statistic,3), ";","p =",round(uu.rhizo.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.rhizo.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_rhizo$tip.label)]
tree_rhizo_herb <- keep.tip(tree_rhizo, sp_herb)

bc.rhizo.mantel = mantel(otu.bc.rhizo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhizo.mantel = mantel(otu.wu.rhizo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhizo.mantel = mantel(otu.uu.rhizo.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.rhizo.mantel$statistic,3), ";","p =",round(bc.rhizo.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.rhizo.mantel$statistic,3), ";","p =",round(wu.rhizo.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.rhizo.mantel$statistic,3), ";","p =",round(uu.rhizo.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.phylosymb.rhizo.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_rhizo$tip.label)]
tree_rhizo_carn <- keep.tip(tree_rhizo, sp_carn)

bc.rhizo.mantel = mantel(otu.bc.rhizo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhizo.mantel = mantel(otu.wu.rhizo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhizo.mantel = mantel(otu.uu.rhizo.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.rhizo.mantel$statistic,3), ";","p =",round(bc.rhizo.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.rhizo.mantel$statistic,3), ";","p =",round(wu.rhizo.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.rhizo.mantel$statistic,3), ";","p =",round(uu.rhizo.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.rhizo.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_rhizo$tip.label)]
tree_rhizo_omni <- keep.tip(tree_rhizo, sp_omni)

bc.rhizo.mantel = mantel(otu.bc.rhizo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.rhizo.mantel = mantel(otu.wu.rhizo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.rhizo.mantel = mantel(otu.uu.rhizo.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.rhizo.mantel$statistic,3), ";","p =",round(bc.rhizo.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.rhizo.mantel$statistic,3), ";","p =",round(wu.rhizo.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.rhizo.mantel$statistic,3), ";","p =",round(uu.rhizo.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

rhizo_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(rhizo_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(rhizo_phylo, file = "rhizo_phylo_results.txt", sep = "\t", quote = F)

# CTU - Vibrionales ----
setwd("../../Vibrionales")
load("core.vibrionales.RData")
load("tree.vibrionales.RData")

# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tri.RData")

bc.vibrionales.mantel = mantel(otu.bc.vibrionales.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.vibrionales.mantel = mantel(otu.wu.vibrionales.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.vibrionales.mantel = mantel(otu.uu.vibrionales.tri, mat_cophenetic_tri , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.vibrionales.mantel$statistic,3), ";","p =",round(bc.vibrionales.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.vibrionales.mantel$statistic,3), ";","p =",round(wu.vibrionales.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.vibrionales.mantel$statistic,3), ";","p =",round(uu.vibrionales.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ **Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.vibrionales.herb.rel.RData")

sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_vibrionales$tip.label)]
tree_vibrionales_herb <- keep.tip(tree_vibrionales, sp_herb)

bc.vibrionales.mantel = mantel(otu.bc.vibrionales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.vibrionales.mantel = mantel(otu.wu.vibrionales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.vibrionales.mantel = mantel(otu.uu.vibrionales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.vibrionales.mantel$statistic,3), ";","p =",round(bc.vibrionales.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.vibrionales.mantel$statistic,3), ";","p =",round(wu.vibrionales.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.vibrionales.mantel$statistic,3), ";","p =",round(uu.vibrionales.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ **Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.RData")
load("core.phylosymb.vibrionales.carn.rel.RData")

sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_vibrionales$tip.label)]
tree_vibrionales_carn <- keep.tip(tree_vibrionales, sp_carn)

bc.vibrionales.mantel = mantel(otu.bc.vibrionales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.vibrionales.mantel = mantel(otu.wu.vibrionales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.vibrionales.mantel = mantel(otu.uu.vibrionales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.vibrionales.mantel$statistic,3), ";","p =",round(bc.vibrionales.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.vibrionales.mantel$statistic,3), ";","p =",round(wu.vibrionales.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.vibrionales.mantel$statistic,3), ";","p =",round(uu.vibrionales.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ **Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.RData")
load("core.phylosymb.vibrionales.omni.rel.RData")

sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_vibrionales$tip.label)]
tree_vibrionales_omni <- keep.tip(tree_vibrionales, sp_omni)

bc.vibrionales.mantel = mantel(otu.bc.vibrionales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.vibrionales.mantel = mantel(otu.wu.vibrionales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.vibrionales.mantel = mantel(otu.uu.vibrionales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.vibrionales.mantel$statistic,3), ";","p =",round(bc.vibrionales.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.vibrionales.mantel$statistic,3), ";","p =",round(wu.vibrionales.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.vibrionales.mantel$statistic,3), ";","p =",round(uu.vibrionales.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

vibrionales_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(vibrionales_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(vibrionales_phylo, file = "vibrionales_phylo_results.txt", sep = "\t", quote = F)

# CTU - Planctomycetes ----
setwd("../Planctomycetes")
load("core.plancto.RData")
load("tree.plancto.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_plancto.RData")
load("core.phylosymb.plancto.tri.rel.RData")
bc.plancto.mantel = mantel(otu.bc.plancto.tri, mat_cophenetic_plancto , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.plancto.mantel = mantel(otu.wu.plancto.tri, mat_cophenetic_plancto , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.plancto.mantel = mantel(otu.uu.plancto.tri, mat_cophenetic_plancto , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.plancto.mantel$statistic,3), ";","p =",round(bc.plancto.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.plancto.mantel$statistic,3), ";","p =",round(wu.plancto.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.plancto.mantel$statistic,3), ";","p =",round(uu.plancto.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.plancto.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_plancto$tip.label)]
tree_plancto_herb <- keep.tip(tree_plancto, sp_herb)

bc.plancto.mantel = mantel(otu.bc.plancto.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.plancto.mantel = mantel(otu.wu.plancto.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.plancto.mantel = mantel(otu.uu.plancto.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.plancto.mantel$statistic,3), ";","p =",round(bc.plancto.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.plancto.mantel$statistic,3), ";","p =",round(wu.plancto.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.plancto.mantel$statistic,3), ";","p =",round(uu.plancto.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.phylosymb.plancto.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_plancto$tip.label)]
tree_plancto_carn <- keep.tip(tree_plancto, sp_carn)

bc.plancto.mantel = mantel(otu.bc.plancto.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.plancto.mantel = mantel(otu.wu.plancto.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.plancto.mantel = mantel(otu.uu.plancto.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.plancto.mantel$statistic,3), ";","p =",round(bc.plancto.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.plancto.mantel$statistic,3), ";","p =",round(wu.plancto.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.plancto.mantel$statistic,3), ";","p =",round(uu.plancto.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.plancto.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_plancto$tip.label)]
tree_plancto_omni <- keep.tip(tree_plancto, sp_omni)

bc.plancto.mantel = mantel(otu.bc.plancto.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.plancto.mantel = mantel(otu.wu.plancto.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.plancto.mantel = mantel(otu.uu.plancto.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.plancto.mantel$statistic,3), ";","p =",round(bc.plancto.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.plancto.mantel$statistic,3), ";","p =",round(wu.plancto.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.plancto.mantel$statistic,3), ";","p =",round(uu.plancto.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

plancto_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(plancto_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(plancto_phylo, file = "plancto_phylo_results.txt", sep = "\t", quote = F)

# CTU - Pirellulales ----
setwd("./Pirellulales")
load("core.pirellulales.RData")
load("tree.pirellulales.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_pirellulales.RData")
load("core.phylosymb.pirellulales.tri.rel.RData")
bc.pirellulales.mantel = mantel(otu.bc.pirellulales.tri, mat_cophenetic_pirellulales , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.pirellulales.mantel = mantel(otu.wu.pirellulales.tri, mat_cophenetic_pirellulales , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.pirellulales.mantel = mantel(otu.uu.pirellulales.tri, mat_cophenetic_pirellulales , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.pirellulales.mantel$statistic,3), ";","p =",round(bc.pirellulales.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.pirellulales.mantel$statistic,3), ";","p =",round(wu.pirellulales.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.pirellulales.mantel$statistic,3), ";","p =",round(uu.pirellulales.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.phylosymb.pirellulales.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_pirellulales$tip.label)]
tree_pirellulales_herb <- keep.tip(tree_pirellulales, sp_herb)

bc.pirellulales.mantel = mantel(otu.bc.pirellulales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.pirellulales.mantel = mantel(otu.wu.pirellulales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.pirellulales.mantel = mantel(otu.uu.pirellulales.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.pirellulales.mantel$statistic,3), ";","p =",round(bc.pirellulales.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.pirellulales.mantel$statistic,3), ";","p =",round(wu.pirellulales.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.pirellulales.mantel$statistic,3), ";","p =",round(uu.pirellulales.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.phylosymb.pirellulales.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_pirellulales$tip.label)]
tree_pirellulales_carn <- keep.tip(tree_pirellulales, sp_carn)

bc.pirellulales.mantel = mantel(otu.bc.pirellulales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.pirellulales.mantel = mantel(otu.wu.pirellulales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.pirellulales.mantel = mantel(otu.uu.pirellulales.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.pirellulales.mantel$statistic,3), ";","p =",round(bc.pirellulales.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.pirellulales.mantel$statistic,3), ";","p =",round(wu.pirellulales.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.pirellulales.mantel$statistic,3), ";","p =",round(uu.pirellulales.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.pirellulales.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_pirellulales$tip.label)]
tree_pirellulales_omni <- keep.tip(tree_pirellulales, sp_omni)

bc.pirellulales.mantel = mantel(otu.bc.pirellulales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.pirellulales.mantel = mantel(otu.wu.pirellulales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.pirellulales.mantel = mantel(otu.uu.pirellulales.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.pirellulales.mantel$statistic,3), ";","p =",round(bc.pirellulales.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.pirellulales.mantel$statistic,3), ";","p =",round(wu.pirellulales.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.pirellulales.mantel$statistic,3), ";","p =",round(uu.pirellulales.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

pirellulales_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(pirellulales_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(pirellulales_phylo, file = "pirellulales_phylo_results.txt", sep = "\t", quote = F)
# CTU - Cyanobacteria ----
setwd("./Cyanobacteria")
load("core.cyanobacteria.RData")
load("tree.cyanobacteria.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_cyanobacteria.RData")
load("core.phylosymb.cyanobacteria.tri.rel.RData")
bc.cyanobacteria.mantel = mantel(otu.bc.cyanobacteria.tri, mat_cophenetic_cyanobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.cyanobacteria.mantel = mantel(otu.wu.cyanobacteria.tri, mat_cophenetic_cyanobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.cyanobacteria.mantel = mantel(otu.uu.cyanobacteria.tri, mat_cophenetic_cyanobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.cyanobacteria.mantel$statistic,3), ";","p =",round(bc.cyanobacteria.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.cyanobacteria.mantel$statistic,3), ";","p =",round(wu.cyanobacteria.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.cyanobacteria.mantel$statistic,3), ";","p =",round(uu.cyanobacteria.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.cyanobacteria.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_cyanobacteria$tip.label)]
tree_cyanobacteria_herb <- keep.tip(tree_cyanobacteria, sp_herb)

bc.cyanobacteria.mantel = mantel(otu.bc.cyanobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.cyanobacteria.mantel = mantel(otu.wu.cyanobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.cyanobacteria.mantel = mantel(otu.uu.cyanobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.cyanobacteria.mantel$statistic,3), ";","p =",round(bc.cyanobacteria.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.cyanobacteria.mantel$statistic,3), ";","p =",round(wu.cyanobacteria.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.cyanobacteria.mantel$statistic,3), ";","p =",round(uu.cyanobacteria.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.cyanobacteria.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_cyanobacteria$tip.label)]
tree_cyanobacteria_carn <- keep.tip(tree_cyanobacteria, sp_carn)

bc.cyanobacteria.mantel = mantel(otu.bc.cyanobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.cyanobacteria.mantel = mantel(otu.wu.cyanobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.cyanobacteria.mantel = mantel(otu.uu.cyanobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.cyanobacteria.mantel$statistic,3), ";","p =",round(bc.cyanobacteria.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.cyanobacteria.mantel$statistic,3), ";","p =",round(wu.cyanobacteria.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.cyanobacteria.mantel$statistic,3), ";","p =",round(uu.cyanobacteria.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.cyanobacteria.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_cyanobacteria$tip.label)]
tree_cyanobacteria_omni <- keep.tip(tree_cyanobacteria, sp_omni)

bc.cyanobacteria.mantel = mantel(otu.bc.cyanobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.cyanobacteria.mantel = mantel(otu.wu.cyanobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.cyanobacteria.mantel = mantel(otu.uu.cyanobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.cyanobacteria.mantel$statistic,3), ";","p =",round(bc.cyanobacteria.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.cyanobacteria.mantel$statistic,3), ";","p =",round(wu.cyanobacteria.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.cyanobacteria.mantel$statistic,3), ";","p =",round(uu.cyanobacteria.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

cyanobacteria_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(cyanobacteria_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(cyanobacteria_phylo, file = "cyanobacteria_phylo_results.txt", sep = "\t", quote = F)

# CTU - Proteobacteria ----
setwd("../Proteobacteria")
load("core.Proteobacteria.RData")
load("tree.Proteobacteria.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_Proteobacteria.RData")
load("core.phylosymb.Proteobacteria.tri.rel.RData")
bc.Proteobacteria.mantel = mantel(otu.bc.Proteobacteria.tri, mat_cophenetic_Proteobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.Proteobacteria.mantel = mantel(otu.wu.Proteobacteria.tri, mat_cophenetic_Proteobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.Proteobacteria.mantel = mantel(otu.uu.Proteobacteria.tri, mat_cophenetic_Proteobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.Proteobacteria.mantel$statistic,3), ";","p =",round(bc.Proteobacteria.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.Proteobacteria.mantel$statistic,3), ";","p =",round(wu.Proteobacteria.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.Proteobacteria.mantel$statistic,3), ";","p =",round(uu.Proteobacteria.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.Proteobacteria.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_Proteobacteria$tip.label)]
tree_Proteobacteria_herb <- keep.tip(tree_Proteobacteria, sp_herb)

bc.Proteobacteria.mantel = mantel(otu.bc.Proteobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.Proteobacteria.mantel = mantel(otu.wu.Proteobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.Proteobacteria.mantel = mantel(otu.uu.Proteobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.Proteobacteria.mantel$statistic,3), ";","p =",round(bc.Proteobacteria.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.Proteobacteria.mantel$statistic,3), ";","p =",round(wu.Proteobacteria.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.Proteobacteria.mantel$statistic,3), ";","p =",round(uu.Proteobacteria.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.Proteobacteria.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_Proteobacteria$tip.label)]
tree_Proteobacteria_carn <- keep.tip(tree_Proteobacteria, sp_carn)

bc.Proteobacteria.mantel = mantel(otu.bc.Proteobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.Proteobacteria.mantel = mantel(otu.wu.Proteobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.Proteobacteria.mantel = mantel(otu.uu.Proteobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.Proteobacteria.mantel$statistic,3), ";","p =",round(bc.Proteobacteria.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.Proteobacteria.mantel$statistic,3), ";","p =",round(wu.Proteobacteria.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.Proteobacteria.mantel$statistic,3), ";","p =",round(uu.Proteobacteria.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.Proteobacteria.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_Proteobacteria$tip.label)]
tree_Proteobacteria_omni <- keep.tip(tree_Proteobacteria, sp_omni)

bc.Proteobacteria.mantel = mantel(otu.bc.Proteobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.Proteobacteria.mantel = mantel(otu.wu.Proteobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.Proteobacteria.mantel = mantel(otu.uu.Proteobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.Proteobacteria.mantel$statistic,3), ";","p =",round(bc.Proteobacteria.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.Proteobacteria.mantel$statistic,3), ";","p =",round(wu.Proteobacteria.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.Proteobacteria.mantel$statistic,3), ";","p =",round(uu.Proteobacteria.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

Proteobacteria_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(Proteobacteria_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(Proteobacteria_phylo, file = "Proteobacteria_phylo_results.txt", sep = "\t", quote = F)

# CTU - Bacteroidetes ----
setwd("../Bacteroidetes")
load("core.bacteroidetes.RData")
load("tree.bacteroidetes.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_bacteroidetes.RData")
load("core.phylosymb.bacteroidetes.tri.rel.RData")
bc.bacteroidetes.mantel = mantel(otu.bc.bacteroidetes.tri, mat_cophenetic_bacteroidetes , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.bacteroidetes.mantel = mantel(otu.wu.bacteroidetes.tri, mat_cophenetic_bacteroidetes , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.bacteroidetes.mantel = mantel(otu.uu.bacteroidetes.tri, mat_cophenetic_bacteroidetes , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.bacteroidetes.mantel$statistic,3), ";","p =",round(bc.bacteroidetes.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.bacteroidetes.mantel$statistic,3), ";","p =",round(wu.bacteroidetes.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.bacteroidetes.mantel$statistic,3), ";","p =",round(uu.bacteroidetes.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.bacteroidetes.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_bacteroidetes$tip.label)]
tree_bacteroidetes_herb <- keep.tip(tree_bacteroidetes, sp_herb)

bc.bacteroidetes.mantel = mantel(otu.bc.bacteroidetes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.bacteroidetes.mantel = mantel(otu.wu.bacteroidetes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.bacteroidetes.mantel = mantel(otu.uu.bacteroidetes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.bacteroidetes.mantel$statistic,3), ";","p =",round(bc.bacteroidetes.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.bacteroidetes.mantel$statistic,3), ";","p =",round(wu.bacteroidetes.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.bacteroidetes.mantel$statistic,3), ";","p =",round(uu.bacteroidetes.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.bacteroidetes.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_bacteroidetes$tip.label)]
tree_bacteroidetes_carn <- keep.tip(tree_bacteroidetes, sp_carn)

bc.bacteroidetes.mantel = mantel(otu.bc.bacteroidetes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.bacteroidetes.mantel = mantel(otu.wu.bacteroidetes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.bacteroidetes.mantel = mantel(otu.uu.bacteroidetes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.bacteroidetes.mantel$statistic,3), ";","p =",round(bc.bacteroidetes.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.bacteroidetes.mantel$statistic,3), ";","p =",round(wu.bacteroidetes.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.bacteroidetes.mantel$statistic,3), ";","p =",round(uu.bacteroidetes.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.bacteroidetes.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_bacteroidetes$tip.label)]
tree_bacteroidetes_omni <- keep.tip(tree_bacteroidetes, sp_omni)

bc.bacteroidetes.mantel = mantel(otu.bc.bacteroidetes.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.bacteroidetes.mantel = mantel(otu.wu.bacteroidetes.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.bacteroidetes.mantel = mantel(otu.uu.bacteroidetes.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.bacteroidetes.mantel$statistic,3), ";","p =",round(bc.bacteroidetes.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.bacteroidetes.mantel$statistic,3), ";","p =",round(wu.bacteroidetes.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.bacteroidetes.mantel$statistic,3), ";","p =",round(uu.bacteroidetes.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

bacteroidetes_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(bacteroidetes_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(bacteroidetes_phylo, file = "bacteroidetes_phylo_results.txt", sep = "\t", quote = F)

# CTU - Fusobacteria ----
setwd("../Fusobacteria")
load("core.fusobacteria.RData")
load("tree.fusobacteria.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_fusobacteria.RData")
load("core.phylosymb.fusobacteria.tri.rel.RData")
bc.fusobacteria.mantel = mantel(otu.bc.fusobacteria.tri, mat_cophenetic_fusobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.fusobacteria.mantel = mantel(otu.wu.fusobacteria.tri, mat_cophenetic_fusobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.fusobacteria.mantel = mantel(otu.uu.fusobacteria.tri, mat_cophenetic_fusobacteria , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.fusobacteria.mantel$statistic,3), ";","p =",round(bc.fusobacteria.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.fusobacteria.mantel$statistic,3), ";","p =",round(wu.fusobacteria.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.fusobacteria.mantel$statistic,3), ";","p =",round(uu.fusobacteria.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.fusobacteria.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_fusobacteria$tip.label)]
tree_fusobacteria_herb <- keep.tip(tree_fusobacteria, sp_herb)

bc.fusobacteria.mantel = mantel(otu.bc.fusobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.fusobacteria.mantel = mantel(otu.wu.fusobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.fusobacteria.mantel = mantel(otu.uu.fusobacteria.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.fusobacteria.mantel$statistic,3), ";","p =",round(bc.fusobacteria.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.fusobacteria.mantel$statistic,3), ";","p =",round(wu.fusobacteria.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.fusobacteria.mantel$statistic,3), ";","p =",round(uu.fusobacteria.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.fusobacteria.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_fusobacteria$tip.label)]
tree_fusobacteria_carn <- keep.tip(tree_fusobacteria, sp_carn)

bc.fusobacteria.mantel = mantel(otu.bc.fusobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.fusobacteria.mantel = mantel(otu.wu.fusobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.fusobacteria.mantel = mantel(otu.uu.fusobacteria.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.fusobacteria.mantel$statistic,3), ";","p =",round(bc.fusobacteria.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.fusobacteria.mantel$statistic,3), ";","p =",round(wu.fusobacteria.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.fusobacteria.mantel$statistic,3), ";","p =",round(uu.fusobacteria.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.fusobacteria.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_fusobacteria$tip.label)]
tree_fusobacteria_omni <- keep.tip(tree_fusobacteria, sp_omni)

bc.fusobacteria.mantel = mantel(otu.bc.fusobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.fusobacteria.mantel = mantel(otu.wu.fusobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.fusobacteria.mantel = mantel(otu.uu.fusobacteria.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.fusobacteria.mantel$statistic,3), ";","p =",round(bc.fusobacteria.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.fusobacteria.mantel$statistic,3), ";","p =",round(wu.fusobacteria.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.fusobacteria.mantel$statistic,3), ";","p =",round(uu.fusobacteria.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

fusobacteria_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(fusobacteria_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(fusobacteria_phylo, file = "fusobacteria_phylo_results.txt", sep = "\t", quote = F)

# CTU - Tenericutes ----
setwd("./Tenericutes")
load("core.tenericutes.RData")
load("tree.tenericutes.RData")
# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_tenericutes.RData")
load("core.phylosymb.tenericutes.tri.rel.RData")
bc.tenericutes.mantel = mantel(otu.bc.tenericutes.tri, mat_cophenetic_tenericutes , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.tenericutes.mantel = mantel(otu.wu.tenericutes.tri, mat_cophenetic_tenericutes , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.tenericutes.mantel = mantel(otu.uu.tenericutes.tri, mat_cophenetic_tenericutes , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.tenericutes.mantel$statistic,3), ";","p =",round(bc.tenericutes.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.tenericutes.mantel$statistic,3), ";","p =",round(wu.tenericutes.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.tenericutes.mantel$statistic,3), ";","p =",round(uu.tenericutes.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.tenericutes.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_tenericutes$tip.label)]
tree_tenericutes_herb <- keep.tip(tree_tenericutes, sp_herb)

bc.tenericutes.mantel = mantel(otu.bc.tenericutes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.tenericutes.mantel = mantel(otu.wu.tenericutes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.tenericutes.mantel = mantel(otu.uu.tenericutes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.tenericutes.mantel$statistic,3), ";","p =",round(bc.tenericutes.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.tenericutes.mantel$statistic,3), ";","p =",round(wu.tenericutes.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.tenericutes.mantel$statistic,3), ";","p =",round(uu.tenericutes.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.tenericutes.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_tenericutes$tip.label)]
tree_tenericutes_carn <- keep.tip(tree_tenericutes, sp_carn)

bc.tenericutes.mantel = mantel(otu.bc.tenericutes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.tenericutes.mantel = mantel(otu.wu.tenericutes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.tenericutes.mantel = mantel(otu.uu.tenericutes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.tenericutes.mantel$statistic,3), ";","p =",round(bc.tenericutes.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.tenericutes.mantel$statistic,3), ";","p =",round(wu.tenericutes.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.tenericutes.mantel$statistic,3), ";","p =",round(uu.tenericutes.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.tenericutes.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_tenericutes$tip.label)]
tree_tenericutes_omni <- keep.tip(tree_tenericutes, sp_omni)

bc.tenericutes.mantel = mantel(otu.bc.tenericutes.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.tenericutes.mantel = mantel(otu.wu.tenericutes.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.tenericutes.mantel = mantel(otu.uu.tenericutes.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.tenericutes.mantel$statistic,3), ";","p =",round(bc.tenericutes.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.tenericutes.mantel$statistic,3), ";","p =",round(wu.tenericutes.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.tenericutes.mantel$statistic,3), ";","p =",round(uu.tenericutes.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

tenericutes_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(tenericutes_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(tenericutes_phylo, file = "tenericutes_phylo_results.txt", sep = "\t", quote = F)

# CTU - Verrucomicrobia ----
setwd("../Verrucomicrobia")
load("core.verrucomicrobia.RData")
load("tree.verrucomicrobia.RData")

# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_verrucomicrobia.RData")
load("core.phylosymb.verrucomicrobia.tri.rel.RData")
bc.verrucomicrobia.mantel = mantel(otu.bc.verrucomicrobia.tri, mat_cophenetic_verrucomicrobia , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.verrucomicrobia.mantel = mantel(otu.wu.verrucomicrobia.tri, mat_cophenetic_verrucomicrobia , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.verrucomicrobia.mantel = mantel(otu.uu.verrucomicrobia.tri, mat_cophenetic_verrucomicrobia , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.verrucomicrobia.mantel$statistic,3), ";","p =",round(bc.verrucomicrobia.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.verrucomicrobia.mantel$statistic,3), ";","p =",round(wu.verrucomicrobia.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.verrucomicrobia.mantel$statistic,3), ";","p =",round(uu.verrucomicrobia.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.verrucomicrobia.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_verrucomicrobia$tip.label)]
tree_verrucomicrobia_herb <- keep.tip(tree_verrucomicrobia, sp_herb)

bc.verrucomicrobia.mantel = mantel(otu.bc.verrucomicrobia.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.verrucomicrobia.mantel = mantel(otu.wu.verrucomicrobia.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.verrucomicrobia.mantel = mantel(otu.uu.verrucomicrobia.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.verrucomicrobia.mantel$statistic,3), ";","p =",round(bc.verrucomicrobia.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.verrucomicrobia.mantel$statistic,3), ";","p =",round(wu.verrucomicrobia.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.verrucomicrobia.mantel$statistic,3), ";","p =",round(uu.verrucomicrobia.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.verrucomicrobia.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_verrucomicrobia$tip.label)]
tree_verrucomicrobia_carn <- keep.tip(tree_verrucomicrobia, sp_carn)

bc.verrucomicrobia.mantel = mantel(otu.bc.verrucomicrobia.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.verrucomicrobia.mantel = mantel(otu.wu.verrucomicrobia.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.verrucomicrobia.mantel = mantel(otu.uu.verrucomicrobia.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.verrucomicrobia.mantel$statistic,3), ";","p =",round(bc.verrucomicrobia.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.verrucomicrobia.mantel$statistic,3), ";","p =",round(wu.verrucomicrobia.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.verrucomicrobia.mantel$statistic,3), ";","p =",round(uu.verrucomicrobia.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

# __ ** Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat_cophenetic_omni.Rdata")
load("core.verrucomicrobia.omni.rel.RData")
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_verrucomicrobia$tip.label)]
tree_verrucomicrobia_omni <- keep.tip(tree_verrucomicrobia, sp_omni)

bc.verrucomicrobia.mantel = mantel(otu.bc.verrucomicrobia.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.verrucomicrobia.mantel = mantel(otu.wu.verrucomicrobia.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.verrucomicrobia.mantel = mantel(otu.uu.verrucomicrobia.omni, mat_cophenetic_omni , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_omni <- paste("r =", round(bc.verrucomicrobia.mantel$statistic,3), ";","p =",round(bc.verrucomicrobia.mantel$signif,3))
phylo_wu_omni <- paste("r =", round(wu.verrucomicrobia.mantel$statistic,3), ";","p =",round(wu.verrucomicrobia.mantel$signif,3))
phylo_uu_omni <- paste("r =", round(uu.verrucomicrobia.mantel$statistic,3), ";","p =",round(uu.verrucomicrobia.mantel$signif,3))
phylo_omni <- cbind("Bray-Curtis" =phylo_bc_omni,"Weighted Unifrac" =phylo_wu_omni,"Unweighted Unifrac" =phylo_uu_omni)

verrucomicrobia_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn,phylo_omni)
rownames(verrucomicrobia_phylo) <- c("Triplicate","Herbivores","Carnivores","Omnivores")

write.table(verrucomicrobia_phylo, file = "verrucomicrobia_phylo_results.txt", sep = "\t", quote = F)

# CTU - Spirochaetes ----
setwd("../Spirochaetes")
load("core.spirochaetes.RData")
load("tree.spirochaetes.RData")

# __ ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat_cophenetic_spirochaetes.RData")
load("core.phylosymb.spirochaetes.tri.rel.RData")
bc.spirochaetes.mantel = mantel(otu.bc.spirochaetes.tri, mat_cophenetic_spirochaetes , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.spirochaetes.mantel = mantel(otu.wu.spirochaetes.tri, mat_cophenetic_spirochaetes , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.spirochaetes.mantel = mantel(otu.uu.spirochaetes.tri, mat_cophenetic_spirochaetes , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_tri <- paste("r =", round(bc.spirochaetes.mantel$statistic,3), ";","p =",round(bc.spirochaetes.mantel$signif,3))
phylo_wu_tri <- paste("r =", round(wu.spirochaetes.mantel$statistic,3), ";","p =",round(wu.spirochaetes.mantel$signif,3))
phylo_uu_tri <- paste("r =", round(uu.spirochaetes.mantel$statistic,3), ";","p =",round(uu.spirochaetes.mantel$signif,3))
phylo_tri <- cbind("Bray-Curtis" =phylo_bc_tri,"Weighted Unifrac" =phylo_wu_tri,"Unweighted Unifrac" =phylo_uu_tri)

# __ ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat_cophenetic_herb.RData")
load("core.spirochaetes.herb.rel.RData")
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_spirochaetes$tip.label)]
tree_spirochaetes_herb <- keep.tip(tree_spirochaetes, sp_herb)

bc.spirochaetes.mantel = mantel(otu.bc.spirochaetes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.spirochaetes.mantel = mantel(otu.wu.spirochaetes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.spirochaetes.mantel = mantel(otu.uu.spirochaetes.herb, mat_cophenetic_herb , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_herb <- paste("r =", round(bc.spirochaetes.mantel$statistic,3), ";","p =",round(bc.spirochaetes.mantel$signif,3))
phylo_wu_herb <- paste("r =", round(wu.spirochaetes.mantel$statistic,3), ";","p =",round(wu.spirochaetes.mantel$signif,3))
phylo_uu_herb <- paste("r =", round(uu.spirochaetes.mantel$statistic,3), ";","p =",round(uu.spirochaetes.mantel$signif,3))
phylo_herb <- cbind("Bray-Curtis" =phylo_bc_herb,"Weighted Unifrac" =phylo_wu_herb,"Unweighted Unifrac" =phylo_uu_herb)

# __ ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat_cophenetic_carn.Rdata")
load("core.spirochaetes.carn.rel.RData")
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_spirochaetes$tip.label)]
tree_spirochaetes_carn <- keep.tip(tree_spirochaetes, sp_carn)

bc.spirochaetes.mantel = mantel(otu.bc.spirochaetes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
wu.spirochaetes.mantel = mantel(otu.wu.spirochaetes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)
uu.spirochaetes.mantel = mantel(otu.uu.spirochaetes.carn, mat_cophenetic_carn , method = "spearman", permutations = 9999, na.rm = TRUE)

phylo_bc_carn <- paste("r =", round(bc.spirochaetes.mantel$statistic,3), ";","p =",round(bc.spirochaetes.mantel$signif,3))
phylo_wu_carn <- paste("r =", round(wu.spirochaetes.mantel$statistic,3), ";","p =",round(wu.spirochaetes.mantel$signif,3))
phylo_uu_carn <- paste("r =", round(uu.spirochaetes.mantel$statistic,3), ";","p =",round(uu.spirochaetes.mantel$signif,3))
phylo_carn <- cbind("Bray-Curtis" =phylo_bc_carn,"Weighted Unifrac" =phylo_wu_carn,"Unweighted Unifrac" =phylo_uu_carn)

verrucomicrobia_phylo <-rbind(phylo_tri,phylo_herb,phylo_carn)
rownames(verrucomicrobia_phylo) <- c("Triplicate","Herbivores","Carnivores")

write.table(verrucomicrobia_phylo, file = "verrucomicrobia_phylo_results.txt", sep = "\t", quote = F)

## Supplementary results table ----






