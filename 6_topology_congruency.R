## Phylosymbiosis (B-divergence) - Congruency tests -----
# Hypothesis: the microbial community relationships dendrograms reflects the host phylogeny relationships.
# 
# Setup ----
setwd("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/phylogeny")
dir.create("./matrices", recursive = T)
dir.create("./Figures", recursive = T)

library("ape")
library("Biostrings")
library("ggplot2")
library(ggtree)
library(phangorn)
library(TreeDist)
library(treeman)
library(dendextend)

load("core.phylosymb.RData")
load("core.phylosymb.tri.rel.RData")
host_tree <- read.tree(file = "fish_tree_modif.tre")
load("host_tree_tri.RData")

# BACTERIOME ----
# ** All species ----
load("matrices/beta_matrices.all.RData")
load("matrices/mat.cophenetic.RData")
tree_jacc_all <- upgma(otu.jacc.phylo, method = "average")
tree_bc_all <- upgma(otu.bc.phylo, method = "average")
tree_wu_all <- upgma(otu.wu.phylo, method = "average")
tree_uu_all <- upgma(otu.uu.phylo, method = "average")

RF_jacc_all <- TreeDistance(host_tree,tree_jacc_all)
RF_bc_all <- TreeDistance(host_tree,tree_bc_all)
RF_wu_all <- TreeDistance(host_tree,tree_wu_all)
RF_uu_all <- TreeDistance(host_tree,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

d_fish_all <- upgma(mat_cophenetic, method = "average") %>% as.dendrogram()
d_bact_jacc_all <- tree_jacc_all %>% as.dendrogram()
d_bact_bc_all <- tree_bc_all %>% as.dendrogram()
d_bact_wu_all <- tree_wu_all %>% as.dendrogram()
d_bact_uu_all <- tree_uu_all %>% as.dendrogram()

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))
plot(rep(1,100),col=(colfunc(100)), pch=19,cex=2)

dendo_all <- dendlist(
  d_fish %>% 
    set("labels_col", value = colfunc(96), k=96) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_all %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value =  "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_all, 
           common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
           margin_inner=12,
           lwd=2)


# ** Triplicate species ----
load("matrices/beta_matrices.tri.RData")
load("matrices/mat.cophenetic.tri.RData")

tree_jacc_tri <- upgma(otu.jacc.phylo.tri, method = "average")
tree_bc_tri <- upgma(otu.bc.phylo.tri, method = "average")
tree_wu_tri <- upgma(otu.wu.phylo.tri, method = "average")
tree_uu_tri <- upgma(otu.uu.phylo.tri, method = "average")

RF_jacc_tri <- TreeDistance(host_tree_tri,tree_jacc_tri)
RF_bc_tri <- TreeDistance(host_tree_tri,tree_bc_tri)
RF_wu_tri <- TreeDistance(host_tree_tri,tree_wu_tri)
RF_uu_tri <- TreeDistance(host_tree_tri,tree_uu_tri)
RF_tri <-  round(c(RF_jacc_tri, RF_bc_tri,RF_wu_tri, RF_uu_tri),2)

d_fish_tri <- upgma(mat_cophenetic_tri, method = "average") %>% as.dendrogram()
d_bact_jacc_tri <- tree_jacc_tri %>% as.dendrogram()
d_bact_bc_tri <- tree_bc_tri %>% as.dendrogram()
d_bact_wu_tri <- tree_wu_tri %>% as.dendrogram()
d_bact_uu_tri <- tree_uu_tri %>% as.dendrogram()

# __ Jacc dendrograms--------
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_all <- dendlist(
  d_fish_tri %>% 
    set("labels_col", value = colfunc(length(tree_jacc_tri$tip.label)), k=length(tree_jacc_tri$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_tri %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_all, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ BC dendrograms--------
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_all <- dendlist(
  d_fish_tri %>% 
    set("labels_col", value = colfunc(length(tree_bc_tri$tip.label)), k=length(tree_bc_tri$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_tri %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_all, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms--------
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_all <- dendlist(
  d_fish_tri %>% 
    set("labels_col", value = colfunc(length(tree_wu_tri$tip.label)), k=length(tree_wu_tri$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_tri %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_all, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)
# __ UU dendrograms--------
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_all <- dendlist(
  d_fish_tri %>% 
    set("labels_col", value = colfunc(length(tree_uu_tri$tip.label)), k=length(tree_uu_tri$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_tri %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_all, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# ** Acanthuridae ----
load("matrices/beta_matrices.acanth.RData")
load("matrices/mat.cophenetic.acanth.RData")

tree_jacc_acanth <- upgma(otu.jacc.phylo.acanth, method = "average")
tree_bc_acanth <- upgma(otu.bc.phylo.acanth, method = "average")
tree_wu_acanth <- upgma(otu.wu.phylo.acanth, method = "average")
tree_uu_acanth <- upgma(otu.uu.phylo.acanth, method = "average")

RF_jacc_acanth <- TreeDistance(host_tree_acanth,tree_jacc_acanth)
RF_bc_acanth <- TreeDistance(host_tree_acanth,tree_bc_acanth)
RF_wu_acanth <- TreeDistance(host_tree_acanth,tree_wu_acanth)
RF_uu_acanth <- TreeDistance(host_tree_acanth,tree_uu_acanth)
RF_acanth <-  round(c(RF_jacc_acanth, RF_bc_acanth,RF_wu_acanth, RF_uu_acanth),2)

d_fish_acanth <- upgma(mat_cophenetic_acanth, method = "average") %>% as.dendrogram()
d_bact_jacc_acanth <- tree_jacc_acanth %>% as.dendrogram()
d_bact_bc_acanth <- tree_bc_acanth %>% as.dendrogram()
d_bact_wu_acanth <- tree_wu_acanth %>% as.dendrogram()
d_bact_uu_acanth <- tree_uu_acanth %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_acanth <- dendlist(
  d_fish_acanth %>% 
    set("labels_col", value = colfunc(length(tree_jacc_acanth$tip.label)), k=length(tree_jacc_acanth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_acanth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_acanth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)
# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_acanth <- dendlist(
  d_fish_acanth %>% 
    set("labels_col", value = colfunc(length(tree_bc_acanth$tip.label)), k=length(tree_bc_acanth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_acanth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_acanth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_acanth <- dendlist(
  d_fish_acanth %>% 
    set("labels_col", value = colfunc(length(tree_wu_acanth$tip.label)), k=length(tree_wu_acanth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_acanth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_acanth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_acanth <- dendlist(
  d_fish_acanth %>% 
    set("labels_col", value = colfunc(length(tree_uu_acanth$tip.label)), k=length(tree_uu_acanth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_acanth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_acanth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# **Holocentridae ----
load("matrices/beta_matrices.holo.RData")
load("matrices/mat.cophenetic.holo.RData")

tree_jacc_holo <- upgma(otu.jacc.phylo.holo, method = "average")
tree_bc_holo <- upgma(otu.bc.phylo.holo, method = "average")
tree_wu_holo <- upgma(otu.wu.phylo.holo, method = "average")
tree_uu_holo <- upgma(otu.uu.phylo.holo, method = "average")

RF_jacc_holo <- TreeDistance(host_tree_holo,tree_jacc_holo)
RF_bc_holo <- TreeDistance(host_tree_holo,tree_bc_holo)
RF_wu_holo <- TreeDistance(host_tree_holo,tree_wu_holo)
RF_uu_holo <- TreeDistance(host_tree_holo,tree_uu_holo)
RF_holo <-  round(c(RF_jacc_holo, RF_bc_holo,RF_wu_holo, RF_uu_holo),2)

d_fish_holo <- upgma(mat_cophenetic_holo, method = "average") %>% as.dendrogram()
d_bact_jacc_holo <- tree_jacc_holo %>% as.dendrogram()
d_bact_bc_holo <- tree_bc_holo %>% as.dendrogram()
d_bact_wu_holo <- tree_wu_holo %>% as.dendrogram()
d_bact_uu_holo <- tree_uu_holo %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_holo <- dendlist(
  d_fish_holo %>% 
    set("labels_col", value = colfunc(length(tree_jacc_holo$tip.label)), k=length(tree_jacc_holo$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_holo %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_holo, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)
# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_holo <- dendlist(
  d_fish_holo %>% 
    set("labels_col", value = colfunc(length(tree_bc_holo$tip.label)), k=length(tree_bc_holo$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_holo %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_holo, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_holo <- dendlist(
  d_fish_holo %>% 
    set("labels_col", value = colfunc(length(tree_wu_holo$tip.label)), k=length(tree_wu_holo$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_holo %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_holo, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_holo <- dendlist(
  d_fish_holo %>% 
    set("labels_col", value = colfunc(length(tree_uu_holo$tip.label)), k=length(tree_uu_holo$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_holo %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_holo, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# ** Labridae ----
load("matrices/beta_matrices.lab.RData")
load("matrices/mat.cophenetic.lab.RData")

tree_jacc_lab <- upgma(otu.jacc.phylo.lab, method = "average")
tree_bc_lab <- upgma(otu.bc.phylo.lab, method = "average")
tree_wu_lab <- upgma(otu.wu.phylo.lab, method = "average")
tree_uu_lab <- upgma(otu.uu.phylo.lab, method = "average")

RF_jacc_lab <- TreeDistance(host_tree_lab,tree_jacc_lab)
RF_bc_lab <- TreeDistance(host_tree_lab,tree_bc_lab)
RF_wu_lab <- TreeDistance(host_tree_lab,tree_wu_lab)
RF_uu_lab <- TreeDistance(host_tree_lab,tree_uu_lab)
RF_lab <-  round(c(RF_jacc_lab, RF_bc_lab,RF_wu_lab, RF_uu_lab),2)

d_fish_lab <- upgma(mat_cophenetic_lab, method = "average") %>% as.dendrogram()
d_bact_jacc_lab <- tree_jacc_lab %>% as.dendrogram()
d_bact_bc_lab <- tree_bc_lab %>% as.dendrogram()
d_bact_wu_lab <- tree_wu_lab %>% as.dendrogram()
d_bact_uu_lab <- tree_uu_lab %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lab <- dendlist(
  d_fish_lab %>% 
    set("labels_col", value = colfunc(length(tree_jacc_lab$tip.label)), k=length(tree_jacc_lab$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_lab %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lab, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lab <- dendlist(
  d_fish_lab %>% 
    set("labels_col", value = colfunc(length(tree_bc_lab$tip.label)), k=length(tree_bc_lab$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_lab %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lab, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)
# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lab <- dendlist(
  d_fish_lab %>% 
    set("labels_col", value = colfunc(length(tree_wu_lab$tip.label)), k=length(tree_wu_lab$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_lab %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lab, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lab <- dendlist(
  d_fish_lab %>% 
    set("labels_col", value = colfunc(length(tree_uu_lab$tip.label)), k=length(tree_uu_lab$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_lab %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lab, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# ** Lethrinidae ----
load("matrices/beta_matrices.leth.RData")
load("matrices/mat.cophenetic.leth.RData")

tree_jacc_leth <- upgma(otu.jacc.phylo.leth, method = "average")
tree_bc_leth <- upgma(otu.bc.phylo.leth, method = "average")
tree_wu_leth <- upgma(otu.wu.phylo.leth, method = "average")
tree_uu_leth <- upgma(otu.uu.phylo.leth, method = "average")

RF_jacc_leth <- TreeDistance(host_tree_leth,tree_jacc_leth)
RF_bc_leth <- TreeDistance(host_tree_leth,tree_bc_leth)
RF_wu_leth <- TreeDistance(host_tree_leth,tree_wu_leth)
RF_uu_leth <- TreeDistance(host_tree_leth,tree_uu_leth)
RF_leth <-  round(c(RF_jacc_leth, RF_bc_leth,RF_wu_leth, RF_uu_leth),2)

d_fish_leth <- upgma(mat_cophenetic_leth, method = "average") %>% as.dendrogram()
d_bact_jacc_leth <- tree_jacc_leth %>% as.dendrogram()
d_bact_bc_leth <- tree_bc_leth %>% as.dendrogram()
d_bact_wu_leth <- tree_wu_leth %>% as.dendrogram()
d_bact_uu_leth <- tree_uu_leth %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_leth <- dendlist(
  d_fish_leth %>% 
    set("labels_col", value = colfunc(length(tree_jacc_leth$tip.label)), k=length(tree_jacc_leth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_leth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_leth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_leth <- dendlist(
  d_fish_leth %>% 
    set("labels_col", value = colfunc(length(tree_bc_leth$tip.label)), k=length(tree_bc_leth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_leth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_leth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_leth <- dendlist(
  d_fish_leth %>% 
    set("labels_col", value = colfunc(length(tree_wu_leth$tip.label)), k=length(tree_wu_leth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_leth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_leth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_leth <- dendlist(
  d_fish_leth %>% 
    set("labels_col", value = colfunc(length(tree_uu_leth$tip.label)), k=length(tree_uu_leth$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_leth %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_leth, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# ** Lutjanidae ----
load("matrices/beta_matrices.lutj.RData")
load("matrices/mat.cophenetic.lutj.RData")

tree_jacc_lutj <- upgma(otu.jacc.phylo.lutj, method = "average")
tree_bc_lutj <- upgma(otu.bc.phylo.lutj, method = "average")
tree_wu_lutj <- upgma(otu.wu.phylo.lutj, method = "average")
tree_uu_lutj <- upgma(otu.uu.phylo.lutj, method = "average")

RF_jacc_lutj <- TreeDistance(host_tree_lutj,tree_jacc_lutj)
RF_bc_lutj <- TreeDistance(host_tree_lutj,tree_bc_lutj)
RF_wu_lutj <- TreeDistance(host_tree_lutj,tree_wu_lutj)
RF_uu_lutj <- TreeDistance(host_tree_lutj,tree_uu_lutj)
RF_lutj <-  round(c(RF_jacc_lutj, RF_bc_lutj,RF_wu_lutj, RF_uu_lutj),2)

d_fish_lutj <- upgma(mat_cophenetic_lutj, method = "average") %>% as.dendrogram()
d_bact_jacc_lutj <- tree_jacc_lutj %>% as.dendrogram()
d_bact_bc_lutj <- tree_bc_lutj %>% as.dendrogram()
d_bact_wu_lutj <- tree_wu_lutj %>% as.dendrogram()
d_bact_uu_lutj <- tree_uu_lutj %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lutj <- dendlist(
  d_fish_lutj %>% 
    set("labels_col", value = colfunc(length(tree_jacc_lutj$tip.label)), k=length(tree_jacc_lutj$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_lutj %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lutj, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lutj <- dendlist(
  d_fish_lutj %>% 
    set("labels_col", value = colfunc(length(tree_bc_lutj$tip.label)), k=length(tree_bc_lutj$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_lutj %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lutj, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lutj <- dendlist(
  d_fish_lutj %>% 
    set("labels_col", value = colfunc(length(tree_wu_lutj$tip.label)), k=length(tree_wu_lutj$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_lutj %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lutj, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_lutj <- dendlist(
  d_fish_lutj %>% 
    set("labels_col", value = colfunc(length(tree_uu_lutj$tip.label)), k=length(tree_uu_lutj$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_lutj %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_lutj, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)



# ** Pomacentridae ----
load("matrices/beta_matrices.poma.RData")
load("matrices/mat.cophenetic.poma.RData")

tree_jacc_poma <- upgma(otu.jacc.phylo.poma, method = "average")
tree_bc_poma <- upgma(otu.bc.phylo.poma, method = "average")
tree_wu_poma <- upgma(otu.wu.phylo.poma, method = "average")
tree_uu_poma <- upgma(otu.uu.phylo.poma, method = "average")

RF_jacc_poma <- TreeDistance(host_tree_poma,tree_jacc_poma)
RF_bc_poma <- TreeDistance(host_tree_poma,tree_bc_poma)
RF_wu_poma <- TreeDistance(host_tree_poma,tree_wu_poma)
RF_uu_poma <- TreeDistance(host_tree_poma,tree_uu_poma)
RF_poma <-  round(c(RF_jacc_poma, RF_bc_poma,RF_wu_poma, RF_uu_poma),2)

d_fish_poma <- upgma(mat_cophenetic_poma, method = "average") %>% as.dendrogram()
d_bact_jacc_poma <- tree_jacc_poma %>% as.dendrogram()
d_bact_bc_poma <- tree_bc_poma %>% as.dendrogram()
d_bact_wu_poma <- tree_wu_poma %>% as.dendrogram()
d_bact_uu_poma <- tree_uu_poma %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_poma <- dendlist(
  d_fish_poma %>% 
    set("labels_col", value = colfunc(length(tree_jacc_poma$tip.label)), k=length(tree_jacc_poma$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_poma %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_poma, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_poma <- dendlist(
  d_fish_poma %>% 
    set("labels_col", value = colfunc(length(tree_bc_poma$tip.label)), k=length(tree_bc_poma$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_poma %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_poma, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_poma <- dendlist(
  d_fish_poma %>% 
    set("labels_col", value = colfunc(length(tree_wu_poma$tip.label)), k=length(tree_wu_poma$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_poma %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_poma, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_poma <- dendlist(
  d_fish_poma %>% 
    set("labels_col", value = colfunc(length(tree_uu_poma$tip.label)), k=length(tree_uu_poma$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_poma %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_poma, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# ** Scaridae ----
load("matrices/beta_matrices.scar.RData")
load("matrices/mat.cophenetic.scar.RData")

tree_jacc_scar <- upgma(otu.jacc.phylo.scar, method = "average")
tree_bc_scar <- upgma(otu.bc.phylo.scar, method = "average")
tree_wu_scar <- upgma(otu.wu.phylo.scar, method = "average")
tree_uu_scar <- upgma(otu.uu.phylo.scar, method = "average")

RF_jacc_scar <- TreeDistance(host_tree_scar,tree_jacc_scar)
RF_bc_scar <- TreeDistance(host_tree_scar,tree_bc_scar)
RF_wu_scar <- TreeDistance(host_tree_scar,tree_wu_scar)
RF_uu_scar <- TreeDistance(host_tree_scar,tree_uu_scar)
RF_scar <-  round(c(RF_jacc_scar, RF_bc_scar,RF_wu_scar, RF_uu_scar),2)

d_fish_scar <- upgma(mat_cophenetic_scar, method = "average") %>% as.dendrogram()
d_bact_jacc_scar <- tree_jacc_scar %>% as.dendrogram()
d_bact_bc_scar <- tree_bc_scar %>% as.dendrogram()
d_bact_wu_scar <- tree_wu_scar %>% as.dendrogram()
d_bact_uu_scar <- tree_uu_scar %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_scar <- dendlist(
  d_fish_scar %>% 
    set("labels_col", value = colfunc(length(tree_jacc_scar$tip.label)), k=length(tree_jacc_scar$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_scar %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_scar, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)
# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_scar <- dendlist(
  d_fish_scar %>% 
    set("labels_col", value = colfunc(length(tree_bc_scar$tip.label)), k=length(tree_bc_scar$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_scar %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_scar, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_scar <- dendlist(
  d_fish_scar %>% 
    set("labels_col", value = colfunc(length(tree_wu_scar$tip.label)), k=length(tree_wu_scar$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_scar %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_scar, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_scar <- dendlist(
  d_fish_scar %>% 
    set("labels_col", value = colfunc(length(tree_uu_scar$tip.label)), k=length(tree_uu_scar$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_scar %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_scar, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# ** Herbivores ----
load("matrices/beta_matrices.herb.RData")
load("matrices/mat.cophenetic.herb.RData")

tree_jacc_herb <- upgma(otu.jacc.phylo.herb, method = "average")
tree_bc_herb <- upgma(otu.bc.phylo.herb, method = "average")
tree_wu_herb <- upgma(otu.wu.phylo.herb, method = "average")
tree_uu_herb <- upgma(otu.uu.phylo.herb, method = "average")

RF_jacc_herb <- TreeDistance(host_tree_herb,tree_jacc_herb)
RF_bc_herb <- TreeDistance(host_tree_herb,tree_bc_herb)
RF_wu_herb <- TreeDistance(host_tree_herb,tree_wu_herb)
RF_uu_herb <- TreeDistance(host_tree_herb,tree_uu_herb)
RF_herb <-  round(c(RF_jacc_herb, RF_bc_herb,RF_wu_herb, RF_uu_herb),2)

d_fish_herb <- upgma(mat_cophenetic_herb, method = "average") %>% as.dendrogram()
d_bact_jacc_herb <- tree_jacc_herb %>% as.dendrogram()
d_bact_bc_herb <- tree_bc_herb %>% as.dendrogram()
d_bact_wu_herb <- tree_wu_herb %>% as.dendrogram()
d_bact_uu_herb <- tree_uu_herb %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_jacc_herb$tip.label)), k=length(tree_jacc_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_bc_herb$tip.label)), k=length(tree_bc_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_wu_herb$tip.label)), k=length(tree_wu_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_uu_herb$tip.label)), k=length(tree_uu_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# ** Carnivores ----
load("matrices/beta_matrices.carn.RData")
load("matrices/mat.cophenetic.carn.RData")

tree_jacc_carn <- upgma(otu.jacc.phylo.carn, method = "average")
tree_bc_carn <- upgma(otu.bc.phylo.carn, method = "average")
tree_wu_carn <- upgma(otu.wu.phylo.carn, method = "average")
tree_uu_carn <- upgma(otu.uu.phylo.carn, method = "average")

RF_jacc_carn <- TreeDistance(host_tree_carn,tree_jacc_carn)
RF_bc_carn <- TreeDistance(host_tree_carn,tree_bc_carn)
RF_wu_carn <- TreeDistance(host_tree_carn,tree_wu_carn)
RF_uu_carn <- TreeDistance(host_tree_carn,tree_uu_carn)
RF_carn <-  round(c(RF_jacc_carn, RF_bc_carn,RF_wu_carn, RF_uu_carn),2)

d_fish_carn <- upgma(mat_cophenetic_carn, method = "average") %>% as.dendrogram()
d_bact_jacc_carn <- tree_jacc_carn %>% as.dendrogram()
d_bact_bc_carn <- tree_bc_carn %>% as.dendrogram()
d_bact_wu_carn <- tree_wu_carn %>% as.dendrogram()
d_bact_uu_carn <- tree_uu_carn %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_carn <- dendlist(
  d_fish_carn %>% 
    set("labels_col", value = colfunc(length(tree_jacc_carn$tip.label)), k=length(tree_jacc_carn$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_carn %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_carn, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_carn <- dendlist(
  d_fish_carn %>% 
    set("labels_col", value = colfunc(length(tree_bc_carn$tip.label)), k=length(tree_bc_carn$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_carn %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_carn, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_carn <- dendlist(
  d_fish_carn %>% 
    set("labels_col", value = colfunc(length(tree_wu_carn$tip.label)), k=length(tree_wu_carn$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_carn %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_carn, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_carn <- dendlist(
  d_fish_carn %>% 
    set("labels_col", value = colfunc(length(tree_uu_carn$tip.label)), k=length(tree_uu_carn$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_carn %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_carn, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# **Omnivores ----
load("matrices/beta_matrices.omni.RData")
load("matrices/mat.cophenetic.omni.RData")

tree_jacc_omni <- upgma(otu.jacc.phylo.omni, method = "average")
tree_bc_omni <- upgma(otu.bc.phylo.omni, method = "average")
tree_wu_omni <- upgma(otu.wu.phylo.omni, method = "average")
tree_uu_omni <- upgma(otu.uu.phylo.omni, method = "average")

RF_jacc_omni <- TreeDistance(host_tree_omni,tree_jacc_omni)
RF_bc_omni <- TreeDistance(host_tree_omni,tree_bc_omni)
RF_wu_omni <- TreeDistance(host_tree_omni,tree_wu_omni)
RF_uu_omni <- TreeDistance(host_tree_omni,tree_uu_omni)
RF_omni <-  round(c(RF_jacc_omni, RF_bc_omni,RF_wu_omni, RF_uu_omni),2)

d_fish_omni <- upgma(mat_cophenetic_omni, method = "average") %>% as.dendrogram()
d_bact_jacc_omni <- tree_jacc_omni %>% as.dendrogram()
d_bact_bc_omni <- tree_bc_omni %>% as.dendrogram()
d_bact_wu_omni <- tree_wu_omni %>% as.dendrogram()
d_bact_uu_omni <- tree_uu_omni %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_omni <- dendlist(
  d_fish_omni %>% 
    set("labels_col", value = colfunc(length(tree_jacc_omni$tip.label)), k=length(tree_jacc_omni$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_omni %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_omni, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_omni <- dendlist(
  d_fish_omni %>% 
    set("labels_col", value = colfunc(length(tree_bc_omni$tip.label)), k=length(tree_bc_omni$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_omni %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_omni, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_omni <- dendlist(
  d_fish_omni %>% 
    set("labels_col", value = colfunc(length(tree_wu_omni$tip.label)), k=length(tree_wu_omni$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_omni %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_omni, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_omni <- dendlist(
  d_fish_omni %>% 
    set("labels_col", value = colfunc(length(tree_uu_omni$tip.label)), k=length(tree_uu_omni$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_omni %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_omni, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# ** Perciformes ----
load("matrices/beta_matrices.perci.RData")
load("matrices/mat.cophenetic.perci.RData")

tree_jacc_perci <- upgma(otu.jacc.phylo.perci, method = "average")
tree_bc_perci <- upgma(otu.bc.phylo.perci, method = "average")
tree_wu_perci <- upgma(otu.wu.phylo.perci, method = "average")
tree_uu_perci <- upgma(otu.uu.phylo.perci, method = "average")

RF_jacc_perci <- TreeDistance(host_tree_perci,tree_jacc_perci)
RF_bc_perci <- TreeDistance(host_tree_perci,tree_bc_perci)
RF_wu_perci <- TreeDistance(host_tree_perci,tree_wu_perci)
RF_uu_perci <- TreeDistance(host_tree_perci,tree_uu_perci)
RF_perci <-  round(c(RF_jacc_perci, RF_bc_perci,RF_wu_perci, RF_uu_perci),2)

d_fish_perci <- upgma(mat_cophenetic_perci, method = "average") %>% as.dendrogram()
d_bact_jacc_perci <- tree_jacc_perci %>% as.dendrogram()
d_bact_bc_perci <- tree_bc_perci %>% as.dendrogram()
d_bact_wu_perci <- tree_wu_perci %>% as.dendrogram()
d_bact_uu_perci <- tree_uu_perci %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_perci <- dendlist(
  d_fish_perci %>% 
    set("labels_col", value = colfunc(length(tree_jacc_perci$tip.label)), k=length(tree_jacc_perci$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_perci %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_perci, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_perci <- dendlist(
  d_fish_perci %>% 
    set("labels_col", value = colfunc(length(tree_bc_perci$tip.label)), k=length(tree_bc_perci$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_perci %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_perci, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_perci <- dendlist(
  d_fish_perci %>% 
    set("labels_col", value = colfunc(length(tree_wu_perci$tip.label)), k=length(tree_wu_perci$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_perci %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_perci, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_perci <- dendlist(
  d_fish_perci %>% 
    set("labels_col", value = colfunc(length(tree_uu_perci$tip.label)), k=length(tree_uu_perci$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_perci %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_perci, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)

# ** Labridae and Scaridae ----
load("matrices/beta_matrices.labsca.RData")
load("matrices/mat.cophenetic.labsca.RData")

tree_jacc_labsca <- upgma(otu.jacc.phylo.labsca, method = "average")
tree_bc_labsca <- upgma(otu.bc.phylo.labsca, method = "average")
tree_wu_labsca <- upgma(otu.wu.phylo.labsca, method = "average")
tree_uu_labsca <- upgma(otu.uu.phylo.labsca, method = "average")

RF_jacc_labsca <- TreeDistance(host_tree_labsca,tree_jacc_labsca)
RF_bc_labsca <- TreeDistance(host_tree_labsca,tree_bc_labsca)
RF_wu_labsca <- TreeDistance(host_tree_labsca,tree_wu_labsca)
RF_uu_labsca <- TreeDistance(host_tree_labsca,tree_uu_labsca)
RF_labsca <-  round(c(RF_jacc_labsca, RF_bc_labsca,RF_wu_labsca, RF_uu_labsca),2)

d_fish_labsca <- upgma(mat_cophenetic_labsca, method = "average") %>% as.dendrogram()
d_bact_jacc_labsca <- tree_jacc_labsca %>% as.dendrogram()
d_bact_bc_labsca <- tree_bc_labsca %>% as.dendrogram()
d_bact_wu_labsca <- tree_wu_labsca %>% as.dendrogram()
d_bact_uu_labsca <- tree_uu_labsca %>% as.dendrogram()

# __ Jacc dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_labsca <- dendlist(
  d_fish_labsca %>% 
    set("labels_col", value = colfunc(length(tree_jacc_labsca$tip.label)), k=length(tree_jacc_labsca$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_jacc_labsca %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_labsca, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_labsca <- dendlist(
  d_fish_labsca %>% 
    set("labels_col", value = colfunc(length(tree_bc_labsca$tip.label)), k=length(tree_bc_labsca$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_bc_labsca %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_labsca, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ WU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_labsca <- dendlist(
  d_fish_labsca %>% 
    set("labels_col", value = colfunc(length(tree_wu_labsca$tip.label)), k=length(tree_wu_labsca$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_wu_labsca %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_labsca, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


# __ UU dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_labsca <- dendlist(
  d_fish_labsca %>% 
    set("labels_col", value = colfunc(length(tree_uu_labsca$tip.label)), k=length(tree_uu_labsca$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_bact_uu_labsca %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_labsca, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=12,
                        lwd=2)


### DISTRIBUTION ------

library("devtools")
install_github("WGS-TB/RFDistributionR")
library("rfdistr")
library("ape")
rfdistr::ntt_polynomial(rtree(5),6)
library("rfdistr")
library("ape")
rfdistr::polynomial(rtree(5),6)



# CTU - Firmicutes ----
dir.create("./Firmicutes/matrices", recursive = T)
setwd("./Firmicutes")
load("../core.phylosymb.tri.RData")
core.firmi <- subset_taxa(core.phylosymb.tri , Phylum == "Firmicutes")
core.firmi <- prune_samples(sample_sums(core.firmi) >0, core.firmi)
save(core.firmi, file = "core.firmi.RData")
host_firmi <- levels(sample_data(core.firmi)$tax2)
tree_firmi <-  keep.tip(host_tree, host_firmi)
save(tree_firmi, file = "tree.firmi.RData")

# __ ** Triplicate species ----
core.phylosymb.firmi.tri.rel <- transform_sample_counts(core.firmi, function(x) x / sum(x))
save(core.phylosymb.firmi.tri.rel, file = "core.phylosymb.firmi.tri.rel.RData")

otu.jacc.firmi.tri <- as.matrix(vegdist(core.phylosymb.firmi.tri.rel@otu_table, method = "jaccard"))
otu.bc.firmi.tri <- as.matrix(vegdist(core.phylosymb.firmi.tri.rel@otu_table, method = "bray"))
otu.wu.firmi.tri <- as.matrix(UniFrac(core.phylosymb.firmi.tri.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.firmi.tri <- as.matrix(UniFrac(core.phylosymb.firmi.tri.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.firmi.tri,otu.bc.firmi.tri,otu.wu.firmi.tri,otu.uu.firmi.tri , file = "matrices/beta_matrices.tri.RData")

tree_jacc_all <- upgma(otu.jacc.firmi.tri, method = "average")
tree_bc_all <- upgma(otu.bc.firmi.tri, method = "average")
tree_wu_all <- upgma(otu.wu.firmi.tri, method = "average")
tree_uu_all <- upgma(otu.uu.firmi.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_firmi,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_firmi,tree_bc_all)
RF_wu_all <- TreeDistance(tree_firmi,tree_wu_all)
RF_uu_all <- TreeDistance(tree_firmi,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

mat_cophenetic_firm <- cophenetic.phylo(tree_firmi)/2
dim(mat_cophenetic_firm) ; dim(otu.bc.firmi.tri)
mat_cophenetic_firm <- mat_cophenetic_firm[match(rownames(otu.jacc.firmi.tri), rownames(mat_cophenetic_firm)),
                                         match(colnames(otu.jacc.firmi.tri), colnames(mat_cophenetic_firm))]

save(mat_cophenetic_firm, file = "matrices/mat_cophenetic_firm")

# __ ** Herbivores ----
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
core.firmi.herb.rel <- subset_samples(core.firmi, tax2 %in% sp_herb)
tree_firmi_herb <- keep.tip(tree_firmi, sp_herb)
mat_cophenetic_herb <- cophenetic.phylo(tree_firmi)/2

otu.jacc.firmi.herb <- as.matrix(vegdist(core.firmi.herb.rel@otu_table, method = "jaccard"))
otu.bc.firmi.herb <- as.matrix(vegdist(core.firmi.herb.rel@otu_table, method = "bray"))
otu.wu.firmi.herb <- as.matrix(UniFrac(core.firmi.herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.firmi.herb <- as.matrix(UniFrac(core.firmi.herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.firmi.herb,otu.bc.firmi.herb,otu.wu.firmi.herb,otu.uu.firmi.herb , file = "matrices/beta_matrices.herb.RData")

tree_jacc_herb <- upgma(otu.jacc.firmi.herb, method = "average")
tree_bc_herb <- upgma(otu.bc.firmi.herb, method = "average")
tree_wu_herb <- upgma(otu.wu.firmi.herb, method = "average")
tree_uu_herb <- upgma(otu.uu.firmi.herb, method = "average")

RF_jacc_herb <- TreeDistance(tree_firmi,tree_jacc_herb)
RF_bc_herb <- TreeDistance(tree_firmi,tree_bc_herb)
RF_wu_herb <- TreeDistance(tree_firmi,tree_wu_herb)
RF_uu_herb <- TreeDistance(tree_firmi,tree_uu_herb)
RF_herb <-  round(c(RF_jacc_herb, RF_bc_herb,RF_wu_herb, RF_uu_herb),2)

dim(mat_cophenetic_herb) ; dim(otu.bc.firmi.herb)
mat_cophenetic_herb <- mat_cophenetic_herb[match(rownames(otu.jacc.firmi.herb), rownames(mat_cophenetic_herb)),
                                           match(colnames(otu.jacc.firmi.herb), colnames(mat_cophenetic_herb))]

save(mat_cophenetic_herb, file = "matrices/mat_cophenetic_herb")

d_fish_herb <- upgma(mat_cophenetic_herb, method = "average") %>% as.dendrogram()
d_firm_jacc_herb <- tree_jacc_herb %>% as.dendrogram()
d_firm_bc_herb <- tree_bc_herb %>% as.dendrogram()
d_firm_wu_herb <- tree_wu_herb %>% as.dendrogram()
d_firm_uu_herb <- tree_uu_herb %>% as.dendrogram()

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_bc_herb$tip.label)), k=length(tree_bc_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_firm_bc_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=15,
                        lwd=1)

pdf(file = "fig_firmicutes_herb.pdf", he = 7, wi = 10)
plot(Fig_trees)
dev.off()

# __ ** Carnivores ----
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_firmi$tip.label)]
core.firmi.carn.rel <- subset_samples(core.firmi, tax2 %in% sp_carn)
tree_firmi_carn <- keep.tip(tree_firmi, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(tree_firmi_carn)/2

otu.jacc.firmi.carn <- as.matrix(vegdist(core.firmi.carn.rel@otu_table, method = "jaccard"))
otu.bc.firmi.carn <- as.matrix(vegdist(core.firmi.carn.rel@otu_table, method = "bray"))
otu.wu.firmi.carn <- as.matrix(UniFrac(core.firmi.carn.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.firmi.carn <- as.matrix(UniFrac(core.firmi.carn.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.firmi.carn,otu.bc.firmi.carn,otu.wu.firmi.carn,otu.uu.firmi.carn , file = "matrices/beta_matrices.carn.RData")

tree_jacc_carn <- upgma(otu.jacc.firmi.carn, method = "average")
tree_bc_carn <- upgma(otu.bc.firmi.carn, method = "average")
tree_wu_carn <- upgma(otu.wu.firmi.carn, method = "average")
tree_uu_carn <- upgma(otu.uu.firmi.carn, method = "average")

RF_jacc_carn <- TreeDistance(tree_firmi_carn,tree_jacc_carn)
RF_bc_carn <- TreeDistance(tree_firmi_carn,tree_bc_carn)
RF_wu_carn <- TreeDistance(tree_firmi_carn,tree_wu_carn)
RF_uu_carn <- TreeDistance(tree_firmi_carn,tree_uu_carn)
RF_carn <-  round(c(RF_jacc_carn, RF_bc_carn,RF_wu_carn, RF_uu_carn),2)

dim(mat_cophenetic_carn) ; dim(otu.bc.firmi.carn)
mat_cophenetic_carn <- mat_cophenetic_carn[match(rownames(otu.jacc.firmi.carn), rownames(mat_cophenetic_carn)),
                                           match(colnames(otu.jacc.firmi.carn), colnames(mat_cophenetic_carn))]

save(mat_cophenetic_carn, file = "matrices/mat_cophenetic_carn")

# __ ** Omnivores ----
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_firmi$tip.label)]
core.firmi.omni.rel <- subset_samples(core.firmi, tax2 %in% sp_omni)
tree_firmi_omni <- keep.tip(tree_firmi, sp_omni)
mat_cophenetic_omni <- cophenetic.phylo(tree_firmi_omni)/2

otu.jacc.firmi.omni <- as.matrix(vegdist(core.firmi.omni.rel@otu_table, method = "jaccard"))
otu.bc.firmi.omni <- as.matrix(vegdist(core.firmi.omni.rel@otu_table, method = "bray"))
otu.wu.firmi.omni <- as.matrix(UniFrac(core.firmi.omni.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.firmi.omni <- as.matrix(UniFrac(core.firmi.omni.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.firmi.omni,otu.bc.firmi.omni,otu.wu.firmi.omni,otu.uu.firmi.omni , file = "matrices/beta_matrices.omni.RData")

tree_jacc_omni <- upgma(otu.jacc.firmi.omni, method = "average")
tree_bc_omni <- upgma(otu.bc.firmi.omni, method = "average")
tree_wu_omni <- upgma(otu.wu.firmi.omni, method = "average")
tree_uu_omni <- upgma(otu.uu.firmi.omni, method = "average")

RF_jacc_omni <- TreeDistance(tree_firmi_omni,tree_jacc_omni)
RF_bc_omni <- TreeDistance(tree_firmi_omni,tree_bc_omni)
RF_wu_omni <- TreeDistance(tree_firmi_omni,tree_wu_omni)
RF_uu_omni <- TreeDistance(tree_firmi_omni,tree_uu_omni)
RF_omni <-  round(c(RF_jacc_omni, RF_bc_omni,RF_wu_omni, RF_uu_omni),2)

dim(mat_cophenetic_omni) ; dim(otu.bc.firmi.omni)
mat_cophenetic_omni <- mat_cophenetic_omni[match(rownames(otu.jacc.firmi.omni), rownames(mat_cophenetic_omni)),
                                           match(colnames(otu.jacc.firmi.omni), colnames(mat_cophenetic_omni))]

save(mat_cophenetic_omni, file = "matrices/mat_cophenetic_omni")

# CTU - Clostridiales ----
dir.create("./Clostridiales/matrices", recursive = T)
setwd("./Clostridiales")
load("../../core.phylosymb.tri.RData")
core.phylosymb.clostri.tri.rel <- subset_taxa(core.phylosymb.tri , Order == "Clostridiales")
core.phylosymb.clostri.tri.rel <- prune_samples(sample_sums(core.phylosymb.clostri.tri.rel) >0, core.clostri)
save(core.phylosymb.clostri.tri.rel, file = "core.phylosymb.clostri.tri.rel.RData")
host_clostri <- levels(sample_data(core.phylosymb.clostri.tri.rel)$tax2)
tree_clostri <-  keep.tip(host_tree, host_clostri)
save(tree_clostri, file = "tree.clostri.RData")

# __ ** Triplicate species ----
otu.jacc.clostri.tri <- as.matrix(vegdist(core.phylosymb.clostri.tri.rel@otu_table, method = "jaccard"))
otu.bc.clostri.tri <- as.matrix(vegdist(core.phylosymb.clostri.tri.rel@otu_table, method = "bray"))
otu.wu.clostri.tri <- as.matrix(UniFrac(core.phylosymb.clostri.tri.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri.tri <- as.matrix(UniFrac(core.phylosymb.clostri.tri.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri.tri,otu.bc.clostri.tri,otu.wu.clostri.tri,otu.uu.clostri.tri , file = "matrices/beta_matrices.tri.RData")

tree_jacc_all <- upgma(otu.jacc.clostri.tri, method = "average")
tree_bc_all <- upgma(otu.bc.clostri.tri, method = "average")
tree_wu_all <- upgma(otu.wu.clostri.tri, method = "average")
tree_uu_all <- upgma(otu.uu.clostri.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_clostri,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_clostri,tree_bc_all)
RF_wu_all <- TreeDistance(tree_clostri,tree_wu_all)
RF_uu_all <- TreeDistance(tree_clostri,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

# __ ** Herbivores ----
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
core.clostri.herb.rel <- subset_samples(core.clostri, tax2 %in% sp_herb)
tree_clostri_herb <- keep.tip(tree_clostri, sp_herb)
mat_cophenetic_herb <- cophenetic.phylo(tree_clostri_herb)/2

otu.jacc.clostri.herb <- as.matrix(vegdist(core.clostri.herb.rel@otu_table, method = "jaccard"))
otu.bc.clostri.herb <- as.matrix(vegdist(core.clostri.herb.rel@otu_table, method = "bray"))
otu.wu.clostri.herb <- as.matrix(UniFrac(core.clostri.herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri.herb <- as.matrix(UniFrac(core.clostri.herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri.herb,otu.bc.clostri.herb,otu.wu.clostri.herb,otu.uu.clostri.herb , file = "matrices/beta_matrices.herb.RData")

tree_jacc_herb <- upgma(otu.jacc.clostri.herb, method = "average")
tree_bc_herb <- upgma(otu.bc.clostri.herb, method = "average")
tree_wu_herb <- upgma(otu.wu.clostri.herb, method = "average")
tree_uu_herb <- upgma(otu.uu.clostri.herb, method = "average")

RF_jacc_herb <- TreeDistance(tree_clostri_herb,tree_jacc_herb)
RF_bc_herb <- TreeDistance(tree_clostri_herb,tree_bc_herb)
RF_wu_herb <- TreeDistance(tree_clostri_herb,tree_wu_herb)
RF_uu_herb <- TreeDistance(tree_clostri_herb,tree_uu_herb)
RF_herb <-  round(c(RF_jacc_herb, RF_bc_herb,RF_wu_herb, RF_uu_herb),2)

dim(mat_cophenetic_herb) ; dim(otu.bc.clostri.herb)
mat_cophenetic_herb <- mat_cophenetic_herb[match(rownames(otu.jacc.clostri.herb), rownames(mat_cophenetic_herb)),
                                           match(colnames(otu.jacc.clostri.herb), colnames(mat_cophenetic_herb))]

save(mat_cophenetic_herb, file = "matrices/mat_cophenetic_herb")

d_fish_herb <- upgma(mat_cophenetic_herb, method = "average") %>% as.dendrogram()
d_firm_jacc_herb <- tree_jacc_herb %>% as.dendrogram()
d_firm_bc_herb <- tree_bc_herb %>% as.dendrogram()
d_firm_wu_herb <- tree_wu_herb %>% as.dendrogram()
d_firm_uu_herb <- tree_uu_herb %>% as.dendrogram()

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_bc_herb$tip.label)), k=length(tree_bc_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_firm_bc_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=15,
                        lwd=1)

pdf(file = "fig_clostridiales_herb.pdf", he = 7, wi = 10)
plot(Fig_trees)
dev.off()

# __ ** Carnivores ----
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_clostri$tip.label)]
core.clostri.carn.rel <- subset_samples(core.clostri, tax2 %in% sp_carn)
tree_clostri_carn <- keep.tip(tree_clostri, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(tree_clostri_carn)/2

otu.jacc.clostri.carn <- as.matrix(vegdist(core.clostri.carn.rel@otu_table, method = "jaccard"))
otu.bc.clostri.carn <- as.matrix(vegdist(core.clostri.carn.rel@otu_table, method = "bray"))
otu.wu.clostri.carn <- as.matrix(UniFrac(core.clostri.carn.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri.carn <- as.matrix(UniFrac(core.clostri.carn.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri.carn,otu.bc.clostri.carn,otu.wu.clostri.carn,otu.uu.clostri.carn , file = "matrices/beta_matrices.carn.RData")

tree_jacc_carn <- upgma(otu.jacc.clostri.carn, method = "average")
tree_bc_carn <- upgma(otu.bc.clostri.carn, method = "average")
tree_wu_carn <- upgma(otu.wu.clostri.carn, method = "average")
tree_uu_carn <- upgma(otu.uu.clostri.carn, method = "average")

RF_jacc_carn <- TreeDistance(tree_clostri_carn,tree_jacc_carn)
RF_bc_carn <- TreeDistance(tree_clostri_carn,tree_bc_carn)
RF_wu_carn <- TreeDistance(tree_clostri_carn,tree_wu_carn)
RF_uu_carn <- TreeDistance(tree_clostri_carn,tree_uu_carn)
RF_carn <-  round(c(RF_jacc_carn, RF_bc_carn,RF_wu_carn, RF_uu_carn),2)

dim(mat_cophenetic_carn) ; dim(otu.bc.clostri.carn)
mat_cophenetic_carn <- mat_cophenetic_carn[match(rownames(otu.jacc.clostri.carn), rownames(mat_cophenetic_carn)),
                                           match(colnames(otu.jacc.clostri.carn), colnames(mat_cophenetic_carn))]

save(mat_cophenetic_carn, file = "matrices/mat_cophenetic_carn")

# __ ** Omnivores ----
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_clostri$tip.label)]
core.clostri.omni.rel <- subset_samples(core.clostri, tax2 %in% sp_omni)
tree_clostri_omni <- keep.tip(tree_clostri, sp_omni)
mat_cophenetic_omni <- cophenetic.phylo(tree_clostri_omni)/2

otu.jacc.clostri.omni <- as.matrix(vegdist(core.clostri.omni.rel@otu_table, method = "jaccard"))
otu.bc.clostri.omni <- as.matrix(vegdist(core.clostri.omni.rel@otu_table, method = "bray"))
otu.wu.clostri.omni <- as.matrix(UniFrac(core.clostri.omni.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri.omni <- as.matrix(UniFrac(core.clostri.omni.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri.omni,otu.bc.clostri.omni,otu.wu.clostri.omni,otu.uu.clostri.omni , file = "matrices/beta_matrices.omni.RData")

tree_jacc_omni <- upgma(otu.jacc.clostri.omni, method = "average")
tree_bc_omni <- upgma(otu.bc.clostri.omni, method = "average")
tree_wu_omni <- upgma(otu.wu.clostri.omni, method = "average")
tree_uu_omni <- upgma(otu.uu.clostri.omni, method = "average")

RF_jacc_omni <- TreeDistance(tree_clostri_omni,tree_jacc_omni)
RF_bc_omni <- TreeDistance(tree_clostri_omni,tree_bc_omni)
RF_wu_omni <- TreeDistance(tree_clostri_omni,tree_wu_omni)
RF_uu_omni <- TreeDistance(tree_clostri_omni,tree_uu_omni)
RF_omni <-  round(c(RF_jacc_omni, RF_bc_omni,RF_wu_omni, RF_uu_omni),2)

dim(mat_cophenetic_omni) ; dim(otu.bc.clostri.omni)
mat_cophenetic_omni <- mat_cophenetic_omni[match(rownames(otu.jacc.clostri.omni), rownames(mat_cophenetic_omni)),
                                           match(colnames(otu.jacc.clostri.omni), colnames(mat_cophenetic_omni))]

save(mat_cophenetic_omni, file = "matrices/mat_cophenetic_omni")

# CTU - Epulopiscium ----
dir.create("../Epulopiscium/matrices", recursive = T)
setwd("../Epulopiscium")
load("../core.phylosymb.tri.RData")
core.epulo <- subset_taxa(core.phylosymb.tri , Genus == "Epulopiscium")
core.epulo <- prune_samples(sample_sums(core.epulo) >0, core.epulo)
save(core.epulo, file = "core.epulo.RData")
host_epulo <- levels(sample_data(core.epulo)$tax2)
tree_epulo <-  keep.tip(host_tree, host_epulo)
save(tree_epulo, file = "tree.epulo.RData")

# __ ** Triplicate species ----
otu.jacc.epulo.tri <- as.matrix(vegdist(core.epulo@otu_table, method = "jaccard"))
otu.bc.epulo.tri <- as.matrix(vegdist(core.epulo@otu_table, method = "bray"))
otu.wu.epulo.tri <- as.matrix(UniFrac(core.epulo, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.epulo.tri <- as.matrix(UniFrac(core.epulo, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.epulo.tri,otu.bc.epulo.tri,otu.wu.epulo.tri,otu.uu.epulo.tri , file = "matrices/beta_matrices.tri.RData")

tree_jacc_all <- upgma(otu.jacc.epulo.tri, method = "average")
tree_bc_all <- upgma(otu.bc.epulo.tri, method = "average")
tree_wu_all <- upgma(otu.wu.epulo.tri, method = "average")
tree_uu_all <- upgma(otu.uu.epulo.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_epulo,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_epulo,tree_bc_all)
RF_wu_all <- TreeDistance(tree_epulo,tree_wu_all)
RF_uu_all <- TreeDistance(tree_epulo,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

# __ ** Herbivores ----
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_epulo$tip.label)]
core.epulo.herb.rel <- subset_samples(core.phylosymb.epulo.tri.rel, tax2 %in% sp_herb)
tree_epulo_herb <- keep.tip(tree_epulo, sp_herb)
mat_cophenetic_herb <- cophenetic.phylo(tree_epulo_herb)/2

otu.jacc.epulo.herb <- as.matrix(vegdist(core.epulo.herb.rel@otu_table, method = "jaccard"))
otu.bc.epulo.herb <- as.matrix(vegdist(core.epulo.herb.rel@otu_table, method = "bray"))
otu.wu.epulo.herb <- as.matrix(UniFrac(core.epulo.herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.epulo.herb <- as.matrix(UniFrac(core.epulo.herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.epulo.herb,otu.bc.epulo.herb,otu.wu.epulo.herb,otu.uu.epulo.herb , file = "matrices/beta_matrices.herb.RData")

tree_jacc_herb <- upgma(otu.jacc.epulo.herb, method = "average")
tree_bc_herb <- upgma(otu.bc.epulo.herb, method = "average")
tree_wu_herb <- upgma(otu.wu.epulo.herb, method = "average")
tree_uu_herb <- upgma(otu.uu.epulo.herb, method = "average")

RF_jacc_herb <- TreeDistance(tree_epulo_herb,tree_jacc_herb)
RF_bc_herb <- TreeDistance(tree_epulo_herb,tree_bc_herb)
RF_wu_herb <- TreeDistance(tree_epulo_herb,tree_wu_herb)
RF_uu_herb <- TreeDistance(tree_epulo_herb,tree_uu_herb)
RF_herb <-  round(c(RF_jacc_herb, RF_bc_herb,RF_wu_herb, RF_uu_herb),2)

dim(mat_cophenetic_herb) ; dim(otu.bc.epulo.herb)
mat_cophenetic_herb <- mat_cophenetic_herb[match(rownames(otu.jacc.epulo.herb), rownames(mat_cophenetic_herb)),
                                           match(colnames(otu.jacc.epulo.herb), colnames(mat_cophenetic_herb))]

save(mat_cophenetic_herb, file = "matrices/mat_cophenetic_herb")

d_fish_herb <- upgma(mat_cophenetic_herb, method = "average") %>% as.dendrogram()
d_epulo_jacc_herb <- tree_jacc_herb %>% as.dendrogram()
d_epulo_bc_herb <- tree_bc_herb %>% as.dendrogram()
d_epulo_wu_herb <- tree_wu_herb %>% as.dendrogram()
d_epulo_uu_herb <- tree_uu_herb %>% as.dendrogram()

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_bc_herb$tip.label)), k=length(tree_bc_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_epulo_bc_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=15,
                        lwd=1)

pdf(file = "fig_epulo_herb.pdf", he = 7, wi = 10)
plot(Fig_trees)
dev.off()

# __ ** Carnivores ----
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_epulo$tip.label)]
core.epulo.carn.rel <- subset_samples(core.phylosymb.epulo.tri.rel, tax2 %in% sp_carn)
tree_epulo_carn <- keep.tip(tree_epulo, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(tree_epulo_carn)/2

otu.jacc.epulo.carn <- as.matrix(vegdist(core.epulo.carn.rel@otu_table, method = "jaccard"))
otu.bc.epulo.carn <- as.matrix(vegdist(core.epulo.carn.rel@otu_table, method = "bray"))
otu.wu.epulo.carn <- as.matrix(UniFrac(core.epulo.carn.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.epulo.carn <- as.matrix(UniFrac(core.epulo.carn.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.epulo.carn,otu.bc.epulo.carn,otu.wu.epulo.carn,otu.uu.epulo.carn , file = "matrices/beta_matrices.carn.RData")

tree_jacc_carn <- upgma(otu.jacc.epulo.carn, method = "average")
tree_bc_carn <- upgma(otu.bc.epulo.carn, method = "average")
tree_wu_carn <- upgma(otu.wu.epulo.carn, method = "average")
tree_uu_carn <- upgma(otu.uu.epulo.carn, method = "average")

RF_jacc_carn <- TreeDistance(tree_epulo_carn,tree_jacc_carn)
RF_bc_carn <- TreeDistance(tree_epulo_carn,tree_bc_carn)
RF_wu_carn <- TreeDistance(tree_epulo_carn,tree_wu_carn)
RF_uu_carn <- TreeDistance(tree_epulo_carn,tree_uu_carn)
RF_carn <-  round(c(RF_jacc_carn, RF_bc_carn,RF_wu_carn, RF_uu_carn),2)

dim(mat_cophenetic_carn) ; dim(otu.bc.epulo.carn)
mat_cophenetic_carn <- mat_cophenetic_carn[match(rownames(otu.jacc.epulo.carn), rownames(mat_cophenetic_carn)),
                                           match(colnames(otu.jacc.epulo.carn), colnames(mat_cophenetic_carn))]

save(mat_cophenetic_carn, file = "matrices/mat_cophenetic_carn")

# __ ** Omnivores ----
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_epulo$tip.label)]
core.epulo.omni.rel <- subset_samples(core.phylosymb.epulo.tri.rel, tax2 %in% sp_omni)
tree_epulo_omni <- keep.tip(tree_epulo, sp_omni)
mat_cophenetic_omni <- cophenetic.phylo(tree_epulo_omni)/2

otu.jacc.epulo.omni <- as.matrix(vegdist(core.epulo.omni.rel@otu_table, method = "jaccard"))
otu.bc.epulo.omni <- as.matrix(vegdist(core.epulo.omni.rel@otu_table, method = "bray"))
otu.wu.epulo.omni <- as.matrix(UniFrac(core.epulo.omni.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.epulo.omni <- as.matrix(UniFrac(core.epulo.omni.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.epulo.omni,otu.bc.epulo.omni,otu.wu.epulo.omni,otu.uu.epulo.omni , file = "matrices/beta_matrices.omni.RData")

tree_jacc_omni <- upgma(otu.jacc.epulo.omni, method = "average")
tree_bc_omni <- upgma(otu.bc.epulo.omni, method = "average")
tree_wu_omni <- upgma(otu.wu.epulo.omni, method = "average")
tree_uu_omni <- upgma(otu.uu.epulo.omni, method = "average")

RF_jacc_omni <- TreeDistance(tree_epulo_omni,tree_jacc_omni)
RF_bc_omni <- TreeDistance(tree_epulo_omni,tree_bc_omni)
RF_wu_omni <- TreeDistance(tree_epulo_omni,tree_wu_omni)
RF_uu_omni <- TreeDistance(tree_epulo_omni,tree_uu_omni)
RF_omni <-  round(c(RF_jacc_omni, RF_bc_omni,RF_wu_omni, RF_uu_omni),2)

dim(mat_cophenetic_omni) ; dim(otu.bc.epulo.omni)
mat_cophenetic_omni <- mat_cophenetic_omni[match(rownames(otu.jacc.epulo.omni), rownames(mat_cophenetic_omni)),
                                           match(colnames(otu.jacc.epulo.omni), colnames(mat_cophenetic_omni))]

save(mat_cophenetic_omni, file = "matrices/mat_cophenetic_omni")

d_fish_omni <- upgma(mat_cophenetic_omni, method = "average") %>% as.dendrogram()
d_epulo_jacc_omni <- tree_jacc_omni %>% as.dendrogram()
d_epulo_bc_omni <- tree_bc_omni %>% as.dendrogram()
d_epulo_wu_omni <- tree_wu_omni %>% as.dendrogram()
d_epulo_uu_omni <- tree_uu_omni %>% as.dendrogram()

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_omni <- dendlist(
  d_fish_omni %>% 
    set("labels_col", value = colfunc(length(tree_bc_omni$tip.label)), k=length(tree_bc_omni$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_epulo_bc_omni %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_omni, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=15,
                        lwd=1)

pdf(file = "fig_epulo_omni.pdf", he = 7, wi = 10)
plot(Fig_trees)
dev.off()

# CTU - Clostridium_sensu_stricto_1 ----
dir.create("../Clostridium_sensu_stricto_1/matrices", recursive = T)
setwd("../Clostridium_sensu_stricto_1")
load("../core.phylosymb.tri.RData")
core.clostri1 <- subset_taxa(core.phylosymb.tri , Genus == "Clostridium_sensu_stricto_1")
save(core.clostri1, file = "core.clostri1.RData")
host_clostri1 <- levels(sample_data(core.clostri1)$tax2)
tree_clostri1 <-  keep.tip(host_tree, host_clostri1)
save(tree_clostri1, file = "tree.clostri1.RData")

# __ ** Triplicate species ----
otu.jacc.clostri1.tri <- as.matrix(vegdist(core.clostri1@otu_table, method = "jaccard"))
otu.bc.clostri1.tri <- as.matrix(vegdist(core.clostri1@otu_table, method = "bray"))
otu.wu.clostri1.tri <- as.matrix(UniFrac(core.clostri1, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri1.tri <- as.matrix(UniFrac(core.clostri1, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri1.tri,otu.bc.clostri1.tri,otu.wu.clostri1.tri,otu.uu.clostri1.tri , file = "matrices/beta_matrices.tri.RData")
mat_cophenetic_tri <- cophenetic.phylo(tree_clostri1)/2
dim(mat_cophenetic_tri) ; dim(otu.bc.clostri1.tri)
mat_cophenetic_tri <- mat_cophenetic_tri[match(rownames(otu.jacc.clostri1.tri), rownames(mat_cophenetic_tri)),
                                         match(colnames(otu.jacc.clostri1.tri), colnames(mat_cophenetic_tri))]
save(mat_cophenetic_tri, file = "matrices/mat_cophenetic_tri.RData")

tree_jacc_all <- upgma(otu.jacc.clostri1.tri, method = "average")
tree_bc_all <- upgma(otu.bc.clostri1.tri, method = "average")
tree_wu_all <- upgma(otu.wu.clostri1.tri, method = "average")
tree_uu_all <- upgma(otu.uu.clostri1.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_clostri1,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_clostri1,tree_bc_all)
RF_wu_all <- TreeDistance(tree_clostri1,tree_wu_all)
RF_uu_all <- TreeDistance(tree_clostri1,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

# __ ** Herbivores ----
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_clostri1$tip.label)]
core.clostri1.herb.rel <- subset_samples(core.clostri1, tax2 %in% sp_herb)
tree_clostri1_herb <- keep.tip(tree_clostri1, sp_herb)
mat_cophenetic_herb <- cophenetic.phylo(tree_clostri1_herb)/2

otu.jacc.clostri1.herb <- as.matrix(vegdist(core.clostri1.herb.rel@otu_table, method = "jaccard"))
otu.bc.clostri1.herb <- as.matrix(vegdist(core.clostri1.herb.rel@otu_table, method = "bray"))
otu.wu.clostri1.herb <- as.matrix(UniFrac(core.clostri1.herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri1.herb <- as.matrix(UniFrac(core.clostri1.herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri1.herb,otu.bc.clostri1.herb,otu.wu.clostri1.herb,otu.uu.clostri1.herb , file = "matrices/beta_matrices.herb.RData")

tree_jacc_herb <- upgma(otu.jacc.clostri1.herb, method = "average")
tree_bc_herb <- upgma(otu.bc.clostri1.herb, method = "average")
tree_wu_herb <- upgma(otu.wu.clostri1.herb, method = "average")
tree_uu_herb <- upgma(otu.uu.clostri1.herb, method = "average")

RF_jacc_herb <- TreeDistance(tree_clostri1_herb,tree_jacc_herb)
RF_bc_herb <- TreeDistance(tree_clostri1_herb,tree_bc_herb)
RF_wu_herb <- TreeDistance(tree_clostri1_herb,tree_wu_herb)
RF_uu_herb <- TreeDistance(tree_clostri1_herb,tree_uu_herb)
RF_herb <-  round(c(RF_jacc_herb, RF_bc_herb,RF_wu_herb, RF_uu_herb),2)

dim(mat_cophenetic_herb) ; dim(otu.bc.clostri1.herb)
mat_cophenetic_herb <- mat_cophenetic_herb[match(rownames(otu.jacc.clostri1.herb), rownames(mat_cophenetic_herb)),
                                           match(colnames(otu.jacc.clostri1.herb), colnames(mat_cophenetic_herb))]

save(mat_cophenetic_herb, file = "matrices/mat_cophenetic_herb.RData")

d_fish_herb <- upgma(mat_cophenetic_herb, method = "average") %>% as.dendrogram()
d_clostri1_jacc_herb <- tree_jacc_herb %>% as.dendrogram()
d_clostri1_bc_herb <- tree_bc_herb %>% as.dendrogram()
d_clostri1_wu_herb <- tree_wu_herb %>% as.dendrogram()
d_clostri1_uu_herb <- tree_uu_herb %>% as.dendrogram()

# __ BC dendrograms----
colfunc<-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

dendo_herb <- dendlist(
  d_fish_herb %>% 
    set("labels_col", value = colfunc(length(tree_bc_herb$tip.label)), k=length(tree_bc_herb$tip.label)) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1),
  d_clostri1_bc_herb %>% 
    set("labels_col", value =  "grey", k=1) %>%
    set("branches_lty", 1) %>%
    set("branches_k_color", value = "grey", k = 1)
)

# Plot them together
Fig_trees <- tanglegram(dendo_herb, 
                        common_subtrees_color_lines = F, highlight_distinct_edges  = T, highlight_branches_lwd=F, 
                        margin_inner=15,
                        lwd=1)

pdf(file = "fig_clostri1_herb.pdf", he = 7, wi = 10)
plot(Fig_trees)
dev.off()

# __ ** Carnivores ----
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_clostri1$tip.label)]
core.clostri1.carn.rel <- subset_samples(core.phylosymb.clostri1.tri.rel, tax2 %in% sp_carn)
tree_clostri1_carn <- keep.tip(tree_clostri1, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(tree_clostri1_carn)/2

otu.jacc.clostri1.carn <- as.matrix(vegdist(core.clostri1.carn.rel@otu_table, method = "jaccard"))
otu.bc.clostri1.carn <- as.matrix(vegdist(core.clostri1.carn.rel@otu_table, method = "bray"))
otu.wu.clostri1.carn <- as.matrix(UniFrac(core.clostri1.carn.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri1.carn <- as.matrix(UniFrac(core.clostri1.carn.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri1.carn,otu.bc.clostri1.carn,otu.wu.clostri1.carn,otu.uu.clostri1.carn , file = "matrices/beta_matrices.carn.RData")

tree_jacc_carn <- upgma(otu.jacc.clostri1.carn, method = "average")
tree_bc_carn <- upgma(otu.bc.clostri1.carn, method = "average")
tree_wu_carn <- upgma(otu.wu.clostri1.carn, method = "average")
tree_uu_carn <- upgma(otu.uu.clostri1.carn, method = "average")

RF_jacc_carn <- TreeDistance(tree_clostri1_carn,tree_jacc_carn)
RF_bc_carn <- TreeDistance(tree_clostri1_carn,tree_bc_carn)
RF_wu_carn <- TreeDistance(tree_clostri1_carn,tree_wu_carn)
RF_uu_carn <- TreeDistance(tree_clostri1_carn,tree_uu_carn)
RF_carn <-  round(c(RF_jacc_carn, RF_bc_carn,RF_wu_carn, RF_uu_carn),2)

dim(mat_cophenetic_carn) ; dim(otu.bc.clostri1.carn)
mat_cophenetic_carn <- mat_cophenetic_carn[match(rownames(otu.jacc.clostri1.carn), rownames(mat_cophenetic_carn)),
                                           match(colnames(otu.jacc.clostri1.carn), colnames(mat_cophenetic_carn))]

save(mat_cophenetic_carn, file = "matrices/mat_cophenetic_carn")

# __ **Omnivores ----
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_clostri1$tip.label)]
core.clostri1.omni.rel <- subset_samples(core.phylosymb.clostri1.tri.rel, tax2 %in% sp_omni)
tree_clostri1_omni <- keep.tip(tree_clostri1, sp_omni)
mat_cophenetic_omni <- cophenetic.phylo(tree_clostri1_omni)/2

otu.jacc.clostri1.omni <- as.matrix(vegdist(core.clostri1.omni.rel@otu_table, method = "jaccard"))
otu.bc.clostri1.omni <- as.matrix(vegdist(core.clostri1.omni.rel@otu_table, method = "bray"))
otu.wu.clostri1.omni <- as.matrix(UniFrac(core.clostri1.omni.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.clostri1.omni <- as.matrix(UniFrac(core.clostri1.omni.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.clostri1.omni,otu.bc.clostri1.omni,otu.wu.clostri1.omni,otu.uu.clostri1.omni , file = "matrices/beta_matrices.omni.RData")

tree_jacc_omni <- upgma(otu.jacc.clostri1.omni, method = "average")
tree_bc_omni <- upgma(otu.bc.clostri1.omni, method = "average")
tree_wu_omni <- upgma(otu.wu.clostri1.omni, method = "average")
tree_uu_omni <- upgma(otu.uu.clostri1.omni, method = "average")

RF_jacc_omni <- TreeDistance(tree_clostri1_omni,tree_jacc_omni)
RF_bc_omni <- TreeDistance(tree_clostri1_omni,tree_bc_omni)
RF_wu_omni <- TreeDistance(tree_clostri1_omni,tree_wu_omni)
RF_uu_omni <- TreeDistance(tree_clostri1_omni,tree_uu_omni)
RF_omni <-  round(c(RF_jacc_omni, RF_bc_omni,RF_wu_omni, RF_uu_omni),2)

dim(mat_cophenetic_omni) ; dim(otu.bc.clostri1.omni)
mat_cophenetic_omni <- mat_cophenetic_omni[match(rownames(otu.jacc.clostri1.omni), rownames(mat_cophenetic_omni)),
                                           match(colnames(otu.jacc.clostri1.omni), colnames(mat_cophenetic_omni))]

save(mat_cophenetic_omni, file = "matrices/mat_cophenetic_omni.RData")

# CTU - Rhodobacterales ----
dir.create("../../Rhodobacterales/matrices", recursive = T)
setwd("../../Rhodobacterales")
load("../core.phylosymb.tri.RData")
core.rhodobacterales <- subset_taxa(core.phylosymb.tri , Order == "Rhodobacterales")
core.rhodobacterales <- prune_samples(sample_sums(core.rhodobacterales) >0, core.rhodobacterales)
save(core.rhodobacterales, file = "core.rhodobacterales.RData")
host_rhodobacterales <- levels(sample_data(core.rhodobacterales)$tax2)
tree_rhodobacterales <-  keep.tip(host_tree, host_rhodobacterales)
save(tree_rhodobacterales, file = "tree.rhodobacterales.RData")

# __ ** Triplicate species ----
otu.jacc.rhodobacterales.tri <- as.matrix(vegdist(core.rhodobacterales@otu_table, method = "jaccard"))
otu.bc.rhodobacterales.tri <- as.matrix(vegdist(core.rhodobacterales@otu_table, method = "bray"))
otu.wu.rhodobacterales.tri <- as.matrix(UniFrac(core.rhodobacterales, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.rhodobacterales.tri <- as.matrix(UniFrac(core.rhodobacterales, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.rhodobacterales.tri,otu.bc.rhodobacterales.tri,otu.wu.rhodobacterales.tri,otu.uu.rhodobacterales.tri , file = "matrices/beta_matrices.tri.RData")

mat_cophenetic_tri <- cophenetic.phylo(tree_rhodobacterales)/2
dim(mat_cophenetic_tri) ; dim(otu.bc.rhodobacterales.tri)
mat_cophenetic_tri <- mat_cophenetic_tri[match(rownames(otu.jacc.rhodobacterales.tri), rownames(mat_cophenetic_tri)),
                                         match(colnames(otu.jacc.rhodobacterales.tri), colnames(mat_cophenetic_tri))]
save(mat_cophenetic_tri, file = "matrices/mat_cophenetic_tri.RData")

tree_jacc_all <- upgma(otu.jacc.rhodobacterales.tri, method = "average")
tree_bc_all <- upgma(otu.bc.rhodobacterales.tri, method = "average")
tree_wu_all <- upgma(otu.wu.rhodobacterales.tri, method = "average")
tree_uu_all <- upgma(otu.uu.rhodobacterales.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_rhodobacterales,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_rhodobacterales,tree_bc_all)
RF_wu_all <- TreeDistance(tree_rhodobacterales,tree_wu_all)
RF_uu_all <- TreeDistance(tree_rhodobacterales,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

# __ **Herbivores ----
sp_herb <- names(table(samp_tri[samp_tri$diet3 == "Herbivorous",]$tax2))
sp_herb <- sp_herb[which(sp_herb %in% tree_rhodobacterales$tip.label)]
core.rhodobacterales.herb.rel <- subset_samples(core.phylosymb.rhodobacterales.tri.rel, tax2 %in% sp_herb)
tree_rhodobacterales_herb <- keep.tip(tree_rhodobacterales, sp_herb)
mat_cophenetic_herb <- cophenetic.phylo(tree_rhodobacterales_herb)/2

otu.jacc.rhodobacterales.herb <- as.matrix(vegdist(core.rhodobacterales.herb.rel@otu_table, method = "jaccard"))
otu.bc.rhodobacterales.herb <- as.matrix(vegdist(core.rhodobacterales.herb.rel@otu_table, method = "bray"))
otu.wu.rhodobacterales.herb <- as.matrix(UniFrac(core.rhodobacterales.herb.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.rhodobacterales.herb <- as.matrix(UniFrac(core.rhodobacterales.herb.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.rhodobacterales.herb,otu.bc.rhodobacterales.herb,otu.wu.rhodobacterales.herb,otu.uu.rhodobacterales.herb , file = "matrices/beta_matrices.herb.RData")

tree_jacc_herb <- upgma(otu.jacc.rhodobacterales.herb, method = "average")
tree_bc_herb <- upgma(otu.bc.rhodobacterales.herb, method = "average")
tree_wu_herb <- upgma(otu.wu.rhodobacterales.herb, method = "average")
tree_uu_herb <- upgma(otu.uu.rhodobacterales.herb, method = "average")

RF_jacc_herb <- TreeDistance(tree_rhodobacterales_herb,tree_jacc_herb)
RF_bc_herb <- TreeDistance(tree_rhodobacterales_herb,tree_bc_herb)
RF_wu_herb <- TreeDistance(tree_rhodobacterales_herb,tree_wu_herb)
RF_uu_herb <- TreeDistance(tree_rhodobacterales_herb,tree_uu_herb)
RF_herb <-  round(c(RF_jacc_herb, RF_bc_herb,RF_wu_herb, RF_uu_herb),2)

dim(mat_cophenetic_herb) ; dim(otu.bc.rhodobacterales.herb)
mat_cophenetic_herb <- mat_cophenetic_herb[match(rownames(otu.jacc.rhodobacterales.herb), rownames(mat_cophenetic_herb)),
                                           match(colnames(otu.jacc.rhodobacterales.herb), colnames(mat_cophenetic_herb))]

save(mat_cophenetic_herb, file = "matrices/mat_cophenetic_herb.RData")

# __ **Carnivores ----
sp_carn <- names(table(samp_tri[samp_tri$diet3 == "Carnivorous",]$tax2))
sp_carn <- sp_carn[which(sp_carn %in% tree_rhodobacterales$tip.label)]
core.rhodobacterales.carn.rel <- subset_samples(core.phylosymb.rhodobacterales.tri.rel, tax2 %in% sp_carn)
tree_rhodobacterales_carn <- keep.tip(tree_rhodobacterales, sp_carn)
mat_cophenetic_carn <- cophenetic.phylo(tree_rhodobacterales_carn)/2

otu.jacc.rhodobacterales.carn <- as.matrix(vegdist(core.rhodobacterales.carn.rel@otu_table, method = "jaccard"))
otu.bc.rhodobacterales.carn <- as.matrix(vegdist(core.rhodobacterales.carn.rel@otu_table, method = "bray"))
otu.wu.rhodobacterales.carn <- as.matrix(UniFrac(core.rhodobacterales.carn.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.rhodobacterales.carn <- as.matrix(UniFrac(core.rhodobacterales.carn.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.rhodobacterales.carn,otu.bc.rhodobacterales.carn,otu.wu.rhodobacterales.carn,otu.uu.rhodobacterales.carn , file = "matrices/beta_matrices.carn.RData")

tree_jacc_carn <- upgma(otu.jacc.rhodobacterales.carn, method = "average")
tree_bc_carn <- upgma(otu.bc.rhodobacterales.carn, method = "average")
tree_wu_carn <- upgma(otu.wu.rhodobacterales.carn, method = "average")
tree_uu_carn <- upgma(otu.uu.rhodobacterales.carn, method = "average")

RF_jacc_carn <- TreeDistance(tree_rhodobacterales_carn,tree_jacc_carn)
RF_bc_carn <- TreeDistance(tree_rhodobacterales_carn,tree_bc_carn)
RF_wu_carn <- TreeDistance(tree_rhodobacterales_carn,tree_wu_carn)
RF_uu_carn <- TreeDistance(tree_rhodobacterales_carn,tree_uu_carn)
RF_carn <-  round(c(RF_jacc_carn, RF_bc_carn,RF_wu_carn, RF_uu_carn),2)

dim(mat_cophenetic_carn) ; dim(otu.bc.rhodobacterales.carn)
mat_cophenetic_carn <- mat_cophenetic_carn[match(rownames(otu.jacc.rhodobacterales.carn), rownames(mat_cophenetic_carn)),
                                           match(colnames(otu.jacc.rhodobacterales.carn), colnames(mat_cophenetic_carn))]

save(mat_cophenetic_carn, file = "matrices/mat_cophenetic_carn.RData")

# __ **Omnivores ----
sp_omni <- names(table(samp_tri[samp_tri$diet3 == "Omnivorous",]$tax2))
sp_omni <- sp_omni[which(sp_omni %in% tree_rhodobacterales$tip.label)]
core.rhodobacterales.omni.rel <- subset_samples(core.rhodobacterales, tax2 %in% sp_omni)
tree_rhodobacterales_omni <- keep.tip(tree_rhodobacterales, sp_omni)
mat_cophenetic_omni <- cophenetic.phylo(tree_rhodobacterales_omni)/2

otu.jacc.rhodobacterales.omni <- as.matrix(vegdist(core.rhodobacterales.omni.rel@otu_table, method = "jaccard"))
otu.bc.rhodobacterales.omni <- as.matrix(vegdist(core.rhodobacterales.omni.rel@otu_table, method = "bray"))
otu.wu.rhodobacterales.omni <- as.matrix(UniFrac(core.rhodobacterales.omni.rel, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.rhodobacterales.omni <- as.matrix(UniFrac(core.rhodobacterales.omni.rel, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.rhodobacterales.omni,otu.bc.rhodobacterales.omni,otu.wu.rhodobacterales.omni,otu.uu.rhodobacterales.omni , file = "matrices/beta_matrices.omni.RData")

tree_jacc_omni <- upgma(otu.jacc.rhodobacterales.omni, method = "average")
tree_bc_omni <- upgma(otu.bc.rhodobacterales.omni, method = "average")
tree_wu_omni <- upgma(otu.wu.rhodobacterales.omni, method = "average")
tree_uu_omni <- upgma(otu.uu.rhodobacterales.omni, method = "average")

RF_jacc_omni <- TreeDistance(tree_rhodobacterales_omni,tree_jacc_omni)
RF_bc_omni <- TreeDistance(tree_rhodobacterales_omni,tree_bc_omni)
RF_wu_omni <- TreeDistance(tree_rhodobacterales_omni,tree_wu_omni)
RF_uu_omni <- TreeDistance(tree_rhodobacterales_omni,tree_uu_omni)
RF_omni <-  round(c(RF_jacc_omni, RF_bc_omni,RF_wu_omni, RF_uu_omni),2)

dim(mat_cophenetic_omni) ; dim(otu.bc.rhodobacterales.omni)
mat_cophenetic_omni <- mat_cophenetic_omni[match(rownames(otu.jacc.rhodobacterales.omni), rownames(mat_cophenetic_omni)),
                                           match(colnames(otu.jacc.rhodobacterales.omni), colnames(mat_cophenetic_omni))]
save(mat_cophenetic_omni, file = "matrices/mat_cophenetic_omni.RData")



# CTU - Rhodobacteraceae ----
dir.create("./rhodobacteraceae/matrices", recursive = T)
setwd("./rhodobacteraceae")
load("../core.phylosymb.tri.RData")
core.rhodobacteraceae <- subset_taxa(core.phylosymb.tri , Family == "Rhodobacteraceae")
core.rhodobacteraceae <- prune_samples(sample_sums(core.rhodobacteraceae) >0, core.rhodobacteraceae)
save(core.rhodobacteraceae, file = "core.rhodobacteraceae.RData")
host_rhodobacteraceae <- levels(sample_data(core.rhodobacteraceae)$tax2)
tree_rhodobacteraceae <-  keep.tip(host_tree, host_rhodobacteraceae)
save(tree_rhodobacteraceae, file = "tree.rhodobacteraceae.RData")

# __ ** Triplicate species ----
otu.jacc.rhodobacteraceae.tri <- as.matrix(vegdist(core.rhodobacteraceae@otu_table, method = "jaccard"))
otu.bc.rhodobacteraceae.tri <- as.matrix(vegdist(core.rhodobacteraceae@otu_table, method = "bray"))
otu.wu.rhodobacteraceae.tri <- as.matrix(UniFrac(core.rhodobacteraceae, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.rhodobacteraceae.tri <- as.matrix(UniFrac(core.rhodobacteraceae, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.rhodobacteraceae.tri,otu.bc.rhodobacteraceae.tri,otu.wu.rhodobacteraceae.tri,otu.uu.rhodobacteraceae.tri , file = "matrices/beta_matrices.tri.RData")

mat_cophenetic_tri <- cophenetic.phylo(tree_rhodobacteraceae)/2
dim(mat_cophenetic_tri) ; dim(otu.bc.rhodobacteraceae.tri)
mat_cophenetic_tri <- mat_cophenetic_tri[match(rownames(otu.jacc.rhodobacteraceae.tri), rownames(mat_cophenetic_tri)),
                                           match(colnames(otu.jacc.rhodobacteraceae.tri), colnames(mat_cophenetic_tri))]
save(mat_cophenetic_tri, file = "matrices/mat_cophenetic_tri.RData")

tree_jacc_all <- upgma(otu.jacc.rhodobacteraceae.tri, method = "average")
tree_bc_all <- upgma(otu.bc.rhodobacteraceae.tri, method = "average")
tree_wu_all <- upgma(otu.wu.rhodobacteraceae.tri, method = "average")
tree_uu_all <- upgma(otu.uu.rhodobacteraceae.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_rhodobacteraceae,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_rhodobacteraceae,tree_bc_all)
RF_wu_all <- TreeDistance(tree_rhodobacteraceae,tree_wu_all)
RF_uu_all <- TreeDistance(tree_rhodobacteraceae,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)

# CTU - Vibrionales ----
dir.create("../../Vibrionales/matrices", recursive = T)
setwd("../../Vibrionales")
load("../core.phylosymb.tri.RData")
core.vibrionales <- subset_taxa(core.phylosymb.tri , Order == "Vibrionales")
core.vibrionales <- prune_samples(sample_sums(core.vibrionales) >0, core.vibrionales)
save(core.vibrionales, file = "core.vibrionales.RData")
host_vibrionales <- levels(sample_data(core.vibrionales)$tax2)
tree_vibrionales <-  keep.tip(host_tree, host_vibrionales)
save(tree_vibrionales, file = "tree.vibrionales.RData")

# __ ** Triplicate species ----
mat_cophenetic_tri <- cophenetic.phylo(tree_vibrionales)/2
dim(mat_cophenetic_tri) ; dim(otu.bc.vibrionales.tri)
mat_cophenetic_tri <- mat_cophenetic_tri[match(rownames(otu.jacc.vibrionales.tri), rownames(mat_cophenetic_tri)),
                                         match(colnames(otu.jacc.vibrionales.tri), colnames(mat_cophenetic_tri))]
save(mat_cophenetic_tri, file = "matrices/mat_cophenetic_tri.RData")

otu.jacc.vibrionales.tri <- as.matrix(vegdist(core.vibrionales@otu_table, method = "jaccard"))
otu.bc.vibrionales.tri <- as.matrix(vegdist(core.vibrionales@otu_table, method = "bray"))
otu.wu.vibrionales.tri <- as.matrix(UniFrac(core.vibrionales, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.vibrionales.tri <- as.matrix(UniFrac(core.vibrionales, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.vibrionales.tri,otu.bc.vibrionales.tri,otu.wu.vibrionales.tri,otu.uu.vibrionales.tri , file = "matrices/beta_matrices.tri.RData")

tree_jacc_all <- upgma(otu.jacc.vibrionales.tri, method = "average")
tree_bc_all <- upgma(otu.bc.vibrionales.tri, method = "average")
tree_wu_all <- upgma(otu.wu.vibrionales.tri, method = "average")
tree_uu_all <- upgma(otu.uu.vibrionales.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_vibrionales,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_vibrionales,tree_bc_all)
RF_wu_all <- TreeDistance(tree_vibrionales,tree_wu_all)
RF_uu_all <- TreeDistance(tree_vibrionales,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)


# CTU - Rhizobiales ----
dir.create("../Rhizobiales/matrices", recursive = T)
setwd("../Rhizobiales")
load("../core.phylosymb.tri.RData")
core.rhizo <- subset_taxa(core.phylosymb.tri , Order == "Rhizobiales")
core.rhizo <- prune_samples(sample_sums(core.rhizo) >0, core.rhizo)
save(core.rhizo, file = "core.rhizo.RData")
host_rhizo <- levels(sample_data(core.rhizo)$tax2)
tree_rhizo <-  keep.tip(host_tree, host_rhizo)
save(tree_rhizo, file = "tree.rhizo.RData")

# __ ** Triplicate species ----
mat_cophenetic_tri <- cophenetic.phylo(tree_rhizo)/2
otu.jacc.rhizo.tri <- as.matrix(vegdist(core.rhizo@otu_table, method = "jaccard"))
otu.bc.rhizo.tri <- as.matrix(vegdist(core.rhizo@otu_table, method = "bray"))
otu.wu.rhizo.tri <- as.matrix(UniFrac(core.rhizo, weighted = T, normalized=F, parallel = F, fast=T))
otu.uu.rhizo.tri <- as.matrix(UniFrac(core.rhizo, weighted = F, normalized=F, parallel = F, fast=T))
save(otu.jacc.rhizo.tri,otu.bc.rhizo.tri,otu.wu.rhizo.tri,otu.uu.rhizo.tri , file = "matrices/beta_matrices.tri.RData")

dim(mat_cophenetic_tri) ; dim(otu.bc.rhizo.tri)
mat_cophenetic_tri <- mat_cophenetic_tri[match(rownames(otu.jacc.rhizo.tri), rownames(mat_cophenetic_tri)),
                                         match(colnames(otu.jacc.rhizo.tri), colnames(mat_cophenetic_tri))]
save(mat_cophenetic_tri, file = "matrices/mat_cophenetic_tri.RData")

tree_jacc_all <- upgma(otu.jacc.rhizo.tri, method = "average")
tree_bc_all <- upgma(otu.bc.rhizo.tri, method = "average")
tree_wu_all <- upgma(otu.wu.rhizo.tri, method = "average")
tree_uu_all <- upgma(otu.uu.rhizo.tri, method = "average")

RF_jacc_all <- TreeDistance(tree_rhizo,tree_jacc_all)
RF_bc_all <- TreeDistance(tree_rhizo,tree_bc_all)
RF_wu_all <- TreeDistance(tree_rhizo,tree_wu_all)
RF_uu_all <- TreeDistance(tree_rhizo,tree_uu_all)
RF_all <-  round(c(RF_jacc_all, RF_bc_all,RF_wu_all, RF_uu_all),2)








