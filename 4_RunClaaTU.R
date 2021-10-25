# ClaaTU ----
# We want to determine a Cladal Taxonomic Unit (CTU) meaning a group of organisms all sharing a
#  common ancestor that manifest a statistical relationship with some ecological covariate.

setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/ClaaTU")
dir.create("./ClaaTU", recursive = T)

setwd("./ClaaTU")
## ** Setup ----
remove(list = ls())
cran_packages   <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny", "microbiome",
                     "tidyverse", "miniUI", "caret", "pls", "e1071", "ggplot2", 
                     "randomForest","entropart", "vegan", "plyr", "dplyr", "here",
                     "ggrepel", "nlme", "R.utils", "gridExtra","grid", "googledrive", 
                     "googlesheets", "phangorn", "devtools", "rmarkdown", "sys",
                     "reshape2", "devtools", "PMA","structSSI","ade4", "ape",
                     "Biostrings", "igraph", "ggnetwork", "intergraph", "ips",
                     "scales", "kableExtra", "pgirmess", "treemap", "knitr","kableExtra",
                     "rstudioapi" ,"data.table","DT","pander","formatR","grDevices","svgPanZoom",
                     "RCurl","plotly","pairwiseAdonis", "stringr")
github_packages <- c("jfukuyama/phyloseqGraphTest")
bioc_packages   <- c("phyloseq", "genefilter", "impute", "dada2", "DECIPHER")
# Install CRAN packages (if not already installed)
#Some packages would be not availbale for your R version
#inst <- cran_packages %in% installed.packages()
#if (any(! inst)) {
#  install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/") }
# 
#inst <- github_packages %in% installed.packages()
#if (any(! inst)) {
#  devtools::install_github(github_packages[!inst]) }
#
#inst <- bioc_packages %in% installed.packages()
#if (any(! inst)) {
#  BiocManager::install(bioc_packages[!inst]) }
# Load libraries
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
sessionInfo()
set.seed(1000)

# We work from the rarefied core.
load("../../Physeq_objects/rff_core.RData")
# 1) We give the newick tree to ClaaTU. It decorates internal nodes with node identifiers that will be 
# used in downstream scripts. 
core_bact_tree <- rff_core@phy_tree
write.tree(core_bact_tree, file= "./core_bact_tree.tre")

#2) We give to ClaaTU the OTU counts in a text format.
#The output ("ctus.txt") is a clade counts with internal node identifiers as columns and sample IDs as rows
otu_table <- t(otu_table(rff_core))
write.table(otu_table, file = "otu_table.txt", sep = "\t", quote = F) #edit the table with #OTU<space>ID<tab> before the entry 

#3) We give the taxonomy file to ClaaTU analysis. A first step is to transform the tax_table() to get
#the appropriate format.To do that:
tax_table <- tax_table(rff_core)
write.table(tax_table, file = "tax_table.txt", sep = "\t", quote = F)
# in the terminal % bin/dada2_tax_convert.pl ClaaTU/tax_table.txt > ClaaTU/tax.txt   will create the good tax.txt with the prefix

#4) We obtain new tree with tax of the nodes (CTUs) and their statistic (if they are overrepresented or not) with a ptest
# We want to get the significance for a group
env <- rff_core %>% sample_data()
env <- cbind("ID" = rownames(env), env)
env <- write.table(env, file ="env.txt", row.names = F, sep ="\t")
env <- read.table(file ="env.txt", header= T, sep ="\t")
groupMap <- env %>% select(ID, tax1, family, diet1, diet3)
write.table(groupMap, file = "groupMap.txt", sep = "\t", quote = F, row.names = F)
#This mapping file should be tab delimited (i.e.,sample_ID tab group_id) and should not contain a header.
write.table(groupMap %>% select(ID , diet3), file = "tgMap.txt", sep = "\t", quote = F, row.names = F,col.names = F)
write.table(groupMap %>% select(ID,tax1), file = "spMap.txt", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(groupMap %>% select(ID,family), file = "famMap.txt", sep = "\t", quote = F, row.names = F, col.names = F)

######### ** End of ClaaTU-----

library(ggtree)
## __ CTUs ----
t = read.tree("./new_prepped_tree.tre")
# Read in the pValue table
pVals = read.table(file = "ptest.txt_stats.txt" , row.names = 1)
# Label the columns
colnames(pVals) = c("obsCore", "meanShuffCore", "sdShuffCore", "Z", "p")
# Select significant nodes
all_sig = rownames(pVals)[which(pVals$p < 0.05)]
N = length(t$tip.label) + length(t$node.label)
#Make a dataframe that we will attach to our tree
d = data.frame(node = 1:N, nodeName = c(t$tip.label, t$node.label),
               pValue = rep (NA, N))
rownames(d) = d$nodeName
d[allSig, "pValue"] = "sig"
#Attach the dataframe to the tree
p = ggtree(t, layout = "circular", branch.length = "none") %<+%  d
p$data
write.table(d, file = "d.all.txt", sep ="\t", row.names = F)
write.table(p$data, file = "p.data.all.txt", sep ="\t", row.names = F)

## There is 379 CTUs that are conserved in all reef fish individuals.

## __ CTus grouped ----

#Read in the groups permutation test
sp_gps = read.table("spPtest_stats.txt", row.names = 1, sep ="\t")
colnames(sp_gps) = c("obsPrev", "meanShuffPrev","sdShuffPrev","Z","p")
fam_gps = read.table("famPtest_stats.txt", row.names = 1, sep ="\t")
colnames(fam_gps) = c("obsPrev", "meanShuffPrev","sdShuffPrev","Z","p")
tg_gps = read.table("tgPtest_stats.txt", row.names = 1, sep ="\t")
colnames(tg_gps) = c("obsPrev", "meanShuffPrev","sdShuffPrev","Z","p")

#Return th group of the rownames
returnGroup = function(v){
  group = strsplit(x = v, split = "_", 2)[[1]][2]
  return(group)
}

# Return the NodeName of the rownames
returnNode = function(v){
  name = strsplit(x = v, split = "_", 2)[[1]][1]
}


# How make a group and node column
sp_gps$node = sapply(rownames(sp_gps), returnNode)
fam_gps$node = sapply(rownames(fam_gps), returnNode)
tg_gps$node = sapply(rownames(tg_gps), returnNode)

sp_gps$group = sapply(rownames(sp_gps), returnGroup)
fam_gps$group = sapply(rownames(fam_gps), returnGroup)
tg_gps$group = sapply(rownames(tg_gps), returnGroup)
write.table(sp_gps, file ="sp_gps.txt", sep = "\t", row.names = F)
write.table(fam_gps, file ="fam_gps.txt", sep = "\t", row.names = F)
write.table(tg_gps, file ="tg_gps.txt", sep = "\t", row.names = F)

# Get the significant nodes for each traits group
sp_gps_sig = sp_gps[which(sp_gps$p < 0.05),] %>% select(node, group)
fam_gps_sig = fam_gps[which(fam_gps$p < 0.05),] %>% select(node, group)
tg_gps_sig = tg_gps[which(tg_gps$p < 0.05),] %>% select(node, group)

CTUs_sp <- levels(as.factor(sp_gps_sig$node))
CTUs_fam <- levels(as.factor(fam_gps_sig$node))
CTUs_tg <- levels(as.factor(tg_gps_sig$node))


# Get off the CTUs associated to the diet 
CTUs_sp <- CTUs_sp[which(!CTUs_sp %in% CTUs_tg)]
CTUs_fam <- CTUs_fam[which(!CTUs_fam %in% CTUs_tg)]

#Make a dataframe that we will attach to our tree
N = length(t$tip.label) + length(t$node.label)
#Make a dataframe that we will attach to our tree
d= data.frame(node = c(t$tip.label, t$node.label),
              pValue = rep (NA, N))
rownames(d) = d$node
d[CTUs_sp, "pValue"] = 1
d[CTUs_fam, "pValue"] = 2

sp_gps_sig2 <- sp_gps_sig %>% filter(node %in% c(CTUs_sp))
fam_gps_sig2 <- fam_gps_sig %>% filter(node %in% c(CTUs_fam),!node %in% c(CTUs_sp))

gps_data <- left_join(d, sp_gps_sig2, by = "node")
gps_data2 <- left_join(gps_data, fam_gps_sig2, by = "node")
colnames(gps_data2)[c(3,4)] <- c("Species", "Family")

data <- gps_data2 %>% filter(pValue >=1)
# Taxonomy of the CTUs ----
taxonomy <- read.table(file = "tax.txt", sep = "\t")
colnames(taxonomy) <- c("OTU", "Kingdom", "Phylum","Class","Order","Family","Genus","Species")
taxonomy <- taxonomy %>% select(Kingdom, Phylum, Class, Order, Family, Genus)

tax_clades <- read.table(file = "tax_clades.txt", sep = "\t")
colnames(tax_clades) <- c("node", "taxonomy")
# Genus
g_clades <- tax_clades %>% filter(str_detect(taxonomy, "g__"), taxonomy != "g__NA") 
colnames(g_clades)[2] <- "Genus"
g_clades <- left_join(g_clades , taxonomy , by ="Genus")
g_clades <- g_clades %>% select(node, Kingdom, Phylum, Class, Order, Family, Genus) %>% unique()
# Family
f_clades <- tax_clades %>% filter(str_detect(taxonomy, "f__"), taxonomy != "f__NA") 
colnames(f_clades)[2] <- "Family"
f_clades <- left_join(f_clades , taxonomy , by ="Family")
f_clades <- f_clades %>% select(node, Kingdom, Phylum, Class, Order, Family) %>% unique()
f_clades <- cbind(f_clades , "Genus" = rep(NA))
# Order
o_clades <- tax_clades %>% filter(str_detect(taxonomy, "o__"), taxonomy != "o__NA") 
colnames(o_clades)[2] <- "Order"
o_clades <- left_join(o_clades , taxonomy , by ="Order")
o_clades <- o_clades %>% select(node, Kingdom, Phylum, Class, Order) %>% unique()
o_clades <- cbind(o_clades , "Family" = rep(NA),"Genus" = rep(NA))
# Class
c_clades <- tax_clades %>% filter(str_detect(taxonomy, "c__"), taxonomy != "c__NA") 
colnames(c_clades)[2] <- "Class"
c_clades <- left_join(c_clades , taxonomy , by ="Class")
c_clades <- c_clades %>% select(node, Kingdom, Phylum, Class) %>% unique()
c_clades <- cbind(c_clades ,"Order"=rep(NA),"Family" = rep(NA),"Genus" = rep(NA))
# Phylum
p_clades <- tax_clades %>% filter(str_detect(taxonomy, "p__"), taxonomy != "p__NA") 
colnames(p_clades)[2] <- "Phylum"
p_clades <- left_join(p_clades , taxonomy , by ="Phylum")
p_clades <- p_clades %>% select(node, Kingdom, Phylum) %>% unique()
p_clades <- cbind(p_clades ,"Class"= rep(NA),"Order"=rep(NA),"Family" = rep(NA),"Genus" = rep(NA))
# Kingdom
k_clades <- tax_clades %>% filter(str_detect(taxonomy, "k__"), taxonomy != "k__NA") 
colnames(k_clades)[2] <- "Kingdom"
k_clades <- left_join(k_clades , taxonomy , by ="Kingdom")
k_clades <- k_clades %>% select(node, Kingdom) %>% unique()
k_clades <- cbind(k_clades,"Phylum"= rep(NA) ,"Class"= rep(NA),"Order"=rep(NA),"Family" = rep(NA),"Genus" = rep(NA))

clades <- rbind(g_clades, f_clades, o_clades,c_clades, p_clades, k_clades)
clades$Genus <- gsub(clades$Genus, pattern = "[a-z]__", replacement = "")
clades$Family <- gsub(clades$Family, pattern = "[a-z]__", replacement = "")
clades$Order <- gsub(clades$Order, pattern = "[a-z]__", replacement = "")
clades$Class <- gsub(clades$Class, pattern = "[a-z]__", replacement = "")
clades$Phylum <- gsub(clades$Phylum, pattern = "[a-z]__", replacement = "")
clades$Kingdom <- gsub(clades$Kingdom, pattern = "[a-z]__", replacement = "")
write.table(clades, file = "clades.txt", sep="\t", row.names = F)

data <- left_join(data, clades, by = "node")
colnames(data)[c(4,10)] <- c("Host_family", "Family")
data <- data %>% select(!taxonomy)
write.table(data, file = "data.txt", sep="\t", row.names = F)

tg_data <- left_join(tg_gps_sig , clades, by ="node")
write.table(tg_data, file = "tg_data.txt", sep="\t", row.names = F)
sort(table(tg_data$Family))

# Heatmap ----
# Prepare the file 
data <- read.table(file = "data.txt", sep ="\t", header= T)
# Keep the more abundant phyla
data$Phylum[is.na(data$Phylum)] <-"z-Other"
data$Phylum <- as.factor(data$Phylum)

levels(data$Phylum)[which(!levels(data$Phylum) %in% 
                                 c("Proteobacteria","Bacteroidetes","Firmicutes","Planctomycetes","Spirochaetes",
                                   "Cyanobacteria", "Verrucomicrobia","Fusobacteria"), T)] <- "z-Other"
# Change the species level orders ----
fish_tree <- read.tree(file = '/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/phylogeny/fish_tree_modif.tre')
env <- read.table(file ="env.txt", sep = "\t", header= T)
colnames(env)[3] <- "Species"
env_sel <- env %>% select(Species, family)

data <- left_join(data, env_sel, by = "Species")

data$Species <- str_replace_all(data$Species , c(" " = "_"))
data$Species <- as.factor(data$Species)
levels(data$Species)[which(!levels(data$Species) %in% fish_tree$tip.label)]
## Add to the 'fish_tree_modif.tre' the missing species
fish_tree <- read.tree(file = 'fish_tree.tre')
levels(data$Species)[which(!levels(data$Species) %in% fish_tree$tip.label)]
host_tree <-  keep.tip(fish_tree, levels(data$Species))
levels_sp <- host_tree$tip.label
data$Species <- factor(data$Species, levels = levels_sp)

levels_fam <- data %>% arrange(Species) %>% select(family) %>% unique()
levels_fam <- levels_fam$family
data$family <- factor(data$family , levels = levels_fam)

# Change the node level orders ----
proteo <- data %>% filter(Phylum == "Proteobacteria") %>% select(node) %>% unique()
bactero <- data %>% filter(Phylum == "Bacteroidetes") %>% select(node) %>% unique()
firmi <- data %>% filter(Phylum == "Firmicutes") %>% select(node) %>% unique()
plancto <- data %>% filter(Phylum == "Planctomycetes") %>% select(node) %>% unique()
spiro <- data %>% filter(Phylum == "Spirochaetes") %>% select(node) %>% unique()
cyano <- data %>% filter(Phylum == "Cyanobacteria") %>% select(node) %>% unique()
verruco <- data %>% filter(Phylum == "Verrucomicrobia") %>% select(node) %>% unique()
fuso <- data %>% filter(Phylum == "Fusobacteria") %>% select(node) %>% unique()
other <-  data %>% filter(Phylum == "z-Other") %>% select(node) %>% unique()

levels_node <- c(proteo$node, bactero$node, firmi$node,plancto$node,spiro$node,cyano$node,verruco$node,fuso$node,other$node)
data$node <- factor(data$node, levels = levels_node)
write.table(data, file ="data_barplot.txt", sep ="\t", row.names= F)


data2 <- data %>% filter(pValue == 1)
data2 <- data %>% filter(Species)

clades <- read.table(file = "clades.txt" , sep ="\t", header= T)
CTUs_data <- left_join(data2 %>% select(node,pValue, Species, family), clades, by = "node")
write.table(CTUs_data, file = "CTUs_data.txt", sep="\t", row.names=F) # delete NA CTUs

data3 <- read.table(file = "CTUs_data.txt", sep="\t", header= T)
data3$Species <- str_replace_all(data3$Species , c("_" = " "))
levels_sp <- str_replace_all(levels_sp , c("_" = " "))
data3$Species <- factor(data3$Species, levels = levels_sp)
data3$family <- factor(data3$family, levels = levels_fam)
data3$node <- factor(data3$node, levels = levels_node)
write.table(data3, file = "data_final.txt", sep="\t", row.names=F) # delete NA CTUs

n <- length(levels(data3$family))
colfunc <-colorRampPalette(c("chocolate4","pink","red","yellow","springgreen","royalblue","darkmagenta"))

p_tile_fam <- ggplot(data3, aes(x = node, y = Species, fill = as.factor(pValue))) + 
  geom_tile(aes(color = family))+
  scale_color_manual(values= colfunc(n))+
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 3, family = "serif", angle = 90,vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 8))

pdf(file = "heatmap_species_col_fam.pdf", he = 12 , wi = 7)
p_tile_fam
dev.off()

p_tile <- ggplot(data3, aes(x = node, y = Species, fill = as.factor(pValue))) + 
  geom_tile()+
  scale_fill_manual(values= "black")+
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 3, family = "serif", angle = 90,vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 8))

pdf(file = "heatmap_species.pdf", he = 12 , wi = 7)
p_tile
dev.off()


tax_CTUs_final <- data3 %>% select(!c(pValue,Species,family)) %>% unique()
write.table(tax_CTUs_final, file ="tax_CTUs_final.txt", sep ="\t", row.names=F)
