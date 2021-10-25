### Core calculation & Descritpion ------------------------------------------
#Setup ----
knitr::opts_chunk$set(eval = FALSE)
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
inst <- cran_packages %in% installed.packages()
if (any(! inst)) {
  install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/") }
# 
inst <- github_packages %in% installed.packages()
if (any(! inst)) {
  devtools::install_github(github_packages[!inst]) }

# Load libraries
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
sessionInfo()
set.seed(1000)

setwd("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis")

#Data preparation ----
load("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/03_assign_taxonomy/ps1_515F-Y.926R.RData")
# Subset the fish and the regions
ps_gut_0 <- subset_samples(ps1 , lineage == "vertebrate")
ps_gut_0 <- prune_taxa(names(which(colSums(ps_gut_0@otu_table)>0)), ps_gut_0)
save(ps_gut_0, file = "ps_gut_0.RData")
# Exctract the organels, the eukaryotes
ps_gut_proka <- subset_taxa(ps_gut_0, Kingdom %in% c("Archaea","Bacteria"))
ps_gut_proka_1 <- subset_taxa(ps_gut_proka, Order != "Chloroplast")
ps_gut_proka_2 <- subset_taxa(ps_gut_proka_1, Family != "Mitochondria")

#Exctract the contaminants (Salter et al., 2019)
cont <- c("Acidovorax","Acinetobacter","Aquabacterium","Bacillus","Beijerinckia","Bradyrhizobium","Brevundimonas","Caulobacter",																																					
          "Chryseobacterium", "Cupriavidus", "Deinococcus", "Delftia","Devosia","Dietzia","Enhydrobacter","Flavobacterium","Hoeflea",																																					
          "Kocuria","Massilia",	"Mesorhizobium","Methylobacterium","Microbacterium","Micrococcus","Novosphingobium","Oxalobacter","Paenibacillus",																																					
          "Paracoccus","Patulibacter","Pedobacter", "Pseudomonas","Pseudoxanthomonas","Psychrobacter","Ralstonia","Roseomonas",	"Sphingobium","Sphingomonas",																																					
          "Stenotrophomonas","Streptococcus", "Undibacterium", "Variovorax")

ps_gut_cleaned <- subset_taxa(ps_gut_proka_2, !Genus %in% cont)
save(ps_gut_cleaned , file = "ps_gut_cleaned.RData")

pdf(file = "gg_gut.pdf", wi = 12, he = 7)
gg_gut <-ggrare(ps_gut_cleaned, step = 500,  label = "tax3", se = FALSE) + 
  theme(axis.title = element_text(family = "serif", size = 20), 
        axis.text = element_text(family = "serif"), axis.text.x = element_text(family = "serif",size = 16), 
        axis.text.y = element_text(family = "serif",size = 16), plot.title = element_text(family = "serif"))
dev.off()

# Core calculation ----
# __ Normalization ----
test_rff <- prune_samples(sample_sums(ps_gut_cleaned) >= 2178 , ps_gut_cleaned)
test_rff <- rarefy_even_depth(ps_gut_cleaned, sample.size = 2178)
cov_gut_rff <- goods(test_rff@otu_table@.Data)
paste0(mean(cov_gut_rff$goods)," +- " ,sd(cov_gut_rff$goods))

ps_gut_nrm <- rarefy_even_depth(ps_gut_cleaned , 2178)
save(ps_gut_nrm, file = "ps_gut_nrm.RData")
# Coverage 99,1% ± 1.6 for a normalization at 2178

# __ Core determination ----
library(labdsv)
otu_table <- ps_gut_nrm@otu_table@.Data
abuoccplot_otu <- abuocc(otu_table)
#sub_objects of abuocc objects
str(abuoccplot_otu)
# transform spc.plt vector into table in order to calculate specific richness
richness_otu <- data.frame(abuoccplot_otu$spc.plt)
# occurence of each OTU
otu_occurence <- data.frame(abuoccplot_otu$plt.spc)
mean.abun_otu <- colSums(otu_table)/otu_occurence
square_otu <- otu_table^2
ss_otu <- data.frame(colSums(square_otu))
# Variance calculation
variance_otu=ss_otu/otu_occurence-mean.abun_otu^2
disp_otu <- (variance_otu/mean.abun_otu)*otu_occurence
# IC calculation for Poisson distribution using Chi square distribution (value and formula within Zar p574)
library(epitools)
poisic_otu = pois.exact(otu_occurence, conf.level = 0.95)
gut_dstat_otu <- cbind(mean.abun_otu, disp_otu, otu_occurence, poisic_otu)
names(gut_dstat_otu) <- c("average","disp", "occurence", "x", "pt", "rate", "lower", "upper", "prob")
# Selection of core ASVs
gut_core_otu <- gut_dstat_otu[gut_dstat_otu$disp > gut_dstat_otu$upper,]
gut_core_otu <- na.exclude(gut_core_otu)
gut_tax <- data.frame(ps_gut_nrm@tax_table@.Data)
gut_core_otu$tax <- gut_tax[rownames(gut_tax) %in% row.names(gut_core_otu),6]
gut_core_otu$phylum <- gut_tax[rownames(gut_tax) %in% row.names(gut_core_otu),2]
save(gut_core_otu, file = "gut_core_otu.RData")
gut_core <- prune_taxa(rownames(gut_core_otu), ps_gut_cleaned)
save(gut_core, file = "gut_core.RData")

# Info ratio
gut_core # 2080 ASVs
sum(sample_sums(gut_core))/sum(sample_sums(ps_gut_cleaned)) #86.1% of the gut microbiome
cov_gut_core <- goods(gut_core@otu_table@.Data)
paste0(mean(cov_gut_core$goods)," +- " ,sd(cov_gut_core$goods)) # 99,9% ± 0.02

# __ Phylogenetic calibration ----
gut_names <- colnames(gut_core@otu_table)
gut_tree <- subset_taxa(gut_core, rownames(gut_core@tax_table) %in% gut_names)
gut_tree
Biostrings::writeXStringSet(gut_core@refseq, file = "gut_core.fasta")
library(ggtree)
gut_tree <- read.tree(file = "core_gut.ntree")
library(picante)
cal1<-makeChronosCalib(gut_tree, node = "root", age.min = 1, age.max = 1, interactive = FALSE, soft.bounds = FALSE) #calibration for ultrametric branchs.
gut_chronogramme<-chronos(gut_tree, lambda=0, model = "discrete", cal=cal1, quiet = FALSE, control=chronos.control(nb.rate.cat=1))
save(gut_chronogramme ,file = "gut_chronogramme.RData")

gut_core <- merge_phyloseq(gut_core, gut_tree)
save(gut_core, file = "gut_core.RData")

table.sample <- sample_data(gut_core)
table(table.sample$order)
write.table(sample_data(gut_core), file = "core.table.txt", sep = "\t" , quote = T)

# __ Figure 1 : Rarefaction curve **----
pdf(file = "gg_core.pdf", he = 7 , wi = 12)
gg_core <- ggrare(gut_core, step = 500,  label = "tax3", se = FALSE) + theme(axis.title = element_text(family = "serif", 
                                                                                                       size = 20), axis.text = element_text(family = "serif"), 
                                                                             axis.text.x = element_text(family = "serif", 
                                                                                                        size = 16), axis.text.y = element_text(family = "serif", 
                                                                                                                                               size = 16), plot.title = element_text(family = "serif"))
dev.off()
# __ Figure 2 : SAD **----
core_physeq_order <- gut_core %>%
  tax_glom(taxrank = "Order") %>%                     # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Order)

# __ Figure 3.A : Treemap ** ----
core_physeq_order$Phylum = as.character(core_physeq_order$Phylum) # Avoid error message with factor for next step
sort(table(factor(core_physeq_order$Phylum)), T)
core_physeq_order[!core_physeq_order$Phylum %in% 
                    c("Proteobacteria","Bacteroidetes","Cyanobacteria","Firmicutes",
                      "Planctomycetes", "Spirochaetes", "Verrucomicrobia","Fusobacteria","Tenericutes"),
                  which(names(core_physeq_order) == "Phylum", T)] <- "Other" #Change phylum to other for those not included in the list

group <-  core_physeq_order$Phylum
subgroup <- core_physeq_order$Order
value <- core_physeq_order$Abundance

gut_core_treemap_data=data.frame(group,subgroup,value)

library(treemap)
pdf(file = "core_treemap.pdf" , he = 7, wi = 7)
gut_core_treemap <- treemap(gut_core_treemap_data,
                            index=c("group","subgroup"), vSize = "value", type = "index",
                            fontcolor.labels=c("white","black"),
                            fontsize.labels=c(12),bg.labels=c("transparent"),
                            fontface.labels=c(2,3),
                            border.col=c("black","white"), border.lwds=c(4,2), 
                            align.labels=list(c("center", "center"),c("left", "bottom")),
                            title="Enteric Core Treemap",fontsize.title=12,
                            fontfamily.title ="serif")

dev.off()

tax_table <- gut_core@tax_table
n_phyla <- length(levels(as.factor(tax_table[,2])))
n_genus <- length(levels(as.factor(tax_table[,6])))
Proteo <- sum(core_physeq_order[core_physeq_order$Phylum == "Proteobacteria",]$Abundance)/378 * 100
Firmicutes <- sum(core_physeq_order[core_physeq_order$Phylum == "Firmicutes",]$Abundance)/378 * 100

# __ Figure 3.B : Barplot ** ----
phylum_colors <- c('Acidobacteria'='lavenderblush4',
                   'Actinobacteria'='darkblue',
                   'Armatimonadetes'='cadetblue3',
                   'Bacteroidetes'='cornflowerblue',
                   'Calditrichaeota'='azure3',
                   'Chloroflexi'='#DCE1D2',
                   'Cyanobacteria'='#DE6554',
                   'Dadabacteria'='brown2',
                   'Deferribacteres'='darkslategray1',
                   'Deinococcus-Thermus'='salmon4',
                   'Dependentiae'='sandybrown',
                   'Epsilonbacteraeota'='darkslateblue',
                   'Elusimicrobia'='plum1',
                   'Euryarchaeota'='hotpink4',
                   'Firmicutes'='brown4',
                   'Fusobacteria'='orange',
                   'Gemmatimonadetes'='darkolivegreen',
                   'Kiritimatiellaeota'='darkkhaki',
                   'Marinimicrobia_(SAR406_clade)'='darkgoldenrod4',
                   'Latescibacteria'='darkseagreen2',
                   'Lentisphaerae'='darkseagreen4',
                   'Patescibacteria'='darkturquoise',
                   'Planctomycetes'='darkslategray',
                   'Proteobacteria'='aquamarine4',
                   'Spirochaetes'='darkolivegreen3',
                   'Tenericutes'='#CA8EA7',
                   'Thaumarchaeota'='gold3',
                   'Verrucomicrobia'='darkgreen',
                   'WPS-2'='thistle2',
                   'Other'= 'black',
                   'Z-Other' = 'black')


core_physeq_order$Phylum <- str_replace_all(core_physeq_order$Phylum,c("Other" = "Z-Other"))

barplot <- core_physeq_order

barplot <- core_physeq_order[!core_physeq_order$tax2 %in% c("Scarus_sp.","Scarus_sp2","Siganus_sp.", "Cantherhines_sp"),]
write.table(barplot, file = "barplot.file.txt", sep= "\t", quote = T)
barplot <- read.table(file = "barplot.file.modif.txt" , sep = "\t", header = T)


host_tree <- read.tree(file = "phylogeny/fish_tree_modif.tre")
tax_names <- levels(factor(barplot$tax2))
missing_names <- tax_names[which(!tax_names %in% host_tree$tip.label)]
plot <- ggplot(barplot, aes(x = tax1 , y = Abundance, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  ylab("Relative Abundance (Order > 2%)") +
  xlab("Species")+
  scale_y_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Phylum Composition of Core Enteric Microbiome") +
  scale_fill_manual(values= phylum_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title = element_text(family = "serif",size = 15), 
        axis.text = element_text(size = 14), 
        axis.text.x = element_text(size = 8, family = "serif", angle = 90,vjust = 0.65, face = "italic"), 
        axis.text.y = element_text(family = "serif", size = 13), 
        plot.title = element_text(family = "serif", size = 15), 
        legend.text = element_text(size = 10,  family = "serif"), 
        legend.title = element_text(size = 12, family = "serif"))

pdf(file = "barplot_compo_by_sp.pdf", he = 7 , wi = 10)
plot
dev.off()
