# Phylum composition -----
library(dplyr)

gut_core_phyla <- gut_core %>%
  tax_glom(taxrank = "Family") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)      

save(gut_core_phyla , file = "gut_core_phyla.RData")
load("gut_core_phyla.RData")

gut_core_phyla$Phylum = as.character(gut_core_phyla$Phylum) # Avoid error message with factor for next step
sort(table(factor(gut_core_phyla$Phylum)), T)
gut_core_phyla[!gut_core_phyla$Phylum %in% c("Proteobacteria","Bacteroidetes","Cyanobacteria","Firmicutes","Planctomycetes", "Spirochaetes", "Verrucomicrobia","Fusobacteria","Tenericutes"),which(names(gut_core_phyla) == "Phylum", T)] <- "z-Other"

samp_data <- sample_data(gut_core)
gut_core_phyla_fam <- subset(gut_core_phyla, family %in% names(which(table(samp_data$family) >= 3)))
gut_core_phyla_fam <- left_join(gut_core_phyla_fam , box3_fam[,c(3,7)], by = "tax3")
gut_core_phyla_fam$family.x <- gut_core_phyla_fam$family.y
gut_core_phyla_fam <- gut_core_phyla_fam[,-29]
colnames(gut_core_phyla_fam)[10] <- "family"
save(gut_core_phyla_fam, file = "gut_core_phyla_fam.RData")
write.table(gut_core_phyla_fam, file = "gut_core_fam.txt", sep = "\t" , quote = F)
gut_core_phyla_fam <- read.table(file= "gut_core_fam.txt", header = T, sep = "\t")

gut_core_phyla_sp <- subset(gut_core_phyla, tax1 %in% names(which(table(samp_data$tax1) >= 3)))
gut_core_phyla_sp <- subset(gut_core_phyla_sp, !tax1 %in% c("Scarus sp2", "Siganus sp."))

phylum_colors <- c('Bacteroidetes'='darkblue',
                   'Cyanobacteria'='#DE6554',
                   'Firmicutes'='brown4',
                   'Fusobacteria'='orange',
                   'Planctomycetes'='darkslategray',
                   'Proteobacteria'='aquamarine4',
                   'Spirochaetes'='darkolivegreen3',
                   'Tenericutes'='#CA8EA7',
                   'Verrucomicrobia'='yellow',
                   'z-Other' = 'black')

pdf(file = "alpha_fam.pdf", he = 7 , wi = 7)
ggplot(gut_core_phyla_fam, aes(x = Abundance , y = family, fill = Phylum)) + 
  geom_bar(stat="identity", position="fill") + 
  ylab("Family") +
  xlab("Relative Abundance (Order > 2%)")+
  scale_x_continuous(expand = c(0,0)) + #remove the space below the 0 of the y axis in the graph
  ggtitle("Phylum Composition of Core Enteric Microbiome") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +
  theme(axis.line = element_line(size = 0), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(), 
        axis.title = element_text(family = "serif",size = 15), 
        axis.text = element_text(size = 14), 
        axis.text.x = element_text(size = 14, family = "serif",vjust = 0.65), 
        axis.text.y = element_text(family = "serif", size = 13, face = "italic"), 
        plot.title = element_text(family = "serif", size = 15), 
        legend.text = element_text(size = 14,  family = "serif"), 
        legend.title = element_text(size = 16, family = "serif"))
dev.off()

### A very shiny archeome ----
setwd( "/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/Composition")
dir.create("./Archaeome", recursive = T)

gut_core <- load("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/Physeq_objects/gut_core.RData")
archaeome <- subset_taxa(gut_core, Kingdom == "Archaea")
archaeome <- prune_samples(sample_sums(archaeome) >0, archaeome)
save(archaeome, file= "archaeome.RData")
tax_table(archaeome) # which taxa?
sort(sample_sums(archaeome)) #How much?
df_arch <- sample_data(archaeome) #which samples?
table(df_arch$family)

write.table(df_arch, file = "core_archaeome_samples.txt", sep = "\t", quote = F)
