# PERMANOVA---- 
#  All dataset -----
core_rel <- transform_sample_counts(gut_core, function(x) x / sum(x) )
samp_data <- data.frame(sample_data(core_rel))
save(core_rel, file = "core_rel.RData")
save(samp_data, file = "samp_data.RData")
# ** Bray-curtis ----
otu.bc.all <- vegdist(core_rel@otu_table, method = "bray")

a.bc = adonis(otu.bc.all ~ family, data = samp_data)
b.bc = adonis(otu.bc.all ~ genus, data = samp_data)
c.bc = adonis(otu.bc.all ~ tax1, data = samp_data)
d.bc = adonis(otu.bc.all ~ diet3, data = samp_data)
e.bc = adonis(otu.bc.all ~ diet1, data = samp_data)
f.bc = adonis(otu.bc.all ~ region, data = samp_data)
#g.bc = adonis(otu.bc.all ~ diet3, data = samp_data, strata = samp_data$family)
#h.bc = adonis(otu.bc.all ~ diet1, data = samp_data, strata = samp_data$family)
#i.bc = adonis(otu.bc.all ~ diet1, data = samp_data, strata = samp_data$diet3)
#j.bc = adonis(otu.bc.all ~ region, data = samp_data, strata = samp_data$family)
#k.bc = adonis(otu.bc.all ~ tax1, data = samp_data, strata = samp_data$region)


results.bc.all <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                    "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                    "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                    "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                    "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                    "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                    #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                    #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                    #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                    #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                    #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.all <- UniFrac(core_rel, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.all ~ family, data = samp_data)
b.wu = adonis(otu.wu.all ~ genus, data = samp_data)
c.wu = adonis(otu.wu.all ~ tax1, data = samp_data)
d.wu = adonis(otu.wu.all ~ diet3, data = samp_data)
e.wu = adonis(otu.wu.all ~ diet1, data = samp_data)
f.wu = adonis(otu.wu.all ~ region, data = samp_data)
#g.wu = adonis(otu.wu.all ~ diet3, data = samp_data, strata = samp_data$family)
#h.wu = adonis(otu.wu.all ~ diet1, data = samp_data, strata = samp_data$family)
#i.wu = adonis(otu.wu.all ~ diet1, data = samp_data, strata = samp_data$diet3)
#j.wu = adonis(otu.wu.all ~ region, data = samp_data, strata = samp_data$family)
#k.wu = adonis(otu.wu.all ~ tax1, data = samp_data, strata = samp_data$region)

results.wu.all <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                        "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                        "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                        "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                        #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)
# ** Unweighted Unifrac ----
otu.uu.all <- UniFrac(core_rel, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.all ~ family, data = samp_data)
b.uu = adonis(otu.uu.all ~ genus, data = samp_data)
c.uu = adonis(otu.uu.all ~ tax1, data = samp_data)
d.uu = adonis(otu.uu.all ~ diet3, data = samp_data)
e.uu = adonis(otu.uu.all ~ diet1, data = samp_data)
f.uu = adonis(otu.uu.all ~ region, data = samp_data)
#g.uu = adonis(otu.uu.all ~ diet3, data = samp_data, strata = samp_data$family)
#h.uu = adonis(otu.uu.all ~ diet1, data = samp_data, strata = samp_data$family)
#i.uu = adonis(otu.uu.all ~ diet1, data = samp_data, strata = samp_data$diet3)
#j.uu = adonis(otu.uu.all ~ region, data = samp_data, strata = samp_data$family)
#k.uu = adonis(otu.uu.all ~ tax1, data = samp_data, strata = samp_data$region)

results.uu.all <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                        "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                        "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                        "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                        #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  > 3 ind/sp -----
names_tri <- rownames(which(table(box3_sp$tax1)>0, T))
gut_core_tri <- subset_samples(gut_core, tax1 %in% names_tri)
core_rel_tri <- transform_sample_counts(gut_core_tri, function(x) x / sum(x) )
samp_data_tri <- data.frame(sample_data(core_rel_tri))

# ** Bray-curtis ----
otu.bc.tri <- vegdist(core_rel_tri@otu_table, method = "bray")

a.bc = adonis(otu.bc.tri ~ family, data = samp_data_tri)
b.bc = adonis(otu.bc.tri ~ genus, data = samp_data_tri)
c.bc = adonis(otu.bc.tri ~ tax1, data = samp_data_tri)
d.bc = adonis(otu.bc.tri ~ diet3, data = samp_data_tri)
e.bc = adonis(otu.bc.tri ~ diet1, data = samp_data_tri)
f.bc = adonis(otu.bc.tri ~ region, data = samp_data_tri)
#g.bc = adonis(otu.bc.tri ~ diet3, data = samp_data_tri, strata = samp_data_tri$family)
#h.bc = adonis(otu.bc.tri ~ diet1, data = samp_data_tri, strata = samp_data_tri$family)
#i.bc = adonis(otu.bc.tri ~ diet1, data = samp_data_tri, strata = samp_data_tri$diet3)
#j.bc = adonis(otu.bc.tri ~ region, data = samp_data_tri, strata = samp_data_tri$family)
#k.bc = adonis(otu.bc.tri ~ tax1, data = samp_data_tri, strata = samp_data_tri$region)

results.bc.tri <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                        "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                        "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                        "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                        #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                        #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.tri <- UniFrac(core_rel_tri, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.tri ~ family, data = samp_data_tri)
b.wu = adonis(otu.wu.tri ~ genus, data = samp_data_tri)
c.wu = adonis(otu.wu.tri ~ tax1, data = samp_data_tri)
d.wu = adonis(otu.wu.tri ~ diet3, data = samp_data_tri)
e.wu = adonis(otu.wu.tri ~ diet1, data = samp_data_tri)
f.wu = adonis(otu.wu.tri ~ region, data = samp_data_tri)
#g.wu = adonis(otu.wu.tri ~ diet3, data = samp_data_tri, strata = samp_data_tri$family)
#h.wu = adonis(otu.wu.tri ~ diet1, data = samp_data_tri, strata = samp_data_tri$family)
#i.wu = adonis(otu.wu.tri ~ diet1, data = samp_data_tri, strata = samp_data_tri$diet3)
#j.wu = adonis(otu.wu.tri ~ region, data = samp_data_tri, strata = samp_data_tri$family)
#k.wu = adonis(otu.wu.tri ~ tax1, data = samp_data_tri, strata = samp_data_tri$region)

results.wu.tri <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                        "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                        "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                        "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                        #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                        #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)
# ** Unweighted Unifrac ----
otu.uu.tri <- UniFrac(core_rel_tri, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.tri ~ family, data = samp_data_tri)
b.uu = adonis(otu.uu.tri ~ genus, data = samp_data_tri)
c.uu = adonis(otu.uu.tri ~ tax1, data = samp_data_tri)
d.uu = adonis(otu.uu.tri ~ diet3, data = samp_data_tri)
e.uu = adonis(otu.uu.tri ~ diet1, data = samp_data_tri)
f.uu = adonis(otu.uu.tri ~ region, data = samp_data_tri)
#g.uu = adonis(otu.uu.tri ~ diet3, data = samp_data_tri, strata = samp_data_tri$family)
#h.uu = adonis(otu.uu.tri ~ diet1, data = samp_data_tri, strata = samp_data_tri$family)
#i.uu = adonis(otu.uu.tri ~ diet1, data = samp_data_tri, strata = samp_data_tri$diet3)
#j.uu = adonis(otu.uu.tri ~ region, data = samp_data_tri, strata = samp_data_tri$family)
#k.uu = adonis(otu.uu.tri ~ tax1, data = samp_data_tri, strata = samp_data_tri$region)

results.uu.tri <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                        "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                        "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                        "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                        "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                        #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                        #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


results_all <- rbind(cbind(results.bc.all, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("All_dataset")),
      cbind(results.wu.all, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("All_dataset")),
      cbind(results.uu.all, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("All_dataset")))
results_tri <- rbind(cbind(results.bc.tri, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("<3 ind/sp")),
                     cbind(results.wu.tri, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("<3 ind/sp")),
                     cbind(results.uu.tri, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("<3 ind/sp")))
results_full_dataset <- rbind(results_all, results_tri)
write.table(results_full_dataset, file = "results_full_dataset.txt", sep = '\t', quote = F)
save(otu.uu.all,otu.bc.all,otu.wu.all, file = "matrices.all.RData")
save(otu.uu.tri,otu.bc.tri,otu.wu.tri, file = "matrices.tri.RData")
save(samp_data, core_rel, file = "samp_data_core_rel.RData")
save(samp_data_tri, core_rel_tri, file = "samp_data_core_rel_tri.RData")

#####################_______ DIET ______________#######------
#  Herbivores -----
gut_core_herb <- subset_samples(gut_core, diet3 == "Herbivorous")
core_rel_herb <- transform_sample_counts(gut_core_herb, function(x) x / sum(x) )
samp_data_herb <- data.frame(sample_data(core_rel_herb))

# ** Bray-curtis ----
otu.bc.herb <- vegdist(core_rel_herb@otu_table, method = "bray")

a.bc = adonis(otu.bc.herb ~ family, data = samp_data_herb)
b.bc = adonis(otu.bc.herb ~ genus, data = samp_data_herb)
c.bc = adonis(otu.bc.herb ~ tax1, data = samp_data_herb)
#d.bc = adonis(otu.bc.herb ~ diet3, data = samp_data_herb)
e.bc = adonis(otu.bc.herb ~ diet1, data = samp_data_herb)
f.bc = adonis(otu.bc.herb ~ region, data = samp_data_herb)
#g.bc = adonis(otu.bc.herb ~ diet3, data = samp_data_herb, strata = samp_data_herb$family)
#h.bc = adonis(otu.bc.herb ~ diet1, data = samp_data_herb, strata = samp_data_herb$family)
#i.bc = adonis(otu.bc.herb ~ diet1, data = samp_data_herb, strata = samp_data_herb$diet3)
#j.bc = adonis(otu.bc.herb ~ region, data = samp_data_herb, strata = samp_data_herb$family)
#k.bc = adonis(otu.bc.herb ~ tax1, data = samp_data_herb, strata = samp_data_herb$region)

results.bc.herb <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.herb <- UniFrac(core_rel_herb, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.herb ~ family, data = samp_data_herb)
b.wu = adonis(otu.wu.herb ~ genus, data = samp_data_herb)
c.wu = adonis(otu.wu.herb ~ tax1, data = samp_data_herb)
#d.wu = adonis(otu.wu.herb ~ diet3, data = samp_data_herb)
e.wu = adonis(otu.wu.herb ~ diet1, data = samp_data_herb)
f.wu = adonis(otu.wu.herb ~ region, data = samp_data_herb)
#g.wu = adonis(otu.wu.herb ~ diet3, data = samp_data_herb, strata = samp_data_herb$family)
#h.wu = adonis(otu.wu.herb ~ diet1, data = samp_data_herb, strata = samp_data_herb$family)
#i.wu = adonis(otu.wu.herb ~ diet1, data = samp_data_herb, strata = samp_data_herb$diet3)
#j.wu = adonis(otu.wu.herb ~ region, data = samp_data_herb, strata = samp_data_herb$family)
#k.wu = adonis(otu.wu.herb ~ tax1, data = samp_data_herb, strata = samp_data_herb$region)

results.wu.herb <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.herb <- UniFrac(core_rel_herb, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.herb ~ family, data = samp_data_herb)
b.uu = adonis(otu.uu.herb ~ genus, data = samp_data_herb)
c.uu = adonis(otu.uu.herb ~ tax1, data = samp_data_herb)
#d.uu = adonis(otu.uu.herb ~ diet3, data = samp_data_herb)
e.uu = adonis(otu.uu.herb ~ diet1, data = samp_data_herb)
f.uu = adonis(otu.uu.herb ~ region, data = samp_data_herb)
#g.uu = adonis(otu.uu.herb ~ diet3, data = samp_data_herb, strata = samp_data_herb$family)
#h.uu = adonis(otu.uu.herb ~ diet1, data = samp_data_herb, strata = samp_data_herb$family)
#i.uu = adonis(otu.uu.herb ~ diet1, data = samp_data_herb, strata = samp_data_herb$diet3)
#j.uu = adonis(otu.uu.herb ~ region, data = samp_data_herb, strata = samp_data_herb$family)
#k.uu = adonis(otu.uu.herb ~ tax1, data = samp_data_herb, strata = samp_data_herb$region)

results.uu.herb <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  HD -----
gut_core_HD <- subset_samples(gut_core, diet1 == "HD")
core_rel_HD <- transform_sample_counts(gut_core_HD, function(x) x / sum(x) )
samp_data_HD <- data.frame(sample_data(core_rel_HD))

# ** Bray-curtis ----
otu.bc.HD <- vegdist(core_rel_HD@otu_table, method = "bray")

a.bc = adonis(otu.bc.HD ~ family, data = samp_data_HD)
b.bc = adonis(otu.bc.HD ~ genus, data = samp_data_HD)
c.bc = adonis(otu.bc.HD ~ tax1, data = samp_data_HD)
#d.bc = adonis(otu.bc.HD ~ diet3, data = samp_data_HD)
#e.bc = adonis(otu.bc.HD ~ diet1, data = samp_data_HD)
f.bc = adonis(otu.bc.HD ~ region, data = samp_data_HD)
#g.bc = adonis(otu.bc.HD ~ diet3, data = samp_data_HD, strata = samp_data_HD$family)
#h.bc = adonis(otu.bc.HD ~ diet1, data = samp_data_HD, strata = samp_data_HD$family)
#i.bc = adonis(otu.bc.HD ~ diet1, data = samp_data_HD, strata = samp_data_HD$diet3)
#j.bc = adonis(otu.bc.HD ~ region, data = samp_data_HD, strata = samp_data_HD$family)
#k.bc = adonis(otu.bc.HD ~ tax1, data = samp_data_HD, strata = samp_data_HD$region)

results.bc.HD <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.HD <- UniFrac(core_rel_HD, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.HD ~ family, data = samp_data_HD)
b.wu = adonis(otu.wu.HD ~ genus, data = samp_data_HD)
c.wu = adonis(otu.wu.HD ~ tax1, data = samp_data_HD)
#d.wu = adonis(otu.wu.HD ~ diet3, data = samp_data_HD)
#e.wu = adonis(otu.wu.HD ~ diet1, data = samp_data_HD)
f.wu = adonis(otu.wu.HD ~ region, data = samp_data_HD)
#g.wu = adonis(otu.wu.HD ~ diet3, data = samp_data_HD, strata = samp_data_HD$family)
#h.wu = adonis(otu.wu.HD ~ diet1, data = samp_data_HD, strata = samp_data_HD$family)
#i.wu = adonis(otu.wu.HD ~ diet1, data = samp_data_HD, strata = samp_data_HD$diet3)
#j.wu = adonis(otu.wu.HD ~ region, data = samp_data_HD, strata = samp_data_HD$family)
#k.wu = adonis(otu.wu.HD ~ tax1, data = samp_data_HD, strata = samp_data_HD$region)

results.wu.HD <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.HD <- UniFrac(core_rel_HD, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.HD ~ family, data = samp_data_HD)
b.uu = adonis(otu.uu.HD ~ genus, data = samp_data_HD)
c.uu = adonis(otu.uu.HD ~ tax1, data = samp_data_HD)
#d.uu = adonis(otu.uu.HD ~ diet3, data = samp_data_HD)
#e.uu = adonis(otu.uu.HD ~ diet1, data = samp_data_HD)
f.uu = adonis(otu.uu.HD ~ region, data = samp_data_HD)
#g.uu = adonis(otu.uu.HD ~ diet3, data = samp_data_HD, strata = samp_data_HD$family)
#h.uu = adonis(otu.uu.HD ~ diet1, data = samp_data_HD, strata = samp_data_HD$family)
#i.uu = adonis(otu.uu.HD ~ diet1, data = samp_data_HD, strata = samp_data_HD$diet3)
#j.uu = adonis(otu.uu.HD ~ region, data = samp_data_HD, strata = samp_data_HD$family)
#k.uu = adonis(otu.uu.HD ~ tax1, data = samp_data_HD, strata = samp_data_HD$region)

results.uu.HD <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

#  H -----
gut_core_H <- subset_samples(gut_core, diet1 == "H")
core_rel_H <- transform_sample_counts(gut_core_H, function(x) x / sum(x) )
samp_data_H <- data.frame(sample_data(core_rel_H))

# ** Bray-curtis ----
otu.bc.H <- vegdist(core_rel_H@otu_table, method = "bray")

a.bc = adonis(otu.bc.H ~ family, data = samp_data_H)
b.bc = adonis(otu.bc.H ~ genus, data = samp_data_H)
c.bc = adonis(otu.bc.H ~ tax1, data = samp_data_H)
#d.bc = adonis(otu.bc.H ~ diet3, data = samp_data_H)
#e.bc = adonis(otu.bc.H ~ diet1, data = samp_data_H)
f.bc = adonis(otu.bc.H ~ region, data = samp_data_H)
#g.bc = adonis(otu.bc.H ~ diet3, data = samp_data_H, strata = samp_data_H$family)
#h.bc = adonis(otu.bc.H ~ diet1, data = samp_data_H, strata = samp_data_H$family)
#i.bc = adonis(otu.bc.H ~ diet1, data = samp_data_H, strata = samp_data_H$diet3)
#j.bc = adonis(otu.bc.H ~ region, data = samp_data_H, strata = samp_data_H$family)
#k.bc = adonis(otu.bc.H ~ tax1, data = samp_data_H, strata = samp_data_H$region)

results.bc.H <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.H <- UniFrac(core_rel_H, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.H ~ family, data = samp_data_H)
b.wu = adonis(otu.wu.H ~ genus, data = samp_data_H)
c.wu = adonis(otu.wu.H ~ tax1, data = samp_data_H)
#d.wu = adonis(otu.wu.H ~ diet3, data = samp_data_H)
#e.wu = adonis(otu.wu.H ~ diet1, data = samp_data_H)
f.wu = adonis(otu.wu.H ~ region, data = samp_data_H)
#g.wu = adonis(otu.wu.H ~ diet3, data = samp_data_H, strata = samp_data_H$family)
#h.wu = adonis(otu.wu.H ~ diet1, data = samp_data_H, strata = samp_data_H$family)
#i.wu = adonis(otu.wu.H ~ diet1, data = samp_data_H, strata = samp_data_H$diet3)
#j.wu = adonis(otu.wu.H ~ region, data = samp_data_H, strata = samp_data_H$family)
#k.wu = adonis(otu.wu.H ~ tax1, data = samp_data_H, strata = samp_data_H$region)

results.wu.H <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.H <- UniFrac(core_rel_H, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.H ~ family, data = samp_data_H)
b.uu = adonis(otu.uu.H ~ genus, data = samp_data_H)
c.uu = adonis(otu.uu.H ~ tax1, data = samp_data_H)
#d.uu = adonis(otu.uu.H ~ diet3, data = samp_data_H)
#e.uu = adonis(otu.uu.H ~ diet1, data = samp_data_H)
f.uu = adonis(otu.uu.H ~ region, data = samp_data_H)
#g.uu = adonis(otu.uu.H ~ diet3, data = samp_data_H, strata = samp_data_H$family)
#h.uu = adonis(otu.uu.H ~ diet1, data = samp_data_H, strata = samp_data_H$family)
#i.uu = adonis(otu.uu.H ~ diet1, data = samp_data_H, strata = samp_data_H$diet3)
#j.uu = adonis(otu.uu.H ~ region, data = samp_data_H, strata = samp_data_H$family)
#k.uu = adonis(otu.uu.H ~ tax1, data = samp_data_H, strata = samp_data_H$region)

results.uu.H <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  Carnivores -----
gut_core_carn <- subset_samples(gut_core, diet3 == "Carnivorous")
core_rel_carn <- transform_sample_counts(gut_core_carn, function(x) x / sum(x) )
samp_data_carn <- data.frame(sample_data(core_rel_carn))

# ** Bray-curtis ----
otu.bc.carn <- vegdist(core_rel_carn@otu_table, method = "bray")

a.bc = adonis(otu.bc.carn ~ family, data = samp_data_carn)
b.bc = adonis(otu.bc.carn ~ genus, data = samp_data_carn)
c.bc = adonis(otu.bc.carn ~ tax1, data = samp_data_carn)
#d.bc = adonis(otu.bc.carn ~ diet3, data = samp_data_carn)
e.bc = adonis(otu.bc.carn ~ diet1, data = samp_data_carn)
f.bc = adonis(otu.bc.carn ~ region, data = samp_data_carn)
#g.bc = adonis(otu.bc.carn ~ diet3, data = samp_data_carn, strata = samp_data_carn$family)
#h.bc = adonis(otu.bc.carn ~ diet1, data = samp_data_carn, strata = samp_data_carn$family)
#i.bc = adonis(otu.bc.carn ~ diet1, data = samp_data_carn, strata = samp_data_carn$diet3)
#j.bc = adonis(otu.bc.carn ~ region, data = samp_data_carn, strata = samp_data_carn$family)
#k.bc = adonis(otu.bc.carn ~ tax1, data = samp_data_carn, strata = samp_data_carn$region)

results.bc.carn <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.carn <- UniFrac(core_rel_carn, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.carn ~ family, data = samp_data_carn)
b.wu = adonis(otu.wu.carn ~ genus, data = samp_data_carn)
c.wu = adonis(otu.wu.carn ~ tax1, data = samp_data_carn)
#d.wu = adonis(otu.wu.carn ~ diet3, data = samp_data_carn)
e.wu = adonis(otu.wu.carn ~ diet1, data = samp_data_carn)
f.wu = adonis(otu.wu.carn ~ region, data = samp_data_carn)
#g.wu = adonis(otu.wu.carn ~ diet3, data = samp_data_carn, strata = samp_data_carn$family)
#h.wu = adonis(otu.wu.carn ~ diet1, data = samp_data_carn, strata = samp_data_carn$family)
#i.wu = adonis(otu.wu.carn ~ diet1, data = samp_data_carn, strata = samp_data_carn$diet3)
#j.wu = adonis(otu.wu.carn ~ region, data = samp_data_carn, strata = samp_data_carn$family)
#k.wu = adonis(otu.wu.carn ~ tax1, data = samp_data_carn, strata = samp_data_carn$region)

results.wu.carn <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.carn <- UniFrac(core_rel_carn, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.carn ~ family, data = samp_data_carn)
b.uu = adonis(otu.uu.carn ~ genus, data = samp_data_carn)
c.uu = adonis(otu.uu.carn ~ tax1, data = samp_data_carn)
#d.uu = adonis(otu.uu.carn ~ diet3, data = samp_data_carn)
e.uu = adonis(otu.uu.carn ~ diet1, data = samp_data_carn)
f.uu = adonis(otu.uu.carn ~ region, data = samp_data_carn)
#g.uu = adonis(otu.uu.carn ~ diet3, data = samp_data_carn, strata = samp_data_carn$family)
#h.uu = adonis(otu.uu.carn ~ diet1, data = samp_data_carn, strata = samp_data_carn$family)
#i.uu = adonis(otu.uu.carn ~ diet1, data = samp_data_carn, strata = samp_data_carn$diet3)
#j.uu = adonis(otu.uu.carn ~ region, data = samp_data_carn, strata = samp_data_carn$family)
#k.uu = adonis(otu.uu.carn ~ tax1, data = samp_data_carn, strata = samp_data_carn$region)

results.uu.carn <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  MI -----
gut_core_MI <- subset_samples(gut_core, diet1 == "MI")
core_rel_MI <- transform_sample_counts(gut_core_MI, function(x) x / sum(x) )
samp_data_MI <- data.frame(sample_data(core_rel_MI))

# ** Bray-curtis ----
otu.bc.MI <- vegdist(core_rel_MI@otu_table, method = "bray")

a.bc = adonis(otu.bc.MI ~ family, data = samp_data_MI)
b.bc = adonis(otu.bc.MI ~ genus, data = samp_data_MI)
c.bc = adonis(otu.bc.MI ~ tax1, data = samp_data_MI)
#d.bc = adonis(otu.bc.MI ~ diet3, data = samp_data_MI)
#e.bc = adonis(otu.bc.MI ~ diet1, data = samp_data_MI)
f.bc = adonis(otu.bc.MI ~ region, data = samp_data_MI)
#g.bc = adonis(otu.bc.MI ~ diet3, data = samp_data_MI, strata = samp_data_MI$family)
#h.bc = adonis(otu.bc.MI ~ diet1, data = samp_data_MI, strata = samp_data_MI$family)
#i.bc = adonis(otu.bc.MI ~ diet1, data = samp_data_MI, strata = samp_data_MI$diet3)
#j.bc = adonis(otu.bc.MI ~ region, data = samp_data_MI, strata = samp_data_MI$family)
#k.bc = adonis(otu.bc.MI ~ tax1, data = samp_data_MI, strata = samp_data_MI$region)

results.bc.MI <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.MI <- UniFrac(core_rel_MI, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.MI ~ family, data = samp_data_MI)
b.wu = adonis(otu.wu.MI ~ genus, data = samp_data_MI)
c.wu = adonis(otu.wu.MI ~ tax1, data = samp_data_MI)
#d.wu = adonis(otu.wu.MI ~ diet3, data = samp_data_MI)
#e.wu = adonis(otu.wu.MI ~ diet1, data = samp_data_MI)
f.wu = adonis(otu.wu.MI ~ region, data = samp_data_MI)
#g.wu = adonis(otu.wu.MI ~ diet3, data = samp_data_MI, strata = samp_data_MI$family)
#h.wu = adonis(otu.wu.MI ~ diet1, data = samp_data_MI, strata = samp_data_MI$family)
#i.wu = adonis(otu.wu.MI ~ diet1, data = samp_data_MI, strata = samp_data_MI$diet3)
#j.wu = adonis(otu.wu.MI ~ region, data = samp_data_MI, strata = samp_data_MI$family)
#k.wu = adonis(otu.wu.MI ~ tax1, data = samp_data_MI, strata = samp_data_MI$region)

results.wu.MI <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.MI <- UniFrac(core_rel_MI, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.MI ~ family, data = samp_data_MI)
b.uu = adonis(otu.uu.MI ~ genus, data = samp_data_MI)
c.uu = adonis(otu.uu.MI ~ tax1, data = samp_data_MI)
#d.uu = adonis(otu.uu.MI ~ diet3, data = samp_data_MI)
#e.uu = adonis(otu.uu.MI ~ diet1, data = samp_data_MI)
f.uu = adonis(otu.uu.MI ~ region, data = samp_data_MI)
#g.uu = adonis(otu.uu.MI ~ diet3, data = samp_data_MI, strata = samp_data_MI$family)
#h.uu = adonis(otu.uu.MI ~ diet1, data = samp_data_MI, strata = samp_data_MI$family)
#i.uu = adonis(otu.uu.MI ~ diet1, data = samp_data_MI, strata = samp_data_MI$diet3)
#j.uu = adonis(otu.uu.MI ~ region, data = samp_data_MI, strata = samp_data_MI$family)
#k.uu = adonis(otu.uu.MI ~ tax1, data = samp_data_MI, strata = samp_data_MI$region)

results.uu.MI <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  SI -----
gut_core_SI <- subset_samples(gut_core, diet1 == "SI")
core_rel_SI <- transform_sample_counts(gut_core_SI, function(x) x / sum(x) )
samp_data_SI <- data.frame(sample_data(core_rel_SI))

# ** Bray-curtis ----
otu.bc.SI <- vegdist(core_rel_SI@otu_table, method = "bray")

a.bc = adonis(otu.bc.SI ~ family, data = samp_data_SI)
b.bc = adonis(otu.bc.SI ~ genus, data = samp_data_SI)
c.bc = adonis(otu.bc.SI ~ tax1, data = samp_data_SI)
#d.bc = adonis(otu.bc.SI ~ diet3, data = samp_data_SI)
#e.bc = adonis(otu.bc.SI ~ diet1, data = samp_data_SI)
f.bc = adonis(otu.bc.SI ~ region, data = samp_data_SI)
#g.bc = adonis(otu.bc.SI ~ diet3, data = samp_data_SI, strata = samp_data_SI$family)
#h.bc = adonis(otu.bc.SI ~ diet1, data = samp_data_SI, strata = samp_data_SI$family)
#i.bc = adonis(otu.bc.SI ~ diet1, data = samp_data_SI, strata = samp_data_SI$diet3)
#j.bc = adonis(otu.bc.SI ~ region, data = samp_data_SI, strata = samp_data_SI$family)
#k.bc = adonis(otu.bc.SI ~ tax1, data = samp_data_SI, strata = samp_data_SI$region)

results.bc.SI <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.SI <- UniFrac(core_rel_SI, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.SI ~ family, data = samp_data_SI)
b.wu = adonis(otu.wu.SI ~ genus, data = samp_data_SI)
c.wu = adonis(otu.wu.SI ~ tax1, data = samp_data_SI)
#d.wu = adonis(otu.wu.SI ~ diet3, data = samp_data_SI)
#e.wu = adonis(otu.wu.SI ~ diet1, data = samp_data_SI)
f.wu = adonis(otu.wu.SI ~ region, data = samp_data_SI)
#g.wu = adonis(otu.wu.SI ~ diet3, data = samp_data_SI, strata = samp_data_SI$family)
#h.wu = adonis(otu.wu.SI ~ diet1, data = samp_data_SI, strata = samp_data_SI$family)
#i.wu = adonis(otu.wu.SI ~ diet1, data = samp_data_SI, strata = samp_data_SI$diet3)
#j.wu = adonis(otu.wu.SI ~ region, data = samp_data_SI, strata = samp_data_SI$family)
#k.wu = adonis(otu.wu.SI ~ tax1, data = samp_data_SI, strata = samp_data_SI$region)

results.wu.SI <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.SI <- UniFrac(core_rel_SI, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.SI ~ family, data = samp_data_SI)
b.uu = adonis(otu.uu.SI ~ genus, data = samp_data_SI)
c.uu = adonis(otu.uu.SI ~ tax1, data = samp_data_SI)
#d.uu = adonis(otu.uu.SI ~ diet3, data = samp_data_SI)
#e.uu = adonis(otu.uu.SI ~ diet1, data = samp_data_SI)
f.uu = adonis(otu.uu.SI ~ region, data = samp_data_SI)
#g.uu = adonis(otu.uu.SI ~ diet3, data = samp_data_SI, strata = samp_data_SI$family)
#h.uu = adonis(otu.uu.SI ~ diet1, data = samp_data_SI, strata = samp_data_SI$family)
#i.uu = adonis(otu.uu.SI ~ diet1, data = samp_data_SI, strata = samp_data_SI$diet3)
#j.uu = adonis(otu.uu.SI ~ region, data = samp_data_SI, strata = samp_data_SI$family)
#k.uu = adonis(otu.uu.SI ~ tax1, data = samp_data_SI, strata = samp_data_SI$region)

results.uu.SI <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  FC -----
gut_core_FC <- subset_samples(gut_core, diet1 == "FC")
core_rel_FC <- transform_sample_counts(gut_core_FC, function(x) x / sum(x) )
samp_data_FC <- data.frame(sample_data(core_rel_FC))

# ** Bray-curtis ----
otu.bc.FC <- vegdist(core_rel_FC@otu_table, method = "bray")

a.bc = adonis(otu.bc.FC ~ family, data = samp_data_FC)
b.bc = adonis(otu.bc.FC ~ genus, data = samp_data_FC)
c.bc = adonis(otu.bc.FC ~ tax1, data = samp_data_FC)
#d.bc = adonis(otu.bc.FC ~ diet3, data = samp_data_FC)
#e.bc = adonis(otu.bc.FC ~ diet1, data = samp_data_FC)
f.bc = adonis(otu.bc.FC ~ region, data = samp_data_FC)
#g.bc = adonis(otu.bc.FC ~ diet3, data = samp_data_FC, strata = samp_data_FC$family)
#h.bc = adonis(otu.bc.FC ~ diet1, data = samp_data_FC, strata = samp_data_FC$family)
#i.bc = adonis(otu.bc.FC ~ diet1, data = samp_data_FC, strata = samp_data_FC$diet3)
#j.bc = adonis(otu.bc.FC ~ region, data = samp_data_FC, strata = samp_data_FC$family)
#k.bc = adonis(otu.bc.FC ~ tax1, data = samp_data_FC, strata = samp_data_FC$region)

results.bc.FC <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.FC <- UniFrac(core_rel_FC, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.FC ~ family, data = samp_data_FC)
b.wu = adonis(otu.wu.FC ~ genus, data = samp_data_FC)
c.wu = adonis(otu.wu.FC ~ tax1, data = samp_data_FC)
#d.wu = adonis(otu.wu.FC ~ diet3, data = samp_data_FC)
#e.wu = adonis(otu.wu.FC ~ diet1, data = samp_data_FC)
f.wu = adonis(otu.wu.FC ~ region, data = samp_data_FC)
#g.wu = adonis(otu.wu.FC ~ diet3, data = samp_data_FC, strata = samp_data_FC$family)
#h.wu = adonis(otu.wu.FC ~ diet1, data = samp_data_FC, strata = samp_data_FC$family)
#i.wu = adonis(otu.wu.FC ~ diet1, data = samp_data_FC, strata = samp_data_FC$diet3)
#j.wu = adonis(otu.wu.FC ~ region, data = samp_data_FC, strata = samp_data_FC$family)
#k.wu = adonis(otu.wu.FC ~ tax1, data = samp_data_FC, strata = samp_data_FC$region)

results.wu.FC <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.FC <- UniFrac(core_rel_FC, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.FC ~ family, data = samp_data_FC)
b.uu = adonis(otu.uu.FC ~ genus, data = samp_data_FC)
c.uu = adonis(otu.uu.FC ~ tax1, data = samp_data_FC)
#d.uu = adonis(otu.uu.FC ~ diet3, data = samp_data_FC)
#e.uu = adonis(otu.uu.FC ~ diet1, data = samp_data_FC)
f.uu = adonis(otu.uu.FC ~ region, data = samp_data_FC)
#g.uu = adonis(otu.uu.FC ~ diet3, data = samp_data_FC, strata = samp_data_FC$family)
#h.uu = adonis(otu.uu.FC ~ diet1, data = samp_data_FC, strata = samp_data_FC$family)
#i.uu = adonis(otu.uu.FC ~ diet1, data = samp_data_FC, strata = samp_data_FC$diet3)
#j.uu = adonis(otu.uu.FC ~ region, data = samp_data_FC, strata = samp_data_FC$family)
#k.uu = adonis(otu.uu.FC ~ tax1, data = samp_data_FC, strata = samp_data_FC$region)

results.uu.FC <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                       "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                       #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                       "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                       #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                       #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

#  Omnivores -----
gut_core_omni <- subset_samples(gut_core, diet3 == "Omnivorous")
core_rel_omni <- transform_sample_counts(gut_core_omni, function(x) x / sum(x) )
samp_data_omni <- data.frame(sample_data(core_rel_omni))

# ** Bray-curtis ----
otu.bc.omni <- vegdist(core_rel_omni@otu_table, method = "bray")

a.bc = adonis(otu.bc.omni ~ family, data = samp_data_omni)
b.bc = adonis(otu.bc.omni ~ genus, data = samp_data_omni)
c.bc = adonis(otu.bc.omni ~ tax1, data = samp_data_omni)
#d.bc = adonis(otu.bc.omni ~ diet3, data = samp_data_omni)
#e.bc = adonis(otu.bc.omni ~ diet1, data = samp_data_omni)
f.bc = adonis(otu.bc.omni ~ region, data = samp_data_omni)
#g.bc = adonis(otu.bc.omni ~ diet3, data = samp_data_omni, strata = samp_data_omni$family)
#h.bc = adonis(otu.bc.omni ~ diet1, data = samp_data_omni, strata = samp_data_omni$family)
#i.bc = adonis(otu.bc.omni ~ diet1, data = samp_data_omni, strata = samp_data_omni$diet3)
#j.bc = adonis(otu.bc.omni ~ region, data = samp_data_omni, strata = samp_data_omni$family)
#k.bc = adonis(otu.bc.omni ~ tax1, data = samp_data_omni, strata = samp_data_omni$region)

results.bc.omni <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.omni <- UniFrac(core_rel_omni, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.omni ~ family, data = samp_data_omni)
b.wu = adonis(otu.wu.omni ~ genus, data = samp_data_omni)
c.wu = adonis(otu.wu.omni ~ tax1, data = samp_data_omni)
#d.wu = adonis(otu.wu.omni ~ diet3, data = samp_data_omni)
#e.wu = adonis(otu.wu.omni ~ diet1, data = samp_data_omni)
f.wu = adonis(otu.wu.omni ~ region, data = samp_data_omni)
#g.wu = adonis(otu.wu.omni ~ diet3, data = samp_data_omni, strata = samp_data_omni$family)
#h.wu = adonis(otu.wu.omni ~ diet1, data = samp_data_omni, strata = samp_data_omni$family)
#i.wu = adonis(otu.wu.omni ~ diet1, data = samp_data_omni, strata = samp_data_omni$diet3)
#j.wu = adonis(otu.wu.omni ~ region, data = samp_data_omni, strata = samp_data_omni$family)
#k.wu = adonis(otu.wu.omni ~ tax1, data = samp_data_omni, strata = samp_data_omni$region)

results.wu.omni <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.omni <- UniFrac(core_rel_omni, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.omni ~ family, data = samp_data_omni)
b.uu = adonis(otu.uu.omni ~ genus, data = samp_data_omni)
c.uu = adonis(otu.uu.omni ~ tax1, data = samp_data_omni)
#d.uu = adonis(otu.uu.omni ~ diet3, data = samp_data_omni)
#e.uu = adonis(otu.uu.omni ~ diet1, data = samp_data_omni)
f.uu = adonis(otu.uu.omni ~ region, data = samp_data_omni)
#g.uu = adonis(otu.uu.omni ~ diet3, data = samp_data_omni, strata = samp_data_omni$family)
#h.uu = adonis(otu.uu.omni ~ diet1, data = samp_data_omni, strata = samp_data_omni$family)
#i.uu = adonis(otu.uu.omni ~ diet1, data = samp_data_omni, strata = samp_data_omni$diet3)
#j.uu = adonis(otu.uu.omni ~ region, data = samp_data_omni, strata = samp_data_omni$family)
#k.uu = adonis(otu.uu.omni ~ tax1, data = samp_data_omni, strata = samp_data_omni$region)

results.uu.omni <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                         #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
                         "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)



save(otu.bc.herb, otu.wu.herb, otu.uu.herb,
     otu.bc.H, otu.wu.H, otu.uu.H, 
     otu.bc.HD, otu.wu.HD, otu.uu.HD, 
     file = "matrices.herb.RData")
save(otu.bc.carn, otu.wu.carn, otu.uu.carn, 
     otu.bc.MI, otu.wu.MI, otu.uu.MI, 
     otu.bc.SI, otu.wu.SI, otu.uu.SI, 
     otu.bc.FC, otu.wu.FC, otu.uu.FC,
     file = "matrices.carn.RData")
save(otu.bc.omni, otu.wu.omni, otu.uu.omni,file= "matrices.omni.RData")
save(samp_data_herb, samp_data_H, samp_data_HD, file = "samp_data_herb.RData")
save(samp_data_carn, samp_data_MI, samp_data_SI,samp_data_FC, file = "samp_data_carn.RData")
save(samp_data_omni, file = "samp_data_omni.RData")
save(core_rel_herb, core_rel_H, core_rel_HD, file = "core_rel_herb.RData")
save(core_rel_carn, core_rel_MI,core_rel_SI,core_rel_FC, file = "core_rel_carn.RData")
save(core_rel_omni, file = "core_rel_omni.RData")

results_herb <- rbind(cbind(results.bc.herb, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Herbivores")),
                     cbind(results.wu.herb, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Herbivores")),
                     cbind(results.uu.herb, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Herbivores")),
                     cbind(results.bc.H, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Strict Herbivores")),
                     cbind(results.wu.H, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Strict Herbivores")),
                     cbind(results.uu.H, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Strict Herbivores")),
                     cbind(results.bc.HD, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Detritivores")),
                     cbind(results.wu.HD, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Detritivores")),
                     cbind(results.uu.HD, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Detritivores")))

results_carn <- rbind(cbind(results.bc.carn, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Carnivores")),
                       cbind(results.wu.carn, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Carnivores")),
                       cbind(results.uu.carn, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Carnivores")),
                       cbind(results.bc.MI, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Mobile Invertivores")),
                       cbind(results.wu.MI, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Mobile Invertivores")),
                       cbind(results.uu.MI, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Mobile Invertivores")),
                       cbind(results.bc.SI, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Sessile Invertivores")),
                       cbind(results.wu.SI, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Sessile Invertivores")),
                       cbind(results.uu.SI, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Sessile Invertivores")),
                      cbind(results.bc.FC, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Piscivorous")),
                      cbind(results.wu.FC, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Piscivorous")),
                      cbind(results.uu.FC, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Piscivorous")))

results_omni <- rbind(cbind(results.bc.omni, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Omnivores")),
                      cbind(results.wu.omni, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Omnivores")),
                      cbind(results.uu.omni, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Omnivores")))

results_diet_dataset <- rbind(results_herb, results_carn, results_omni)
write.table(results_diet_dataset, file = "results_diet_dataset.txt", sep = '\t', quote = F)

#####################_______ REGION ______________#######------
# Martinique -----
gut_core_marti <- subset_samples(gut_core, region == "Martinique")
core_rel_marti <- transform_sample_counts(gut_core_marti, function(x) x / sum(x) )
samp_data_marti <- data.frame(sample_data(core_rel_marti))

# ** Bray-curtis ----
otu.bc.marti <- vegdist(core_rel_marti@otu_table, method = "bray")

a.bc = adonis(otu.bc.marti ~ family, data = samp_data_marti)
b.bc = adonis(otu.bc.marti ~ genus, data = samp_data_marti)
c.bc = adonis(otu.bc.marti ~ tax1, data = samp_data_marti)
d.bc = adonis(otu.bc.marti ~ diet3, data = samp_data_marti)
e.bc = adonis(otu.bc.marti ~ diet1, data = samp_data_marti)
#f.bc = adonis(otu.bc.marti ~ region, data = samp_data_marti)
#g.bc = adonis(otu.bc.marti ~ diet3, data = samp_data_marti, strata = samp_data_marti$family)
#h.bc = adonis(otu.bc.marti ~ diet1, data = samp_data_marti, strata = samp_data_marti$family)
#i.bc = adonis(otu.bc.marti ~ diet1, data = samp_data_marti, strata = samp_data_marti$diet3)
#j.bc = adonis(otu.bc.marti ~ region, data = samp_data_marti, strata = samp_data_marti$family)
#k.bc = adonis(otu.bc.marti ~ tax1, data = samp_data_marti, strata = samp_data_marti$region)

results.bc.marti <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1])
                         #"Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
                         #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.marti <- UniFrac(core_rel_marti, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.marti ~ family, data = samp_data_marti)
b.wu = adonis(otu.wu.marti ~ genus, data = samp_data_marti)
c.wu = adonis(otu.wu.marti ~ tax1, data = samp_data_marti)
d.wu = adonis(otu.wu.marti ~ diet3, data = samp_data_marti)
e.wu = adonis(otu.wu.marti ~ diet1, data = samp_data_marti)
#f.wu = adonis(otu.wu.marti ~ region, data = samp_data_marti)
#g.wu = adonis(otu.wu.marti ~ diet3, data = samp_data_marti, strata = samp_data_marti$family)
#h.wu = adonis(otu.wu.marti ~ diet1, data = samp_data_marti, strata = samp_data_marti$family)
#i.wu = adonis(otu.wu.marti ~ diet1, data = samp_data_marti, strata = samp_data_marti$diet3)
#j.wu = adonis(otu.wu.marti ~ region, data = samp_data_marti, strata = samp_data_marti$family)
#k.wu = adonis(otu.wu.marti ~ tax1, data = samp_data_marti, strata = samp_data_marti$region)

results.wu.marti <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1])
                         #"Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.marti <- UniFrac(core_rel_marti, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.marti ~ family, data = samp_data_marti)
b.uu = adonis(otu.uu.marti ~ genus, data = samp_data_marti)
c.uu = adonis(otu.uu.marti ~ tax1, data = samp_data_marti)
d.uu = adonis(otu.uu.marti ~ diet3, data = samp_data_marti)
e.uu = adonis(otu.uu.marti ~ diet1, data = samp_data_marti)
#f.uu = adonis(otu.uu.marti ~ region, data = samp_data_marti)
#g.uu = adonis(otu.uu.marti ~ diet3, data = samp_data_marti, strata = samp_data_marti$family)
#h.uu = adonis(otu.uu.marti ~ diet1, data = samp_data_marti, strata = samp_data_marti$family)
#i.uu = adonis(otu.uu.marti ~ diet1, data = samp_data_marti, strata = samp_data_marti$diet3)
#j.uu = adonis(otu.uu.marti ~ region, data = samp_data_marti, strata = samp_data_marti$family)
#k.uu = adonis(otu.uu.marti ~ tax1, data = samp_data_marti, strata = samp_data_marti$region)

results.uu.marti <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                         "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                         "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1])
                         #"Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
                         #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                         #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

# Seychelles -----
gut_core_sey <- subset_samples(gut_core, region == "Seychelles")
core_rel_sey <- transform_sample_counts(gut_core_sey, function(x) x / sum(x) )
samp_data_sey <- data.frame(sample_data(core_rel_sey))

# ** Bray-curtis ----
otu.bc.sey <- vegdist(core_rel_sey@otu_table, method = "bray")

a.bc = adonis(otu.bc.sey ~ family, data = samp_data_sey)
b.bc = adonis(otu.bc.sey ~ genus, data = samp_data_sey)
c.bc = adonis(otu.bc.sey ~ tax1, data = samp_data_sey)
d.bc = adonis(otu.bc.sey ~ diet3, data = samp_data_sey)
e.bc = adonis(otu.bc.sey ~ diet1, data = samp_data_sey)
#f.bc = adonis(otu.bc.sey ~ region, data = samp_data_sey)
#g.bc = adonis(otu.bc.sey ~ diet3, data = samp_data_sey, strata = samp_data_sey$family)
#h.bc = adonis(otu.bc.sey ~ diet1, data = samp_data_sey, strata = samp_data_sey$family)
#i.bc = adonis(otu.bc.sey ~ diet1, data = samp_data_sey, strata = samp_data_sey$diet3)
#j.bc = adonis(otu.bc.sey ~ region, data = samp_data_sey, strata = samp_data_sey$family)
#k.bc = adonis(otu.bc.sey ~ tax1, data = samp_data_sey, strata = samp_data_sey$region)

results.bc.sey <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]) 
                          #"Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.sey <- UniFrac(core_rel_sey, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.sey ~ family, data = samp_data_sey)
b.wu = adonis(otu.wu.sey ~ genus, data = samp_data_sey)
c.wu = adonis(otu.wu.sey ~ tax1, data = samp_data_sey)
d.wu = adonis(otu.wu.sey ~ diet3, data = samp_data_sey)
e.wu = adonis(otu.wu.sey ~ diet1, data = samp_data_sey)
#f.wu = adonis(otu.wu.sey ~ region, data = samp_data_sey)
#g.wu = adonis(otu.wu.sey ~ diet3, data = samp_data_sey, strata = samp_data_sey$family)
#h.wu = adonis(otu.wu.sey ~ diet1, data = samp_data_sey, strata = samp_data_sey$family)
#i.wu = adonis(otu.wu.sey ~ diet1, data = samp_data_sey, strata = samp_data_sey$diet3)
#j.wu = adonis(otu.wu.sey ~ region, data = samp_data_sey, strata = samp_data_sey$family)
#k.wu = adonis(otu.wu.sey ~ tax1, data = samp_data_sey, strata = samp_data_sey$region)

results.wu.sey <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1])
                          #"Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.sey <- UniFrac(core_rel_sey, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.sey ~ family, data = samp_data_sey)
b.uu = adonis(otu.uu.sey ~ genus, data = samp_data_sey)
c.uu = adonis(otu.uu.sey ~ tax1, data = samp_data_sey)
d.uu = adonis(otu.uu.sey ~ diet3, data = samp_data_sey)
e.uu = adonis(otu.uu.sey ~ diet1, data = samp_data_sey)
#f.uu = adonis(otu.uu.sey ~ region, data = samp_data_sey)
#g.uu = adonis(otu.uu.sey ~ diet3, data = samp_data_sey, strata = samp_data_sey$family)
#h.uu = adonis(otu.uu.sey ~ diet1, data = samp_data_sey, strata = samp_data_sey$family)
#i.uu = adonis(otu.uu.sey ~ diet1, data = samp_data_sey, strata = samp_data_sey$diet3)
#j.uu = adonis(otu.uu.sey ~ region, data = samp_data_sey, strata = samp_data_sey$family)
#k.uu = adonis(otu.uu.sey ~ tax1, data = samp_data_sey, strata = samp_data_sey$region)

results.uu.sey <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1])
                          #"Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

# Europa -----
gut_core_europa <- subset_samples(gut_core, region == "Europa")
core_rel_europa <- transform_sample_counts(gut_core_europa, function(x) x / sum(x) )
samp_data_europa <- data.frame(sample_data(core_rel_europa))

# ** Bray-curtis ----
otu.bc.europa <- vegdist(core_rel_europa@otu_table, method = "bray")

a.bc = adonis(otu.bc.europa ~ family, data = samp_data_europa)
b.bc = adonis(otu.bc.europa ~ genus, data = samp_data_europa)
c.bc = adonis(otu.bc.europa ~ tax1, data = samp_data_europa)
d.bc = adonis(otu.bc.europa ~ diet3, data = samp_data_europa)
e.bc = adonis(otu.bc.europa ~ diet1, data = samp_data_europa)
#f.bc = adonis(otu.bc.europa ~ region, data = samp_data_europa)
#g.bc = adonis(otu.bc.europa ~ diet3, data = samp_data_europa, strata = samp_data_europa$family)
#h.bc = adonis(otu.bc.europa ~ diet1, data = samp_data_europa, strata = samp_data_europa$family)
#i.bc = adonis(otu.bc.europa ~ diet1, data = samp_data_europa, strata = samp_data_europa$diet3)
#j.bc = adonis(otu.bc.europa ~ region, data = samp_data_europa, strata = samp_data_europa$family)
#k.bc = adonis(otu.bc.europa ~ tax1, data = samp_data_europa, strata = samp_data_europa$region)

results.bc.europa <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1])
                          #"Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.europa <- UniFrac(core_rel_europa, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.europa ~ family, data = samp_data_europa)
b.wu = adonis(otu.wu.europa ~ genus, data = samp_data_europa)
c.wu = adonis(otu.wu.europa ~ tax1, data = samp_data_europa)
d.wu = adonis(otu.wu.europa ~ diet3, data = samp_data_europa)
e.wu = adonis(otu.wu.europa ~ diet1, data = samp_data_europa)
#f.wu = adonis(otu.wu.europa ~ region, data = samp_data_europa)
#g.wu = adonis(otu.wu.europa ~ diet3, data = samp_data_europa, strata = samp_data_europa$family)
#h.wu = adonis(otu.wu.europa ~ diet1, data = samp_data_europa, strata = samp_data_europa$family)
#i.wu = adonis(otu.wu.europa ~ diet1, data = samp_data_europa, strata = samp_data_europa$diet3)
#j.wu = adonis(otu.wu.europa ~ region, data = samp_data_europa, strata = samp_data_europa$family)
#k.wu = adonis(otu.wu.europa ~ tax1, data = samp_data_europa, strata = samp_data_europa$region)

results.wu.europa <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1])
                          #"Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
                          #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.europa <- UniFrac(core_rel_europa, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.europa ~ family, data = samp_data_europa)
b.uu = adonis(otu.uu.europa ~ genus, data = samp_data_europa)
c.uu = adonis(otu.uu.europa ~ tax1, data = samp_data_europa)
d.uu = adonis(otu.uu.europa ~ diet3, data = samp_data_europa)
e.uu = adonis(otu.uu.europa ~ diet1, data = samp_data_europa)
#f.uu = adonis(otu.uu.europa ~ region, data = samp_data_europa)
#g.uu = adonis(otu.uu.europa ~ diet3, data = samp_data_europa, strata = samp_data_europa$family)
#h.uu = adonis(otu.uu.europa ~ diet1, data = samp_data_europa, strata = samp_data_europa$family)
#i.uu = adonis(otu.uu.europa ~ diet1, data = samp_data_europa, strata = samp_data_europa$diet3)
#j.uu = adonis(otu.uu.europa ~ region, data = samp_data_europa, strata = samp_data_europa$family)
#k.uu = adonis(otu.uu.europa ~ tax1, data = samp_data_europa, strata = samp_data_europa$region)

results.uu.europa <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1])
                          #"Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

# Juan de Nova -----
gut_core_jdn <- subset_samples(gut_core, region == "Juan_de_nova")
core_rel_jdn <- transform_sample_counts(gut_core_jdn, function(x) x / sum(x) )
samp_data_jdn <- data.frame(sample_data(core_rel_jdn))

# ** Bray-curtis ----
otu.bc.jdn <- vegdist(core_rel_jdn@otu_table, method = "bray")

a.bc = adonis(otu.bc.jdn ~ family, data = samp_data_jdn)
b.bc = adonis(otu.bc.jdn ~ genus, data = samp_data_jdn)
c.bc = adonis(otu.bc.jdn ~ tax1, data = samp_data_jdn)
d.bc = adonis(otu.bc.jdn ~ diet3, data = samp_data_jdn)
e.bc = adonis(otu.bc.jdn ~ diet1, data = samp_data_jdn)
#f.bc = adonis(otu.bc.jdn ~ region, data = samp_data_jdn)
#g.bc = adonis(otu.bc.jdn ~ diet3, data = samp_data_jdn, strata = samp_data_jdn$family)
#h.bc = adonis(otu.bc.jdn ~ diet1, data = samp_data_jdn, strata = samp_data_jdn$family)
#i.bc = adonis(otu.bc.jdn ~ diet1, data = samp_data_jdn, strata = samp_data_jdn$diet3)
#j.bc = adonis(otu.bc.jdn ~ region, data = samp_data_jdn, strata = samp_data_jdn$family)
#k.bc = adonis(otu.bc.jdn ~ tax1, data = samp_data_jdn, strata = samp_data_jdn$region)

results.bc.jdn <- rbind("Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]) 
                          #"Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.jdn <- UniFrac(core_rel_jdn, weighted = T, normalized=F, parallel = F, fast=T)

a.wu = adonis(otu.wu.jdn ~ family, data = samp_data_jdn)
b.wu = adonis(otu.wu.jdn ~ genus, data = samp_data_jdn)
c.wu = adonis(otu.wu.jdn ~ tax1, data = samp_data_jdn)
d.wu = adonis(otu.wu.jdn ~ diet3, data = samp_data_jdn)
e.wu = adonis(otu.wu.jdn ~ diet1, data = samp_data_jdn)
#f.wu = adonis(otu.wu.jdn ~ region, data = samp_data_jdn)
#g.wu = adonis(otu.wu.jdn ~ diet3, data = samp_data_jdn, strata = samp_data_jdn$family)
#h.wu = adonis(otu.wu.jdn ~ diet1, data = samp_data_jdn, strata = samp_data_jdn$family)
#i.wu = adonis(otu.wu.jdn ~ diet1, data = samp_data_jdn, strata = samp_data_jdn$diet3)
#j.wu = adonis(otu.wu.jdn ~ region, data = samp_data_jdn, strata = samp_data_jdn$family)
#k.wu = adonis(otu.wu.jdn ~ tax1, data = samp_data_jdn, strata = samp_data_jdn$region)

results.wu.jdn <- rbind("Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]) 
                          #"Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)


# ** Unweighted Unifrac ----
otu.uu.jdn <- UniFrac(core_rel_jdn, weighted = F, normalized=F, parallel = F, fast=T)

a.uu = adonis(otu.uu.jdn ~ family, data = samp_data_jdn)
b.uu = adonis(otu.uu.jdn ~ genus, data = samp_data_jdn)
c.uu = adonis(otu.uu.jdn ~ tax1, data = samp_data_jdn)
d.uu = adonis(otu.uu.jdn ~ diet3, data = samp_data_jdn)
e.uu = adonis(otu.uu.jdn ~ diet1, data = samp_data_jdn)
#f.uu = adonis(otu.uu.jdn ~ region, data = samp_data_jdn)
#g.uu = adonis(otu.uu.jdn ~ diet3, data = samp_data_jdn, strata = samp_data_jdn$family)
#h.uu = adonis(otu.uu.jdn ~ diet1, data = samp_data_jdn, strata = samp_data_jdn$family)
#i.uu = adonis(otu.uu.jdn ~ diet1, data = samp_data_jdn, strata = samp_data_jdn$diet3)
#j.uu = adonis(otu.uu.jdn ~ region, data = samp_data_jdn, strata = samp_data_jdn$family)
#k.uu = adonis(otu.uu.jdn ~ tax1, data = samp_data_jdn, strata = samp_data_jdn$region)

results.uu.jdn <- rbind("Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1])
                          #"Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
                          #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

save(otu.bc.marti, otu.wu.marti,otu.uu.marti, samp_data_marti, core_rel_marti, file = "beta_marti.RData")
save(otu.bc.sey, otu.wu.sey,otu.uu.sey, samp_data_sey, core_rel_sey, file = "beta_seychelles.RData")
save(otu.bc.europa, otu.wu.europa,otu.uu.europa, samp_data_europa, core_rel_europa, file = "beta_europa.RData")
save(otu.bc.jdn, otu.wu.jdn,otu.uu.jdn, samp_data_jdn, core_rel_jdn, file = "beta_jdn.RData")

results_region <- rbind(cbind(results.bc.marti, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Martinique")),
                      cbind(results.wu.marti, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Martinique")),
                      cbind(results.uu.marti, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Martinique")),
                      cbind(results.bc.sey, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Seychelles")),
                      cbind(results.wu.sey, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Seychelles")),
                      cbind(results.uu.sey, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Seychelles")),
                      cbind(results.bc.europa, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Europa")),
                      cbind(results.wu.europa, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Europa")),
                      cbind(results.uu.europa, "Distance" = rep("Unweighted Unifrac"), "Sampling" =rep("Europa")),
                      cbind(results.bc.jdn, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Juan de Nova")),
                      cbind(results.wu.jdn, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Juan de Nova")),
                      cbind(results.uu.jdn, "Distance" = rep("Unweighted Unifrac"), "Sampling" =rep("Juan de Nova")))

write.table(results_region, file = "results_region.txt", sep = "\t", quote = F)

#####################_______ TAXONOMY ______________#######------
# Acanthuridae -----
gut_core_acanth <- subset_samples(gut_core, family == "Acanthuridae")
core_rel_acanth <- transform_sample_counts(gut_core_acanth, function(x) x / sum(x) )
samp_data_acanth <- data.frame(sample_data(core_rel_acanth))

gut_core_acanth_sp <- subset_samples(gut_core_acanth, tax1 %in% rownames(which(table(samp_data_acanth$tax1)>=3, T)))
core_rel_acanth_sp <- transform_sample_counts(gut_core_acanth_sp, function(x) x / sum(x) )
samp_data_acanth_sp <- data.frame(sample_data(core_rel_acanth_sp))

gut_core_acanth_region <- subset_samples(gut_core_acanth, tax1 %in% c("Acanthurus nigrofuscus" , "Ctenochaetus striatus"))
core_rel_acanth_region <- transform_sample_counts(gut_core_acanth_region, function(x) x / sum(x) )
samp_data_acanth_region <- data.frame(sample_data(core_rel_acanth_region))

# ** Bray-curtis ----
otu.bc.acanth <- vegdist(core_rel_acanth@otu_table, method = "bray")
otu.bc.acanth_sp <- vegdist(core_rel_acanth_sp@otu_table, method = "bray")
otu.bc.acanth_region <- vegdist(core_rel_acanth_region@otu_table, method = "bray")

#a.bc = adonis(otu.bc.acanth ~ family, data = samp_data_acanth)
b.bc = adonis(otu.bc.acanth ~ genus, data = samp_data_acanth)
c.bc = adonis(otu.bc.acanth_sp ~ tax1, data = samp_data_acanth_sp) # We keep only triplicate species
#d.bc = adonis(otu.bc.acanth ~ diet3, data = samp_data_acanth) They are all herbivores
e.bc = adonis(otu.bc.acanth ~ diet1, data = samp_data_acanth)
f.bc = adonis(otu.bc.acanth ~ region, data = samp_data_acanth)
#g.bc = adonis(otu.bc.acanth ~ diet3, data = samp_data_acanth, strata = samp_data_acanth$family)
#h.bc = adonis(otu.bc.acanth ~ diet1, data = samp_data_acanth, strata = samp_data_acanth$family)
#i.bc = adonis(otu.bc.acanth ~ diet1, data = samp_data_acanth, strata = samp_data_acanth$diet3)
#j.bc = adonis(otu.bc.acanth ~ region, data = samp_data_acanth, strata = samp_data_acanth$family)
k.bc = adonis(otu.bc.acanth_region ~ region, data = samp_data_acanth_region, strata = samp_data_acanth_region$tax1)

results.bc.acanth <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
                          "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
                          #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
                          "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
                          "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
                          #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
                          #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
                          "Species : Region" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)


# ** Weighted Unifrac ----
otu.wu.acanth <- UniFrac(core_rel_acanth, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.acanth_sp <- UniFrac(core_rel_acanth_sp, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.acanth_region <- UniFrac(core_rel_acanth_region, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.acanth ~ family, data = samp_data_acanth)
b.wu = adonis(otu.wu.acanth ~ genus, data = samp_data_acanth)
c.wu = adonis(otu.wu.acanth_sp ~ tax1, data = samp_data_acanth_sp) # We keep only triplicate species
#d.wu = adonis(otu.wu.acanth ~ diet3, data = samp_data_acanth) They are all herbivores
e.wu = adonis(otu.wu.acanth ~ diet1, data = samp_data_acanth)
f.wu = adonis(otu.wu.acanth ~ region, data = samp_data_acanth)
#g.wu = adonis(otu.wu.acanth ~ diet3, data = samp_data_acanth, strata = samp_data_acanth$family)
#h.wu = adonis(otu.wu.acanth ~ diet1, data = samp_data_acanth, strata = samp_data_acanth$family)
#i.wu = adonis(otu.wu.acanth ~ diet1, data = samp_data_acanth, strata = samp_data_acanth$diet3)
#j.wu = adonis(otu.wu.acanth ~ region, data = samp_data_acanth, strata = samp_data_acanth$family)
k.wu = adonis(otu.wu.acanth_region ~ region, data = samp_data_acanth_region, strata = samp_data_acanth_region$tax1)

results.wu.acanth <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.acanth <- UniFrac(core_rel_acanth, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.acanth_sp <- UniFrac(core_rel_acanth_sp, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.acanth_region <- UniFrac(core_rel_acanth_region, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.acanth ~ family, data = samp_data_acanth)
b.uu = adonis(otu.uu.acanth ~ genus, data = samp_data_acanth)
c.uu = adonis(otu.uu.acanth_sp ~ tax1, data = samp_data_acanth_sp) # We keep only triplicate species
#d.uu = adonis(otu.uu.acanth ~ diet3, data = samp_data_acanth) They are all herbivores
e.uu = adonis(otu.uu.acanth ~ diet1, data = samp_data_acanth)
f.uu = adonis(otu.uu.acanth ~ region, data = samp_data_acanth)
#g.uu = adonis(otu.uu.acanth ~ diet3, data = samp_data_acanth, strata = samp_data_acanth$family)
#h.uu = adonis(otu.uu.acanth ~ diet1, data = samp_data_acanth, strata = samp_data_acanth$family)
#i.uu = adonis(otu.uu.acanth ~ diet1, data = samp_data_acanth, strata = samp_data_acanth$diet3)
#j.uu = adonis(otu.uu.acanth ~ region, data = samp_data_acanth, strata = samp_data_acanth$family)
k.uu = adonis(otu.uu.acanth_region ~ region, data = samp_data_acanth_region, strata = samp_data_acanth_region$tax1)

results.uu.acanth <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


# Holocentridae -----
gut_core_holo <- subset_samples(gut_core, family == "Holocentridae")
core_rel_holo <- transform_sample_counts(gut_core_holo, function(x) x / sum(x) )
samp_data_holo <- data.frame(sample_data(core_rel_holo))

gut_core_holo_sp <- subset_samples(gut_core_holo, tax1 %in% rownames(which(table(samp_data_holo$tax1)>=3, T)))
core_rel_holo_sp <- transform_sample_counts(gut_core_holo_sp, function(x) x / sum(x) )
samp_data_holo_sp <- data.frame(sample_data(core_rel_holo_sp))
# ** Bray-curtis ----
otu.bc.holo <- vegdist(core_rel_holo@otu_table, method = "bray")
otu.bc.holo_sp <- vegdist(core_rel_holo_sp@otu_table, method = "bray")

#a.bc = adonis(otu.bc.holo ~ family, data = samp_data_holo)
b.bc = adonis(otu.bc.holo ~ genus, data = samp_data_holo)
c.bc = adonis(otu.bc.holo_sp ~ tax1, data = samp_data_holo_sp) # We keep only triplicate species
d.bc = adonis(otu.bc.holo ~ diet3, data = samp_data_holo)
e.bc = adonis(otu.bc.holo ~ diet1, data = samp_data_holo) # d and e are the same 
f.bc = adonis(otu.bc.holo ~ region, data = samp_data_holo)
#g.bc = adonis(otu.bc.holo ~ diet3, data = samp_data_holo, strata = samp_data_holo$family)
#h.bc = adonis(otu.bc.holo ~ diet1, data = samp_data_holo, strata = samp_data_holo$family)
#i.bc = adonis(otu.bc.holo ~ diet1, data = samp_data_holo, strata = samp_data_holo$diet3) 
#j.bc = adonis(otu.bc.holo ~ region, data = samp_data_holo, strata = samp_data_holo$family)
#k.bc = adonis(otu.bc.holo ~ tax1, data = samp_data_holo, strata = samp_data_holo$region)

results.bc.holo <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
  "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.holo <- UniFrac(core_rel_holo, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.holo_sp <- UniFrac(core_rel_holo_sp, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.holo ~ family, data = samp_data_holo)
b.wu = adonis(otu.wu.holo ~ genus, data = samp_data_holo)
c.wu = adonis(otu.wu.holo_sp ~ tax1, data = samp_data_holo_sp) # We keep only triplicate species
d.wu = adonis(otu.wu.holo ~ diet3, data = samp_data_holo)
e.wu = adonis(otu.wu.holo ~ diet1, data = samp_data_holo) 
f.wu = adonis(otu.wu.holo ~ region, data = samp_data_holo)
#g.wu = adonis(otu.wu.holo ~ diet3, data = samp_data_holo, strata = samp_data_holo$family)
#h.wu = adonis(otu.wu.holo ~ diet1, data = samp_data_holo, strata = samp_data_holo$family)
#i.wu = adonis(otu.wu.holo ~ diet1, data = samp_data_holo, strata = samp_data_holo$diet3) 
#j.wu = adonis(otu.wu.holo ~ region, data = samp_data_holo, strata = samp_data_holo$family)
#k.wu = adonis(otu.wu.holo ~ tax1, data = samp_data_holo, strata = samp_data_holo$region)

results.wu.holo <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.holo <- UniFrac(core_rel_holo, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.holo_sp <- UniFrac(core_rel_holo_sp, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.holo ~ family, data = samp_data_holo)
b.uu = adonis(otu.uu.holo ~ genus, data = samp_data_holo)
c.uu = adonis(otu.uu.holo_sp ~ tax1, data = samp_data_holo_sp) # We keep only triplicate species
d.uu = adonis(otu.uu.holo ~ diet3, data = samp_data_holo)
e.uu = adonis(otu.uu.holo ~ diet1, data = samp_data_holo) # d and e are the same 
f.uu = adonis(otu.uu.holo ~ region, data = samp_data_holo)
#g.uu = adonis(otu.uu.holo ~ diet3, data = samp_data_holo, strata = samp_data_holo$family)
#h.uu = adonis(otu.uu.holo ~ diet1, data = samp_data_holo, strata = samp_data_holo$family)
#i.uu = adonis(otu.uu.holo ~ diet1, data = samp_data_holo, strata = samp_data_holo$diet3) 
#j.uu = adonis(otu.uu.holo ~ region, data = samp_data_holo, strata = samp_data_holo$family)
#k.uu = adonis(otu.uu.holo ~ tax1, data = samp_data_holo, strata = samp_data_holo$region)

results.uu.holo <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


# Labridae -----
gut_core_lab <- subset_samples(gut_core, family == "Labridae")
core_rel_lab <- transform_sample_counts(gut_core_lab, function(x) x / sum(x) )
samp_data_lab <- data.frame(sample_data(core_rel_lab))

gut_core_lab_sp <- subset_samples(gut_core_lab, tax1 %in% rownames(which(table(samp_data_lab$tax1)>=3, T)))
core_rel_lab_sp <- transform_sample_counts(gut_core_lab_sp, function(x) x / sum(x) )
samp_data_lab_sp <- data.frame(sample_data(core_rel_lab_sp))

gut_core_lab_region <- subset_samples(gut_core_lab, tax1 == "Halichoeres hortulanus")
core_rel_lab_region <- transform_sample_counts(gut_core_lab_region, function(x) x / sum(x) )
samp_data_lab_region <- data.frame(sample_data(core_rel_lab_region))

# ** Bray-curtis ----
otu.bc.lab <- vegdist(core_rel_lab@otu_table, method = "bray")
otu.bc.lab_sp <- vegdist(core_rel_lab_sp@otu_table, method = "bray")
otu.bc.lab_region <- vegdist(core_rel_lab_region@otu_table, method = "bray")

#a.bc = adonis(otu.bc.lab ~ family, data = samp_data_lab)
b.bc = adonis(otu.bc.lab_sp ~ genus, data = samp_data_lab_sp) # We keep only triplicate genus (which is actually the same that sp)
c.bc = adonis(otu.bc.lab_sp ~ tax1, data = samp_data_lab_sp) # We keep only triplicate species
#d.bc = adonis(otu.bc.lab ~ diet3, data = samp_data_lab) # They are all carnivorous
#e.bc = adonis(otu.bc.lab ~ diet1, data = samp_data_lab) # They are all MI (only 1 FC)
f.bc = adonis(otu.bc.lab ~ region, data = samp_data_lab)
#g.bc = adonis(otu.bc.lab ~ diet3, data = samp_data_lab, strata = samp_data_lab$family)
#h.bc = adonis(otu.bc.lab ~ diet1, data = samp_data_lab, strata = samp_data_lab$family)
#i.bc = adonis(otu.bc.lab ~ diet1, data = samp_data_lab, strata = samp_data_lab$diet3) 
#j.bc = adonis(otu.bc.lab ~ region, data = samp_data_lab, strata = samp_data_lab$family)
k.bc = adonis(otu.bc.lab_region ~ region, data = samp_data_lab_region) # Only on Halichoeres hortulanus

results.bc.lab <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.lab <- UniFrac(core_rel_lab, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.lab_sp <- UniFrac(core_rel_lab_sp, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.lab_region <- UniFrac(core_rel_lab_region, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.lab ~ family, data = samp_data_lab)
b.wu = adonis(otu.wu.lab ~ genus, data = samp_data_lab)
c.wu = adonis(otu.wu.lab_sp ~ tax1, data = samp_data_lab_sp) # We keep only triplicate species
#d.wu = adonis(otu.wu.lab ~ diet3, data = samp_data_lab)
#e.wu = adonis(otu.wu.lab ~ diet1, data = samp_data_lab)
f.wu = adonis(otu.wu.lab ~ region, data = samp_data_lab)
#g.wu = adonis(otu.wu.lab ~ diet3, data = samp_data_lab, strata = samp_data_lab$family)
#h.wu = adonis(otu.wu.lab ~ diet1, data = samp_data_lab, strata = samp_data_lab$family)
#i.wu = adonis(otu.wu.lab ~ diet1, data = samp_data_lab, strata = samp_data_lab$diet3) 
#j.wu = adonis(otu.wu.lab ~ region, data = samp_data_lab, strata = samp_data_lab$family)
k.wu = adonis(otu.wu.lab_region ~ region, data = samp_data_lab_region)

results.wu.lab <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.lab <- UniFrac(core_rel_lab, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.lab_sp <- UniFrac(core_rel_lab_sp, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.lab_region <- UniFrac(core_rel_lab_region, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.lab ~ family, data = samp_data_lab)
b.uu = adonis(otu.uu.lab ~ genus, data = samp_data_lab)
c.uu = adonis(otu.uu.lab_sp ~ tax1, data = samp_data_lab_sp) # We keep only triplicate species
#d.uu = adonis(otu.uu.lab ~ diet3, data = samp_data_lab)
#e.uu = adonis(otu.uu.lab ~ diet1, data = samp_data_lab) # d and e are the same 
f.uu = adonis(otu.uu.lab ~ region, data = samp_data_lab)
#g.uu = adonis(otu.uu.lab ~ diet3, data = samp_data_lab, strata = samp_data_lab$family)
#h.uu = adonis(otu.uu.lab ~ diet1, data = samp_data_lab, strata = samp_data_lab$family)
#i.uu = adonis(otu.uu.lab ~ diet1, data = samp_data_lab, strata = samp_data_lab$diet3) 
#j.uu = adonis(otu.uu.lab ~ region, data = samp_data_lab, strata = samp_data_lab$family)
k.uu = adonis(otu.uu.lab_region ~ region, data = samp_data_lab_region)

results.uu.lab <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

#  Lethrinidae -----
gut_core_leth <- subset_samples(gut_core, family == "Lethrinidae")
core_rel_leth <- transform_sample_counts(gut_core_leth, function(x) x / sum(x) )
samp_data_leth <- data.frame(sample_data(core_rel_leth))

gut_core_leth_sp <- subset_samples(gut_core_leth, tax1 %in% rownames(which(table(samp_data_leth$tax1)>=3, T)))
core_rel_leth_sp <- transform_sample_counts(gut_core_leth_sp, function(x) x / sum(x) )
samp_data_leth_sp <- data.frame(sample_data(core_rel_leth_sp))

gut_core_leth_region <- subset_samples(gut_core_leth, tax1 == "Monotaxis grandoculis")
core_rel_leth_region <- transform_sample_counts(gut_core_leth_region, function(x) x / sum(x) )
samp_data_leth_region <- data.frame(sample_data(core_rel_leth_region))

# ** Bray-curtis ----
otu.bc.leth <- vegdist(core_rel_leth@otu_table, method = "bray")
otu.bc.leth_sp <- vegdist(core_rel_leth_sp@otu_table, method = "bray")
otu.bc.leth_region <- vegdist(core_rel_leth_region@otu_table, method = "bray")

#a.bc = adonis(otu.bc.leth ~ family, data = samp_data_leth)
b.bc = adonis(otu.bc.leth ~ genus, data = samp_data_leth) # We keep only triplicate genus (which is actually the same that sp)
c.bc = adonis(otu.bc.leth_sp ~ tax1, data = samp_data_leth_sp) # We keep only triplicate species
#d.bc = adonis(otu.bc.leth ~ diet3, data = samp_data_leth) # They are all carnivorous
#e.bc = adonis(otu.bc.leth ~ diet1, data = samp_data_leth) # They are all MI (only 1 FC)
f.bc = adonis(otu.bc.leth ~ region, data = samp_data_leth)
#g.bc = adonis(otu.bc.leth ~ diet3, data = samp_data_leth, strata = samp_data_leth$family)
#h.bc = adonis(otu.bc.leth ~ diet1, data = samp_data_leth, strata = samp_data_leth$family)
#i.bc = adonis(otu.bc.leth ~ diet1, data = samp_data_leth, strata = samp_data_leth$diet3) 
#j.bc = adonis(otu.bc.leth ~ region, data = samp_data_leth, strata = samp_data_leth$family)
k.bc = adonis(otu.bc.leth_region ~ region, data = samp_data_leth_region) #Only on M. grandoculis

results.bc.leth <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.leth <- UniFrac(core_rel_leth, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.leth_sp <- UniFrac(core_rel_leth_sp, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.leth_region <- UniFrac(core_rel_leth_region, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.leth ~ family, data = samp_data_leth)
b.wu = adonis(otu.wu.leth ~ genus, data = samp_data_leth)
c.wu = adonis(otu.wu.leth_sp ~ tax1, data = samp_data_leth_sp) # We keep only triplicate species
#d.wu = adonis(otu.wu.leth ~ diet3, data = samp_data_leth)
#e.wu = adonis(otu.wu.leth ~ diet1, data = samp_data_leth) # d and e are the same 
f.wu = adonis(otu.wu.leth ~ region, data = samp_data_leth)
#g.wu = adonis(otu.wu.leth ~ diet3, data = samp_data_leth, strata = samp_data_leth$family)
#h.wu = adonis(otu.wu.leth ~ diet1, data = samp_data_leth, strata = samp_data_leth$family)
#i.wu = adonis(otu.wu.leth ~ diet1, data = samp_data_leth, strata = samp_data_leth$diet3) 
#j.wu = adonis(otu.wu.leth ~ region, data = samp_data_leth, strata = samp_data_leth$family)
k.wu = adonis(otu.wu.leth_region ~ region, data = samp_data_leth_region) #Only on M. grandoculis

results.wu.leth <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.leth <- UniFrac(core_rel_leth, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.leth_sp <- UniFrac(core_rel_leth_sp, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.leth_region <- UniFrac(core_rel_leth_region, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.leth ~ family, data = samp_data_leth)
b.uu = adonis(otu.uu.leth ~ genus, data = samp_data_leth)
c.uu = adonis(otu.uu.leth_sp ~ tax1, data = samp_data_leth_sp) # We keep only triplicate species
#d.uu = adonis(otu.uu.leth ~ diet3, data = samp_data_leth)
#e.uu = adonis(otu.uu.leth ~ diet1, data = samp_data_leth) # d and e are the same 
f.uu = adonis(otu.uu.leth ~ region, data = samp_data_leth)
#g.uu = adonis(otu.uu.leth ~ diet3, data = samp_data_leth, strata = samp_data_leth$family)
#h.uu = adonis(otu.uu.leth ~ diet1, data = samp_data_leth, strata = samp_data_leth$family)
#i.uu = adonis(otu.uu.leth ~ diet1, data = samp_data_leth, strata = samp_data_leth$diet3) 
#j.uu = adonis(otu.uu.leth ~ region, data = samp_data_leth, strata = samp_data_leth$family)
k.uu = adonis(otu.uu.leth_region ~ region, data = samp_data_leth_region) # Only on Monotaxis grandoculis

results.uu.leth <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)


#  Lutjanidae -----
gut_core_lutj <- subset_samples(gut_core, family == "Lutjanidae")
core_rel_lutj <- transform_sample_counts(gut_core_lutj, function(x) x / sum(x) )
samp_data_lutj <- data.frame(sample_data(core_rel_lutj))

gut_core_lutj_sp <- subset_samples(gut_core_lutj, tax1 %in% rownames(which(table(samp_data_lutj$tax1)>=3, T)))
core_rel_lutj_sp <- transform_sample_counts(gut_core_lutj_sp, function(x) x / sum(x) )
samp_data_lutj_sp <- data.frame(sample_data(core_rel_lutj_sp))

gut_core_lutj_region <- subset_samples(gut_core_lutj, tax1 %in% c("Lutjanus bohar", "Lutjanus kasmira"))
core_rel_lutj_region <- transform_sample_counts(gut_core_lutj_region, function(x) x / sum(x) )
samp_data_lutj_region <- data.frame(sample_data(core_rel_lutj_region))

# ** Bray-curtis ----
otu.bc.lutj <- vegdist(core_rel_lutj@otu_table, method = "bray")
otu.bc.lutj_sp <- vegdist(core_rel_lutj_sp@otu_table, method = "bray")
otu.bc.lutj_region <- vegdist(core_rel_lutj_region@otu_table, method = "bray")

#a.bc = adonis(otu.bc.lutj ~ family, data = samp_data_lutj)
b.bc = adonis(otu.bc.lutj ~ genus, data = samp_data_lutj) # We keep only triplicate genus (which is actually the same that sp)
c.bc = adonis(otu.bc.lutj_sp ~ tax1, data = samp_data_lutj_sp) # We keep only triplicate species
#d.bc = adonis(otu.bc.lutj ~ diet3, data = samp_data_lutj) # They are all carnivorous
e.bc = adonis(otu.bc.lutj ~ diet1, data = samp_data_lutj) 
f.bc = adonis(otu.bc.lutj_sp ~ region, data = samp_data_lutj_sp) # To avoid the only sample from Martinique
#g.bc = adonis(otu.bc.lutj ~ diet3, data = samp_data_lutj, strata = samp_data_lutj$family)
#h.bc = adonis(otu.bc.lutj ~ diet1, data = samp_data_lutj, strata = samp_data_lutj$family)
#i.bc = adonis(otu.bc.lutj ~ diet1, data = samp_data_lutj, strata = samp_data_lutj$diet3) 
#j.bc = adonis(otu.bc.lutj ~ region, data = samp_data_lutj, strata = samp_data_lutj$family)
k.bc = adonis(otu.bc.lutj_region ~ region, data = samp_data_lutj_region, strata = samp_data_lutj_region$tax1) # only on 2 species

results.bc.lutj <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.lutj <- UniFrac(core_rel_lutj, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.lutj_sp <- UniFrac(core_rel_lutj_sp, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.lutj_region <- UniFrac(core_rel_lutj_region, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.lutj ~ family, data = samp_data_lutj)
b.wu = adonis(otu.wu.lutj ~ genus, data = samp_data_lutj)
c.wu = adonis(otu.wu.lutj_sp ~ tax1, data = samp_data_lutj_sp) # We keep only triplicate species
#d.wu = adonis(otu.wu.lutj ~ diet3, data = samp_data_lutj) # They are all carnivores
e.wu = adonis(otu.wu.lutj ~ diet1, data = samp_data_lutj) 
f.wu = adonis(otu.wu.lutj_sp ~ region, data = samp_data_lutj_sp)
#g.wu = adonis(otu.wu.lutj ~ diet3, data = samp_data_lutj, strata = samp_data_lutj$family)
#h.wu = adonis(otu.wu.lutj ~ diet1, data = samp_data_lutj, strata = samp_data_lutj$family)
#i.wu = adonis(otu.wu.lutj ~ diet1, data = samp_data_lutj, strata = samp_data_lutj$diet3) 
#j.wu = adonis(otu.wu.lutj ~ region, data = samp_data_lutj, strata = samp_data_lutj$family)
k.wu = adonis(otu.wu.lutj_region ~ region, data = samp_data_lutj_region, strata = samp_data_lutj_region$tax1)

results.wu.lutj <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.lutj <- UniFrac(core_rel_lutj, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.lutj_sp <- UniFrac(core_rel_lutj_sp, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.lutj_region <- UniFrac(core_rel_lutj_region, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.lutj ~ family, data = samp_data_lutj)
b.uu = adonis(otu.uu.lutj ~ genus, data = samp_data_lutj)
c.uu = adonis(otu.uu.lutj_sp ~ tax1, data = samp_data_lutj_sp) # We keep only triplicate species
#d.uu = adonis(otu.uu.lutj ~ diet3, data = samp_data_lutj)
e.uu = adonis(otu.uu.lutj ~ diet1, data = samp_data_lutj) 
f.uu = adonis(otu.uu.lutj_sp ~ region, data = samp_data_lutj_sp)
#g.uu = adonis(otu.uu.lutj ~ diet3, data = samp_data_lutj, strata = samp_data_lutj$family)
#h.uu = adonis(otu.uu.lutj ~ diet1, data = samp_data_lutj, strata = samp_data_lutj$family)
#i.uu = adonis(otu.uu.lutj ~ diet1, data = samp_data_lutj, strata = samp_data_lutj$diet3) 
#j.uu = adonis(otu.uu.lutj ~ region, data = samp_data_lutj, strata = samp_data_lutj$family)
k.uu = adonis(otu.uu.lutj_region ~ region, data = samp_data_lutj_region, strata = samp_data_lutj_region$tax1)

results.uu.lutj <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  "Species : Region" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

#  Pomacentridae -----
gut_core_poma <- subset_samples(gut_core, family == "Pomacentridae")
core_rel_poma <- transform_sample_counts(gut_core_poma, function(x) x / sum(x) )
samp_data_poma <- data.frame(sample_data(core_rel_poma))

gut_core_poma_sp <- subset_samples(gut_core_poma, tax1 %in% rownames(which(table(samp_data_poma$tax1)>=3, T)))
core_rel_poma_sp <- transform_sample_counts(gut_core_poma_sp, function(x) x / sum(x) )
samp_data_poma_sp <- data.frame(sample_data(core_rel_poma_sp))

gut_core_poma_gen <- subset_samples(gut_core_poma, !tax1 %in% c("Stegastes partitus", "Lepidozygus tapeinosoma", "Microspathodon chrysurus"))
core_rel_poma_gen <- transform_sample_counts(gut_core_poma_gen, function(x) x / sum(x) )
samp_data_poma_gen <- data.frame(sample_data(core_rel_poma_gen))

# ** Bray-curtis ----
otu.bc.poma <- vegdist(core_rel_poma@otu_table, method = "bray")
otu.bc.poma_sp <- vegdist(core_rel_poma_sp@otu_table, method = "bray")
otu.bc.poma_gen <- vegdist(core_rel_poma_gen@otu_table, method = "bray")

#a.bc = adonis(otu.bc.poma ~ family, data = samp_data_poma)
b.bc = adonis(otu.bc.poma_gen ~ genus, data = samp_data_poma_gen) # We keep only triplicate genus 
c.bc = adonis(otu.bc.poma_sp ~ tax1, data = samp_data_poma_sp) # We keep only triplicate species
d.bc = adonis(otu.bc.poma ~ diet3, data = samp_data_poma) 
#e.bc = adonis(otu.bc.poma ~ diet1, data = samp_data_poma) # Not enough effectives
f.bc = adonis(otu.bc.poma_sp ~ region, data = samp_data_poma_sp) # To avoid the only sample from Martinique
#g.bc = adonis(otu.bc.poma ~ diet3, data = samp_data_poma, strata = samp_data_poma$family)
#h.bc = adonis(otu.bc.poma ~ diet1, data = samp_data_poma, strata = samp_data_poma$family)
#i.bc = adonis(otu.bc.poma ~ diet1, data = samp_data_poma, strata = samp_data_poma$diet3) 
#j.bc = adonis(otu.bc.poma ~ region, data = samp_data_poma, strata = samp_data_poma$family)
#k.bc = adonis(otu.bc.poma_region ~ tax1, data = samp_data_poma_region, strata = samp_data_poma_region$region) # only on 2 species

results.bc.poma <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
  "Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.poma <- UniFrac(core_rel_poma, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.poma_sp <- UniFrac(core_rel_poma_sp, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.poma_gen <- UniFrac(core_rel_poma_gen, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.poma ~ family, data = samp_data_poma)
b.wu = adonis(otu.wu.poma_gen ~ genus, data = samp_data_poma_gen) # We keep only triplicate genus
c.wu = adonis(otu.wu.poma_sp ~ tax1, data = samp_data_poma_sp) # We keep only triplicate species
d.wu = adonis(otu.wu.poma ~ diet3, data = samp_data_poma) 
#e.wu = adonis(otu.wu.poma ~ diet1, data = samp_data_poma) 
f.wu = adonis(otu.wu.poma ~ region, data = samp_data_poma)
#g.wu = adonis(otu.wu.poma ~ diet3, data = samp_data_poma, strata = samp_data_poma$family)
#h.wu = adonis(otu.wu.poma ~ diet1, data = samp_data_poma, strata = samp_data_poma$family)
#i.wu = adonis(otu.wu.poma ~ diet1, data = samp_data_poma, strata = samp_data_poma$diet3) 
#j.wu = adonis(otu.wu.poma ~ region, data = samp_data_poma, strata = samp_data_poma$family)
#k.wu = adonis(otu.wu.poma_region ~ tax1, data = samp_data_poma_region, strata = samp_data_poma_region$region)

results.wu.poma <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.poma <- UniFrac(core_rel_poma, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.poma_sp <- UniFrac(core_rel_poma_sp, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.poma_gen <- UniFrac(core_rel_poma_gen, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.poma ~ family, data = samp_data_poma)
b.uu = adonis(otu.uu.poma_gen ~ genus, data = samp_data_poma_gen)
#c.uu = adonis(otu.uu.poma_sp ~ tax1, data = samp_data_poma_sp) # We keep only triplicate species
d.uu = adonis(otu.uu.poma ~ diet3, data = samp_data_poma)
e.uu = adonis(otu.uu.poma ~ diet1, data = samp_data_poma) 
f.uu = adonis(otu.uu.poma ~ region, data = samp_data_poma)
#g.uu = adonis(otu.uu.poma ~ diet3, data = samp_data_poma, strata = samp_data_poma$family)
#h.uu = adonis(otu.uu.poma ~ diet1, data = samp_data_poma, strata = samp_data_poma$family)
#i.uu = adonis(otu.uu.poma ~ diet1, data = samp_data_poma, strata = samp_data_poma$diet3) 
#j.uu = adonis(otu.uu.poma ~ region, data = samp_data_poma, strata = samp_data_poma$family)
#k.uu = adonis(otu.uu.poma_region ~ tax1, data = samp_data_poma_region, strata = samp_data_poma_region$region)

results.uu.poma <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)

#  Scaridae -----
gut_core_scar <- subset_samples(gut_core, family == "Scaridae")
core_rel_scar <- transform_sample_counts(gut_core_scar, function(x) x / sum(x) )
samp_data_scar <- data.frame(sample_data(core_rel_scar))

gut_core_scar_sp <- subset_samples(gut_core_scar, tax1 %in% rownames(which(table(samp_data_scar$tax1)>=3, T)))
core_rel_scar_sp <- transform_sample_counts(gut_core_scar_sp, function(x) x / sum(x) )
samp_data_scar_sp <- data.frame(sample_data(core_rel_scar_sp))

# ** Bray-curtis ----
otu.bc.scar <- vegdist(core_rel_scar@otu_table, method = "bray")
otu.bc.scar_sp <- vegdist(core_rel_scar_sp@otu_table, method = "bray")

#a.bc = adonis(otu.bc.scar ~ family, data = samp_data_scar)
b.bc = adonis(otu.bc.scar ~ genus, data = samp_data_scar) # We keep only triplicate genus 
c.bc = adonis(otu.bc.scar_sp ~ tax1, data = samp_data_scar_sp) # We keep only triplicate species
#d.bc = adonis(otu.bc.scar ~ diet3, data = samp_data_scar) 
e.bc = adonis(otu.bc.scar ~ diet1, data = samp_data_scar) # Not enough effectives
f.bc = adonis(otu.bc.scar_sp ~ region, data = samp_data_scar_sp) # To avoid the only sample from Martinique
#g.bc = adonis(otu.bc.scar ~ diet3, data = samp_data_scar, strata = samp_data_scar$family)
#h.bc = adonis(otu.bc.scar ~ diet1, data = samp_data_scar, strata = samp_data_scar$family)
#i.bc = adonis(otu.bc.scar ~ diet1, data = samp_data_scar, strata = samp_data_scar$diet3) 
#j.bc = adonis(otu.bc.scar ~ region, data = samp_data_scar, strata = samp_data_scar$family)
#k.bc = adonis(otu.bc.scar_region ~ tax1, data = samp_data_scar_region, strata = samp_data_scar_region$region) # only on 2 species

results.bc.scar <- rbind(#"Family" = paste("R2 =",a.bc$aov.tab$R2[1],";","p =",(a.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.bc$aov.tab$R2[1],";","p =",(b.bc$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.bc$aov.tab$R2[1],";","p =",(c.bc$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.bc$aov.tab$R2[1],";","p =",(d.bc$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.bc$aov.tab$R2[1],";","p =",(e.bc$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.bc$aov.tab$R2[1],";","p =",(f.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.bc$aov.tab$R2[1],";","p =",(g.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.bc$aov.tab$R2[1],";","p =",(h.bc$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.bc$aov.tab$R2[1],";","p =",(i.bc$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.bc$aov.tab$R2[1],";","p =",(j.bc$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.bc$aov.tab$R2[1],";","p =",(k.bc$aov.tab$`Pr(>F)`)[1])
)

# ** Weighted Unifrac ----
otu.wu.scar <- UniFrac(core_rel_scar, weighted = T, normalized=F, parallel = F, fast=T)
otu.wu.scar_sp <- UniFrac(core_rel_scar_sp, weighted = T, normalized=F, parallel = F, fast=T)

#a.wu = adonis(otu.wu.scar ~ family, data = samp_data_scar)
b.wu = adonis(otu.wu.scar ~ genus, data = samp_data_scar) # We keep only triplicate genus
c.wu = adonis(otu.wu.scar_sp ~ tax1, data = samp_data_scar_sp) # We keep only triplicate species
#d.wu = adonis(otu.wu.scar ~ diet3, data = samp_data_scar) 
e.wu = adonis(otu.wu.scar ~ diet1, data = samp_data_scar) 
f.wu = adonis(otu.wu.scar ~ region, data = samp_data_scar)
#g.wu = adonis(otu.wu.scar ~ diet3, data = samp_data_scar, strata = samp_data_scar$family)
#h.wu = adonis(otu.wu.scar ~ diet1, data = samp_data_scar, strata = samp_data_scar$family)
#i.wu = adonis(otu.wu.scar ~ diet1, data = samp_data_scar, strata = samp_data_scar$diet3) 
#j.wu = adonis(otu.wu.scar ~ region, data = samp_data_scar, strata = samp_data_scar$family)
#k.wu = adonis(otu.wu.scar_region ~ tax1, data = samp_data_scar_region, strata = samp_data_scar_region$region)

results.wu.scar <- rbind(#"Family" = paste("R2 =",a.wu$aov.tab$R2[1],";","p =",(a.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.wu$aov.tab$R2[1],";","p =",(b.wu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.wu$aov.tab$R2[1],";","p =",(c.wu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.wu$aov.tab$R2[1],";","p =",(d.wu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.wu$aov.tab$R2[1],";","p =",(e.wu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.wu$aov.tab$R2[1],";","p =",(f.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.wu$aov.tab$R2[1],";","p =",(g.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.wu$aov.tab$R2[1],";","p =",(h.wu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.wu$aov.tab$R2[1],";","p =",(i.wu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.wu$aov.tab$R2[1],";","p =",(j.wu$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.wu$aov.tab$R2[1],";","p =",(k.wu$aov.tab$`Pr(>F)`)[1])
)

# ** Unweighted Unifrac ----
otu.uu.scar <- UniFrac(core_rel_scar, weighted = F, normalized=F, parallel = F, fast=T)
otu.uu.scar_sp <- UniFrac(core_rel_scar_sp, weighted = F, normalized=F, parallel = F, fast=T)

#a.uu = adonis(otu.uu.scar ~ family, data = samp_data_scar)
b.uu = adonis(otu.uu.scar ~ genus, data = samp_data_scar)
c.uu = adonis(otu.uu.scar_sp ~ tax1, data = samp_data_scar_sp) # We keep only triplicate species
#d.uu = adonis(otu.uu.scar ~ diet3, data = samp_data_scar)
e.uu = adonis(otu.uu.scar ~ diet1, data = samp_data_scar) 
f.uu = adonis(otu.uu.scar ~ region, data = samp_data_scar)
#g.uu = adonis(otu.uu.scar ~ diet3, data = samp_data_scar, strata = samp_data_scar$family)
#h.uu = adonis(otu.uu.scar ~ diet1, data = samp_data_scar, strata = samp_data_scar$family)
#i.uu = adonis(otu.uu.scar ~ diet1, data = samp_data_scar, strata = samp_data_scar$diet3) 
#j.uu = adonis(otu.uu.scar ~ region, data = samp_data_scar, strata = samp_data_scar$family)
#k.uu = adonis(otu.uu.scar_region ~ tax1, data = samp_data_scar_region, strata = samp_data_scar_region$region)

results.uu.scar <- rbind(#"Family" = paste("R2 =",a.uu$aov.tab$R2[1],";","p =",(a.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Genus" = paste("R2 =",b.uu$aov.tab$R2[1],";","p =",(b.uu$aov.tab$`Pr(>F)`)[1]) ,
  "Species"= paste("R2 =",c.uu$aov.tab$R2[1],";","p =",(c.uu$aov.tab$`Pr(>F)`)[1]) , 
  #"Diet3" = paste("R2 =",d.uu$aov.tab$R2[1],";","p =",(d.uu$aov.tab$`Pr(>F)`)[1]) , 
  "Diet1" = paste("R2 =",e.uu$aov.tab$R2[1],";","p =",(e.uu$aov.tab$`Pr(>F)`)[1]), 
  "Region" = paste("R2 =",f.uu$aov.tab$R2[1],";","p =",(f.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Diet3" = paste("R2 =",g.uu$aov.tab$R2[1],";","p =",(g.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Diet1" = paste("R2 =",h.uu$aov.tab$R2[1],";","p =",(h.uu$aov.tab$`Pr(>F)`)[1]),
  #"Family : Region" = paste("R2 =",i.uu$aov.tab$R2[1],";","p =",(i.uu$aov.tab$`Pr(>F)`)[1])
  #"Family : Species" = paste("R2 =",j.uu$aov.tab$R2[1],";","p =",(j.uu$aov.tab$`Pr(>F)`)[1]),
  #"Region : Species" = paste("R2 =",k.uu$aov.tab$R2[1],";","p =",(k.uu$aov.tab$`Pr(>F)`)[1])
)



save(samp_data_acanth,samp_data_acanth_sp,samp_data_acanth_region,core_rel_acanth, core_rel_acanth_sp, core_rel_acanth_region,file = "data_acanth.RData")
save(samp_data_holo,samp_data_holo_sp,core_rel_holo, core_rel_holo_sp,file = "data_holo.RData")
save(samp_data_lab,samp_data_lab_sp,samp_data_lab_region,core_rel_lab, core_rel_lab_sp,core_rel_lab_region,file = "data_lab.RData")
save(samp_data_leth,samp_data_leth_sp,samp_data_leth_region,core_rel_leth, core_rel_leth_sp,core_rel_leth_region,file = "data_leth.RData")
save(samp_data_lutj,samp_data_lutj_sp,samp_data_lutj_region,core_rel_lutj, core_rel_lutj_sp,core_rel_lutj_region,file = "data_lutj.RData")
save(samp_data_poma,samp_data_poma_sp,samp_data_poma_gen,core_rel_poma, core_rel_poma_sp,core_rel_poma_gen,file = "data_poma.RData")
save(samp_data_scar,samp_data_scar_sp,core_rel_scar, core_rel_scar_sp,file = "data_scar.RData")

save(otu.bc.acanth, otu.bc.acanth_region, otu.bc.acanth_sp, 
     otu.wu.acanth, otu.wu.acanth_region,otu.wu.acanth_sp,
     otu.uu.acanth,otu.uu.acanth_region,otu.uu.acanth_sp, file = "acanth_matrices.RData")
save(otu.bc.lab, otu.bc.lab_region, otu.bc.lab_sp, 
     otu.wu.lab, otu.wu.lab_region,otu.wu.lab_sp,
     otu.uu.lab,otu.uu.lab_region,otu.uu.lab_sp, file = "lab_matrices.RData")
save(otu.bc.leth, otu.bc.leth_region, otu.bc.leth_sp, 
     otu.wu.leth, otu.wu.leth_region,otu.wu.leth_sp,
     otu.uu.leth,otu.uu.leth_region,otu.uu.leth_sp, file = "leth_matrices.RData")
save(otu.bc.lutj, otu.bc.lutj_region, otu.bc.lutj_sp, 
     otu.wu.lutj, otu.wu.lutj_region,otu.wu.lutj_sp,
     otu.uu.lutj,otu.uu.lutj_region,otu.uu.lutj_sp, file = "lutj_matrices.RData")
save(otu.bc.holo,  otu.bc.holo_sp, 
     otu.wu.holo, otu.wu.holo_sp,
     otu.uu.holo,otu.uu.holo_sp, file = "holo_matrices.RData")
save(otu.bc.poma, otu.bc.poma_gen, otu.bc.poma_sp, 
     otu.wu.poma, otu.wu.poma_gen,otu.wu.poma_sp,
     otu.uu.poma,otu.uu.poma_gen,otu.uu.poma_sp, file = "poma_matrices.RData")
save(otu.bc.scar,  otu.bc.scar_sp, 
     otu.wu.scar, otu.wu.scar_sp,
     otu.uu.scar,otu.uu.scar_sp, file = "scar_matrices.RData")


results_tax <- rbind(cbind(results.bc.acanth, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Acanthuridae")),
                        cbind(results.wu.acanth, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Acanthuridae")),
                        cbind(results.uu.acanth, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Acanthuridae")),
                        cbind(results.bc.holo, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Holocentridae")),
                        cbind(results.wu.holo, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Holocentridae")),
                        cbind(results.uu.holo, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Holocentridae")),
                        cbind(results.bc.lab, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Labridae")),
                        cbind(results.wu.lab, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Labridae")),
                        cbind(results.uu.lab, "Distance" = rep("Unweighted Unifrac"), "Sampling" =rep("Labridae")),
                        cbind(results.bc.leth, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Lethrinidae")),
                        cbind(results.wu.leth, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Lethrinidae")),
                        cbind(results.uu.leth, "Distance" = rep("Unweighted Unifrac"), "Sampling" =rep("Lethrinidae")),
                     cbind(results.bc.lutj, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Lutjanidae")),
                     cbind(results.wu.lutj, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Lutjanidae")),
                     cbind(results.uu.lutj, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Lutjanidae")),
                     cbind(results.bc.poma, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Pomacentridae")),
                     cbind(results.wu.poma, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Pomacentridae")),
                     cbind(results.uu.poma, "Distance" = rep("Unweighted Unifrac"),"Sampling" =rep("Pomacentridae")),
                     cbind(results.bc.scar, "Distance" = rep("Bray-Curtis"),"Sampling" =rep("Scaridae")),
                     cbind(results.wu.scar, "Distance" = rep("Weighted Unifrac"), "Sampling" =rep("Scaridae")),
                     cbind(results.uu.scar, "Distance" = rep("Unweighted Unifrac"), "Sampling" =rep("Scaridae")))

write.table(results_tax, file = "results_tax2.txt", sep = "\t", quote = F)



#####################_______ REGION ______________#######------
dir.create("./region_on_host")
setwd("./region_on_host")

load("../../Physeq_objects/gut_core.RData")
samp_data <- gut_core %>% sample_data() %>% as.data.frame()

sp_jdn <- names(which(table(samp_data[samp_data$region == "Juan_de_nova",]$tax1) >=3))
sp_euro <- names(which(table(samp_data[samp_data$region == "Europa",]$tax1) >=3))
sp_sey <- names(which(table(samp_data[samp_data$region == "Seychelles",]$tax1) >=3))

sp_jdn[which(sp_jdn %in% sp_euro)]
#"Acanthurus nigrofuscus", "Ctenochaetus striatus","Halichoeres hortulanus",
#"Lutjanus bohar","Monotaxis grandoculis","Zanclus cornutus"
sp_jdn[which(sp_jdn %in% sp_sey)] # "Lutjanus kasmira"
sp_euro[which(sp_euro %in% sp_sey)]

## Acanthurus nigrofuscus ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./Ac_nigrofuscus")
setwd("./Ac_nigrofuscus")

ac_nigro <- subset_samples(gut_core, tax1 == "Acanthurus nigrofuscus")
ac_nigro <- prune_taxa(names(which(colSums(ac_nigro@otu_table)>0)), ac_nigro)
ac_nigro_rel <- transform_sample_counts(ac_nigro, function(x) x / sum(x) )
ac_nigro_dat <- data.frame(sample_data(ac_nigro_rel))

# ** Bray-curtis ----
ac_nigro_bc <- vegdist(ac_nigro_rel@otu_table, method = "bray")
ac_nigro_bc_region = adonis2(ac_nigro_bc ~ region, data = ac_nigro_dat)

# ** Weighted Unifrac ----
ac_nigro_wu <- UniFrac(ac_nigro_rel , weighted = T, normalized=F, parallel = F, fast=T)
ac_nigro_wu_region = adonis2(ac_nigro_wu ~ region, data = ac_nigro_dat)

# ** Unweighted Unifrac ----
ac_nigro_uu <- UniFrac(ac_nigro_rel, weighted = F, normalized=F, parallel = F, fast=T)
ac_nigro_uu_region = adonis2(ac_nigro_uu ~ region, data = ac_nigro_dat)

results.ac_nigro <- cbind("Bray-Curtis" = paste("R2 =",round(ac_nigro_bc_region$R2[1],3),";","p =",(ac_nigro_bc_region$`Pr(>F)`)[1]) ,
                          "W-Uni" = paste("R2 =",round(ac_nigro_wu_region$R2[1],3),";","p =",(ac_nigro_wu_region$`Pr(>F)`)[1]) ,
                          "Un-Uni" = paste("R2 =",round(ac_nigro_uu_region$R2[1],3),";","p =",(ac_nigro_uu_region$`Pr(>F)`)[1]))

save(ac_nigro, ac_nigro_rel, ac_nigro_dat, file ="ac_nigro_objects.RData")
save(ac_nigro_bc, ac_nigro_uu, ac_nigro_wu,file = "ac_nigro_matrices.RData")

# ** __ PCoA ----


## Ctenochaetus striatus ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./Ct_striatus")
setwd("./Ct_striatus")

ct_stri <- subset_samples(gut_core, tax1 == "Ctenochaetus striatus")
ct_stri <- prune_taxa(names(which(colSums(ct_stri@otu_table)>0)), ct_stri)
ct_stri_rel <- transform_sample_counts(ct_stri, function(x) x / sum(x) )
ct_stri_dat <- data.frame(sample_data(ct_stri_rel))

# ** Bray-curtis ----
ct_stri_bc <- vegdist(ct_stri_rel@otu_table, method = "bray")
ct_stri_bc_region = adonis2(ct_stri_bc ~ region, data = ct_stri_dat)

# ** Weighted Unifrac ----
ct_stri_wu <- UniFrac(ct_stri_rel , weighted = T, normalized=F, parallel = F, fast=T)
ct_stri_wu_region = adonis2(ct_stri_wu ~ region, data = ct_stri_dat)

# ** Unweighted Unifrac ----
ct_stri_uu <- UniFrac(ct_stri_rel, weighted = F, normalized=F, parallel = F, fast=T)
ct_stri_uu_region = adonis2(ct_stri_uu ~ region, data = ct_stri_dat)

results.ct_stri <- cbind("Bray-Curtis" = paste("R2 =",round(ct_stri_bc_region$R2[1],3),";","p =",(ct_stri_bc_region$`Pr(>F)`)[1]) ,
                          "W-Uni" = paste("R2 =",round(ct_stri_wu_region$R2[1],3),";","p =",(ct_stri_wu_region$`Pr(>F)`)[1]) ,
                          "Un-Uni" = paste("R2 =",round(ct_stri_uu_region$R2[1],3),";","p =",(ct_stri_uu_region$`Pr(>F)`)[1]))

save(ct_stri, ct_stri_rel, ct_stri_dat, file ="ct_stri_objects.RData")
save(ct_stri_bc, ct_stri_uu, ct_stri_wu,file = "ct_stri_matrices.RData")


## Halichoeres hortulanus ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./Ha_hortulanus")
setwd("./Ha_hortulanus")

ha_hort <- subset_samples(gut_core, tax1 == "Halichoeres hortulanus")
ha_hort <- prune_taxa(names(which(colSums(ha_hort@otu_table)>0)), ha_hort)
ha_hort_rel <- transform_sample_counts(ha_hort, function(x) x / sum(x) )
ha_hort_dat <- data.frame(sample_data(ha_hort_rel))

# ** Bray-curtis ----
ha_hort_bc <- vegdist(ha_hort_rel@otu_table, method = "bray")
ha_hort_bc_region = adonis2(ha_hort_bc ~ region, data = ha_hort_dat)

# ** Weighted Unifrac ----
ha_hort_wu <- UniFrac(ha_hort_rel , weighted = T, normalized=F, parallel = F, fast=T)
ha_hort_wu_region = adonis2(ha_hort_wu ~ region, data = ha_hort_dat)

# ** Unweighted Unifrac ----
ha_hort_uu <- UniFrac(ha_hort_rel, weighted = F, normalized=F, parallel = F, fast=T)
ha_hort_uu_region = adonis2(ha_hort_uu ~ region, data = ha_hort_dat)

results.ha_hort <- cbind("Bray-Curtis" = paste("R2 =",round(ha_hort_bc_region$R2[1],3),";","p =",(ha_hort_bc_region$`Pr(>F)`)[1]) ,
                         "W-Uni" = paste("R2 =",round(ha_hort_wu_region$R2[1],3),";","p =",(ha_hort_wu_region$`Pr(>F)`)[1]) ,
                         "Un-Uni" = paste("R2 =",round(ha_hort_uu_region$R2[1],3),";","p =",(ha_hort_uu_region$`Pr(>F)`)[1]))

save(ha_hort, ha_hort_rel, ha_hort_dat, file ="ha_hort_objects.RData")
save(ha_hort_bc, ha_hort_uu, ha_hort_wu,file = "ha_hort_matrices.RData")


## Lutjanus bohar ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./Lu_bohar")
setwd("./Lu_bohar")

lu_boh <- subset_samples(gut_core, tax1 == "Lutjanus bohar")
lu_boh <- prune_taxa(names(which(colSums(lu_boh@otu_table)>0)), lu_boh)
lu_boh_rel <- transform_sample_counts(lu_boh, function(x) x / sum(x) )
lu_boh_dat <- data.frame(sample_data(lu_boh_rel))

# ** Bray-curtis ----
lu_boh_bc <- vegdist(lu_boh_rel@otu_table, method = "bray")
lu_boh_bc_region = adonis2(lu_boh_bc ~ region, data = lu_boh_dat)

# ** Weighted Unifrac ----
lu_boh_wu <- UniFrac(lu_boh_rel , weighted = T, normalized=F, parallel = F, fast=T)
lu_boh_wu_region = adonis2(lu_boh_wu ~ region, data = lu_boh_dat)

# ** Unweighted Unifrac ----
lu_boh_uu <- UniFrac(lu_boh_rel, weighted = F, normalized=F, parallel = F, fast=T)
lu_boh_uu_region = adonis2(lu_boh_uu ~ region, data = lu_boh_dat)

results.lu_boh <- cbind("Bray-Curtis" = paste("R2 =",round(lu_boh_bc_region$R2[1],3),";","p =",(lu_boh_bc_region$`Pr(>F)`)[1]) ,
                         "W-Uni" = paste("R2 =",round(lu_boh_wu_region$R2[1],3),";","p =",(lu_boh_wu_region$`Pr(>F)`)[1]) ,
                         "Un-Uni" = paste("R2 =",round(lu_boh_uu_region$R2[1],3),";","p =",(lu_boh_uu_region$`Pr(>F)`)[1]))

save(lu_boh, lu_boh_rel, lu_boh_dat, file ="lu_boh_objects.RData")
save(lu_boh_bc, lu_boh_uu, lu_boh_wu,file = "lu_boh_matrices.RData")

## Monotaxis grandoculis ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./Mo_grandoculis", recursive = T)
setwd("./Mo_grandoculis")

mo_gra <- subset_samples(gut_core, tax1 == "Monotaxis grandoculis")
mo_gra <- prune_taxa(names(which(colSums(mo_gra@otu_table)>0)), mo_gra)
mo_gra_rel <- transform_sample_counts(mo_gra, function(x) x / sum(x) )
mo_gra_dat <- data.frame(sample_data(mo_gra_rel))

# ** Bray-curtis ----
mo_gra_bc <- vegdist(mo_gra_rel@otu_table, method = "bray")
mo_gra_bc_region = adonis2(mo_gra_bc ~ region, data = mo_gra_dat)

# ** Weighted Unifrac ----
mo_gra_wu <- UniFrac(mo_gra_rel , weighted = T, normalized=F, parallel = F, fast=T)
mo_gra_wu_region = adonis2(mo_gra_wu ~ region, data = mo_gra_dat)

# ** Unweighted Unifrac ----
mo_gra_uu <- UniFrac(mo_gra_rel, weighted = F, normalized=F, parallel = F, fast=T)
mo_gra_uu_region = adonis2(mo_gra_uu ~ region, data = mo_gra_dat)

results.mo_gra <- cbind("Bray-Curtis" = paste("R2 =",round(mo_gra_bc_region$R2[1],3),";","p =",(mo_gra_bc_region$`Pr(>F)`)[1]) ,
                         "W-Uni" = paste("R2 =",round(mo_gra_wu_region$R2[1],3),";","p =",(mo_gra_wu_region$`Pr(>F)`)[1]) ,
                         "Un-Uni" = paste("R2 =",round(mo_gra_uu_region$R2[1],3),";","p =",(mo_gra_uu_region$`Pr(>F)`)[1]))

save(mo_gra, mo_gra_rel, mo_gra_dat, file ="mo_gra_objects.RData")
save(mo_gra_bc, mo_gra_uu, mo_gra_wu,file = "mo_gra_matrices.RData")

## Zanclus cornutus ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./za_cornutus")
setwd("./za_cornutus")

za_cornu <- subset_samples(gut_core, tax1 == "Zanclus cornutus")
za_cornu <- prune_taxa(names(which(colSums(za_cornu@otu_table)>0)), za_cornu)
za_cornu_rel <- transform_sample_counts(za_cornu, function(x) x / sum(x) )
za_cornu_dat <- data.frame(sample_data(za_cornu_rel))

# ** Bray-curtis ----
za_cornu_bc <- vegdist(za_cornu_rel@otu_table, method = "bray")
za_cornu_bc_region = adonis2(za_cornu_bc ~ region, data = za_cornu_dat)

# ** Weighted Unifrac ----
za_cornu_wu <- UniFrac(za_cornu_rel , weighted = T, normalized=F, parallel = F, fast=T)
za_cornu_wu_region = adonis2(za_cornu_wu ~ region, data = za_cornu_dat)

# ** Unweighted Unifrac ----
za_cornu_uu <- UniFrac(za_cornu_rel, weighted = F, normalized=F, parallel = F, fast=T)
za_cornu_uu_region = adonis2(za_cornu_uu ~ region, data = za_cornu_dat)

results.za_cornu <- cbind("Bray-Curtis" = paste("R2 =",round(za_cornu_bc_region$R2[1],3),";","p =",(za_cornu_bc_region$`Pr(>F)`)[1]) ,
                         "W-Uni" = paste("R2 =",round(za_cornu_wu_region$R2[1],3),";","p =",(za_cornu_wu_region$`Pr(>F)`)[1]) ,
                         "Un-Uni" = paste("R2 =",round(za_cornu_uu_region$R2[1],3),";","p =",(za_cornu_uu_region$`Pr(>F)`)[1]))

save(za_cornu, za_cornu_rel, za_cornu_dat, file ="za_cornu_objects.RData")
save(za_cornu_bc, za_cornu_uu, za_cornu_wu,file = "za_cornu_matrices.RData")

## Lutjanus kasmira ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")
dir.create("./Lu_kasmira")
setwd("./Lu_kasmira")

lu_kas <- subset_samples(gut_core, tax1 == "Lutjanus kasmira")
lu_kas <- prune_taxa(names(which(colSums(lu_kas@otu_table)>0)), lu_kas)
lu_kas_rel <- transform_sample_counts(lu_kas, function(x) x / sum(x) )
lu_kas_dat <- data.frame(sample_data(lu_kas_rel))

# ** Bray-curtis ----
lu_kas_bc <- vegdist(lu_kas_rel@otu_table, method = "bray")
lu_kas_bc_region = adonis2(lu_kas_bc ~ region, data = lu_kas_dat)

# ** Weighted Unifrac ----
lu_kas_wu <- UniFrac(lu_kas_rel , weighted = T, normalized=F, parallel = F, fast=T)
lu_kas_wu_region = adonis2(lu_kas_wu ~ region, data = lu_kas_dat)

# ** Unweighted Unifrac ----
lu_kas_uu <- UniFrac(lu_kas_rel, weighted = F, normalized=F, parallel = F, fast=T)
lu_kas_uu_region = adonis2(lu_kas_uu ~ region, data = lu_kas_dat)

results.lu_kas <- cbind("Bray-Curtis" = paste("R2 =",round(lu_kas_bc_region$R2[1],3),";","p =",(lu_kas_bc_region$`Pr(>F)`)[1]) ,
                         "W-Uni" = paste("R2 =",round(lu_kas_wu_region$R2[1],3),";","p =",(lu_kas_wu_region$`Pr(>F)`)[1]) ,
                         "Un-Uni" = paste("R2 =",round(lu_kas_uu_region$R2[1],3),";","p =",(lu_kas_uu_region$`Pr(>F)`)[1]))

save(lu_kas, lu_kas_rel, lu_kas_dat, file ="lu_kas_objects.RData")
save(lu_kas_bc, lu_kas_uu, lu_kas_wu,file = "lu_kas_matrices.RData")

## RESULTS ----
setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/04_data_cleaning/Meta_analisis/beta_div/region_on_host/")

results <- as.data.frame(cbind("Species" = c("A.nigrofuscus","C.striatus","H.hortulanus",
                                                 "L.bohar","M.grandoculis","Z.cornutus", "L.kasmira"), 
                               rbind(results.ac_nigro,results.ct_stri,
                                     results.ha_hort,results.lu_boh, 
                                     results.mo_gra,results.za_cornu,results.lu_kas)))

write.table(results, file = "results.txt", sep= "\t", row.names = F,quote = F)


## Plot ----
pdf(file="PCOA.diet.pdf", wi = 7, he = 7)
pcoa.carn <- ggplot() +
  geom_hline(yintercept=0, colour="lightgrey", linetype = 2) +
  geom_vline(xintercept=0, colour="lightgrey", linetype = 2) +
  geom_point(data = hull, aes(x=Axis.1, y=Axis.2, color = diet3 ), alpha = 0.7, size = 4, shape = 17) +
  scale_color_manual(values =c("Darkred","darkgreen", "Darkblue")) +
  xlab(paste("PCo1 (", round(pcoa.sub$values$Relative_eig[1]*100, 1), "%)")) +
  ylab(paste("PCo2 (", round(pcoa.sub$values$Relative_eig[2]*100, 1), "%)"))  +
  theme_bw() +
  coord_equal() +
  theme(axis.title.x = element_text(size=16, family = "serif"), # remove x-axis labels
        axis.title.y = element_text(size=16,family = "serif"), # remove y-axis labels
        axis.text = element_text(family = "serif", size = 14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank(),
        axis.title = element_text(family = "serif", size = 16), 
        legend.text = element_text(size = 14, family = "serif"),
        legend.title = element_text(size = 14,family = "serif")) +
  facet_wrap(~ region)+
  labs(colour = "Diet")

pcoa.bray
dev.off()




















