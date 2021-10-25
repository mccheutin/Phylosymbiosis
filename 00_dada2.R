################################################################################
#  __  __ _____ _____  ___  _____  ___  _   _ __   _  ___
# |  \/  | ____|  ___|/ _ \|  ___|/ _ \| | | |  \ | |/ _ \
# | |\/| |  _| | | __  |_| | |_  | |_| | | | |   \| | |_| |
# | |  | | |___| |_| | | | |  _| | | | | |_| | |\   | | | |
# |_|  |_|_____|_____|_| |_|_|   |_| |_|_____|_| \__|_| |_|

################################################################################
##                               HOW TO PROCEED TO YOUR SEQUENCING ANALISIS   ##
################################################################################

## Preparation of your working space with the library packaging ####
 
rm(list=ls())
cran_packages   <- c("knitr", "phyloseqGraphTest", "phyloseq", "shiny", "microbiome",
"tidyverse", "miniUI", "caret", "pls", "e1071", "ggplot2",
"randomForest","entropart", "vegan", "plyr", "dplyr","here", "ggrepel", "nlme", "R.utils", "gridExtra", "googledrive",
"googlesheets", "phangorn", "devtools", "rmarkdown", "sys",
"reshape2", "devtools", "PMA","structSSI","ade4", "ape",
"Biostrings", "igraph", "ggnetwork", "intergraph", "ips",
"scales", "kableExtra", "pgirmess", "treemap")
github_packages <- c("jfukuyama/phyloseqGraphTest")
bioc_packages   <- c("phyloseq", "genefilter", "impute", "dada2", "DECIPHER")
inst <- cran_packages %in% installed.packages()
if (any(! inst)) {
    install.packages(cran_packages[!inst], repos = "http://cran.rstudio.com/") }
inst <- github_packages %in% installed.packages()
if (any(! inst)) {
    devtools::install_github(github_packages[!inst]) }
inst <- bioc_packages %in% installed.packages()
if (any(! inst)) {
    source("http://bioconductor.org/biocLite.R")
    BiocManager::install(bioc_packages[!inst]) }
sapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
library(dada2)
library(stringr)

################################################################################
###                             setup                         ####
################################################################################
# Create the following folders
path <- setwd("/Users/marie-charlottecheutin/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/")
dir_analyses <- paste0(path,"/analyses/")
dir_sample_selection <- paste0(dir_analyses, "01_select_samples/")
dir_seq_processing   <- paste0(dir_analyses, "02_process_sequences/")
dir_taxa_assign      <- paste0(dir_analyses, "03_assign_taxonomy/")
dir_data_cleaning    <- paste0(dir_analyses, "04_data_cleaning/")
dir_data_source <- paste0(path, "/data_sources/")
dir_primers   <- paste0(dir_data_source, "primers_sequences/")
dir_refdb   <- paste0(dir_data_source, "reference_databases/")

# ------------------------------------------------------------------------------
# path to the folder containing sequences
# names of and paths to the different runs
dir_fastq_source <- paste0(dir_data_source, "sequences/")
dir.create(dir_fastq_source) # paste your sequences in the folder
nms_seq_runs     <- list.files(dir_fastq_source)
# path to the folder containing reference databases
# names of and paths to the different reference databases
dir_refdb   <- paste0(dir_data_source, "reference_databases/")
dir.create(dir_refdb)
nms_refdb   <- list.files(dir_refdb)

# path to the primer sequences
# names of and paths to the different primer sequences
dir_primers   <- paste0(dir_data_source, "primers_sequences/")
dir.create(dir_primers)
nms_primers   <- list.files(dir_primers)

dir.create(dir_scripts)
dir.create(dir_analyses)
dir.create(dir_sample_selection)
dir.create(dir_seq_processing)
dir.create(dir_taxa_assign)
dir.create(dir_data_cleaning)
dir.create(dir_data_source)

# ------------------------------------------------------------------------------


################################################################################
###                             dada2_process                 ####
################################################################################
# to proceed to your sequences

## Put you selected sequences in the right folder "WP5.3_data/analyses/01_select_samples
# (first test with a dozen of R1R2 to test some parameters)
# - Get the list of your samples
fns <- sort(list.files(dir_sample_selection, full.names = TRUE))
fns <- fns[str_detect(basename(fns), ".fastq")]
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]
if(length(fns_R1) != length(fns_R2)) stop("Forward and reverse files do not match.")
# - Exctract the sample names
sample_names <- sapply(strsplit(basename(fns_R1), "_"), `[`, 1)
sample_names

# - Remove primers sequences to your data ####
# Note: put a text.file for each primer sequence (primer_name) in your dir_primers folder
primer_set_fwd <- read_lines(paste0(dir_primers, "primer_", "515F-Y" , ".txt"))
primer_set_rev <- read_lines(paste0(dir_primers, "primer_", "926R", ".txt"))
primer_length_fwd <- str_length(primer_set_fwd[1])
primer_length_rev <- str_length(primer_set_rev[1])
# - See your quality profiles of reads
plotQualityProfile(fns_R1[1])
plotQualityProfile(fns_R2[1])
# - Filter and Trim
dir_filtered <- paste0(dir_seq_processing, "fastq_filtered/")
filt_R1 <- str_c(dir_filtered, sample_names, "_R1_filt.fastq")
filt_R2 <- str_c(dir_filtered, sample_names, "_R2_filt.fastq")
names(filt_R1) <- sample_names
names(filt_R2) <- sample_names
set.seed(1000)

# - Filter the forwards and reverse reads ####
out <- filterAndTrim(fns_R1, filt_R1, fns_R2, filt_R2, truncLen=c(240,240),
maxN=0, maxEE=c(2,2), truncQ=2, trimLeft=c(primer_length_fwd,primer_length_rev) , rm.phix=TRUE,
compress=TRUE, multithread=TRUE)
head(out)

sample_names <- sapply(strsplit(basename(filt_R1), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample_namesR <- sapply(strsplit(basename(filt_R2), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample_names, sample_namesR)) stop("Forward and reverse files do not match.")
names(filt_R1) <- sample_names
names(filt_R2) <- sample_names
set.seed(1000)
# - Learn the error rates
errF <- learnErrors(filt_R1,nbases=1e8, multithread=TRUE)
errR <- learnErrors(filt_R2, nbases=1e8, multithread=TRUE)
plotErrors_F <- plotErrors(errF, nominalQ=TRUE)
plotErrors_R <- plotErrors(errR, nominalQ=TRUE)

ggsave(plotErrors_F, file = paste0(dir_seq_processing, "/quality_pdf/plots/" , "plotErrors_F.png"))
ggsave(plotErrors_R, file = paste0(dir_seq_processing, "/quality_pdf/plots/" , "plotErrors_R.png"))
### If you have big data, please proceed directly at the BIG DATA script
# - Dereplication ####
```{r Dereplicate the filtered fastq files, echo=F, message=F, warning=FALSE, include=F}
derepFs <- derepFastq(filt_R1, verbose=TRUE)
derepRs <- derepFastq(filt_R2, verbose=TRUE)
names(derepFs) <- sample_names
names(derepRs) <- sample_names

# - Sample inference ####
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = T)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = T)
dadaFs[[1]]
dadaRs[[1]]

# - Merge paire reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# - Construct sequence table
seqtab <- makeSequenceTable(mergers)
row.names(seqtab) = sample_names
dim(seqtab)
table(nchar(getSequences(seqtab)))

# - Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, paste0(dir_seq_processing,"/seqtab.nochim.rds"))
# - Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample_names
head(track)
save(track, file = "track_A515F-Y_926R.RData")

################################## THIS IS FOR BIG DATA ##########################################
# Inference and merge pair-end reads for BIG DATA
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names
for(sam in sample_names) {
    cat("Processing:", sam, "\n")
    derepFs <- derepFastq(filt_R1[[sam]])
    ddF <- dada(derepFs, err=errF, multithread=TRUE)
    derepRs <- derepFastq(filt_R2[[sam]])
    ddR <- dada(derepRs, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepFs, ddR, derepRs)
    mergers[[sam]] <- merger
}
rm(derepFs); rm(derepRs)


seqtab <- makeSequenceTable(mergers)
row.names(seqtab) = sample_names
dim(seqtab)
table(nchar(getSequences(seqtab)))

saveRDS(seqtab, paste0(dir_seq_processing,"/seqtab.rds"))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "tabled")
rownames(track) <- sample_names
head(track)
save(track, file = "table_track_A515F-Y_926R.RData")

# Merge multiple runs (if necessary)
st1 <- readRDS(paste0(dir_seq_processing,"seqtab_SEY.rds"))
st2 <- readRDS(paste0(dir_seq_processing,"seqtab_YCA.rds"))
st3 <- readRDS(paste0(dir_seq_processing,"seqtab_MOI.rds"))

seqtab <- mergeSequenceTables(st1, st2, st3)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
saveRDS(seqtab.nochim, paste0(dir_seq_processing,"seqtab.nochim.rds"))
################################## END OF BIG DATA ##########################################

# Taxonomy assignment ####
load(dir_seq_processing , "seqtab.nochim.rds")
dir_refdb   <- paste0(dir_data_source, "reference_databases/") #where you put all your reference database

#get the entire path to put in the function
path_reference_db <- paste0(dir_refdb, "silva_nr_v132_train_set.fa.gz")
path_species_db   <- paste0(dir_refdb, "silva_species_assignment_v132.fa.gz")

taxaRC <- assignTaxonomy(seqtab.nochim, path_reference_db, tryRC=TRUE)
taxaSp <- addSpecies(taxaRC, path_species_db)

taxa.print <- taxaSp # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

saveRDS(seqtab.nochim,taxaRC, taxaSp, file= paste0(dir_seq_processing,"dada2_files.rds"))
################################################################################
################     phyloseq_process  ################
################################################################################
# Empty and prepare the workspace
rm(list=ls()) # in order to free r stack memory
load("~/Drive/Thesis/MEGAFAUNA_GUT/WP5_seq/WP5.3_data/analyses/03_assign_taxonomy/dada2_files.rds")
load(paste0(dir_taxa_assign,"seqtab.nochim_515F-Y.926R.RData")

## Phyloseq object creation
library(phyloseq)
library(Biostrings)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
tax_table(taxaRC))
# Change name of the sequences in ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
# First results
ntaxa(ps)
nsamples(ps)
sample_names(ps)[1:5]
rank_names(ps)
otu_table(ps)[1:5, 1:5]
tax_table(ps)[1:5, 1:6]
#saving the files
save(ps, taxaRC,seqtab.nochim, seqtab, file=paste0(dir_taxa_assign ,"ps_515F-Y.926R.RData"))

################################################################################
##                             Merge your environmental table            #####
################################################################################
# Merge new data with current phyloseq object:
load(paste0(dir_taxa_assign,"ps_515F-Y.926R.RData"))
envdata <- read.table(paste0(dir_refdb,"env.txt") ,sep="\t",row.names = 1,fill = T, header=T)
DAT <- sample_data(envdata)
DAT
sort(sample_names(DAT)) == sort(sample_names(ps))

ps1<- merge_phyloseq(ps, DAT)
ps1
save(ps1, envdata,file=paste0(dir_taxa_assign ,"ps1_515F-Y.926R.RData"))

# Sample Minimum OTU
min(rowSums(ps1@otu_table@.Data))
readsumsdf = data.frame(nreads = sort(taxa_sums(ps1),TRUE), sorted = 1:ntaxa(ps1), type = "ASVs")
readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ps1), TRUE), sorted = 1:nsamples(ps1), type = "Samples"))

title = "Total number of reads"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")

# Sample selection ####
load(paste0(dir_taxa_assign, "ps1_515F-Y.926R.RData"))
sample_variables(ps1)

# Fasta importation
refseq <- ps1@refseq
Biostrings::writeXStringSet(ps1, "ps1.fasta")
################################################################################
##                             END OF THE PHYLOSEQ CREATION                 ####
################################################################################

