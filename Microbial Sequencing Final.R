###Microbiome Sequencing: Final Trial to get data correct for publication
library (dada2); packageVersion("dada2") 
library(devtools)
devtools::install_github("benjjneb/dada2", ref="v1.16") # change the ref argument to get other versions
library(phangorn)
library(Biostrings)
library(phyloseq)
library(ggplot2)
library(DECIPHER); packageVersion("DECIPHER")
library(tidyverse)
library(ggtext)
library(vegan)
library(dplyr)
library(abdiv)
library(picante)
library(RColorBrewer)

##Tutorial from package designers: https://benjjneb.github.io/dada2/tutorial.html
##Youtube tutorial walkthrough with tree construction: https://www.youtube.com/watch?v=wV5_z7rR6yw&t=2574s

#****set working directory/folder to the one that contains the fastq files*******
setwd("C:/BaseSpace/MINI720_Kyle-Emerson_42188-348319972/FASTQ_Generation_2022-06-16_14_53_56Z-574196980")
path <- "C:/BaseSpace/MINI720_Kyle-Emerson_42188-348319972/FASTQ_Generation_2022-06-16_14_53_56Z-574196980" #****change to the directory containing the fastq files*****
list.files(path) #list the files within that directory/folder within the path

# Forward and reverse fastq filenames have format: 
# Forward = SAMPLENAME_R1_001.fastq and 
# Reverse = SAMPLENAME_R2_001.fastq
#This pulls each FASTQ based on naming pattern and places as Fwd or Rvs
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#Removes extra values besides the "B1A-2" & "FA-B1A"
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Visualizing quality profiles 
#This will generate a plot for FORWARD READS (R1)
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnFs[3:4])
plotQualityProfile(fnFs[5:6])
plotQualityProfile(fnFs[7:8])
#Changing the numbers will show you the quality profile of different sequences

#In gray-scale is a heat map of the frequency of each quality score at each 
#base position. 
#The mean quality score at each position is shown by the green line, and 
#the quartiles of the quality score distribution by the orange lines. 
#The red line shows the scaled proportion of reads that extend to at least 
#that position (this is more useful for other sequencing technologies, 
#as Illumina reads are typically all the same length, hence the flat red line).

#Generally advise to trim the last few nucleotides to avoid less well-controlled errors. 

#Visualize the quality profile of the reverse reads (R2) within a plot:
plotQualityProfile(fnRs[1:2])
plotQualityProfile(fnRs[3:4])
plotQualityProfile(fnRs[5:6])
plotQualityProfile(fnRs[7:8])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
#Using the assigned path above, filtered sequences are named either 
# _F_filt or _R_filt whether forward or reserve reads. 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = c(19, 20),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE, matchIDs = TRUE) # On Windows set multithread=FALSE
head(out)
# Trim left is my trimming of primers (19 is forward primer (515F) length)
# 20 is my length of the 806R primer length
# Based on quality profile plots, the only trimming required is that of my forward and reverse primers
# On Windows set multithread=FALSE
# The truncQ = 2 code, truncates the reads at the first instance of a quality 
# score less than or equal to 2. 
# The “matchIDs = true” code is used for the paired-read filtering, to enforce 
# matching between the id-line sequence identifiers of the forward and reverse 
# fastq files. Since this is true, ONLY paired reads that share sequence ID 
# fields are shown in the output. 

#Learn the Error Rates

# The DADA2 algorithm makes use of a parametric error model (err) and every 
#amplicon dataset has a different set of error rates. 
#The learnErrors method learns this error model from the data, 
#by alternating estimation of the error rates and inference of sample 
#composition until they converge on a jointly consistent solution. 
#As in many machine-learning problems, the algorithm must begin with an initial guess, 
#for which the maximum possible error rates in this data are used 
#(the error rates if only the most abundant sequence is correct and all the rest are errors).

#Error model, A<>C<>T<>G changes of nucleic bases
errF <- learnErrors(filtFs, multithread=FALSE) 
#100274141 total bases in 748509 reads from 8 samples will be used for learning the error rates.
errR <- learnErrors(filtRs, multithread=FALSE)
#113964090 total bases in 857157 reads from 9 samples will be used for learning the error rates.
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#The error rates for each possible transition (A→C, A→G, …) are shown. 
#Points are the observed error rates for each consensus quality score. 
#The black line shows the estimated error rates after convergence of the 
#machine-learning algorithm. 
#The red line shows the error rates expected under the 
#nominal definition of the Q-score. 
#Here the estimated error rates (black line) are a good fit to the 
#observed rates (points), 
#and the error rates drop with increased quality as expected. 
#Everything looks reasonable and we proceed with confidence.

######################################################
#Sample Inference

#We are now ready to apply the core sample inference algorithm to the filtered 
#and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=FALSE)
dadaRs <- dada(filtRs, err=errR, multithread=FALSE)
#Inspecting the returned dada-class object:
dadaFs[[1]]
print(dadaFs[[1]])
print(dadaFs)
#dada-class: object describing DADA2 denoising results
#Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40
#, BAND_SIZE = 16
#The DADA2 algorithm inferred 260 true sequence variants from the 18354 
#unique sequences in the first sample. There is much more to the dada-class 
#return object than this (see help("dada-class") for some info).

#####################################
##Merge paired reads
#We now merge the forward and reverse reads together to obtain the full denoised sequences. 
#Merging is performed by aligning the denoised forward reads with the 
#reverse-complement of the corresponding denoised reverse reads, and then 
#constructing the merged “contig” sequences. 
#By default, merged sequences are only output if the forward and reverse reads 
#overlap by at least 12 bases, and are identical to each other in the overlap 
#region (but these conditions can be changed via function arguments).

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
#Appears to have trimmed correctly!

#Constructing a sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab) 
# 24 samples, 6567 unique ASVs
# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))
# Considerations for your own data: Sequences that are much longer or shorter 
# than expected may be the result of non-specific priming. 
# You can remove non-target-length sequences from your sequence table 
# (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). 
# This is analogous to “cutting a band” in-silico to get amplicons of the targeted length
# Because I have some amplicons that are a bit short, maybe because of primer issues, 
# I am going to trim the desired length to 250-255
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:255]
dim(seqtab2)
#24, 6547
table(nchar(getSequences(seqtab2)))
#250  251  252   253    254   255 
# 6    5   678   5589   258   11 

#Remove chimeras
#The core dada method corrects substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of sequence variants after denoising makes identifying 
#chimeric ASVs simpler than when dealing with fuzzy OTUs. 
#Chimeric sequences are identified if they can be exactly reconstructed by 
#combining a left-segment and a right-segment from two more abundant “parent” sequences.
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", 
                                    multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)
#The frequency of chimeric sequences varies substantially from dataset to dataset, 
#and depends on on factors including experimental procedures and sample complexity. 
#Here chimeras make up about 75% (6547-1627) of the merged sequence variants, 
#but when we account for the abundances of those variants we see they account 
#for only about 8-9% of the merged sequence reads.
#Thus, this shows that ~92% of my reads are not chimeric. It is not uncommon for a large
#amount of sequence variants be removed (4920), but the majority of your reads should remain
#as mine does here, losing only about 8.5%. 

####Removal of Singleton reads
##According to tutorial and video, singleton reads are automatically filtered out
##and are not classified as unique sequences
##Our seqtab.nochim is our current ASV table (unrarefied). Currently, I see
##A few sequences at the end of the table with 1 read. Unsure if they are classified as singletons,
##But the link below is from the author and is what I will use to remove these sequences
## https://github.com/benjjneb/dada2/issues/1519

is1 <- colSums(seqtab.nochim) <= 1
seqtab.nochim1 <- seqtab.nochim[,!is1]

#####Track reads through the pipeline

#As a final check of our progress, we’ll look at the number of reads that 
#made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim1))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
list(track)
#Tracking looks good, the majority of our raw reads have successfully
#merged and have not been filtered out through any of our processing steps
#We are good to continue onward!

##### Assign Taxonomy

#It is common at this point, especially in 16S/18S/ITS amplicon sequencing, 
#to assign taxonomy to the sequence variants. The DADA2 package provides a native 
#implementation of the naive Bayesian classifier method for this purpose. 
#The assignTaxonomy function takes as input a set of sequences to be classified 
#and a training set of reference sequences with known taxonomy, 
#and outputs taxonomic assignments with at least minBoot bootstrap confidence.
#We maintain formatted training fastas for the RDP training set, 
#GreenGenes clustered at 97% identity, and the Silva reference database, 
#and additional trainings fastas suitable for protists and certain specific 
#environments have been contributed. For fungal taxonomy, the General Fasta 
#release files from the UNITE ITS database can be used as is. 
#To follow along, download the silva_nr_v132_train_set.fa.gz file, 
#These files can be downloaded here: https://benjjneb.github.io/dada2/training.html
#and place it in the directory with the fastq files.

taxa <- assignTaxonomy(seqtab.nochim1, "C:/BaseSpace/MINI720_Kyle-Emerson_42188-348319972/FASTQ_Generation_2022-06-16_14_53_56Z-574196980/silva_nr99_v138.1_train_set.fa.gz", multithread=FALSE)
#Assigning the taxonomy and adding the species level will take a considerable
#amount of time. Best to run the code and check back later. 

write.table(taxa, file = "EM1_taxa.txt", sep="\t")
#Creates a txt file in the path listed of the taxa information. Note:
#this goes into your designated R folder (for me: Emerson-Microbial-1)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#### Species Level Assignment (not doing for this project) 
##Secondary file to assign species will need downloaded
# taxaSpecies <- addSpecies(taxa, "C:/BaseSpace/MINI720_Kyle-Emerson_42188-348319972/FASTQ_Generation_2022-06-16_14_53_56Z-574196980/silva_species_assignment_v138.1.fa.gz")
# write.table(taxaSpecies, "C:/BaseSpace/MINI720_Kyle-Emerson_42188-348319972/FASTQ_Generation_2022-06-16_14_53_556Z-574196980", sep="\t")
# 
# taxaSpecies.print <- taxaSpecies # Removing sequence rownames for display only
# rownames(taxaSpecies.print) <- NULL
# list(taxaSpecies.print)
# ##Species level!!


#####Making a phylogenetic tree
## Note: this is strictly for phylogenetic analysis
## I.e. Faiths. Will not impact ASVs, Shannon, Bray Curtis, etc
## Information found in the youtube tutorial listed at beginning of doc

# Align Sequences for Phylogeny
# Extract sequences from DADA2 output
sequences <- getSequences(seqtab.nochim1)
names(sequences)<-sequences
# Run sequence alignment (MSA) using DECIPHER
alignment <- AlignSeqs(DNAStringSet(sequences), anchor = NA)


phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# Change sequence alignment output into a phyDat structure

dm <- dist.ml(phang.align)
# Making a distance matrix

treeNJ <- NJ(dm)
# Making a neighbor joining tree

fit = pml(treeNJ, data = phang.align)
# Fitting for internal maximum likelihood

fitGTR <- update(fit, k =4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

### Have my tables (unrarefied) and my tree
theme_set(theme_bw())
file.choose()
meta3 <- read.csv("C:\\R\\Emerson-Microbial-1\\Excel Files\\EM1_MicrobialMap.csv")
#This is my map file. This file is necessary to give R a point of reference
#to merge files downstream. The most important consistency is including the
#sample names that are identical across all files

row.names(meta3) <- meta3$Sample_ID
#ID needs to be row names for correct merging of files

meta3$Microbial_Trtmt = factor(meta3$Microbial_Trtmt)
meta3$Replicate = factor(meta3$Replicate)

ps <- phyloseq(otu_table(seqtab.nochim1, taxa_are_rows=FALSE), 
               tax_table(taxa), sample_data(meta3),
               phy_tree(fitGTR$tree))

###### Rooting my phylogenetic tree

set.seed(711)
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(ps))

##### Creating my phyloseq object for analysis

ps <- merge_phyloseq(ps, map)
#Merge ps object with map
#Creates a phyloseq object that uses sequences with chimeras removed from the 
#sequencing run, data from the meta table, and information on taxa from the downloaded
#Silva file 
#had to change sample data to sample names

taxa_names(ps)

wholetax <- do.call(paste, c(as.data.frame(tax_table(ps))
                             [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")],
                             sep = "__"))  # to distinguish from "_" within tax ranks
#generate a vector containing the full taxonomy path for all ASVs

otu_export <- as.data.frame(otu_table(ps))
tmp <- names(otu_export)
# turn the otu_table into a data.frame

for(i in 1:length(tmp)){
  names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}
# paste wholetax and OTU_ids together

names(otu_export) <- names(tmp)
head(otu_export)[5]
#overwrite old names

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
list(ps)
#1605 ASVs with all singletons and chimeric sequences removed!

refseqPS <-refseq(ps)
#return non-empty slot names of phyloseq object

getslots.phyloseq(ps)
#we have "otu_table" "tax_table" "sam_data" "refseq" "phy_tree"
#our completed phyloseq object

########### Cleaning up my data in R (https://www.youtube.com/watch?v=e3rKYipvdJo)
write.table(otu_table(ps), file = "EM1_otutable.csv", sep = ",")
write.table(tax_table(ps), file = "EM1_taxatable.csv", sep = ",")
write.table(refseq(ps), file = "EM1_refseqtable.csv", sep = ",")

ps_ra = transform_sample_counts(ps,function(x){x/sum(x)})
write.table(otu_table(ps_ra), file = "EM1_relabundtable.csv", sep = ",")
#This is our relative abundance table. 

file.choose()
otu_counts <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 2 - Emerson Microbial 1\\Emersion Microbial Experiment 1 (2021)\\Microbiome Analysis\\Microbiome Analysis Final - Outlier included\\EM1_otutable.csv") %>%
  pivot_longer(-Sample_ID, names_to="ASV", values_to="count")

taxonomy <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 2 - Emerson Microbial 1\\Emersion Microbial Experiment 1 (2021)\\Microbiome Analysis\\Microbiome Analysis Final - Outlier included\\EM1_taxatable.csv", fileEncoding="UTF-8-BOM") %>%
  filter(Kingdom != 'Archaea')

#This filter step allowed me to remove any archaea from the df
#fileEncoding="UTF-8-BOM" can be used when Excel -> CSV gives weird symbols in column headers

#Now, I want to join my dataframes together
#meta3
#otu_counts
#taxonomy
################ Relative Abundance
# https://www.statology.org/dplyr-remove-rows/#:~:text=You%20can%20use%20the%20following%20basic%20syntax%20to,4%29%29%205%205.%20Remove%20rows%20based%20on%20condition

otu_rel_abund <- inner_join(meta3, otu_counts, by = "Sample_ID") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Sample_ID) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "ASV"),
               names_to="level",
               values_to= "taxon") %>%
  filter(level != 'ASV')


##Now, our data frame only includes bacteria

############### Relative Abundance data frame: Phylum
phyla_abundance <- otu_rel_abund %>%
  select(-Microbial_Trtmt, -Replicate) %>%
  filter(level =="Phylum") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  pivot_wider(names_from = "taxon", values_from = "rel_abund")

write.table(phyla_abundance, file = "EM1 Final phyla_abundance.csv", sep = "," )  

############### Relative Abundance data frame: Genus
genus_abundance <- otu_rel_abund %>%
  select(-Microbial_Trtmt, -Replicate) %>%
  filter(level =="Genus") %>%
  group_by(Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund)) %>%
  pivot_wider(names_from = "taxon", values_from = "rel_abund")

write.table(genus_abundance, file = "EM1 Final genus_abundance.csv", sep = "," )

#Check to make sure our relative abundances add up to 1
otu_rel_abund %>%
  group_by(Sample_ID) %>%
  summarize(total = sum(rel_abund))
#Success, says 6 but that is because each taxa is counted 6 times due to 
#that being the depth of taxonomic assignment

## https://www.youtube.com/watch?v=NVym44SdcaE&t=309s 
## Making stacked barcharts for relative abundance

Pond_Water <- c("Natural", "Autoclaved")

otu_rel_abund %>%
  filter(Sample_ID != 'Kyle024-PCR1Blank') %>%
  filter(level =="Phylum") %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  ggplot(aes(x = Microbial_Trtmt, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_discrete(name = NULL) +
  labs(x = "Pond Water", y = "Mean Relative Abundance (%)") +
  scale_x_discrete(labels = Pond_Water) +
  theme_classic() +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund <- otu_rel_abund %>%
  filter(Sample_ID != 'Kyle024-PCR1Blank') %>%
  filter(level =="Phylum") %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool <- taxon_rel_abund %>%
  group_by(taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop") 

inner_join(taxon_rel_abund, taxon_pool, by = "taxon") %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = sum(mean_rel_abund),
            mean = min(mean), .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = Microbial_Trtmt, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Actinobacteriota", "Bacteroidota",
                               "Firmicutes", "Fusobacteriota",
                               "Proteobacteria", "Other"),
                    values = c(brewer.pal(5, "Accent"), "gray")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Pond Water", y = "Mean Relative Abundance (%)") +
  scale_x_discrete(labels = Pond_Water) +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(face = "bold", size = 14))  +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

#Final Mean Relative Abundance chart

### Now, we still want mean relative abundance charts for individuals within each treatment
## Natural Treatment

otu_rel_abund %>%
  filter(Sample_ID != 'Kyle024-PCR1Blank') %>%
  filter(level =="Phylum") %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop") %>%
  ggplot(aes(x = Microbial_Trtmt, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_discrete(name = NULL) +
  labs(x = "Pond Water", y = "Mean Relative Abundance (%)") +
  scale_x_discrete(labels = Pond_Water) +
  theme_classic() +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(10, "pt"))

taxon_rel_abund_natural <- otu_rel_abund %>%
  filter(Sample_ID != 'Kyle024-PCR1Blank') %>%
  filter(level =="Phylum") %>%
  group_by(Microbial_Trtmt, Sample_ID, taxon) %>%
  summarize(rel_abund = sum(rel_abund), .groups="drop") %>%
  group_by(Sample_ID, Microbial_Trtmt, taxon) %>%
  summarize(mean_rel_abund = 100*mean(rel_abund), .groups="drop")

taxon_pool_natural <- taxon_rel_abund_natural %>%
  group_by(Sample_ID, Microbial_Trtmt, taxon) %>%
  summarize(pool = max(mean_rel_abund) < 1, 
            mean = mean(mean_rel_abund),
            .groups = "drop") 

taxon_pool_natural %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Sample_ID, taxon) %>%
  filter(Microbial_Trtmt == 1) %>%
  summarize(mean_rel_abund = sum(mean),
            .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean_rel_abund, .desc = TRUE),
         taxon = fct_shift(taxon, n = 1)) %>%
  ggplot(aes(x = Sample_ID, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Actinobacteriota", "Bacteroidota", "Chloroflexi",
                               "Cyanobacteria" , "Firmicutes",
                               "Fusobacteriota", "Myxococcota", 
                               "Planctomycetota", "Proteobacteria",
                               "Verrucomicrobiota", "NA", "Other"),
                    values = c(brewer.pal(12, "Spectral"), "gray")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Natural Treatment Bins", y = "Mean Relative Abundance (%)", las = 2) +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(face = "bold", size = 14))  +
  scale_x_discrete(labels = c('1', '2', '3', 
                              '4', '5', '6',
                              '7', '8', '9',
                              '10', '11', '12')) +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))


brewer.pal(n = 11, name = "Spectral")

natural_rel_abund <- taxon_pool_natural %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Sample_ID, taxon) %>%
  filter(Microbial_Trtmt == 1) %>%
  summarize(mean_rel_abund = sum(mean),
            .groups = "drop")

write.table(natural_rel_abund, file = "EM1 Final Natural Bar Chart Abundance.csv", sep = "," )

####Autoclaved

taxon_pool_natural %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Sample_ID, taxon) %>%
  filter(Microbial_Trtmt == 2) %>%
  summarize(mean_rel_abund = sum(mean),
            .groups = "drop") %>%
  mutate(taxon = factor(taxon),
         taxon = fct_reorder(taxon, mean_rel_abund, .desc = TRUE),
         taxon = fct_shift(taxon, n = 2)) %>%
  ggplot(aes(x = Sample_ID, y=mean_rel_abund, fill=taxon)) +
  geom_col() +
  scale_fill_manual(name = NULL,
                    breaks = c("Actinobacteriota", "Bacteroidota",
                               "Bdellovibrionota",  "Cyanobacteria" , 
                               "Firmicutes", "Fusobacteriota", 
                               "Proteobacteria", "Other"),
                    values = c("#9E0142", "#D53E4F", "#3288BD",
                               "#FDAE61", "#FEE08B", "#FFFFBF", 
                               "#66C2A5", "gray"))+
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Autoclaved Treatment Bins", y = "Mean Relative Abundance (%)") +
  theme_classic() +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text = element_text(face = "bold", size = 14))  +
  scale_x_discrete(labels = c('1', '2', '3', 
                              '4', '5', '6',
                              '7', '8', '9',
                              '10', '12')) +
  theme(axis.title = element_text(face = "bold", size = 16)) +
  theme(legend.text = element_text(face="italic"),
        legend.key.size = unit(15, "pt"))

natural_rel_abund2 <- taxon_pool_natural %>%
  mutate(taxon = if_else(pool, "Other", taxon)) %>%
  group_by(Sample_ID, taxon) %>%
  filter(Microbial_Trtmt == 2) %>%
  summarize(mean_rel_abund = sum(mean),
            .groups = "drop")

write.table(natural_rel_abund2, file = "EM1 Final Autoclaved Bar Chart Abundance.csv", sep = "," )

#######Rarefaction

rarefy_asv <- inner_join(meta3, otu_counts, by = "Sample_ID") %>%
  inner_join(., taxonomy, by = "ASV") %>%
  group_by(Sample_ID) %>%
  select(Sample_ID, Microbial_Trtmt, Replicate, ASV, count)
# This above command turned our otu_counts table into a new table
# that has archaea removed and is strictly bacteria

sampling_coverage <- rarefy_asv %>%
  group_by(Sample_ID) %>%
  summarize(n_seqs = sum(count)) 
#gives us the number of sequences in each sample
#Lowest number in an actual sample is 26142 reads

sampling_coverage %>%
  ggplot(aes(x = n_seqs)) +
  geom_histogram(binwidth = 5000) +
  coord_cartesian(xlim = c(0, 90000))
# Helps visualize if there are any breaks in my data that
# will allow me to rarefy my sequences

sampling_coverage %>%
  ggplot(aes(x = 1, y = n_seqs)) +
  geom_jitter()

sampling_coverage %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1: nrow (.), y = n_seqs)) +
  geom_line()

sampling_coverage %>%
  arrange(n_seqs) %>%
  print(n = 24)
# Based on this, my candidate threshold is about 70,000 sequences where there is
# A nice break in the data
# But, SF rarefied to the lowest sequence number (~14K)
# So I will rarefy to 26142, which will cause me to drop off the blank 
# while still keeping the ASVs found in the kit (and subsequently in any sample)
# in the analysis

rarefied_asv_table <- rarefy_asv %>%
  group_by(Sample_ID) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 26142) %>%
  select(-n)  %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  column_to_rownames("Sample_ID")
# Our new, rarefied (26142reads) ASV counts table without archaea

#write.table(rarefied_asv_table, file = "EM1 rarefied asv table.csv", sep = "," )

##### Alpha Diversity Metrics
## No. Observed ASVs

observed_asvs <- rarefied_asv_table %>%
  select(-Microbial_Trtmt, -Replicate) %>%
  rarefy(sample = 26142) %>%
  as_tibble(rownames = "Sample_ID") %>%
  select(Sample_ID, Obs_ASVs=value)

#write.table(observed_asvs, file = "EM1 no. observed asvs.csv", sep = "," )

## Shannon Diversity Index

rarefied_asv_table %>%
  select(-Microbial_Trtmt, -Replicate) %>%
  rrarefy(sample = 26142) %>%
  diversity()

#Now, we have shannon diversity values for our samples
#But, we want to repeat this to get the average
#So we are going to repeat this process, average that output for our rarefied shannon

shannon_iteration <- function(){
  
  rarefied_asv_table %>%
    select(-Microbial_Trtmt, -Replicate) %>%
    rrarefy(sample = 26142) %>%
    diversity()
  
}

rarefied_shannon <- replicate (100, shannon_iteration()) %>% 
  as_tibble(rownames = "Sample_ID", .name_repair = "unique" ) %>%
  pivot_longer(-Sample_ID) %>%
  group_by(Sample_ID) %>%
  summarize(shannon = mean(value))

write.table(rarefied_shannon, file = "EM1 shannon diversity.csv", sep = "," )

## Faiths Phylogenetic Diversity
# #https://search.r-project.org/CRAN/refmans/abdiv/html/faith_pd.html

faiths.df <- rarefied_asv_table %>%
  select(-Microbial_Trtmt, -Replicate)

faiths.df <- as.matrix(faiths.df)
#Turned our rarefied asv table that we have been using into a matrix

Faith_PD <- pd(faiths.df, phy_tree(ps), include.root = TRUE) 
#These are our faiths phylogenetic diversity values!
#SR represents species richness values, but we will keep our other ones

#write.table(Faith_PD, file = "EM1 Final FaithsPD.csv", sep = "," )

####### Beta Diversity
##Bray Curtis
#https://www.youtube.com/watch?v=G5Qckqq5Erw create PCoA plot
#https://www.youtube.com/watch?v=xyufizOpc5I nmds plot

distance_matrix <- rarefy_asv %>%
  group_by(Sample_ID) %>%
  select(-Microbial_Trtmt, -Replicate) %>%
  mutate(n = sum(count)) %>%
  ungroup() %>%
  filter(n >= 26142) %>%
  ungroup() %>%
  group_by(ASV) %>%
  select(-n) %>%
  mutate(total = sum(count)) %>%
  filter(total !=0) %>%
  ungroup() %>%
  select(-total) %>%
  pivot_wider(names_from = "ASV", values_from = "count") %>%
  as.data.frame()

rownames(distance_matrix) <- distance_matrix$Sample_ID
distance_matrix <- distance_matrix[,-1]
distance_matrix <- as.matrix(distance_matrix)

set.seed(19950406)
dist <- avgdist(distance_matrix, dmethod = "bray", sample = 26142)
set.seed(17)
nmds <- metaMDS(dist)

metadata_nmds <- nmds$points %>%
  as_tibble(rownames = "Sample_ID") %>%
  inner_join(., meta3, by = "Sample_ID") 

metadata_nmds$Microbial_Trtmt = factor(metadata_nmds$Microbial_Trtmt)

Pond_Water <- c("Natural", "Autoclaved")

metadata_nmds %>%
  ggplot(aes(x=MDS1, y=MDS2, color = Microbial_Trtmt)) +
  geom_point(aes(shape = Microbial_Trtmt, size = 3)) +
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 14))  +
  theme(axis.title = element_text(face = "bold",size = 16)) +
  labs(x = "Bray-Curtis NMDS 1", y = "Bray-Curtis NMDS 2") +
  scale_color_manual(labels = c("Natural", "Autoclaved"),
                     values = c("seagreen3", "skyblue3"),
                     name = "Pond Water") +
  scale_shape_manual(values = c(16,17)) +
  theme(legend.text = element_text(size = 10),
        legend.title = element_text( size = 12),
        legend.position = "none")

##Plot of our Bray Curtis NMDS (rarefied)
## Bray Curtis Analysis: https://www.youtube.com/watch?v=oLf0EpMJ4yA
adonis2(dist~metadata_nmds$Microbial_Trtmt)
#Essentially our PERMANOVA. Significant!
bd <- betadisper(dist, metadata_nmds$Microbial_Trtmt)
anova(bd)
#Our Permdisp, nonsignificant