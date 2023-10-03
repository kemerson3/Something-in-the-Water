#Microbial Sequencing Final - MaAsLin Analysis
library(lme4)
library(car)
library(readr)
library(moments)
library(psych)
library(pastecs)
library(ggplot2)
library(Maaslin2)

##import metadata and phyla/genus abundance data

file.choose()
Maslin_meta <- read.table("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 2 - Emerson Microbial 1\\Emersion Microbial Experiment 1 (2021)\\Text Files\\EmersonMicrobialExp1_Batch1.txt", sep = '\t',
                   header=T)
Maslin_meta$Microbial_Trtmt = factor(Maslin_meta$Microbial_Trtmt)
Maslin_meta$Replicate = factor(Maslin_meta$Replicate)

Genera <- read.table("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 2 - Emerson Microbial 1\\Emersion Microbial Experiment 1 (2021)\\Text Files\\EM1 Final genus_abundance.txt", sep = '\t',
                     header=T)
Phyla <- read.table("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 2 - Emerson Microbial 1\\Emersion Microbial Experiment 1 (2021)\\Text Files\\EM1 Final phyla_abundance.txt", sep = '\t',
                    header=T)

##create tmp files to work around bug: "https://github.com/biobakery/Maaslin2/issues/1"
library(tidyverse)
tmp_data_phyla = tempfile(pattern = "data")
write_delim(Phyla, tmp_data_phyla, delim = "\t")

tmp_data_gen = tempfile(pattern = "data")
write_delim(Genera, tmp_data_gen, delim = "\t")

tmp_metadata = tempfile(pattern = "metadata")
write_delim(Maslin_meta, tmp_metadata, delim = "\t")

##Next step will not work if the Sample_IDs do not match
##Need to go into excel files for phyla and genera, make sure Sample_IDs
##Match in all 3 files. Save new txt files if need be.


##Run at phyla level
Maaslin_phyla_brains <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_treatment_output', 
  transform = "AST",
  fixed_effects = c('Microbial_Trtmt'),
  standardize = FALSE)

##Run at genus level
Maaslin_gen_brains <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_treatment_output', 
  transform = "AST",
  fixed_effects = c('Microbial_Trtmt'),
  standardize = FALSE)

#Still more to do regarding MaAsLin ~ Brain and Behavior measurements
#This above analysis included the outlier, at the recommendation of KDK
#Want to double back and make sure that this is correct before advancing, as it can impact results & interpretations


####Relative Brain Mass
##Run at phyla level (Mass)
Maaslin_phyla_brainshape1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_brainmass', 
  transform = "AST",
  fixed_effects = c('MA_BrainMass'),
  standardize = FALSE)

##Run at genus level (Mass)
Maaslin_genus_brains <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_brainmass', 
  transform = "AST",
  fixed_effects = c('MA_BrainMass'),
  standardize = FALSE)


####Relative Brain Shape
##Run at phyla level (PC1)
Maaslin_phyla_brainshape1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_brainshape_PC1', 
  transform = "AST",
  fixed_effects = c('RC1'),
  standardize = FALSE)

##Run at genus level (PC1)
Maaslin_genus_brainshape1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_brainshape_PC1', 
  transform = "AST",
  fixed_effects = c('RC1'),
  standardize = FALSE)

##Run at phyla level (PC2)
Maaslin_phyla_brainshape2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_brainshape_PC2', 
  transform = "AST",
  fixed_effects = c('RC2'),
  standardize = FALSE)

##Run at genus level (PC2)
Maaslin_genus_brainshape2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_brainshape_PC2', 
  transform = "AST",
  fixed_effects = c('RC2'),
  standardize = FALSE)

##Run at phyla level (PC3)
Maaslin_phyla_brainshape3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_brainshape_PC3', 
  transform = "AST",
  fixed_effects = c('RC3'),
  standardize = FALSE)

##Run at genus level (PC3)
Maaslin_genus_brainshape3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_brainshape_PC3', 
  transform = "AST",
  fixed_effects = c('RC3'),
  standardize = FALSE)

####Behavior - Visual Empty
##Run at phyla level (PC1)
Maaslin_phyla_visempty1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_visempty_PC1', 
  transform = "AST",
  fixed_effects = c('VisEmpty_RC1'),
  standardize = FALSE)

##Run at genus level (PC1)
Maaslin_genus_visempty1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_visempty_PC1', 
  transform = "AST",
  fixed_effects = c('VisEmpty_RC1'),
  standardize = FALSE)

##Run at phyla level (PC2)
Maaslin_phyla_visempty2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_visempty_PC2', 
  transform = "AST",
  fixed_effects = c('VisEmpty_RC2'),
  standardize = FALSE)

##Run at genus level (PC2)
Maaslin_genus_visempty2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_visempty_PC2', 
  transform = "AST",
  fixed_effects = c('VisEmpty_RC2'),
  standardize = FALSE)

##Run at phyla level (PC3)
Maaslin_phyla_visempty3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_visempty_PC3', 
  transform = "AST",
  fixed_effects = c('VisEmpty_RC3'),
  standardize = FALSE)

##Run at genus level (PC3)
Maaslin_genus_visempty3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_visempty_PC3', 
  transform = "AST",
  fixed_effects = c('VisEmpty_RC3'),
  standardize = FALSE)


####Behavior - Visual Food
##Run at phyla level (PC1)
Maaslin_phyla_visfood1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_visfood_PC1', 
  transform = "AST",
  fixed_effects = c('VisFood_RC1'),
  standardize = FALSE)

##Run at genus level (PC1)
Maaslin_genus_visfood1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_visfood_PC1', 
  transform = "AST",
  fixed_effects = c('VisFood_RC1'),
  standardize = FALSE)

##Run at phyla level (PC2)
Maaslin_phyla_visfood2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_visfood_PC2', 
  transform = "AST",
  fixed_effects = c('VisFood_RC2'),
  standardize = FALSE)

##Run at genus level (PC2)
Maaslin_genus_visfood2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_visfood_PC2', 
  transform = "AST",
  fixed_effects = c('VisFood_RC2'),
  standardize = FALSE)

##Run at phyla level (PC3)
Maaslin_phyla_visfood3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_visfood_PC3', 
  transform = "AST",
  fixed_effects = c('VisFood_RC3'),
  standardize = FALSE)

##Run at genus level (PC3)
Maaslin_genus_visfood3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_visfood_PC3', 
  transform = "AST",
  fixed_effects = c('VisFood_RC3'),
  standardize = FALSE)

####Behavior - Olfacoty
##Run at phyla level (PC1)
Maaslin_phyla_olf1 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_olfactory_PC1', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC1'),
  standardize = FALSE)

##Run at genus level (PC1)
Maaslin_genus_olf1 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_olfactory_PC1', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC1'),
  standardize = FALSE)

##Run at phyla level (PC2)
Maaslin_phyla_olf2 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_olfactory_PC2', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC2'),
  standardize = FALSE)

##Run at genus level (PC2)
Maaslin_genus_olf2 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_olfactory_PC2', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC2'),
  standardize = FALSE)

##Run at phyla level (PC3)
Maaslin_phyla_olf3 <- Maaslin2(
  tmp_data_phyla, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_phyla_olfactory_PC3', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC3'),
  standardize = FALSE)

##Run at genus level (PC3)
Maaslin_genus_olf3 <- Maaslin2(
  tmp_data_gen, tmp_metadata, 'C:/R/Emerson-Microbial-1/EM1_Final_Maaslin_genus_olfactory_PC3', 
  transform = "AST",
  fixed_effects = c('Olfactory_RC3'),
  standardize = FALSE)
