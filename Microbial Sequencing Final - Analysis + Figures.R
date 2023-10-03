###Microbial Sequencing Final - Analysis and Figures
library(devtools)
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
library(lme4)
library(car)

##### Alpha Diversity
## ASVs

file.choose()
df <- read.csv("C:\\Users\\kjeme\\OneDrive\\Desktop\\Woodley Lab\\Aim 2 - Emerson Microbial 1\\Emersion Microbial Experiment 1 (2021)\\EmersonMicrobialExp1_Batch1.csv")
df$Microbial_Trtmt = factor(df$Microbial_Trtmt)
df$Replicate = factor(df$Replicate)

asvs_glmm <- glm(ASV~Microbial_Trtmt, data = df, family = "gaussian")
Anova(asvs_glmm)
# p = 0.014

Pond_Water <- c("Natural", "Autoclaved")

##Boxplot help: https://www.youtube.com/watch?v=mbvU8fF4eGM
ggplot(df, aes(x= Microbial_Trtmt, y = ASV, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  annotate(x = 1, xend = 2, y = 282, yend = 282, geom = "segment") +
  annotate(x=1.5, y= 290, label = "p < 0.05", geom = "text", size = 5.2) +
  theme_classic() +
  stat_boxplot(geom = "errorbar", width = .35) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Pond Water", y = "No. Observed ASVs") +
  scale_x_discrete(labels = Pond_Water) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)

df.outlier <- df[-c(7),] 
#Outlier with 830 ASVs detected. leaving in for stats
#because no evidence of contamination, but can remove from figure

ggplot(df.outlier, aes(x= Microbial_Trtmt, y = ASV, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  annotate(x = 1, xend = 2, y = 282, yend = 282, geom = "segment") +
  annotate(x=1.5, y= 290, label = "p = 0.014", geom = "text", size = 5.2) +
  coord_cartesian(ylim = c(0, 300)) +
  theme_classic() +
  stat_boxplot(geom = "errorbar", width = .35) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Pond Water", y = "No. Observed ASVs") +
  scale_x_discrete(labels = Pond_Water) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)


## Shannon Diversity Index

shannon_glmm <- glm(Shannon~Microbial_Trtmt, data = df, family = "gaussian")
Anova(shannon_glmm)
# p = 0.0014

ggplot(df, aes(x= Microbial_Trtmt, y = Shannon, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  coord_cartesian(ylim = c(0, 4)) +
  annotate(x = 1, xend = 2, y = 3.9, yend = 3.9, geom = "segment") +
  annotate(x=1.5, y= 4, label = "p = 0.0014", geom = "text", size = 5.2) +
  theme_classic() +
  stat_boxplot(geom = "errorbar", width = .35) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Pond Water", y = "Shannon Diversity Index") +
  scale_x_discrete(labels = Pond_Water) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)

ggplot(df.outlier, aes(x= Microbial_Trtmt, y = Shannon, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  coord_cartesian(ylim = c(0, 4)) +
  annotate(x = 1, xend = 2, y = 3.9, yend = 3.9, geom = "segment") +
  annotate(x=1.5, y= 4, label = "p = 0.0014", geom = "text", size = 5.2) +
  theme_classic() +
  stat_boxplot(geom = "errorbar", width = .35) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Pond Water", y = "Shannon Diversity Index") +
  scale_x_discrete(labels = Pond_Water) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)


## Faiths Phylogenetic Diversity

faith_glmm <- glm(Faiths~Microbial_Trtmt, data = df, family = "gaussian")
Anova(faith_glmm)
# p = 0.0014

ggplot(df, aes(x= Microbial_Trtmt, y = Faiths, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  annotate(x = 1, xend = 2, y = 59, yend = 59, geom = "segment") +
  annotate(x=1.5, y= 60.5, label = "p < 0.001", geom = "text", size = 5.2) +
  theme_classic() +
  stat_boxplot(geom = "errorbar", width = .35) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Pond Water", y = "Faith's Phylogenetic Diversity") +
  scale_x_discrete(labels = Pond_Water) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)

ggplot(df.outlier, aes(x= Microbial_Trtmt, y = Faiths, fill = Microbial_Trtmt)) +
  geom_boxplot(width = 0.5, outlier.colour = "transparent", alpha =0.35) +
  geom_jitter(width = .2, size = 2.5, shape = 21, color = "black") +
  annotate(x = 1, xend = 2, y = 62, yend = 62, geom = "segment") +
  annotate(x=1.5, y= 63.5, label = "p = 0.0014", geom = "text", size = 5.2) +
  coord_cartesian(ylim = c(0, 65)) +
  theme_classic() +
  stat_boxplot(geom = "errorbar", width = .35) +
  scale_fill_manual(values = c("seagreen3", "lightcyan"),
                    name = "Pond Water",
                    labels = c("Natural", "Autoclaved")) +
  labs(x = "Pond Water", y = "Faith's Phylogenetic Diversity") +
  scale_x_discrete(labels = Pond_Water) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1)

######### Now, we want to evaluate if our metrics of alpha diversity are significant 
## predictors of brain shape. 

## Alpha Diversity ~ Relative Brain Mass

alpha_brainmass_glmm <- glm(Relative_BrainMass~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_brainmass_glmm)

asv_mass_glmm <- glm(Relative_BrainMass~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_mass_glmm)

shannon_mass_glmm <- glm(Relative_BrainMass~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_mass_glmm)

Faiths_mass_glmm <- glm(Relative_BrainMass~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_mass_glmm)

####Alpha Diversity ~ Brain Shape (PC1)

alpha_PC1_glmm <- glm(RC1~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_PC1_glmm)

asv_PC1_glmm <- glm(RC1~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_PC1_glmm)

shannon_PC1_glmm <- glm(RC1~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_PC1_glmm)

Faiths_PC1_glmm <- glm(RC1~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_PC1_glmm)

####Alpha Diversity ~ Brain Shape (PC2)

alpha_PC2_glmm <- glm(RC2~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_PC2_glmm)

asv_PC2_glmm <- glm(RC2~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_PC2_glmm)

shannon_PC2_glmm <- glm(RC2~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_PC2_glmm)

Faiths_PC2_glmm <- glm(RC2~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_PC2_glmm)

####Alpha Diversity ~ Brain Shape (PC3)

alpha_PC3_glmm <- glm(RC3~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_PC3_glmm)

asv_PC3_glmm <- glm(RC3~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_PC3_glmm)

shannon_PC3_glmm <- glm(RC3~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_PC3_glmm)

Faiths_PC3_glmm <- glm(RC3~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_PC3_glmm)

###Brain Shape ~ Alpha Diversity figures
#RC3 ~ Shannon

ggplot(df, aes(x= Shannon, y = RC3, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  theme_classic() +
  labs(x = "Shannon Diversity Index", y = "Medulla Width (PC-3)") +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural", "Autoclaved")) +
  scale_shape_manual(values = c(16,2)) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))

###### Now, we can do the same thing but with behavioral parameters

####Visual Empty PC1 ~ Alpha Diversity
alpha_visemptyPC1_glmm <- glm(VisEmpty_RC1~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_visemptyPC1_glmm)

asv_visemptyPC1_glmm <- glm(VisEmpty_RC1~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_visemptyPC1_glmm)

shannon_visemptyPC1_glmm <- glm(VisEmpty_RC1~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_visemptyPC1_glmm)

Faiths_visemptyPC1_glmm <- glm(VisEmpty_RC1~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_visemptyPC1_glmm)

####Visual Empty PC2 ~ Alpha Diversity

alpha_visemptyPC2_glmm <- glm(VisEmpty_RC2~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_visemptyPC2_glmm)

asv_visemptyPC2_glmm <- glm(VisEmpty_RC2~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_visemptyPC2_glmm)

shannon_visemptyPC2_glmm <- glm(VisEmpty_RC2~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_visemptyPC2_glmm)

Faiths_visemptyPC2_glmm <- glm(VisEmpty_RC2~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_visemptyPC2_glmm)

####Visual Empty PC3 ~ Alpha Diversity

alpha_visemptyPC3_glmm <- glm(VisEmpty_RC3~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_visemptyPC3_glmm)

asv_visemptyPC3_glmm <- glm(VisEmpty_RC3~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_visemptyPC3_glmm)

shannon_visemptyPC3_glmm <- glm(VisEmpty_RC3~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_visemptyPC3_glmm)

Faiths_visemptyPC3_glmm <- glm(VisEmpty_RC3~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_visemptyPC3_glmm)

####Visual Food PC1 ~ Alpha Diversity

alpha_visfoodPC1_glmm <- glm(VisFood_RC1~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_visfoodPC1_glmm)

asv_visfoodPC1_glmm <- glm(VisFood_RC1~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_visfoodPC1_glmm)

shannon_visfoodPC1_glmm <- glm(VisFood_RC1~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_visfoodPC1_glmm)

Faiths_visfoodPC1_glmm <- glm(VisFood_RC1~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_visfoodPC1_glmm)

####Visual Food PC2 ~ Alpha Diversity

alpha_visfoodPC2_glmm <- glm(VisFood_RC2~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_visfoodPC2_glmm)

asv_visfoodPC2_glmm <- glm(VisFood_RC2~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_visfoodPC2_glmm)

shannon_visfoodPC2_glmm <- glm(VisFood_RC2~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_visfoodPC2_glmm)

Faiths_visfoodPC2_glmm <- glm(VisFood_RC2~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_visfoodPC2_glmm)

####Visual Food PC3 ~ Alpha Diversity

alpha_visfoodPC3_glmm <- glm(VisFood_RC3~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_visfoodPC3_glmm)

asv_visfoodPC3_glmm <- glm(VisFood_RC3~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_visfoodPC3_glmm)

shannon_visfoodPC3_glmm <- glm(VisFood_RC3~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_visfoodPC3_glmm)

Faiths_visfoodPC3_glmm <- glm(VisFood_RC3~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_visfoodPC3_glmm)

####Olfactory PC1 ~ Alpha Diversity

alpha_olfactoryPC1_glmm <- glm(Olfactory_RC1~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_olfactoryPC1_glmm)

asv_olfactoryPC1_glmm <- glm(Olfactory_RC1~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_olfactoryPC1_glmm)

shannon_olfactoryPC1_glmm <- glm(Olfactory_RC1~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_olfactoryPC1_glmm)

Faiths_olfactoryPC1_glmm <- glm(Olfactory_RC1~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_olfactoryPC1_glmm)

####Olfactory PC2 ~ Alpha Diversity

alpha_olfactoryPC2_glmm <- glm(Olfactory_RC2~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_olfactoryPC2_glmm)

asv_olfactoryPC2_glmm <- glm(Olfactory_RC2~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_olfactoryPC2_glmm)

shannon_olfactoryPC2_glmm <- glm(Olfactory_RC2~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_olfactoryPC2_glmm)

Faiths_olfactoryPC2_glmm <- glm(Olfactory_RC2~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_olfactoryPC2_glmm)

####Olfactory PC3 ~ Alpha Diversity

alpha_olfactoryPC3_glmm <- glm(Olfactory_RC3~ASV*Shannon*Faiths, data = df, family = "gaussian")
Anova(alpha_olfactoryPC3_glmm)

asv_olfactoryPC3_glmm <- glm(Olfactory_RC3~Microbial_Trtmt*ASV, data = df, family = "gaussian")
Anova(asv_olfactoryPC3_glmm)

shannon_olfactoryPC3_glmm <- glm(Olfactory_RC3~Microbial_Trtmt*Shannon, data = df, family = "gaussian")
Anova(shannon_olfactoryPC3_glmm)

Faiths_olfactoryPC3_glmm <- glm(Olfactory_RC3~Microbial_Trtmt*Faiths, data = df, family = "gaussian")
Anova(Faiths_olfactoryPC3_glmm)


#Shannon and Brain Mass correlation plot

ggplot(df, aes(x= Shannon, y = Relative_BrainMass, color = Microbial_Trtmt)) +
  geom_point (size = 4)+
  theme_classic() +
  labs(x = "Shannon", y = "Relative Brain Mass") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = c(1.2, .9)) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))



#Shannon and Brain Shape correlation plot
#Fitting a linear regression model to the data

model <- lm(RC3~Shannon, data=df)
summary(model)

ggplot(df, aes(x= Shannon, y = RC3, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  theme_classic() +
  geom_smooth(method=lm)+
  labs(x = "Shannon Diversity Index", y = "Medulla Width (PC-3)") +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural", "Autoclaved")) +
  scale_shape_manual(values = c(16,17)) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))

#Shannon and visual food locomotion correlation plot

ggplot(df, aes(x= Shannon, y = VisFood_RC2, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  theme_classic() +
  geom_smooth(method=lm)+
  labs(x = "Shannon Diversity Index", y = "Locomotory Activity (Vis Food - PC2)") +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural", "Autoclaved")) +
  scale_shape_manual(values = c(16,2)) +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))

#ASV and visual food locomotion correlation plot
#Fitting a linear regression model to the data

model <- lm(VisFood_RC1~ASV, data=df)
summary(model)

ggplot(df, aes(x= ASV, y = VisFood_RC1, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  geom_smooth(method=lm)+
  coord_cartesian(xlim = c(0, 350)) +
  theme_classic() +
  labs(x = "No. Observed ASVs", y = "Change in Locomotory Activity (PC-1): Vial with Food") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 17.5)) +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural" , "Autoclaved")) +
  scale_shape_manual(values = c(16,17)) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))


#Faith and visual food locomotion correlation plot

ggplot(df, aes(x= Faiths, y = VisFood_RC1, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  theme_classic() +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "Faith's Phylogenetic Diversity", y = "Change in Locomotory Activity (PC1)") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural" , "Autoclaved")) +
  scale_shape_manual(values = c(16,2)) +
  theme(legend.position = c("none")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14)) 

#ASV and olfactory locomotion correlation plot 

ggplot(df, aes(x= ASV, y = Olfactory_RC2, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  geom_smooth(method=lm)+
  coord_cartesian(xlim = c(0, 350)) +
  theme_classic() +
  labs(x = "No. Observed ASVs", y = "Change in Locomotory Activity (Olfactory - PC2)") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural" , "Autoclaved")) +
  scale_shape_manual(values = c(16,2)) +
  theme(legend.position = "none") +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14))

# Faith and olfactory locomotion correlation plot

ggplot(df, aes(x= Faiths, y = Olfactory_RC2, color = Microbial_Trtmt)) +
  geom_point (aes(shape = Microbial_Trtmt, size = 3))+
  theme_classic() +
  coord_cartesian(xlim = c(0, 60)) +
  labs(x = "Faith's Phylogenetic Diversity", y = "Change in Locomotory Activity (Olfactory - PC2)") +
  theme(aspect.ratio = 1) +
  theme(axis.text = element_text(face = "bold", size = 16))  +
  theme(axis.title = element_text(face = "bold", size = 18)) +
  scale_color_manual(values = c("seagreen3", "skyblue3"),
                     labels = c("Natural" , "Autoclaved")) +
  scale_shape_manual(values = c(16,2)) +
  theme(legend.position = c("none")) +
  theme(legend.title = element_blank()) +
  theme(legend.text = element_text(size = 14)) 
