#####Libraries#####
library(tidyr)
library(dplyr)
library(readxl)
library(ggpubr)
library(ggplot2)
library(grid)
library(gridExtra)

####
#### this script utilizes amino acid specific metrics to model the relationship between the average amino acid frequencies
#### of the dataset and the optimal amino acid frequencies corresponding to minimal loss (differential retention). Also looked at 
#### was the relationship between the average amino acid frequencies and equilibrium frequencies (decent with modification)
#### 
#### calculated metrics:
####    1. average amino acid frequencies. Calculated with the average domain lengths of the Pfam's instances as weights
####
#### figures:
####    1. average amino acid frequencies predicted by amino acid frequencies for minimal loss (differential retention)
####    2. average amino acid frequencies predicted by equilibrium frequencies (decent with modification)
####
#### from the average amino acid frequencies ~ frequencies for minimal loss figure, it can be observed that differential retention is
#### a candidate explanation for why those averages are the averages
####

#### 
#### UPDATE SCRIPT WHEN GET 7 MISSING PFAMS - recalculate average aa frequencies with new length data
####

# suppress Rplots.pdf generation 
pdf(NULL)

# create results folder if doesnt exist
dir.create("results", showWarnings = FALSE)


##### setup
#####
# read in amino acid specific frequencies
amino_acid_specific_freq <-  read.csv("data/amino_acid_specific_metrics.csv", row.names = 1)

# read in pfam metrics
pfam_metrics <- read.csv("data/Pfam_metrics.csv")

#####



##### Pfams that we don't have length data for...but know their composition
##### these Pfams did not contribute to the weighted average amino acid frequency
#####
# get rows that dont have length data
which(is.na(pfam_metrics[, c("MeanDomainLength")]))

# print PfamUID for rows without length data
pfam_metrics$PfamUID[2191]
pfam_metrics$PfamUID[4111]
pfam_metrics$PfamUID[4127]
pfam_metrics$PfamUID[4440]
pfam_metrics$PfamUID[4780]
pfam_metrics$PfamUID[5008]
pfam_metrics$PfamUID[6639]

#####



##### compute weighted average amino acid frequency across entire data set for each amino acid using average domain lengths as weights
##### 20 values
#####
# omit rows with no length data
pfam_metrics_sub <- pfam_metrics[!is.na(pfam_metrics$MeanDomainLength), ]

# vector of one letter AA
amino_acids <- colnames(pfam_metrics_sub)[2:21]

# initialize vector to hold avg AA comp across data set
avg_AAfreq_acrossDataSet_vect <- vector()

# get average AA comp across entire data set weighted by avg domain lengths
# computes avg in same order that the minimum amino acid freq are in, so their columns line up
for (i in amino_acids) {
  weighted_avg <- (sum(pfam_metrics_sub[[i]] * pfam_metrics_sub$MeanDomainLength)) / sum(pfam_metrics_sub$MeanDomainLength)
  print(weighted_avg)
  avg_AAfreq_acrossDataSet_vect <- c(avg_AAfreq_acrossDataSet_vect, weighted_avg)
}

# store average amino acid freq in data with rows as amino acids
# columns?
avg_AAfreq_acrossDataSet <- data.frame(cbind(amino_acids, avg_AAfreq_acrossDataSet_vect))
colnames(avg_AAfreq_acrossDataSet) <- c("amino_acids", "weighted_avgAA_freq")

#####



##### build models and plots
#####
# generate linear model of avgAAfreq ~ minLossAAfreq
# this model is looking at how differential retention is related to the weighted average amino acid frequencies
minAAfreq_avgAAfreq_model <- lm(AverageAminoAcidFrequency ~ MinLossAminoAcidFrequency, amino_acid_specific_freq)
minAAfreq_avgAAfreq_coeff <- summary(minAAfreq_avgAAfreq_model)
minAAfreq_avgAAfreq_pvalue <- minAAfreq_avgAAfreq_coeff$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
plot_pvalue1 <- signif(minAAfreq_avgAAfreq_pvalue, digits = 1)
plot_rsquared1 <- signif(minAAfreq_avgAAfreq_coeff$adj.r.squared, digits = 2)

# generate linear model of avgAAfreq ~ equilAAfreq
# this model is looking at how descent with modification is related to the weighted average amino acid frequencies
equilAAfreq_avgAAfreq_model <- lm(AverageAminoAcidFrequency ~ EquilibriumAminoAcidFrequency, amino_acid_specific_freq)
equilAAfreq_avgAAfreq_coeff <- summary(equilAAfreq_avgAAfreq_model)
equilAAfreq_avgAAfreq_pvalue <- equilAAfreq_avgAAfreq_coeff$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
plot_pvalue2 <- signif(equilAAfreq_avgAAfreq_pvalue, digits = 1)
plot_rsquared2 <- signif(equilAAfreq_avgAAfreq_coeff$adj.r.squared, digits = 2)

# weighted average amino acid frequency predicted by amino acid frequenciees correesponding to minimal loss. Looking at 
# how much differential retention explains the average frequencies. 
plot_minAAfreq_avgAAfreq <- ggplot() +
  geom_smooth(data = amino_acid_specific_freq, method = "lm", aes(x = MinLossAminoAcidFrequency, y = AverageAminoAcidFrequency), color = 'midnightblue') +
  scale_x_continuous(breaks = c(0, 0.025, 0.050, 0.075, 0.10), limits = c(0, 0.10)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.050, 0.075, 0.10)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.25) +
  geom_text(data = amino_acid_specific_freq, aes(x = MinLossAminoAcidFrequency, y = AverageAminoAcidFrequency, label = amino_acids), size = 30/.pt) +
  annotate(geom = "text", label = paste("R^2 == ", plot_rsquared1), x = 0.0175, y = 0.10, size = 30/.pt, parse = TRUE) +
  annotate(geom = "text", label = bquote("P = " ~ 10^-09), x = 0.018, y = 0.09, size = 30/.pt) +
  theme_bw() +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(aspect.ratio = 1) +
  labs(x = "Amino acid frequency for minimum loss rate") +
  labs(y = "Average amino acid frequency") +
  theme(axis.title = element_text(size = 38, face = "bold")) +
  theme(axis.text = element_text(size = 30))

# size = 15 x 9 for annotations to be in correct positions
pdf("results/minAAfreq_avgAAfreq_plot.pdf", width = 15, height = 9)
plot_minAAfreq_avgAAfreq
dev.off()

# weighted average amino acid frequency predicted by amino acid frequencies corresponding to descent with modification
plot_equilAAfreq_avgAAfreq <- ggplot() +
  geom_smooth(data = amino_acid_specific_freq, method = "lm", aes(x = EquilibriumAminoAcidFrequency , y = AverageAminoAcidFrequency), color = 'midnightblue') +
  scale_x_continuous(breaks = c(0, 0.025, 0.050, 0.075, 0.10), limits = c(0, 0.10)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.050, 0.075, 0.10), limits = c(0, 0.11)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 0.25) +
  geom_text(data = amino_acid_specific_freq, aes(x = EquilibriumAminoAcidFrequency , y = AverageAminoAcidFrequency, label = amino_acids), size = 30/.pt) +
  annotate(geom = "text", label = paste("R^2 == ", plot_rsquared2), x = 0.0175, y = 0.10, size = 30/.pt, parse = TRUE) +
  annotate(geom = "text", label = bquote("P = " ~ 10^-12), x = 0.018, y = 0.09, size = 30/.pt) +
  theme_bw() +
  theme_linedraw() +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(aspect.ratio = 1) +
  labs(x = "Equilibrium amino acid frequency") +
  labs(y = "Average amino acid frequency") +
  theme(axis.title = element_text(size = 38, face = "bold")) +
  theme(axis.text = element_text(size = 30))

# size = 15 x 9 for annotations to be in correct positions
# plot_equilAAfreq_avgAAfreq


#####
