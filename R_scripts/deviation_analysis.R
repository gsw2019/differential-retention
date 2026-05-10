#####Libraries#####
library(dplyr)
library(ggplot2)
library(geoR)
library(ggpubr)
library(grid)
library(gridExtra)
library(scales)
library(ggallin)
library(ggrepel)
library(tidyverse)


####
#### this script utilizes Pfam amino acid frequencies, Pfam domain lengths, and amino acid specific metrics from the dataset to
#### calculate the absolute deviation of a Pfam from the average amino acid frequencies. The deviation is then used in a model to
#### predict loss rate. 
####
#### the issue of smaller Pfams deviating more just by chance, and therefore appearing as lost more often, is also investigated
####
#### sigma stuff...
####
#### finally larger domain driving our loess curve away from the sigma expectation is explored
####
#### calculated metrics:
####    1. absolute deviation of each Pfam from the average amino acid frequencies
####    2. log10 absolute deviation of each Pfam from the average amino acid frequencies
####    3. 
####
#### figures:
####    1. lm and loess of absolute deviation from average amino acid frequency VS loss rate
####    2. lm and loess of log10(absolute deviation from average amino acid frequency) VS loss rate
####    3. loess of log10(median domain length) VS loss rate
####    4. loess sigma plot, (deviation*sqrt(median domain length)) ~ log10(medianDomainLength)
####    5. 20 loess figures, a plot of log10(median domain length) VS frequency for each amino acid
####    6. 20 loess lines on one plot of log10(median domain length) VS cumulative frequencies
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

# read in Pfam metrics
pfam_metrics <- read.csv("data/Pfam_metrics.csv")

amino_acids <- colnames(pfam_metrics[2:21])

# normalizing loss rate and getting lambda values
Loss.box = boxcoxfit(pfam_metrics$MLLossRate, lambda2=TRUE)
Loss.box
lambda = Loss.box$lambda[1]       # 0.0536
lambda2 = Loss.box$lambda[2]      # 6.1e-6
if(lambda==0){loss_norm = log(pfam_metrics$MLLossRate + lambda2)}
if(lambda!=0){loss_norm = ((pfam_metrics$MLLossRate + lambda2) ^ lambda - 1) / lambda}

# required functions for plotting
bc.transform <- function(x,L){(x^L-1)/L}
bc.backtransform <- function(y,L){(L*y+1)^(1/L)}

#####



##### compute absolute deviation and log10(absolute deviation) from average amino acid frequencies
#####
# initialize vectors to store each Pfams deviation from average amino acid frequencies
deviation_from_avg_AAfreq <- vector()

# compute absolute deviation for each Pfam
for (i in 1:nrow(pfam_metrics)) {
  # variable to store Pfam's sum of squared differences
  pfam_sum_sd <- 0
  
  for (j in amino_acids) {
    # calculate each squared difference (observed - average)^2 and add to total
    temp_sd <- (pfam_metrics[i, j] - amino_acid_specific_freq[j, "AverageAminoAcidFrequency"])^2
    pfam_sum_sd <- pfam_sum_sd + temp_sd
  }
  
  # take the square root of the sum of squared differences and store absolute deviation from average amino acid frequency
  pfam_abs_dev <- sqrt(pfam_sum_sd)
  deviation_from_avg_AAfreq <- c(deviation_from_avg_AAfreq, pfam_abs_dev)
}

log10_deviation_from_avg_AAfreq <- log10(deviation_from_avg_AAfreq)

#####



##### build plots of each deviation metric vs. loss rate 
##### -
##### 1. absolute deviation from meanAAfreq VS loss rate
#####
# linear model
LossNorm_AbsDeviation_lm <- lm(LossNormalized ~ DeviationFromAvgAminoAcidFrequency, pfam_metrics)
LossNorm_AbsDeviation_lm_coeff <- summary(LossNorm_AbsDeviation_lm)


# lm plot
plot_LossNorm_AbsDeviation_lm <- ggplot()+
  # geom_violin(data = master_df, aes(group=cut_interval((DeviationFromAvgAAcomp), n = 20), x = DeviationFromAvgAAcomp, y = Loss_norm),
  #             fill= "azure2", colour = "darkslategray4", position = position_nudge(-0.04)) +
  geom_jitter(data = pfam_metrics, aes(x = DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color="black", size=0.4, alpha=0.3) +
  geom_smooth(data = pfam_metrics, method = "lm", aes(x = DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color = 'midnightblue') +
  scale_y_continuous(breaks = c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), labels = signif(bc.backtransform(c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), lambda), digits = 1)) +
  scale_x_continuous(limits = c(0.025, 0.3), breaks = seq(0.025, 0.3, 0.055)) +
  # x-axis omits one value less than 0.025 deviation from meanAAfreq and cuts off some sparse larger values
  labs(x = "Deviation From Mean Amino Acid Frequency") +
  labs(y = "Losses/million years") +
  geom_text(label = paste('R^2 =', signif(LossNorm_AbsDeviation_lm_coeff$adj.r.squared, digits = 3)), aes(x = 0.06, y = 0), color = "blue", size = 30/.pt) +
  geom_text(label = paste('p =', signif(LossNorm_AbsDeviation_lm_coeff$coefficients[2,4], digits = 2)), aes(x = 0.06, y = -0.50), color = "blue", size = 30/.pt) +
  
  theme_bw() +
  theme_linedraw() +
  theme(text = element_text(size=30)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

# plot_LossNorm_AbsDeviation_lm


# loess plot
plot_meanAAfreqAbsDeviation_loess <- ggplot()+
  # geom_violin(data = master_df, aes(group=cut_interval((DeviationFromAvgAAcomp), n = 20), x = DeviationFromAvgAAcomp, y = Loss_norm),
  #             fill= "azure2", colour = "darkslategray4", position = position_nudge(-0.04)) +
  geom_jitter(data = pfam_metrics, aes(x = DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color="black", size=0.4, alpha=0.3) +
  geom_smooth(data = pfam_metrics, method = "loess", aes(x = DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color = 'midnightblue') +
  scale_y_continuous(breaks = c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), labels = signif(bc.backtransform(c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), lambda), digits = 1)) +
  scale_x_continuous(limits = c(0.025, 0.3), breaks = seq(0.025, 0.3, 0.055)) +
  # x-axis omits one value less than 0.025 deviation from meanAAfreq and cuts off some sparse larger values
  labs(x = "Deviation from mean amino acid frequency") +
  labs(y = "Losses/million years") +
  theme_bw() +
  theme_linedraw() +
  theme(text = element_text(size=30)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

# plot_meanAAfreqAbsDeviation_loess

#####
##### 2. log10(absolute deviation from meanAAfreq) VS loss rate with 10% and 90% quantiles
#####
# linear model
LossNorm_log10AbsDeviation_lm <- lm(LossNormalized ~ Log10DeviationFromAvgAminoAcidFrequency, pfam_metrics)
LossNorm_log10AbsDeviation_lm_coeff <- summary(LossNorm_log10AbsDeviation_lm)

# getting full pvalue from model
LossNorm_log10AbsDeviation_pvalue <- LossNorm_log10AbsDeviation_lm_coeff$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
plot_pvalue <- signif(LossNorm_log10AbsDeviation_pvalue, digits = 1)

# getting rsquared from model
plot_rsquared <- signif(LossNorm_log10AbsDeviation_lm_coeff$adj.r.squared, digits = 2)

# x-axis for graph
x_axis <- seq(-1.60, -0.25, 0.45)

# back transforming x-axis values
log10_x_axis <- signif(10^x_axis, digits = 2)

# get quantiles, predict loss rate from quantiles, calculate fold difference, and ensure equal data on both sides
# Calculate the 10% quantile
LossNorm_log10AbsDeviation_quantile10 <- quantile(pfam_metrics$Log10DeviationFromAvgAminoAcidFrequency, 0.1)

# Calculate the 90% quantile
LossNorm_log10AbsDeviation_quantile90 <- quantile(pfam_metrics$Log10DeviationFromAvgAminoAcidFrequency, 0.9)

# check the amount of data is the same on the left of the 10% quantile and right of 90% quantile
length(which(pfam_metrics$Log10DeviationFromAvgAminoAcidFrequency <= LossNorm_log10AbsDeviation_quantile10))
length(which(pfam_metrics$Log10DeviationFromAvgAminoAcidFrequency >= LossNorm_log10AbsDeviation_quantile90))

# make data frame with 10% quantile of deviations and 90% quantile of deviations
LossNorm_log10AbsDeviation_quantiles <- data.frame(Log10DeviationFromAvgAminoAcidFrequency = c(LossNorm_log10AbsDeviation_quantile10, LossNorm_log10AbsDeviation_quantile90))

# use the linear model to predict the y-axis values (loss rate) from the deviation quantile values
loss_quantiles1 <- predict(LossNorm_log10AbsDeviation_lm, LossNorm_log10AbsDeviation_quantiles)

# transform the loss rates
loss_quantiles_transformed1 <- bc.backtransform(loss_quantiles1, lambda)

# add the loss rates and transformed loss rates to the data frame
LossNorm_log10AbsDeviation_quantiles <- cbind(LossNorm_log10AbsDeviation_quantiles, loss_quantiles1, loss_quantiles_transformed1)

# compute fold difference between 90% quantile and 10% quantile loss rate, 4.86
(LossNorm_log10AbsDeviation_quantiles$loss_quantiles_transformed1[2])/(LossNorm_log10AbsDeviation_quantiles$loss_quantiles_transformed1[1])

# lm plot
plot_LossNorm_log10AbsDeviation_lm <- ggplot() +
  geom_jitter(data = pfam_metrics, aes(x = Log10DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color="black", size=0.4, alpha=0.3) +
  geom_smooth(data = pfam_metrics, method = "lm", aes(x = Log10DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color = 'midnightblue') +
  scale_y_continuous(breaks = c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), labels = signif(bc.backtransform(c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), lambda), digits = 1)) +
  scale_x_continuous(limits = c(-1.6, -0.25), breaks = x_axis, labels = log10_x_axis) +
  # scale_x_continuous(limits = c(-1.75, -0.25), breaks = x_axis, labels = log10_x_axis) +
  geom_point(aes(x = LossNorm_log10AbsDeviation_quantiles[1,1], y = LossNorm_log10AbsDeviation_quantiles[1, 2]), color = "red", size = 6) +  # 10% quantile
  geom_point(aes(x = LossNorm_log10AbsDeviation_quantiles[2,1], y = LossNorm_log10AbsDeviation_quantiles[2, 2]), color = "red", size = 6) +  # 90% quantile
  labs(x = "Deviation from expected amino acid frequency") +
  labs(y = "Losses/million years") +
  annotate(geom = "text", label = paste('R^2 == ', plot_rsquared), x = -1.4, y = -0.75, size = 40/.pt, parse = TRUE) + 
  annotate(geom = "text", label = bquote("P <" ~ 10^-16), x = -1.4, y = -1.75, size = 40/.pt) +
  theme_bw() + 
  theme_linedraw() +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(axis.text = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

# size = 14 x 8 for annotations in correct positions
pdf("results/loss_norm_log10_abs_deviation_plot.pdf", width = 14, height = 8)
plot_LossNorm_log10AbsDeviation_lm
dev.off()

# loess plot
plot_LossNorm_log10AbsDeviation_loess <- ggplot()+
  # geom_violin(data = master_df, aes(group=cut_interval((DeviationFromAvgAAcomp), n = 20), x = DeviationFromAvgAAcomp, y = Loss_norm),
  #             fill= "azure2", colour = "darkslategray4", position = position_nudge(-0.04)) +
  geom_jitter(data = pfam_metrics, aes(x = Log10DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color="black", size=0.4, alpha=0.3) +
  geom_smooth(data = pfam_metrics, method = "loess", aes(x = Log10DeviationFromAvgAminoAcidFrequency, y = LossNormalized), color = 'midnightblue') +
  scale_y_continuous(breaks = c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), labels = signif(bc.backtransform(c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), lambda), digits = 1)) +
  scale_x_continuous(limits = c(-1.6, -0.25), breaks = x_axis, labels = log10_x_axis) +
  labs(x = "Deviation from expected amino acid frequency") +
  labs(y = "Losses/million years") +
  theme_bw() + 
  theme_linedraw() +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(axis.text = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

# plot_LossNorm_log10AbsDeviation_loess

#####



##### looking at issue of shorter Pfams deviating more by chance and therefore appearing as being lost more often
##### -
##### a plot of log10(median domain length) predicting loss rate with 10% and 90% quantiles
##### -
##### also computed are model summaries showing how log10(median domain length) and log10(absolute deviation) 
##### predict loss rate and that their independent R^2s almost add up to their combined prediction R^2. We are
##### using this as evidence that their effects on loss rate are mostly independent. 
#####
# linear and loess model of log10(median domain length) vs loss rate
LossNorm_log10MedianDomainLength_lm <- lm(LossNormalized ~ Log10MedianDomainLength, pfam_metrics)
LossNorm_log10MedianDomainLength_loess <- loess(LossNormalized ~ Log10MedianDomainLength, pfam_metrics)
LossNorm_log10MedianDomainLength_lm_coeff <- summary(LossNorm_log10MedianDomainLength_lm)

# getting full pvalue from linear model
LossNorm_log10MedianDomainLength_pvalue <- LossNorm_log10MedianDomainLength_lm_coeff$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
plot_pvalue2 <- signif(LossNorm_log10MedianDomainLength_pvalue, digits = 1)

# getting rsquared from model
plot_rsquared2 <- signif(LossNorm_log10MedianDomainLength_lm_coeff$adj.r.squared, digits = 2)

# create data frame to pass to predict and find y-axis (loss rate) coordinates for quantiles
# calculate the 10% quantile 
LossNorm_log10MedianDomainLength_quantile10 <- quantile(pfam_metrics$Log10MedianDomainLength, 0.10, na.rm = T)

# calculate the 90% quantile
LossNorm_log10MedianDomainLength_quantile90 <- quantile(pfam_metrics$Log10MedianDomainLength, 0.9, na.rm = T)

# check the amount of data is the same or very close to the same on the left of the 10% quantile and right of 90% quantile
length(which(pfam_metrics$Log10MedianDomainLength <= LossNorm_log10MedianDomainLength_quantile10))
length(which(pfam_metrics$Log10MedianDomainLength >= LossNorm_log10MedianDomainLength_quantile90))

# make data frame with 10% quantile of x-axis and 90% quantile of x-axis
LossNorm_log10MedianDomainLength_quantiles <- data.frame(Log10MedianDomainLength = c(LossNorm_log10MedianDomainLength_quantile10, LossNorm_log10MedianDomainLength_quantile90))

# use the loess model to predict the y-axis values from the x-axis values
loss_quantiles2 <- predict(LossNorm_log10MedianDomainLength_loess, LossNorm_log10MedianDomainLength_quantiles)

# transform the loss rates
loss_quantiles_transformed2 <- bc.backtransform(loss_quantiles2, lambda)

# add the loss rates and transformed loss rates to the data frame
LossNorm_log10MedianDomainLength_quantiles <- cbind(LossNorm_log10MedianDomainLength_quantiles, loss_quantiles2, loss_quantiles_transformed2)

# calculate fold difference b/w 90% quantile and 10% quantile
# 1.70x difference in loss rate between Pfams with a 10% quantile length and Pfams with a 90% quantile length
(LossNorm_log10MedianDomainLength_quantiles$loss_quantiles_transformed2[1])/(LossNorm_log10MedianDomainLength_quantiles$loss_quantiles_transformed2[2])

# plot
plot_LossNorm_log10MedianDomainLength_loess <- ggplot() +
  geom_jitter(data = pfam_metrics, aes(x = Log10MedianDomainLength, y = LossNormalized), color="black", size=0.4, alpha=0.3) +
  geom_smooth(data = pfam_metrics, method = "loess", aes(x = Log10MedianDomainLength, y = LossNormalized), color = 'midnightblue') +
  scale_y_continuous(breaks = c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), labels = signif(bc.backtransform(c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), lambda), digits = 1)) +
  #scale_x_continuous(limits = c(0, 750), breaks = seq(0, 750, 250)) +
  scale_x_continuous(limits = c(0.8, 3.2), breaks = seq(0.8, 3.2, 0.4), labels = signif(10^seq(0.8, 3.2, 0.4), digits = 2)) +
  # scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  geom_point(aes(x = LossNorm_log10MedianDomainLength_quantiles[1,1], y = LossNorm_log10MedianDomainLength_quantiles[1, 2]), color = "red", size = 6) +  # 10% quantile
  geom_point(aes(x = LossNorm_log10MedianDomainLength_quantiles[2,1], y = LossNorm_log10MedianDomainLength_quantiles[2, 2]), color = "red", size = 6) +  # 90% quantile
  labs(y = "Losses/million years") +
  labs(x = "Median Domain length") +
  annotate(geom = "text", label = paste('R^2 == ', plot_rsquared2), x = 2.8, y = -1, size = 40/.pt, parse = TRUE) + 
  annotate(geom = "text", label = paste('P = ',  plot_pvalue2), x = 2.8, y = -2, size = 40/.pt) +
  theme_bw() + 
  theme_linedraw() +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(axis.text = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())
  
# 14 x 8 for annotations in correct position
pdf("results/loss_norm_log10_median_domain_length_plot.pdf", width = 14, height = 8)
plot_LossNorm_log10MedianDomainLength_loess
dev.off()

# lm of LossNormalized predicted by log10(median domain length) and log10(Deviation from mean AA freq)
# shorter domain length predicts high loss rate + higher aa deviation
# length and deviation affect loss mostly independently, based on R^2
# the two factors independently predicting loss rate have R^2s that almost add up to their combined prediction
summary(lm(pfam_metrics$LossNormalized ~ pfam_metrics$Log10MedianDomainLength))
summary(lm(pfam_metrics$LossNormalized ~ pfam_metrics$Log10MedianDomainLength + pfam_metrics$Log10DeviationFromAvgAminoAcidFrequency))
summary(lm(pfam_metrics$LossNormalized ~ pfam_metrics$Log10DeviationFromAvgAminoAcidFrequency))

#####



##### calculating sigma
##### sigma represents...
##### current value of sigma = 0.9797
#####
# create data frame with 20 rows and 20 columns (amino acids). Each row has 19 0's and one 1. The 1 is in a 
# different column for each row, representing a peptide of length 1 with respective amino acid. 
observed_len1 <- data.frame(matrix(ncol = 20, nrow = 20))
colnames(observed_len1) <- amino_acids

# add 1s
for (i in 1:nrow(observed_len1)) {
  observed_len1[i, i] <- 1
}

# add 0s to empty cells
na_cells <- which(is.na(observed_len1), arr.ind = TRUE)
observed_len1[na_cells] <- 0

# calculating expected sigma value
# returns a data frame with each column as an amino acid and the rows as squared differences 
seq_len1_dev <- sapply(1:20, function(aa) {(observed_len1[, aa] - amino_acid_specific_freq[aa, "AverageAminoAcidFrequency"])^2})
colnames(seq_len1_dev) <- amino_acids

# sum the squared differences of each amino acid, 20 values
seq_len1_dev <- sapply(1:20, function(aa) {(sum(seq_len1_dev[, aa]))})
names(seq_len1_dev) <- amino_acids

# multiply each sum of squared differences by the average frequency, sum them, then square root
sigma <- sqrt(sum(seq_len1_dev*amino_acid_specific_freq$AverageAminoAcidFrequency))

#####



##### utilizing sigma and building plots that help shed light on confounding factor (domain length) to our 
##### loss rate ~ deviation result
#####
# # deviation ~ medianDomainLength linear model
# length_deviation_lm <- lm(DeviationFromAvgAminoAcidFrequency ~ MedianDomainLength, pfam_metrics)
# length_deviation_lm_coeff <- summary(length_deviation_lm)
# 
# # getting slope of model
# length_deviation_slope <- signif(length_deviation_lm_coeff$coefficients[2, 1], digits = 2)
# 
# # getting full pvalue from model
# length_deviation_pvalue <- length_deviation_lm_coeff$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
# length_deviation_pvalue <- signif(length_deviation_pvalue, digits = 1)
# 
# # getting rsquared from model
# length_deviation_rsquared <- signif(length_deviation_lm_coeff$adj.r.squared, digits = 2)



# # log10(deviation) ~ log10(medianDomainLength) linear model
# length_deviation_lm2 <- lm(Log10DeviationFromAvgAminoAcidFrequency ~ Log10MedianDomainLength, pfam_metrics)
# length_deviation_lm_coeff2 <- summary(length_deviation_lm2)
# 
# # getting slope of model
# length_deviation_slope2 <- signif(length_deviation_lm_coeff2$coefficients[2, 1], digits = 2)
# 
# # getting full pvalue from model
# length_deviation_pvalue2 <- length_deviation_lm_coeff2$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
# length_deviation_pvalue2 <- signif(length_deviation_pvalue2, digits = 1)
# 
# # getting rsquared from model
# length_deviation_rsquared2 <- signif(length_deviation_lm_coeff2$adj.r.squared, digits = 2)



# # deviation ~ 1/(medianDomainLength)^1/3) linear model
# temp_domainLength <- 1/((pfam_metrics$MedianDomainLength)^(1/3))
# length_deviation_lm3 <- lm(pfam_metrics$DeviationFromAvgAminoAcidFrequency ~ temp_domainLength, pfam_metrics)
# length_deviation_lm_coeff3 <- summary(length_deviation_lm3)
# 
# # fix ticks
# temp_domainLength_ticks <- signif(1/(seq(0.085, 0.5, 0.1)^3), digits = 2)
# 
# # getting slope of model
# length_deviation_slope3 <- signif(length_deviation_lm_coeff3$coefficients[2, 1], digits = 2)
# 
# # getting full pvalue from model
# length_deviation_pvalue3 <- length_deviation_lm_coeff3$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
# length_deviation_pvalue3 <- signif(length_deviation_pvalue3, digits = 1)
# 
# # getting rsquared from model
# length_deviation_rsquared3 <- signif(length_deviation_lm_coeff3$adj.r.squared, digits = 2)



# # deviation ~ 1/(medianDomainLength)^1/2) linear model
# temp_domainLength2 <- 1/((pfam_metrics$MedianDomainLength)^(1/2))
# length_deviation_lm4 <- lm(DeviationFromAvgAminoAcidFrequency ~ temp_domainLength2, pfam_metrics)
# length_deviation_lm_coeff4 <- summary(length_deviation_lm4)
# 
# # fix ticks
# temp_domainLength2_ticks <- signif(1/(seq(0.025, 0.4, 0.085)^2), digits = 2)
# 
# # getting slope of model
# length_deviation_slope4 <- signif(length_deviation_lm_coeff4$coefficients[2, 1], digits = 2)
# 
# # getting full pvalue from model
# length_deviation_pvalue4 <- length_deviation_lm_coeff4$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
# length_deviation_pvalue4 <- signif(length_deviation_pvalue4, digits = 1)
# 
# # getting rsquared from model
# length_deviation_rsquared4 <- signif(length_deviation_lm_coeff4$adj.r.squared, digits = 2)



# # log10(deviaton) ~ medianDomainLength linear model
# length_deviation_lm5 <- lm(Log10DeviationFromAvgAminoAcidFrequency ~ MedianDomainLength, pfam_metrics)
# length_deviation_lm_coeff5 <- summary(length_deviation_lm5)
# 
# # getting slope of model
# length_deviation_slope5 <- signif(length_deviation_lm_coeff5$coefficients[2, 1], digits = 2)
# 
# # getting full pvalue from model
# length_deviation_pvalue5 <- summary(length_deviation_lm5)$fstatistic %>% {unname(pf(.[1],.[2],.[3],lower.tail=F))}
# length_deviation_pvalue5 <- signif(length_deviation_pvalue5, digits = 1)
# 
# # getting rsquared from model
# length_deviation_rsquared5 <- signif(length_deviation_lm_coeff4$adj.r.squared, digits = 2)



# (deviation*sqrt(median domain length)) ~ log10(medianDomainLength)
sigma_model_y_axis = pfam_metrics$DeviationFromAvgAminoAcidFrequency*sqrt(pfam_metrics$MedianDomainLength)    # * for SE of mean, / for SE of sums

# temp = pfam_deviations2_vect/sqrt(master_df$MedianDomainLength)    # was for attempted counts
sigma_model <- lm(sigma_model_y_axis ~ pfam_metrics$Log10MedianDomainLength)
summary(sigma_model)


# plots for all the above models
#
plot_sigma_model <- ggplot() +
  # deviation ~ medianDomainLength
  #
  # geom_jitter(data = master_df, aes(x = MedianDomainLength, y = DeviationFromAvgAAcomp), color="black", size=0.4, alpha=0.3) +
  # geom_smooth(data = master_df, method = "loess", aes(x = MedianDomainLength, y = DeviationFromAvgAAcomp), color = 'midnightblue') +
  # scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  # annotate(geom = "text", label = paste('R^2 == ', length_deviation_rsquared), x = 1000, y = 0.51, size = 40/.pt, parse = TRUE) +
  # # annotate(geom = "text", label = paste("slope = ", length_deviation_slope),  x = 1000, y = 0.45, size = 40/.pt) +
  # annotate(geom = "text", label = bquote("P <" ~ 10^-16), x = 1000, y = 0.39, size = 40/.pt) +
  # labs(y = "Deviation from mean amino acid frequency") +
  # labs(x = "Median Domain length") +
  
  # log10(deviation) ~ log10(medianDomainLength)
  # 
  # geom_jitter(data = master_df, aes(x = log10MedianDomainLength, y = log10DeviationFromAvgAAcomp), color="black", size=0.4, alpha=0.3) +
  # geom_smooth(data = master_df, method = "loess", aes(x = log10MedianDomainLength, y = log10DeviationFromAvgAAcomp), color = 'midnightblue') +
  # geom_hline(yintercept = -0.652, color = "red", linewidth = 0.25) +
  # scale_x_continuous(limits = c(0.8, 3.2), breaks = seq(0.8, 3.2, 0.4), labels = signif(10^seq(0.8, 3.2, 0.4), digits = 2)) +
  # scale_y_continuous(limits = c(-1.6, -0.1), breaks = seq(-1.6, -0.1, 0.3), labels = signif(10^seq(-1.6, -0.1, 0.3), digits = 2)) +
  # annotate(geom = "text", label = paste('R^2 == ', length_deviation_rsquared2), x = 2.8, y = -0.15, size = 40/.pt, parse = TRUE) +
  # annotate(geom = "text", label = paste("slope = ", length_deviation_slope2),  x = 2.8, y = -0.27, size = 40/.pt) +
  # annotate(geom = "text", label = bquote("P <" ~ 10^-16),  x = 2.8, y = -0.39, size = 40/.pt) +
  # labs(x = "Median Domain length (log10(x) transformed)") +
  # labs(y = "Deviation from mean amino acid \n frequency (log10(x) transformed)") +
  
  # deviation ~ 1/((medianDomainLength)^1/3)
  #
  # geom_jitter(data = master_df, aes(x = temp_domainLength, y = DeviationFromAvgAAcomp), color="black", size=0.4, alpha=0.3) +
  # geom_smooth(data = master_df, method = "loess", aes(x = temp_domainLength, y = DeviationFromAvgAAcomp), color = 'midnightblue') +
  # scale_x_continuous(limits = c(0.085, 0.5), breaks = seq(0.085, 0.5, 0.1), labels = temp_domainLength_ticks) +
  # scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.6, 0.1)) +
  # annotate(geom = "text", label = paste('R^2 == ', length_deviation_rsquared3), x = 0.2, y = 0.51, size = 40/.pt, parse = TRUE) +
  # annotate(geom = "text", label = paste("slope = ", length_deviation_slope3),  x = 0.22, y = 0.45, size = 40/.pt) +
  # annotate(geom = "text", label = bquote("P <" ~ 10^-16),  x = 0.22, y = 0.39, size = 40/.pt) +
  # labs(x = "Median Domain length (1/x^(1/3) transformed)") +
  # labs(y = "Deviation from mean amino acid frequency") +
  
  # deviation ~ 1/(medianDomainLength)^1/2) linear model
  #
  # geom_jitter(data = master_df, aes(y = (DeviationFromAvgAAcomp*sqrt(MedianDomainLength)), x = log10MedianDomainLength), color="black", size=0.4, alpha=0.3) +
  # geom_smooth(data = master_df, method = "loess", aes(y = (DeviationFromAvgAAcomp*sqrt(MedianDomainLength)), x = log10MedianDomainLength), color = 'midnightblue') +
  # geom_hline(yintercept = 0.233, color = "red", linewidth = 0.25) +
  # scale_y_continuous(limits = c(0, 5.5), breaks = seq(0, 5.5, 1)) +
  # scale_x_continuous(limits = c(0, 3.2), breaks = seq(0, 3.2, 0.8)) +
  # annotate(geom = "text", label = paste('R^2 == ', length_deviation_rsquared4), x = 0.8, y = 5, size = 40/.pt, parse = TRUE) +
  # annotate(geom = "text", label = paste("slope = ", length_deviation_slope4),  x = 0.8, y = 4.5 , size = 40/.pt) +
  # annotate(geom = "text", label = bquote("P <" ~ 10^16),  x = 0.8, y = 4, size = 40/.pt) +
  # labs(x = "log(medianLength)") +
  # labs(y = "deviation*sqrt(medianLength)") +
  
  # log10(deviation) ~ medianDomainLength
  #
  # geom_jitter(data = master_df, aes(x = MedianDomainLength, y = log10DeviationFromAvgAAcomp), color="black", size=0.4, alpha=0.3) +
  # geom_smooth(data = master_df, method = "loess", aes(x = MedianDomainLength, y = log10DeviationFromAvgAAcomp), color = 'midnightblue') +
  # scale_x_continuous(limits = c(0, 1600), breaks = seq(0, 1600, 400)) +
  # scale_y_continuous(limits = c(-1.6, -0.25), breaks = deviation_axis, labels = log10_deviation_axis) +
  # annotate(geom = "text", label = paste('R^2 == ', length_deviation_rsquared5), x = 1200, y = -0.35, size = 40/.pt, parse = TRUE) +
  # annotate(geom = "text", label = paste("slope = ", length_deviation_slope5),  x = 1200, y = -0.48, size = 40/.pt) +
  # annotate(geom = "text", label = bquote("P <" ~ 10^-16),  x = 1200, y = -0.6, size = 40/.pt) +
  # labs(y = "Deviation from mean amino acid frequency \n (log10(x) transformed)") +
  # labs(x = "Median Domain length") +

  # (deviation*sqrt(median domain length)) ~ log10(medianDomainLength)
  #
  geom_jitter(data = pfam_metrics, aes(x = Log10MedianDomainLength, y = sigma_model_y_axis), color="black", size=0.4, alpha=0.3) +
  geom_smooth(data = pfam_metrics, method = "loess", aes(x = Log10MedianDomainLength, y = sigma_model_y_axis), color = 'midnightblue') +
  geom_hline(yintercept = sigma, color = "red", linewidth = 0.25) +
  scale_x_continuous(limits = c(0.8, 3.2), breaks = seq(0.8, 3.2, 0.4), labels = signif(10^seq(0.8, 3.2, 0.4), digits = 2)) +
  scale_y_continuous(limits = c(0.25, 5.5), breaks = seq(0.25, 5.5, 1.25), labels = signif(seq(0.25, 5.5, 1.25), digits = 2)) +
  # annotate(geom = "text", label = paste('R^2 == ', 0.090), x = 1.2, y = 5.2, size = 40/.pt, parse = TRUE) +
  # annotate(geom = "text", label = paste("slope = ", 0.448),  x = 2.8, y = -0.27, size = 40/.pt) +
  # annotate(geom = "text", label = bquote("P <" ~ 10^-16),  x = 1.2, y = 4.5, size = 40/.pt) +
  labs(x = "Median domain length") +
  labs(y = "Per amino acid contribution \n to the deviation") +
  # labs(y = "Deviation from mean amino acid \n frequency * sqrt(median domain length)") 
  
  # constants for all plots
  theme_bw() + 
  theme_linedraw() +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(axis.text = element_text(size = 20)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank())

# size 10x8 
pdf("results/sigma_model_plot.pdf", width = 10, height = 8)
plot_sigma_model
dev.off()

#####



##### hypothesis: larger domains have more hydrophobic AA in their core and this accounts for the push of the LOESS above 
##### sigma in above plot
##### - 
##### plot each amino acid freq as a function of log10(median domain length)
#####
# make 20 LOESS plots of freq ~ log10(median_domain_length) for visualization
plot_list <- list()
index <- 1

## cut off x-axis at length of 40
for (i in amino_acids) {
  print(i)
  model <- loess(pfam_metrics[[i]] ~ pfam_metrics$Log10MedianDomainLength)
  
  plot <- ggplot() +
    geom_smooth(data = pfam_metrics, method = "loess", aes(x = Log10MedianDomainLength, y = .data[[i]])) +
    geom_jitter(data = pfam_metrics, aes(x = Log10MedianDomainLength, y = .data[[i]]), color="black", size=0.5, alpha=0.4) +
    scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05)) +
    # scale_x_continuous(limits = c(0.8, 3.2), breaks = seq(0.8, 3.2, 0.4), labels = signif(10^seq(0.8, 3.2, 0.4), digits = 2)) +
    scale_x_continuous(limits = c(1.6, 3.2), breaks = seq(1.6, 3.2, 0.4), labels = signif(10^seq(1.6, 3.2, 0.4), digits = 2)) +
    labs(x = "Median domain length") +
    labs(y = "Amino acid frequency") +
    geom_text(aes(x = 3.1, y = 0.125), label = i, size = 30/.pt) +
    theme(axis.title = element_text(size = 30, face = "bold")) +
    theme(axis.text = element_text(size = 20)) 
  
  # print(plot)
  
  plot_list[[index]] <- ggplotGrob(plot)
  index <- index + 1
}

# arrange all plots on one page
# size 36x36
pdf("results/20_loess_log10_median_domain_length_plot.pdf", width = 36, height = 36)
grid.arrange(grobs = plot_list, nrow = 5, ncol = 4)
dev.off()

# hydrophobic AA
# figure with order by greater solvent accessibility, increasing from bottom to top
# from more soluble surface aa to less soluble core aa

# sort the solubility values in increasing order, higher values more hydrophobic and found in cores
sorted_mean_RSA <- amino_acid_specific_freq$MeanRelativeSolventAccessibility
names(sorted_mean_RSA) <- amino_acids
sorted_mean_RSA <-  sort(sorted_mean_RSA, decreasing = T)

# make a sub dataframe of all the amino acid frequencies
pfam_metrics_sub <- pfam_metrics[2:21]

# match the sub dataframes names to the names of the solubility vector, thus putting the sub dataframe in increasing order
pfam_metrics_sub <- pfam_metrics_sub[,match(names(sorted_mean_RSA), colnames(pfam_metrics_sub))]
                     
# cumulatively sum the columns
pfam_metrics_matrix <- sapply(2:20, function(col) {rowSums(pfam_metrics_sub[, 1:col])})

# coerce matrix to dataframe
pfam_metrics_df <- as.data.frame(pfam_metrics_matrix)

# add back in first column that wasn't summed
pfam_metrics_df <- cbind(pfam_metrics_sub$K, pfam_metrics_df)

# change dataframe column names
colnames(pfam_metrics_df) <- paste0('cumulative_', colnames(pfam_metrics_sub)[1:20])

# add PfamUID column and log10(medianDomainLength) column
pfam_metrics_cummulative_df <- cbind(pfam_metrics$PfamUID, pfam_metrics$Log10MedianDomainLength, pfam_metrics_df)
pfam_metrics_cummulative_df <- pfam_metrics_cummulative_df %>% rename("PfamUID" = `pfam_metrics$PfamUID`, "Log10MedianDomainLength" = `pfam_metrics$Log10MedianDomainLength`)

pfam_metrics_cummulative_df <- pfam_metrics_cummulative_df %>% 
  pivot_longer(cols = -c(PfamUID, Log10MedianDomainLength), 
               names_to = "cumulative_aa", 
               values_to = "cumulative_freq")

# decreasing solvent accessibility from bottom to top
legend_order <- c("cumulative_K", "cumulative_E","cumulative_D",
                  "cumulative_Q", "cumulative_R", "cumulative_P",
                  "cumulative_N", "cumulative_S", "cumulative_T",
                  "cumulative_H", "cumulative_G", "cumulative_A",
                  "cumulative_Y", "cumulative_W", "cumulative_M",
                  "cumulative_L", "cumulative_V", "cumulative_F",
                  "cumulative_I", "cumulative_C")

# edit for plot legend (decreasing order of solvent accessibility)
pfam_metrics_cummulative_df$cumulative_aa <- factor(pfam_metrics_cummulative_df$cumulative_aa, levels = legend_order)

##
## cut off x-axis at 40 (1.6)
## fix legend to match order of lines
##
solventAccessability_plot <- ggplot(pfam_metrics_cummulative_df, aes(x = Log10MedianDomainLength, y = cumulative_freq, color = cumulative_aa, 
                                              group = cumulative_aa)) +
  geom_smooth(method = "loess") +
  scale_color_manual(values = c("cumulative_K" = "blue","cumulative_E" = "red","cumulative_D" = "green",
                                "cumulative_Q" = "orange", "cumulative_R" = "purple", "cumulative_N" = "yellow",
                                "cumulative_P" = "cyan", "cumulative_S" = "magenta", "cumulative_T" = "brown",
                                "cumulative_G" = "gray", "cumulative_H" = "pink", "cumulative_A" = "lightblue",
                                "cumulative_Y" = "lightgreen", "cumulative_W" = "lightcoral", "cumulative_M" = "peachpuff",
                                "cumulative_L" = "orchid", "cumulative_V" = "lightyellow", "cumulative_F" = "lightcyan",
                                "cumulative_I" = "darkgray", "cumulative_C" = "black")) +
  # scale_x_continuous(limits = c(0.8, 3.2), breaks = seq(0.8, 3.2, 0.4), labels = signif(10^seq(0.8, 3.2, 0.4), digits = 2)) +
  scale_x_continuous(limits = c(1.6, 3.2), breaks = seq(1.6, 3.2, 0.4), labels = signif(10^seq(1.6, 3.2, 0.4), digits = 2)) +
  scale_y_continuous(limits = c(0, 1.01), breaks = seq(0,1,0.25), labels = seq(0,1,0.25)) +
  labs(x = "Median domain length") +
  labs(y = "Cumulative frequencies") +
  guides(colour = guide_legend(title = "cumulative AA")) +
  theme(axis.title = element_text(size = 30, face = "bold")) +
  theme(axis.text = element_text(size = 20)) 

# size 10x8
pdf("results/solvent_accessability_plot.pdf", width = 10, height = 8)
solventAccessability_plot
dev.off()

#####
  
