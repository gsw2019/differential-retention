#######Libraries#######
library(ggplot2)
library(geoR)
library(grid)
library(gridExtra)
library(cmna)
library(dplyr)
library(tidyverse)
library(cutpointr)
library(forcats)


####
#### this script utilizes amino acid frequency and maximum likelihood loss rates of 6700 Pfams extracted from the Masel 
#### Lab MySQL database to model loss rate predicted by amino acid frequency
####
#### calculated metrics:
####    1. normalized loss rates
####    2. the amino acid frequency corresponding to minimum loss or, the optimal amino acid frequency, for each amino acid
####		3. the loss rate at the optimal amino acid frequencies
####
#### figures: 
####    1. 20 LOESS models figure, where each model is a prediction of loss rate by amino acid frequency. The red point on the
####       plots represent the optimal amino acid frequency. Can also add boxplots of clan age data, must uncomment in build chunk. 
####
#### from this figure it can be observed that there is an optimal intermediate amino acid frequency for each of the 20 amino acids.
#### The optimal amino acid frequencies are key indicators of differential retention. This frequency data is showing that there is
#### variable retention dependent on a peptides amino acid frequency
####


# suppress Rplots.pdf generation 
pdf(NULL)

# create results folder if doesnt exist
dir.create("results", showWarnings = FALSE)


##### setup
#####

# read in pfam data
pfam_metrics <- read.csv("data/Pfam_metrics.csv")

# read in age data from Sawsan
barplotData <- read.csv("data/Clan_AgesandAAC.csv", header = T)

# normalizing loss rate
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



##### build individual plots, models, and a data frame that holds minimum amino acid frequency and minimum loss for each amino acid
#####

# create new column to group data by
for (i in 1:nrow(barplotData))(
  if (barplotData$Clan_ancestor[i] == 'preLUCA') {barplotData$PhylostrataAge[i] <- 1}
  else if (barplotData$Clan_ancestor[i] == 'LUCA') {barplotData$PhylostrataAge[i] <- 2}
  else if (barplotData$Clan_ancestor[i] == 'postLUCA') {barplotData$PhylostrataAge[i] <- 3}
  else if (barplotData$Clan_ancestor[i] == 'modern') {barplotData$PhylostrataAge[i] <- 4})


# vector of all amino acid one letter codes
amino_acids <- colnames(pfam_metrics)[2:21]

# initialize df and vectors to store data each iteration
minAAfreq_minLoss_df <- data.frame(amino_acids)
minLossAAfreq <- vector()
minLoss <- vector()

# make list of plots to arrange onto one page
plot_list <- list()
index <- 1

# plotting the ML loss rate predicted by amino acid frequency for each amino acid
for (i in amino_acids) {
    print(i)
    model <- loess(pfam_metrics$LossNormalized ~ pfam_metrics[[i]])
    # print(summary(model))
    
    # each amino acid's unique features for its figure
    if (i == "A") {
      optRange = c(0, 0.125)
      ylabel = "Losses/million years"
      ytext = element_text()
      yticks = element_line()
      name <- "Alanine"
    }     
    if (i == "R") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Arginine"
    }     
    if (i == "N") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Asparagine"
    }    
    if (i == "D") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Aspartic Acid"
    }     
    if (i == "C") {
      optRange = c(0, 0.125)
      ylabel = "Losses/million years"
      ytext = element_text()
      yticks = element_line()
      name <- "Cysteine"
    }     
    if (i == "E") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Glutamic Acid"
    }     
    if (i == "Q") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Glutamine"
    }    
    if (i == "G") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Glycine"
    } 
    if (i == "H") {
      optRange = c(0, 0.125)
      ylabel = "Losses/million years"
      ytext = element_text()
      yticks = element_line()
      name <- "Histidine"
    } 
    if (i == "I") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Isoleucine"
    }     
    if (i == "L") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Leucine"
    }    
    if (i == "K") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Lysine"
    }     
    if (i == "M") {
      optRange = c(0, 0.125)
      ylabel = "Losses/million years"
      ytext = element_text()
      yticks = element_line()
      name <- "Methionine"
    }     
    if (i == "F") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Phenylalinine"
    }     
    if (i == "P") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Proline"
    }     
    if (i == "S") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Serine"
    }     
    if (i == "T") {
      optRange = c(0, 0.125)
      ylabel = "Losses/million years"
      ytext = element_text()
      yticks = element_line()
      name <- "Threonine"
    }    
    if (i == "W") {
      optRange = c(0, 0.025)      # if c(0, 0.125) visually incorrect minimum?
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Tryptophan"
    }     
    if (i == "Y") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Tyrosine"
    }    
    if (i == "V") {
      optRange = c(0, 0.125)
      ylabel = ""
      ytext = element_blank()
      yticks = element_blank()
      name <- "Valine"
    }
    
    xBreaks = seq(0, 0.125, 0.025)
    xEnd = 0.125
    
    # function used in optimize()
    f <- function(x) {
      predict(model, x)
    }
    
    # optimize function
    optimization <- optimize(f, optRange)
    
    # get min x (AA comp)
    min_x <- optimization$minimum
    
    # get min y (loss rate)
    min_y <- optimization$objective
    
    # get min y (back transformed loss rate)
    corrected_min_y <- signif(bc.backtransform(min_y, lambda), digits = 1)
    
    # boxplot stuff
    AA_median <- barplotData  %>%
      group_by(PhylostrataAge) %>%
      summarize(median = median(.data[[i]])) %>%
      ungroup()
    
    plot <- ggplot()+
      #### base plot
      geom_jitter(data = pfam_metrics, aes(x=.data[[i]], y=LossNormalized), color="black", size=0.4, alpha=0.3) +          # optional comment out for boxplots
      geom_smooth(data = pfam_metrics, method = "loess", aes(x=.data[[i]], y=LossNormalized), color = 'midnightblue') +
      scale_y_continuous(breaks = c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), labels = signif(bc.backtransform(c(0, -1, -2, -3, -4, -5, -6, -7, -8, -9), lambda), digits = 1), limits = c(-9, 0)) +
      scale_x_continuous(limits = c(0, xEnd), breaks = xBreaks) +           
    	geom_point(aes(x = min_x, y = min_y), color = "red", size = 6) +
      labs(x = paste(name, "frequency", sep = " ")) +
      labs(y = ylabel) +
      
      #### for adding age boxplots to base plot
      #### will add some extreme and sparse frequency data, curves wont look the same as without boxplots. Needed to show whisker data
      # scale_x_continuous(limits = c(0, 0.25), breaks = seq(0, 0.25, 0.05)) +
      # geom_jitter(aes(x = 0.25, y = 0), color = "black", size = 0.1, alpha = 0.3) +
      # geom_jitter(aes(x = 0.25, y = -9), color = "black", size = 0.1, alpha = 0.3) +
      # geom_boxplot(data = barplotData, aes(y = -4, x = .data[[i]], group = as.factor(PhylostrataAge), color = as.factor(PhylostrataAge)), width = 3, show.legend = FALSE) +     # show.legend = TRUE makes legend on each plot
      # scale_color_manual(labels = c("preLUCA","LUCA", "postLUCA", "modern"), values = c("blue", "red", "green", "orange")) +
      # guides(colour = guide_legend(reverse = TRUE, title = "Phylostrata Age")) +

      #### theme settings
      theme_bw() + 
      theme_linedraw() + 
      theme(text = element_text(size = 30), axis.text = element_text(size = 25)) +
      theme(panel.grid.minor = element_blank()) +
      theme(panel.grid.major = element_blank()) +
    	theme(axis.text.y = ytext, axis.ticks.y = yticks)

    # print(plot)
    plot_list[[index]] <- ggplotGrob(plot)
    index <- index + 1
    
    # add each min_x to minAA vector and each min_y to minLoss vector
    minLossAAfreq <- c(minLossAAfreq, min_x)
    minLoss <- c(minLoss, min_y)
} 

#####


##### organize data frame and final figure
#####

# data frame with minimum amino acid freq and minimum loss from each loess plot
minLoss_bac <- signif(bc.backtransform(minLoss, lambda), digits = 1)
minLossAAfreq_minLoss_df <- cbind(minAAfreq_minLoss_df, minLossAAfreq, minLoss, minLoss_bac)

# arrange all plots on one page
# 36x36 for PDF
# 3600x3600 for SVG
graphics.off()
svg("results/Pfam_metrics_analysis_plot.svg", width = 36, height = 36)
grid.arrange(grobs = plot_list, nrow = 5, ncol = 4)
dev.off()

#####
