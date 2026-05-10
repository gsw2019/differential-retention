#####Libraries#####
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(car)


####
#### this script builds linear models predicting loss rate from fldpnn ISD, dn/ds, and fldpnn ISD + dn/ds to investigate
#### their effects on each other
####
#### calculated metrics:
####    None
####
#### figures:
####    None
####
#### from these models, it can be observed that fldpnn ISD and dn/ds have mostly independent affects, based on R^2 values
####


##### setup
#####

# read in all Pfam data
pfam_metrics <- read_csv("data/Pfam_metrics.csv")

#####


##### models
#####

# lm predicting Loss_norm with ISD and dN/ds and ISD + dN/ds
summary(lm(LossNormalized ~ AvgFldpnnDisordAcross10Species, data = pfam_metrics))                  # adjusted R^2 = 0.05671
summary(lm(LossNormalized ~ dN.dS, data = pfam_metrics))                                           # adjusted R^2 = 0.08259
summary(lm(LossNormalized ~ AvgFldpnnDisordAcross10Species + dN.dS, data = pfam_metrics))          # adjusted R^2 = 0.1191

#####
