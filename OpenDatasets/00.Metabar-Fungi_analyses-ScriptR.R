##########################################################################
#                                                                        #
#                                                                        #
#          The distribution of mycorrhizal fungi on tree barks           #
#          correlates with the host preference of the tropical           #
#         epiphytic orchid Bulbophyllum variegatum on La Réunion         #
#                                                                        #
#                        Petrolli et al. 2026                            #
#                                                                        #
#                                                                        #
##########################################################################

## ReadMe ##
# R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
# Copyright (C) 2021 The R Foundation for Statistical Computing
# Platform: x86_64-apple-darwin17.0 (64-bit)
#
# Script from Remi Petrolli (remi.petrolli@laposte.net)
# Date : 29/04/2024
#
# Note : each section can be run separately once the tables are loaded
## ##  ## ## 

rm(list = ls())
par(mfrow = c(1, 1))
par(ask=FALSE)

# Please indicate path to print plots :
SortiePath <- ""
date <- ""

#### Packages ###############################################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(vegan) 
library(reshape2)
library(adespatial)

# DDA
library(ANCOMBC)
library(tidyverse)
library(phyloseq)
library(ALDEx2)
##

library(car) # Anova
library(wesanderson)

#### Functions ##############################################################

# CompareTwoLists : detects differences in lists and indicates the source (position or difference in length)
CompareTwoLists <- function(A, B){
  A <- get(A)
  B <- get(B)
  if (length(A)!=length(B)){
    stop("Difference in length")}
  for (k in 1:length(A)){
    if (A[k]!=B[k]){
      stop(paste("Different lists, position:", k))}}
  paste("Identical lists")}

# getTaxoForAllOTUs : return taxonomy given OTU_ID
getTaxoForAllOTUs <- function(x){ #x = OTU_ID
  taxtxt <- TaxoTable$Taxonomy[TaxoTable$OTU_ID%in%x]
  GroupofOTU <- WhichFam(taxtxt)
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichOrder(taxtxt)}
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichClass(taxtxt)}
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichPhyllum(taxtxt)}
  return(GroupofOTU)
}
getTaxoForAllOTUsFromGenus2 <- function(x){ 
  GroupofOTU <- WhichGenera(x)
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichFam(x)}
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichOrder(x)}
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichClass(x)}
  if (GroupofOTU == "unidentified"){
    GroupofOTU <- WhichPhyllum(x)}
  return(GroupofOTU)
}
WhichSpecies <- function(x){
  debut <- str_locate(x, ";s__")
  str_sub(x, debut[2]+1, -1) # Renvoi juste le nom de la famille
}
WhichGenera <- function(x){
  fin <- str_locate(x, ";s__")
  debut <- str_locate(x, ";g__")
  str_sub(x, debut[2]+1, fin[1]-1) # Renvoi juste le nom de la famille
}
WhichFam <- function(x){
  fin <- str_locate(x, ";g__")
  debut <- str_locate(x, ";f__")
  str_sub(x, debut[2]+1, fin[1]-1) # Renvoi juste le nom de la famille
}
WhichOrder <- function(x){
  fin <- str_locate(x, ";f__")
  debut <- str_locate(x, ";o__")
  str_sub(x, debut[2]+1, fin[1]-1) # Renvoi juste le nom de l'rdre
}
WhichClass <- function(x){
  fin <- str_locate(x, ";o__")
  debut <- str_locate(x, ";c__")
  str_sub(x, debut[2]+1, fin[1]-1) # Renvoi juste le nom de l'rdre
}
WhichPhyllum <- function(x){
  fin <- str_locate(x, ";c__")
  debut <- str_locate(x, ";p__")
  str_sub(x, debut[2]+1, fin[1]-1) 
}

# https://gist.github.com/mcgoodman/58c9d1257fd1625954a4ffa1c3301939
pairwise_permanova <- function(sp_matrix, group_var, dist, adj = "fdr", perm = 999) {

  require(vegan)
  
  ## list contrasts
  group_var <- as.character(group_var)
  groups <- as.data.frame(t(combn(unique(group_var), m = 2)))
  
  contrasts <- data.frame(
    group1 = groups$V1, group2 = groups$V2,
    R2 = NA, F_value = NA, df1 = NA, df2 = NA, p_value = NA,
    dfb = NA, R2b = NA, Fb = NA, p_valueb = NA
  )
  
  for (i in seq(nrow(contrasts))) {

    sp_subset <- group_var == contrasts$group1[i] | group_var == contrasts$group2[i] 
    contrast_matrix <- sp_matrix[sp_subset,]
    
    ## fit contrast using adonis
    fit <- vegan::adonis2(
      contrast_matrix ~ group_var[sp_subset],
      method = dist, 
      perm = perm
    )
    
    contrast_veg_matrix <- vegdist(contrast_matrix, method = dist)
    betad <- anova( betadisper(contrast_veg_matrix, factor(group_var[sp_subset])) )
    
    contrasts$R2[i] <- round(fit$R2[1], digits = 3)
    contrasts$F_value[i] <- round(fit[["F"]][1], digits = 3)
    contrasts$df1[i] <- fit$Df[1]
    contrasts$df2[i] <- fit$Df[2]
    contrasts$p_value[i] <- fit$`Pr(>F)`[1]
    
    contrasts$dfb[i] <- betad$Df[1]
    contrasts$R2b[i] <- round( betad$`Sum Sq`[1]/(betad$`Sum Sq`[1]+betad$`Sum Sq`[2]) , digits = 3)
    contrasts$Fb[i] <- round( betad$`F value`[1] , digits = 3)
    contrasts$p_valueb[i] <- betad$`Pr(>F)`[1]
    
  }
  
  ## adjust p-values for multiple comparisons
  contrasts$p_value <- round(p.adjust(contrasts$p_value, method = adj), digits = 3)
  contrasts$p_valueb <- round(p.adjust(contrasts$p_valueb, method = adj), digits = 3)
  
  return(list(
    contrasts = contrasts, 
    "p-value adjustment" = adj, 
    permutations = perm
  ))
}

#### 1. Load tables ########################################################
#### 1.1. Taxonomy ####

TaxoTable <- readRDS("Path_to/24_05_07_TaxoTable.rds")

#### 1.2. Metadata ####

Metadata_Bark <- readRDS("Path_to/24_05_07_Metadata_Bark.rds")

Metadata_Roots <- readRDS("Path_to/24_05_07_Metadata_Roots.rds")

#### 1.3. OTU Tables ####

# Table A:
# Relative abundance data with primers ITS86-F/ITS4 only
TableA <- read.table("Path_to/23_09_08_Table_A_OTU97_Full_Filtered_AbsAbundances_WithThreshold-0.8.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Table B:
#Presence/absence data (i.e., binary) with primers ITS86-F/ITS4 and 5.8S-OF/ITS4-Tul
TableB <- read.table("Path_to/23_09_08_Table_B_OTU97_Full_Filtered_Binary_WithThreshold-0.8.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Table C:
# Same as Table B but with orchid mycorrhizal fungi only, as defined in Dearnaley et al., (2012)
# "Orchid mycorrhizas: molecular ecology, physiology, evolution and conservation aspects"
TableC <- read.table("Path_to/23_09_08_Table_C_OTU97_Full_Filtered_Binary_WithThreshold-0.8.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Table D: 
# Relative abundance data with primers 5.8S-OF/ITS4-Tul only
TableD <- read.table("Path_to/23_09_08_Table_D_OTU97_Full_Filtered_AbsAbundances_WithThreshold-0.8.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#### 2. Beta diversity analyses - NMDS ##################################### 
#### 2.1. NMDS - Fig. S6-B  - Jaccard - Bark ################################

Table <- TableB
rownames(Table) <- Table$OTU_ID
Table <- Table[,c(9:ncol(Table))]
Table <- Table[,which(colnames(Table)%in%Metadata_Bark$Ech_ID)]  
Table <- Table[apply(Table, 1, sum)>0,]

###
Table <- Table[which(apply(Table, 1, sum)>0),] # keep OTUs in bark only
Table <- Table[c(1:50),] 
Table <- Table[, which(apply(Table, 2, sum)>0)] 
Table <- as.data.frame(t(Table))
Table.dist <- vegdist(Table, method = "jaccard", binary = TRUE)
set.seed(920)
selec_NMDS <- metaMDS(Table.dist,k=2,trymax=100, halfchange = FALSE)

x <- selec_NMDS$points[,1]
y <- selec_NMDS$points[,2]
NMDS_PlotTab <- data.frame(x=x, y=y)

NMDS_PlotTab$Host <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID%in%x)]}))
NMDS_PlotTab$Host2 <- unlist(lapply(NMDS_PlotTab$Host, function(x){
  ifelse(x%in%c("AgaSal", "AntBor", "LabCal", "NuxVer"), x, "Other")
}))

ColorLegend <- c("AgaSal"="#8f2d56", "LabCal"="#ec8931", "NuxVer"="#f5b940", "AntBor"="#9bc9e3", 
                 "CyaBor" = "#0f2f45", "DorApe" = "#165577",
                 "SyzCim" = "#49a2c3","SyzBor"="#49a2c3", "HomPan" = "#cde5f2", "Other" = "#d6d6d6")

ggplot(data = NMDS_PlotTab, aes(x=x, y=y)) + 
  geom_point(aes(fill = factor(Host2)), color = "black", size = 5, alpha = .90, shape = 21) + 
  theme_classic() + 
  scale_fill_manual(values = ColorLegend)+
  scale_color_manual(values = ColorLegend)+
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dotted", color = "darkgrey") +
  xlab("NMDS1") + ylab("NMDS2") + #ylim(-.3, .55) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        legend.title = element_text(size=0)) +
  stat_ellipse(type = "t", geom = "polygon", aes(group = factor(Host2), fill = factor(Host2), col = factor(Host2)), alpha = .05)+
  annotate(geom="text", x=-0.35, y=-0.25, label=paste("Stress = ", round(selec_NMDS$stress, 4), sep =""),
         col = "black", size = 4.8)
#ggsave(paste(SortiePath,date, "NMDS_TableB_Total-Species-Ellipse.pdf",  sep=""), width = 8, height = 6, plot = last_plot())

## Significance

listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

NMDS_PlotTab$Host2 <- factor(NMDS_PlotTab$Host2)

tab_adonis <- adonis2(formula = Table ~ Host2, data = NMDS_PlotTab, method = "jaccard", nperm = 999) 
tab_adonis

# adonis2(formula = Table ~ Host2, data = NMDS_PlotTab, method = "jaccard", nperm = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# Host2      4   19.816 0.14036 12.695  0.001 ***
# Residual 311  121.367 0.85964                  
# Total    315  141.183 1.00000                


mod <- betadisper(Table.dist, factor(NMDS_PlotTab$Host2))
anova(mod)

# Response: Distances
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups      4 0.39254 0.098135  10.548 5.012e-08 ***
#   Residuals 311 2.89336 0.009303
0.39254/(0.39254+2.89336)

## Pairwise significance

listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

pperm_res <- pairwise_permanova(sp_matrix = Table, group_var = NMDS_PlotTab$Host2, dist = "jaccard")
pperm_res

#### 2.2. NMDS - Fig. S6-C  - Bray - Bark ###################################

Table <- TableA

rownames(Table) <- Table$OTU_ID
Table <- Table[,c(9:ncol(Table))]
Table <- subset(Table, select = -c(E287)) # alters NMDS calculation (only one OTU in this sample)
Table <- Table[,which(colnames(Table)%in%Metadata_Bark$Ech_ID)]  
Table <- Table[which(apply(Table, 1, sum)>0),] # keep OTUs in bark only

Table <- Table[, which(apply(Table, 2, sum)>0)] 
Table <- as.data.frame(t(Table))
Table.dist <- vegdist(Table, method = "bray", binary = FALSE)
set.seed(920)
selec_NMDS <- metaMDS(Table.dist,k=2,trymax=100, halfchange = FALSE) 

x <- selec_NMDS$points[,1]
y <- selec_NMDS$points[,2]
NMDS_PlotTab <- data.frame(x=x, y=y)

NMDS_PlotTab$Host <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID%in%x)]}))
NMDS_PlotTab$Host2 <- unlist(lapply(NMDS_PlotTab$Host, function(x){
  ifelse(x%in%c("AgaSal", "AntBor", "LabCal", "NuxVer"), x, "Other")
}))

ColorLegend <- c("AgaSal"="#8f2d56", "LabCal"="#ec8931", "NuxVer"="#f5b940", "AntBor"="#9bc9e3", 
                 "CyaBor" = "#0f2f45", "DorApe" = "#165577",
                 "SyzCim" = "#49a2c3","SyzBor"="#49a2c3", "HomPan" = "#cde5f2", "Other" = "#d6d6d6")

ggplot(data = NMDS_PlotTab, aes(x=x, y=y)) + 
  geom_point(aes(fill = factor(Host2)), color = "black", size = 5, alpha = .90, shape = 21) + 
  theme_classic() + 
  scale_fill_manual(values = ColorLegend)+
  scale_color_manual(values = ColorLegend)+
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dotted", color = "darkgrey") +
  xlab("NMDS1") + ylab("NMDS2") + #ylim(-.3, .55) +
  xlim(-0.4, .35)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size=28)) +
  stat_ellipse(type = "t", geom = "polygon", aes(group = factor(Host2), fill = factor(Host2), col = factor(Host2)), alpha = .05)+
  annotate(geom="text", x=-0.28, y=-0.23, label=paste("Stress = ", round(selec_NMDS$stress, 4), sep =""),
           col = "black", size = 4.8)
#ggsave(paste(SortiePath,date, "NMDS_TableA_Total-Species-Ellipse.pdf",  sep=""), width = 8, height = 6, plot = last_plot())

## Significance

listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

NMDS_PlotTab$Host2 <- factor(NMDS_PlotTab$Host2)

tab_adonis <- adonis2(formula = Table ~ Host2, data = NMDS_PlotTab, method = "bray", nperm = 999) 
tab_adonis

# adonis2(formula = Table ~ Host2, data = NMDS_PlotTab, method = "bray", nperm = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# Host2      4   14.231 0.07892 7.9041  0.001 ***
#   Residual 369  166.095 0.92108                  
# Total    373  180.326 1.00000               


mod <- betadisper(Table.dist, factor(NMDS_PlotTab$Host2))
amod <- anova(mod)
amod
amod$`Sum Sq`[1]/(amod$`Sum Sq`[1]+amod$`Sum Sq`[2])

# Response: Distances
# Df  Sum Sq   Mean Sq F value    Pr(>F)    
# Groups      4 0.07265 0.0181630  14.134 9.652e-11 ***
#   Residuals 369 0.47420 0.0012851

## Pairwise significance

listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

pperm_res <- pairwise_permanova(sp_matrix = Table, group_var = NMDS_PlotTab$Host2, dist = "bray")
pperm_res

#### 2.3. NMDS - Fig. 3-C - Jaccard - Total ###################################

Table <- TableB
rownames(Table) <- Table$OTU_ID
Table <- Table[,c(9:ncol(Table))]
Table <- subset(Table, select = -c(E287, E12)) # alters NMDS calculation (only one OTU in E287, almost one in E12)
Table <- Table[which(apply(Table, 1, sum)>0),] 
Table <- Table[, which(apply(Table, 2, sum)>0)] 

Table <- as.data.frame(t(Table))
Table.dist <- vegdist(Table, method = "jaccard", binary = TRUE)
set.seed(122)
selec_NMDS <- metaMDS(Table.dist,k=2,trymax=100, halfchange = FALSE)

x <- selec_NMDS$points[,1]
y <- selec_NMDS$points[,2]
NMDS_PlotTab <- data.frame(x=x, y=y)

# Just rename some samples
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp1")] <- "RBV53b(P)_1"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp2")] <- "RBV53b(P)_2"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp3")] <- "RBV53b(P)_3"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp4")] <- "RBV53b(P)_4"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp5")] <- "RBV53b(P)_5"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV94_1.1")] <- "RBV94_1"
  
  
NMDS_PlotTab$SampType <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  ValEch <- NULL
  if (x %in% Metadata_Bark$Ech_ID){
    if (Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID==x)]=="AgaSal"){ValEch <- "Bark: AgaSal"}
    else if (Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID==x)]!="AgaSal"){ValEch <- "Bark: Other"}
    else (print(x))
  }
  else if (x %in% Metadata_Roots$Sample_ID){
    ValEch <- "Root: BulVar"
  }
  else (print(x))
  ValEch
  }))

ggplot(data = NMDS_PlotTabr, aes(x=x, y=y)) + 
  geom_point(aes(fill = factor(SampType)), color = "black", size = 5, alpha = .80, shape = 21) + 
  theme_classic() + 
  scale_fill_manual(values = c("Bark: AgaSal" = "#bc3250", "Bark: Other" = "#02c39a",
                               "Root: BulVar" = "#fbb13c"))+
  scale_color_manual(values = c("Bark: AgaSal" = "#bc3250", "Bark: Other" = "#02c39a",
                                "Root: BulVar" = "#fbb13c"))+
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dotted", color = "darkgrey") +
  xlab("NMDS1") + ylab("NMDS2") +
  xlim(-0.4, .4)+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(.5, "cm"),
        legend.position = c(.15,.12),
        legend.background = element_rect(size = .25, linetype = "solid", colour = "black")) +
  stat_ellipse(type = "t", geom = "polygon", aes(group = factor(SampType), fill = factor(SampType), col = factor(SampType)), alpha = .2)+
  annotate(geom="text", x=-0.335, y=-0.075, label=paste("Stress = ", round(selec_NMDS$stress, 4), sep =""),
           col = "black", size = 4.8)
#ggsave(paste(SortiePath,date, "NMDS_FULLDATA_Total-Species-Ellipse2.pdf",  sep=""), width = 8, height = 6, plot = last_plot())

## Significance

listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

NMDS_PlotTab$SampType <- factor(NMDS_PlotTab$SampType)

tab_adonis <- adonis2(formula = Table ~ SampType, data = NMDS_PlotTab, method = "jaccard", nperm = 999) 
tab_adonis

# adonis2(formula = Table ~ SampType, data = NMDS_PlotTab, method = "jaccard", nperm = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# SampType   2    18.09 0.05166 20.265  0.001 ***
#   Residual 744   332.02 0.94834                  
# Total    746   350.11 1.00000

mod <- betadisper(Table.dist, factor(NMDS_PlotTab$SampType))
anova(mod)

# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups      2 0.23744 0.118721  91.806 < 2.2e-16 ***
#   Residuals 744 0.96213 0.001293

#### 2.4. NMDS - Fig. S6-A  - Bray - Total ###################################

Table <- TableA

rownames(Table) <- Table$OTU_ID
Table <- Table[,c(9:ncol(Table))]
Table <- subset(Table, select = -c(E287, E12)) # alters NMDS calculation (only one OTU in E287, almost one in E12)
Table <- Table[which(apply(Table, 1, sum)>0),] 
Table <- Table[, which(apply(Table, 2, sum)>0)] 

Table <- as.data.frame(t(Table))
Table.dist <- vegdist(Table, method = "bray", binary = FALSE)
set.seed(122)
selec_NMDS <- metaMDS(Table.dist,k=2,trymax=100, halfchange = FALSE)

x <- selec_NMDS$points[,1]
y <- selec_NMDS$points[,2]
NMDS_PlotTab <- data.frame(x=x, y=y)

rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp1")] <- "RBV53b(P)_1"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp2")] <- "RBV53b(P)_2"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp3")] <- "RBV53b(P)_3"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp4")] <- "RBV53b(P)_4"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV53bp5")] <- "RBV53b(P)_5"
rownames(selec_NMDS$points)[which(rownames(selec_NMDS$points)=="RBV94_1.1")] <- "RBV94_1"


NMDS_PlotTab$SampType <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  ValEch <- NULL
  if (x %in% Metadata_Bark$Ech_ID){
    if (Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID==x)]=="AgaSal"){ValEch <- "Bark: AgaSal"}
    else if (Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID==x)]!="AgaSal"){ValEch <- "Bark: Other"}
    else (print(x))
  }
  else if (x %in% Metadata_Roots$Sample_ID){
    ValEch <- "Root: BulVar"
  }
  else (print(x))
  ValEch
}))
unique(NMDS_PlotTab$SampType)


ggplot(data = NMDS_PlotTab, aes(x=x, y=y)) + 
  geom_point(aes(fill = factor(SampType)), color = "black", size = 5, alpha = .80, shape = 21) + 
  theme_classic() + 
  scale_fill_manual(values = c("Bark: AgaSal" = "#bc3250", "Bark: Other" = "#02c39a",
                               "Root: BulVar" = "#fbb13c"))+
  scale_color_manual(values = c("Bark: AgaSal" = "#bc3250", "Bark: Other" = "#02c39a",
                                "Root: BulVar" = "#fbb13c"))+
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dotted", color = "darkgrey") +
  xlab("NMDS1") + ylab("NMDS2") +
  #xlim(-0.4, .4)+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=20),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(.5, "cm"),
        legend.position = c(.83,.88),
        legend.background = element_rect(size = .25, linetype = "solid", colour = "black")) +
  stat_ellipse(type = "t", geom = "polygon", aes(group = factor(SampType), fill = factor(SampType), col = factor(SampType)), alpha = .2)+
  annotate(geom="text", x=-0.45, y=-0.22, label=paste("Stress = ", round(selec_NMDS$stress, 4), sep =""),
           col = "black", size = 4.8)
#ggsave(paste(SortiePath,date, "NMDS_FULLDATA_Total-Species-Ellipse2-TableA.pdf",  sep=""), width = 8, height = 6, plot = last_plot())


## Significance

listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

NMDS_PlotTab$SampType <- factor(NMDS_PlotTab$SampType)

tab_adonis <- adonis2(formula = Table ~ SampType, data = NMDS_PlotTab, method = "bray", nperm = 999) 
tab_adonis

# adonis2(formula = Table ~ SampType, data = NMDS_PlotTab, method = "bray", nperm = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# SampType   2    17.10 0.04802 18.765  0.001 ***
#   Residual 744   339.02 0.95198                  
# Total    746   356.12 1.00000                  


mod <- betadisper(Table.dist, factor(NMDS_PlotTab$SampType))
anova(mod)

# Response: Distances
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups      2 0.17179 0.085895  43.373 < 2.2e-16 ***
#   Residuals 744 1.47339 0.001980
0.17179/(0.17179+1.47339)

#### 2.5. NMDS - Fig. S6-D  - Jaccard - Roots ###############################

Table <- TableB
Table <- Table[,c(9:ncol(Table))]

SampleToRemove <- c("RAM1", "RAM2", "RAM3", "RBV32_3", "Eden4", "RBV99_2", "BB4", "D14", "RBV98_3")
# RAM are not B. variegatum and others alter NMDS calculation (i.e. convergence)

Table <- Table[,which(!colnames(Table)%in%SampleToRemove)]
Table <- Table[,which(colnames(Table)%in%Metadata_Roots$Sample_ID)]
Table <- Table[apply(Table, 1, sum)>0,] 
Table <- Table[c(1:100),]

Table <- Table[, which(apply(Table, 2, sum)>0)]
Table <- as.data.frame(t(Table))
set.seed(122)
Table.dist <- vegdist(Table, method = "jaccard", binary = TRUE)
selec_NMDS <- metaMDS(Table.dist,k=2,trymax=100, halfchange = FALSE)

x <- selec_NMDS$points[,1]
y <- selec_NMDS$points[,2]
NMDS_PlotTab <- data.frame(x=x, y=y)

NMDS_PlotTab$SampleSite <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Roots$Site[which(Metadata_Roots$Sample_ID%in%x)]}))
NMDS_PlotTab$Host <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Roots$Species[which(Metadata_Roots$Sample_ID%in%x)]}))
NMDS_PlotTab$Year <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Roots$Year[which(Metadata_Roots$Sample_ID%in%x)]}))
NMDS_PlotTab$SYH <- paste(NMDS_PlotTab$SampleSite, NMDS_PlotTab$Year, NMDS_PlotTab$Host, sep ="-")

ggplot(data = NMDS_PlotTab, aes(x=x, y=y)) + 
  geom_point(aes(fill = factor(SYH)), color = "black", size = 5, alpha = .90, shape = 21) + 
  theme_classic() + 
  scale_fill_manual(values = c("ML-2021-AgaSal"="white", 
                               "ML-2006--" = "white",
                               "ML-2021-SyzBor"="#3498db", "ML-2021-SyzCim"="#85c1e9", "ML-2021-LabCal" = "#ec8931",
                               "BB-2006--" = "#b2babb", "BB-2021-AgaSal" = "#b2babb",
                               "Ed-2006--" = "#5f6a6a", "Ed-2021-AgaSal" = "#5f6a6a",
                               "BV-2021-AgaSal" = "#34495e",
                               "GE-2021-EucRob"="#04a896", 
                               "Tr-2021-PanSyl"="#f1c40f"))+
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dotted", color = "darkgrey") +
  xlab("NMDS1") + ylab("NMDS2") + #ylim(-.3, .55) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size=28),
        legend.position = "none") +
  annotate(geom="text", x=-0.22, y=-.32, label=paste("Stress = ", round(selec_NMDS$stress, 4), sep =""),
           col = "black", size = 4.8)
#ggsave(paste(SortiePath,date, "NMDS_TableB_Roots-Tot.pdf",  sep=""), height = 6, width = 7, plot = last_plot())


listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

tab_adonis <- adonis2(formula = Table ~ SampleSite, data = NMDS_PlotTab, method = "jaccard", nperm = 999) 
tab_adonis

# adonis2(formula = Table ~ SampleSite, data = NMDS_PlotTab, method = "jaccard", nperm = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# SampleSite   5   12.576 0.09094 7.0624  0.001 ***
#   Residual   353  125.722 0.90906                  
# Total      358  138.298 1.00000          

mod <- betadisper(Table.dist, factor(NMDS_PlotTab$SampleSite))
amod <- anova(mod)
amod
amod$`Sum Sq`[1]/(amod$`Sum Sq`[1]+amod$`Sum Sq`[2])
# Response: Distances
# Df  Sum Sq  Mean Sq F value    Pr(>F)    
# Groups      5 0.55189 0.110377  21.738 < 2.2e-16 ***
#   Residuals 353 1.79239 0.005078


#### 2.6. NMDS - Fig. S6-E  - Bray - Roots #############################################

Table <- TableA
Table <- Table[,c(9:ncol(Table))]
SampleToRemove <- c("RAM1", "RAM2", "RAM3", "RBV94_1", "Eden4")
Table <- Table[,which(!colnames(Table)%in%SampleToRemove)]
Table <- Table[,which(colnames(Table)%in%Metadata_Roots$Sample_ID)]

Table <- Table[, which(apply(Table, 2, sum)>0)]
Table <- Table[which(apply(Table, 1, sum)>0),]
Table <- as.data.frame(t(Table))
set.seed(122)
Table.dist <- vegdist(Table, method = "bray", binary = FALSE)
selec_NMDS <- metaMDS(Table.dist,k=2,trymax=100, halfchange = FALSE)

x <- selec_NMDS$points[,1]
y <- selec_NMDS$points[,2]
NMDS_PlotTab <- data.frame(x=x, y=y)

NMDS_PlotTab$SampleSite <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Roots$Site[which(Metadata_Roots$Sample_ID%in%x)]}))

NMDS_PlotTab$Host <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Roots$Species[which(Metadata_Roots$Sample_ID%in%x)]}))
NMDS_PlotTab$Year <- unlist(lapply(rownames(selec_NMDS$points), function(x){
  Metadata_Roots$Year[which(Metadata_Roots$Sample_ID%in%x)]}))
NMDS_PlotTab$SY <- paste(NMDS_PlotTab$SampleSite, NMDS_PlotTab$Year, sep ="-")
NMDS_PlotTab$SYH <- paste(NMDS_PlotTab$SampleSite, NMDS_PlotTab$Year, NMDS_PlotTab$Host, sep ="-")

ggplot(data = NMDS_PlotTab, aes(x=x, y=y)) + 
  geom_point(aes(fill = factor(SYH)), color = "black", size = 5, alpha = .90, shape = 21) + # shape = factor(Year)
  theme_classic() + 
  scale_fill_manual(values = c("ML-2021-AgaSal"="white", 
                               "ML-2006--" = "white",
                               "ML-2021-SyzBor"="#3498db", "ML-2021-SyzCim"="#85c1e9", "ML-2021-LabCal" = "#ec8931",
                               "BB-2006--" = "#b2babb", "BB-2021-AgaSal" = "#b2babb",
                               "Ed-2006--" = "#5f6a6a", "Ed-2021-AgaSal" = "#5f6a6a",
                               "BV-2021-AgaSal" = "#34495e",
                               "GE-2021-EucRob"="#04a896", 
                               "Tr-2021-PanSyl"="#f1c40f"))+
  geom_hline(yintercept=0, linetype="dotted", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dotted", color = "darkgrey") +
  xlab("NMDS1") + ylab("NMDS2") + ylim(-.3, .30) + xlim(-.22, .6)+
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        plot.title = element_text(size=28),
        legend.position = "none") +
  annotate(geom="text", x=.5, y=-.30, label=paste("Stress = ", round(selec_NMDS$stress, 4), sep =""),
           col = "black", size = 4.8)
#ggsave(paste(SortiePath,date, "NMDS_TableA_Roots-Tot.pdf",  sep=""), height = 6, width = 7, plot = last_plot())


listA <- rownames(Table)
listB <- rownames(NMDS_PlotTab)
CompareTwoLists("listA", "listB")

tab_adonis <- adonis2(formula = Table ~ SampleSite, data = NMDS_PlotTab, method = "bray", nperm = 999) 
tab_adonis

# adonis2(formula = Table ~ SampleSite, data = NMDS_PlotTab, method = "bray", nperm = 999)
# Df SumOfSqs      R2      F Pr(>F)    
# SampleSite   5   11.324 0.07091 5.4496  0.001 ***
#   Residual   357  148.372 0.92909                  
# Total      362  159.697 1.00000           

mod <- betadisper(Table.dist, factor(NMDS_PlotTab$SampleSite))
amod <- anova(mod)
amod
amod$`Sum Sq`[1]/(amod$`Sum Sq`[1]+amod$`Sum Sq`[2])
# Response: Distances
# Df Sum Sq  Mean Sq F value    Pr(>F)    
# Groups      5 1.0115 0.202300   16.53 1.101e-14 ***
#   Residuals 357 4.3691 0.012238

#### 3. Differencial analyses (DA) #########################################
#### 3.1. ANCOM-BC ###########################################################
# https://github.com/FrederickHuangLin/ANCOMBC

ps86f <- readRDS("Path_to/23_09_21_PhyloSeq-Table_A-ANCOM-BC.rds")
pstul <- readRDS("Path_to/23_09_21_PhyloSeq-Table_D-ANCOM-BC.rds")


out86f <- ancombc(phyloseq = ps86f, formula = "Tree_species2", p_adj_method = "holm",
               zero_cut = .90, lib_cut = 0, group = "Tree_species2", struc_zero = TRUE, neg_lb = TRUE, conserve = TRUE,
               alpha = 0.05, global = TRUE)
outtul <- ancombc(phyloseq = pstul, formula = "Tree_species2", p_adj_method = "holm",
                  zero_cut = .90, lib_cut = 0, group = "Tree_species2", struc_zero = TRUE, neg_lb = TRUE, conserve = TRUE,
                  alpha = 0.05, global = TRUE)

ancom_res_df86f <- data.frame(
  Species = row.names(out86f$res$beta),
  beta = unlist(out86f$res$beta),
  se = unlist(out86f$res$se),
  W = unlist(out86f$res$W),
  p_val = unlist(out86f$res$p_val),
  q_val = unlist(out86f$res$q_val),
  diff_abn = unlist(out86f$res$diff_abn))

holm_ancom <- ancom_res_df86f %>%
  dplyr::filter(q_val < 0.05)
holm_ancom

ancom_res_dftul <- data.frame(
  Species = row.names(outtul$res$beta),
  beta = unlist(outtul$res$beta),
  se = unlist(outtul$res$se),
  W = unlist(outtul$res$W),
  p_val = unlist(outtul$res$p_val),
  q_val = unlist(outtul$res$q_val),
  diff_abn = unlist(outtul$res$diff_abn))

holm_ancom <- ancom_res_dftul %>%
  dplyr::filter(q_val < 0.05)
holm_ancom

###
df_fig1.86f = data.frame(Species = ancom_res_df86f$Species,
                         beta = ancom_res_df86f$beta * ancom_res_df86f$diff_abn)
df_fig2.86f = data.frame(Species = ancom_res_df86f$Species,
                         SE = ancom_res_df86f$se * ancom_res_df86f$diff_abn)

df_fig1.tul = data.frame(Species = ancom_res_dftul$Species,
                         beta = ancom_res_dftul$beta * ancom_res_dftul$diff_abn)
df_fig2.tul = data.frame(Species = ancom_res_dftul$Species,
                         SE = ancom_res_dftul$se * ancom_res_dftul$diff_abn)

# Combine two primers
df_fig1 <- rbind(df_fig1.86f, df_fig1.tul)
df_fig2 <- rbind(df_fig2.86f, df_fig2.tul)


df_fig = df_fig1 %>% left_join(df_fig2, by = "Species") %>%
  transmute(Species, beta, SE) %>%
  filter(beta != 0) %>% arrange(beta) %>%
  mutate(group = ifelse(beta > 0, "Others", "AgaSal"))
df_fig$Species = factor(df_fig$Species, levels = df_fig$Species)

# Inverse the barchart: AgaSal > 0
df_fig$beta <- -df_fig$beta
df_fig$OTU_ID <- sapply(df_fig$Species, function(x){unlist(str_split(x, "-"))[1]})
df_fig$OTU_ID = factor(df_fig$OTU_ID, levels = df_fig$OTU_ID)


ggplot(data = df_fig, 
            aes(x = OTU_ID, y = beta)) + 
  geom_hline(yintercept=c(-0.2,.2), linetype = "dashed", col = "grey", linewidth = .5)+
  geom_bar(aes(fill = group), stat = "identity", width = 0.7,
           position = position_dodge(width = 0.4)) +
  scale_fill_manual(values = c("AgaSal"= "#bc3250", "Others"="#04a896"), name = "Tree species")+
  geom_errorbar(aes(ymin = beta - SE, ymax = beta + SE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "'AgaSal' vs. 'other tree species' effect size \n(log fold change)") + 
  theme_classic() + 
  theme(legend.position = c(.2, .15),
        axis.text = element_text(size = 15, angle = 45, hjust = 1.1),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(.5, "cm"),
        legend.background = element_rect(size = .25, linetype = "solid", colour = "black"))
#ggsave(paste(SortiePath, date, "ANCOM-BC-OUT.pdf",  sep=""), width = 8, height = 8, plot = last_plot())


#### 3.2. ALDEx2 ###########################################################

# require the homemade function "getTaxoForAllOTUs"
# the package ggrepel is also needed for plot

for (tbl in c("TableA", "TableD")){
  Table <- get(tbl)
  Table$OTU_ID <- paste(Table$OTU_ID, unlist(lapply(Table$OTU_ID, getTaxoForAllOTUs)), sep = "-")
  
  rownames(Table) <- Table$OTU_ID
  
  Table <- Table[,c(9:ncol(Table))]
  Table <- Table[,which(colnames(Table)%in%Metadata_Bark$Ech_ID)]
  
  Conds <- sapply(colnames(Table), function(x){ifelse(Metadata_Bark$Tree_species[which(Metadata_Bark$Ech_ID==x)]=="AgaSal", "AgaSal", "Other")})

  print(paste("Processing:", tbl))
  x <- aldex.clr(round(Table*10000), Conds) # transformation in centered log ratio = clr

  x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE) # Stat test x 2 (Wilcox+Welch's)
  x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
  aldex_out <- data.frame(x_tt, x_effect)
  par(mfrow = c(1, 2))
  aldex.plot(aldex_out,
             type = "MA",
             test = "welch",
             xlab = "Log-ratio abundance",
             ylab = "Difference",
             cutoff = 0.05)
  aldex.plot(aldex_out,
             type = "MW",
             test = "welch",
             xlab = "Dispersion",
             ylab = "Difference",
             cutoff = 0.05)
  par(mfrow = c(1, 1))

  aldex_out %>%
    rownames_to_column(var = "OTUs") %>%
    # here we choose the wilcoxon output rather than t-test output
    filter(wi.eBH <= 0.1)  %>%
    dplyr::select(OTUs, we.eBH, wi.eBH, effect, overlap) %>%
    knitr::kable()
  # Rq : wi.eBH = Expected Benjamini-Hochberg corrected P value of Wilcoxon test
  
  aldex_out$wi.eBH_neg <- -log10(aldex_out$wi.eBH)
  
  groupCol <- c()
  for (k in 1:nrow(aldex_out)){
    fl <- NULL
    if (aldex_out$wi.eBH_neg[k]>1 & aldex_out$diff.btw[k]<0){fl <- "AgaSal"}
    else if (aldex_out$wi.eBH_neg[k]>1 & aldex_out$diff.btw[k]>0){fl <- "Others"}
    else(fl <- "None")
    groupCol <- c(groupCol, fl)
  }
  aldex_out$group <- groupCol
  aldex_out$OTU_ID <- sapply(rownames(aldex_out), function(x){unlist(str_split(x, "-"))[1]})
  aldex_out$OTU_ID[which(aldex_out$group=="None")] <- NA
  
  if (tbl == "TableA"){aldex_out$OTU_ID[which(rownames(aldex_out) =="OTU4-Sebacinales")] <- "OTU4"}
  assign(paste("aldex_out_", tbl, sep =""), aldex_out)
  
}

# Combine the two primers
aldex_out <- rbind(aldex_out_TableA, aldex_out_TableD)

ggplot(aldex_out, aes(x = diff.btw, y = wi.eBH_neg)) +
  geom_vline(xintercept=0, linetype = "dashed", col = "grey", linewidth = .3)+
  geom_vline(xintercept=c(-2, 2), linetype = "dashed", col = "black", linewidth = .5)+
  geom_hline(yintercept=c(1,2,3), linetype = "dashed", col = "grey", linewidth = .5)+
  geom_point(aes(x = diff.btw, y = wi.eBH_neg, col = group), size = 3, alpha = .75) +
  scale_color_manual(values=c("AgaSal" = "#bc3250", "Others" = "#04a896" , "None"="black"), name = "Tree species")+
  scale_y_continuous(trans = "log10") +
  theme_classic() +
  xlab("log2 (fold change)")+ylab("-log10(corrected p-val)")+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(.5, "cm"),
        legend.position = c(.80, .2),
        legend.background = element_rect(size = .25, linetype = "solid", colour = "black"))+
  geom_text_repel(aes(label = OTU_ID), box.padding = 0.6, min.segment.length = Inf, size = 5)+
  annotate(geom = "rect", fill = "white", xmin = 9.4, xmax = 14.4, ymin = .9, ymax = 1.1)+
  annotate(geom="text", x = 12, y = 1.01, label = "p-value = 0.1", col = "darkgrey", size = 5)+
  annotate(geom = "rect", fill = "white", xmin = 9, xmax = 14.4, ymin = 1.9, ymax = 2.1)+
  annotate(geom="text", x = 11.8, y = 2.02, label = "p-value = 0.01", col = "darkgrey", size = 5)+
  annotate(geom = "rect", fill = "white", xmin = 8.4, xmax = 14.4, ymin = 2.9, ymax = 3.1)+
  annotate(geom="text", x = 11.6, y = 3.02, label = "p-value = 0.001", col = "darkgrey", size = 5)
#ggsave(paste(SortiePath, date, "ALDEx2-OUT.pdf",  sep=""), width = 8, height = 8, plot = last_plot())



#### 4. Diversity profiles #################################################
#### 4.1. Fig. 3-B  ###################################

Table <- TableB

Table$Taxo <- unlist(lapply(Table$taxonomy, WhichOrder))
Table$Taxo2 <- unlist(lapply(Table$taxonomy, WhichPhyllum))
Table <- Table %>% relocate(Taxo, 1)
Table$Taxo <- paste(Table$Taxo, Table$Taxo2, sep ="-")
Table$Taxo2 <- NULL
rownames(Table) <- Table$OTU_ID
Table <- Table[,c(1, 10:ncol(Table))] # Keep only Taxo and samples

# Choose here bark or roots <== required
#KeptDataset <- which(colnames(Table)%in%Metadata_Bark$Ech_ID)
KeptDataset <- which(colnames(Table)%in%Metadata_Roots$Sample_ID)

Table <- Table[, c(1, KeptDataset)]

# If roots : remove 2006 sample <== required
Ech2006Remove <- Metadata_Roots$Sample_ID[which(Metadata_Roots$Year=="2006")]
Table <- Table[,!colnames(Table)%in%Ech2006Remove]
Table <- Table[,!colnames(Table)%in%c("RAM1", "RAM2", "RAM3")]

Table <- Table %>% group_by(Taxo) %>%
  mutate_at(vars(-group_cols()), .funs = sum)
if (sum(duplicated(Table)) > 0){Table <- as.data.frame(Table[-which(duplicated(Table)),])}
rownames(Table) <- Table$Taxo
Table <- subset(Table, select = -c(Taxo))
df <- data.frame(tax = rownames(Table),
                 val = apply(Table, 1, sum))
df <- df[c(order(df$val, decreasing = TRUE)),]

df$Taxo <- unlist(lapply(df$tax, function(x){unlist(str_split(x, "-"))[1]}))
df$Taxo2 <- unlist(lapply(df$tax, function(x){unlist(str_split(x, "-"))[2]}))
df$tax <- NULL

df$valP <- decostand(df$val, method = "total", MARGIN = 2)
df$valP <- round(df$valP, 4)*100

df <- df[c(order(df$Taxo2)),]

df2 <- data.frame(val = sum(subset(df, Taxo2 == "Ascomycota" & valP <= 3)$val),
                  Taxo = "Other Ascomycota",
                  Taxo2 = "Ascomycota",
                  valP = sum(subset(df, Taxo2 == "Ascomycota" & valP <= 3)$valP))



df3 <- data.frame(val = sum(subset(df, Taxo2 == "Basidiomycota" & valP <= 3)$val),
                  Taxo = "Other Basidiomycota",
                  Taxo2 = "Basidiomycota",
                  valP = sum(subset(df, Taxo2 == "Basidiomycota" & valP <= 3)$valP))

df4 <- data.frame(val = sum(subset(df, !Taxo2 %in% c("Ascomycota", "Basidiomycota") & valP <= 3)$val),
                  Taxo = "Others",
                  Taxo2 = "Others",
                  valP = sum(subset(df, !Taxo2 %in% c("Ascomycota", "Basidiomycota") & valP <= 3)$valP))

df <- subset(df, valP > 3)
df <- rbind(df, df2, df3, df4)
df <- df[c(order(df$Taxo2)),]



df$Taxo[which(df$Taxo%in%"unidentified" & df$Taxo2%in%"Ascomycota")] <- "Unidentified Asco."
df$Taxo[which(df$Taxo%in%"unidentified" & df$Taxo2%in%"Basidiomycota")] <- "Unidentified Basi."

df$Taxo2[which(df$Taxo2%in%"unidentified")] <- "Unidentified"
df$Taxo[which(df$Taxo%in%"unidentified")] <- "Unidentified"


df$Ordre <- seq(1, nrow(df))

df2 <- data.frame(Taxo2 = c("Ascomycota", "Basidiomycota", "Unidentified"),
                  valP = c(sum(subset(df, Taxo2=="Ascomycota")$valP), 
                           sum(subset(df, Taxo2=="Basidiomycota")$valP),
                           sum(subset(df, Taxo2=="unidentified")$valP)),
                  Ordre = c(1:3))
ggplot(data = df) + 
  geom_bar(aes(x=3, y=valP, fill=reorder(Taxo, -Ordre)), width = 1, color="black", stat = "identity") +
  geom_bar(aes(x=2.4, y=valP, fill=reorder(Taxo2, -Ordre)), width = .6, color="black", stat = "identity") +
  coord_polar("y", start = 0, direction = 1) +
  theme_minimal() +
  xlim(0.5, 4.2)+
  scale_fill_manual(values = c("Ascomycota" = "#FFDDF8", 
                               "Helotiales" = "#8f2d56", "Hypocreales" = "#C34F80", "Capnodiales" = "#E55F98", 
                               "Unidentified Asco." = "#FBACD1", "Other Ascomycota" = "#FBD7E7", 
                               "Basidiomycota" = "#d6eaf8",
                               "Cantharellales"="#218380", "Sebacinales" = "#5EBDBB", 
                               "Unidentified Basi." = "#A2DDD3",
                               "Agaricales" = "#06d6a0", "Hymenochaetales" = "#9bc9e3",
                               "Auriculariales" = "#02c39a", "Other Basidiomycota" = "#BAF7E8",
                               "Unidentified" = "#bdc3c7",
                               "Others" = "#ecf0f1"))+
  geom_text(aes(label = paste0(round(valP, 1), "%"), x=4, y = valP), position = position_stack(vjust = .5), size = 6) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        panel.grid=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks = element_blank())


#### 4.2. Fig. 4-A  ###############################################

# reshape2 package is needed for this section

listeOTUs <- c("OTU1", "OTU3", "OTU32", "OTU20", "OTU22", "OTU79", "OTU7", "OTU4",
               "OTU351", "OTU240", "OTU89", "OTU53", "OTU81", "OTU70", 
               "OTU47", "OTU2") 
Table <- TableB
Table <- Table[Table$OTU_ID%in%listeOTUs,] # <== == == == 
rownames(Table) <- paste(Table$OTU_ID, Table$taxonomy, sep =":")

Table <- Table[,c(9:ncol(Table))]

SampleToRemove <- c("RAM1", "RAM2", "RAM3", "RBV53bp1", "RBV53bp2", "RBV53bp3", "RBV53bp4", "RBV53bp5")
Table <- Table[,which(!colnames(Table)%in%SampleToRemove)]
# Retirer aussi les echantillons de 2006
Ech2006Remove <- Metadata_Roots$Sample_ID[which(Metadata_Roots$Year=="2006")]
Table <- Table[,!colnames(Table)%in%Ech2006Remove]
# Adjust order to avoid OTUs detected in RAM (e.g. TUL39) or in 2006 only

Table <- Table[which(apply(Table, 1, sum)>0),]
OrderTable <- apply(Table, 1, sum)
Table <- Table[order(OrderTable, decreasing = TRUE),]


###
if (nrow(Table) > 30){Table <- Table[c(1:30),]} 
tTable <- as.data.frame(t(Table))
tTable$Cat <- sapply(rownames(tTable), function(x){
  if (x %in% Metadata_Roots$Sample_ID){Metadata_Roots$Site[which(Metadata_Roots$Sample_ID==x)]}
  else if (x %in% Metadata_Bark$Ech_ID){
    if (x %in% Metadata_Bark$Ech_ID[Metadata_Bark$Tree_species=="AgaSal"]){"AgaSal"}
    else ("Other")}
  else (print(x))
  
})
tTable <- tTable %>% relocate(Cat, 1)
tTable$Tot <- 1

#Merge
tTable <- tTable %>% group_by(Cat) %>%
  mutate_at(vars(-group_cols()),.funs = sum)
if (sum(duplicated(tTable)) > 0){tTable <- as.data.frame(tTable[-which(duplicated(tTable)),])}
rownames(tTable) <- tTable$Cat
tTable$Cat <- NULL
rownames(tTable) <- paste(rownames(tTable), " (", tTable$Tot, ")", sep ="")
TotalRoots <- as.data.frame(apply(tTable[1:6,], 2, sum)) # /!\ Peut changer selon les categories choisies
colnames(TotalRoots) <- "Total"
tTable <- rbind(tTable, t(TotalRoots))
for (k in 1:nrow(tTable)){
  for (i in 1:ncol(tTable)){
    tTable[k, i] <- tTable[k, i]/tTable$Tot[k]
  }
}
tTable$Tot <- NULL
Table <- as.data.frame(t(tTable))
Table <- round(Table, 2)*100 # round values
OTUProv <- unlist(lapply(rownames(Table), function(x){unlist(str_split(x, ":"))[1]}))
TaxoProv <- sapply(OTUProv, getTaxoForAllOTUs)
Table$OTU_IDTaxo <- paste(OTUProv, TaxoProv, sep = ":")


# Construct graph

# 1. OMF families
TableOMF <- Table[which(OTUProv %in% TaxoTable$OTU_ID[which(TaxoTable$OMF==1)]),]
dfOMF <- melt(TableOMF, id = "OTU_IDTaxo")
dfOMF$OTU_IDTaxo <- factor(dfOMF$OTU_IDTaxo, 
                           levels=(rev(TableOMF$OTU_IDTaxo)))
dfOMF$variable <- factor(dfOMF$variable,
                         levels = c("ML (141)", "Ed (46)", "BV (51)", "BB (55)", "GE (28)", "Tr (4)", "Total",
                                    "AgaSal (97)", "Other (279)"))
colVal <- unlist(lapply(dfOMF$value, function(x){ifelse(x>0 & x < 80, "#212f3c", "white")}))
plot1 <- ggplot(dfOMF, aes(x = variable, y = OTU_IDTaxo, fill = value)) +
  geom_tile(color = "white", linetype = 1, lwd = 1.5) + theme_classic()+
  scale_fill_gradient(low = "white", high = "#04a896", name = "Frequency")+
  geom_text(aes(variable, OTU_IDTaxo, label = value), color = colVal, size = 3)+
  theme(axis.text.x = element_text(angle =45, vjust = .6),
        axis.title = element_text(size = 0))
plot1

# 2. non-OMF families

TableOth <- Table[which(!OTUProv %in% TaxoTable$OTU_ID[which(TaxoTable$OMF==1)]),]
TableOth$OTU_IDTaxo[TableOth$OTU_IDTaxo=="OTU79:Hymenochaetales_fam_Incertae_sedis"] <- "OTU79:Hymenochaetales"
dfOth <- melt(TableOth, id = "OTU_IDTaxo")
dfOth$OTU_IDTaxo <- factor(dfOth$OTU_IDTaxo, 
                           levels=(rev(TableOth$OTU_IDTaxo)))
dfOth$variable <- factor(dfOth$variable,
                         levels = c("ML (141)", "Ed (46)", "BV (51)", "BB (55)", "GE (28)", "Tr (4)", "Total",
                                    "AgaSal (97)", "Other (279)"))

colVal <- unlist(lapply(dfOth$value, function(x){ifelse(x>0 & x < 80, "#212f3c", "white")}))

plot2 <- ggplot(dfOth, aes(x = variable, y = OTU_IDTaxo, fill = value)) +
  geom_tile(color = "white", linetype = 1, lwd = 1.5) + theme_classic()+
  scale_fill_gradient(low = "white", high = "#bc3250", name = "Frequency") +
  geom_text(aes(variable, OTU_IDTaxo, label = value), color = colVal, size = 3)+
  theme(axis.text.x = element_text(angle =45, vjust = .6),
        axis.title = element_text(size = 0))
plot2
plot_grid(plot1, plot2, nrow = 2, rel_heights = c(1.1, 1), align = "hv") # for listOTUs
#ggsave(paste(SortiePath,date, "CompoTableau-OTUSelects.pdf",  sep=""), plot = last_plot())

#### 4.3. Fig. 4-C  ##############################################

# library "car" is needed for this section (Anova)

dff <- data.frame(Abd = NA, Size = NA, OTU = NA)
ListOTU <- c("OTU1", "OTU3", "OTU4", "OTU7")
for (otu in ListOTU){
  tax <- TaxoTable$Taxonomy[TaxoTable$OTU_ID==otu]
  if (!str_detect(tax, "Tulasne")){Table <- TableA}
  else
    {Table <- TableD}
  
  rownames(Table) <- Table$OTU_ID
  SampleToRemove <- c("RAM1", "RAM2", "RAM3", "RBV53bp1", "RBV53bp2", "RBV53bp3", "RBV53bp4", "RBV53bp5") 
  Table <- Table[,which(!colnames(Table)%in%SampleToRemove)]
  Ech2006Remove <- Metadata_Roots$Sample_ID[which(Metadata_Roots$Year=="2006")]
  Table <- Table[,!colnames(Table)%in%Ech2006Remove]
  EchBarkRemove <- colnames(Table)[!colnames(Table)%in%Metadata_Roots$Sample_ID]
  Table <- Table[,!colnames(Table)%in%EchBarkRemove]
  Table <- Table[,c(9:ncol(Table))]
  
  
  L <- sapply(colnames(Table), function(x){
    Metadata_Roots$Longueur[which(Metadata_Roots$Sample_ID==x)]
  })
  L <- as.numeric(L)
  
  df <- data.frame(Abd = unlist(c(Table[otu,])),
                   Size = L)
  df <- df[!is.na(df$Size),]
  
  if(Anova(lm(df, formula = Abd ~ Size))["Size", "Pr(>F)"]< 0.05){
    df$OTU <- otu
    dff <- rbind(dff, df)
    ListOTU <- c(ListOTU, otu)

  }else{ 
    df$OTU <- otu
    dff <- rbind(dff, df)
    ListOTU <- c(ListOTU, otu)
  }
  
}

ListOTU 
dff <- dff[2:nrow(dff),]
ggplot(dff, aes(x = Size, y = Abd, group = OTU, col = OTU, fill = OTU)) + theme_classic() +
  geom_smooth(formula = y~x, method="lm", se=TRUE, level=0.95) +
  scale_color_manual(values = c("OTU1" = "#00808f", "OTU3" = "#04a896", "OTU4"="#02c39a", "OTU7" = "#bc3250"))+
  scale_fill_manual(values = c("OTU1" = "#00808f", "OTU3" = "#04a896", "OTU4"="#02c39a", "OTU7" = "#bc3250")) +
  annotate(geom = "text", x = c(25, 33, 18, 7), y = c(.75, .3, .24, .13), 
           col = c("#00808f","#04a896","#02c39a", "#bc3250"),
           label = c("OTU1***", "OTU3***", "OTU4***", "OTU7 (n.s)"), size = 5) +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        legend.position = "none") + xlab("Leaf length") + ylab("Abundance")

#ggsave(paste(SortiePath, date, "OTUs-EvolWithLeaf.pdf",  sep=""), width = 8, height = 8, plot = last_plot())


#### 4.4. Fig. S13  ##################################

NonAgaSalRoots <- c(paste0("RBV26_", 1:5), paste0("RBV27_", 1:3),
                    paste0("RBV106_", 1:3), paste0("RBV107_", 1:3),
                    paste0("RBV28_", 1:3), 
                    paste0("RBV29_", 1:3), paste0("RBV30_", 1:2),
                    paste0("RBV71_", 1:3), paste0("RBV72_", 1:3), paste0("RBV73_", 1:3), paste0("RBV74_", 1:3), paste0("RBV75_", 1:3),
                    paste0("RBV76_", 1:3), paste0("RBV77_", 1:3), paste0("RBV78_", 1:3), paste0("RBV79_", 1:3), paste0("RBV80_", 1:3),
                    "TR1_1BV", "TR1_2BV", "TR2_2BV", "TRBV1_2bis")

Table <- TableD # <== choose here between TableA (ITS86-f/ITS4) or TableD (5.8S-OF/ITS4-Tul)
rownames(Table) <- Table$OTU_ID
Table <- Table[,which(colnames(Table)%in%NonAgaSalRoots)]
Table <- Table[  which(apply(Table, MARGIN = 1, FUN = sum)>0),  ]
Table$taxonomy <- sapply(rownames(Table), getTaxoForAllOTUs)
Table$taxonomy[which(Table$taxonomy=="Helotiaceae")] <- "Helotiales"

for (k in 1:nrow(Table)){
  otu <- rownames(Table)[k]
  taxo <- Table$taxonomy[k]
  if (all(otu%in%c("OTU1", "OTU3", "OTU4", "OTU7"))==FALSE){next}
  else( Table$taxonomy[k] <- paste(taxo, otu, sep = ":")  )
}


Table$taxonomy[which(!Table$taxonomy%in%c("Helotiaceae", "Helotiales", "Helotiales:OTU7",
                                          "Sebacinales",
                                          "Sebacinales:OTU4", "Tulasnellaceae",
                                          "Tulasnellaceae:OTU1", "Tulasnellaceae:OTU3"))] <- "Others"

Table$Sum <- NULL
Table <- Table %>% group_by(taxonomy) %>%
  mutate_at(vars(-group_cols()),.funs = sum)
if (sum(duplicated(Table)) > 0){Table <- as.data.frame(Table[-which(duplicated(Table)),])}
mTable <- melt(Table, id = "taxonomy")
mTable$variable <- factor(mTable$variable, levels = NonAgaSalRoots[NonAgaSalRoots%in%colnames(Table)]) # Reorder


ggplot(data=mTable, aes(x=variable, y=value, fill=as.factor(taxonomy))) +
  geom_bar(stat="identity", col = "white")+
  theme_classic()+
  ylab("Cumulative percentage \n of sequences (%)") + #xlab("Orchid individual") +
  #facet_grid(. ~ Category_f, scale = "free_x", space = "free_x")
  scale_fill_manual(values = c("Sebacinales" = "#06d6a0",
                               "Sebacinales:OTU4" = "#04a896",
                               "Helotiales" = "#DCA0A4",
                               "Helotiales:OTU7" = "#8f2d56",
                               "Tulasnellaceae" = "#0f2f45",
                               "Tulasnellaceae:OTU1" = "#fca311",
                               "Tulasnellaceae:OTU3" = "#9bc9e3", "Others" = "#e5e5e5"))+
  theme(axis.text.x = element_text(size = 15, angle = 90, vjust = .5),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        plot.title = element_text(size=0),
        legend.title = element_text(size = 0),
        legend.background = element_rect(size = .25, linetype = "solid", colour = "black"),
        legend.text = element_text(size = 16),
        legend.spacing.x = unit(.5, "cm"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),)

#ggsave(paste(SortiePath, date, "BarPlotCompo_BV_Roots_on_Non_AgaSalTrees.pdf",  sep=""), width = 10, height = 8, plot = last_plot())
#ggsave(paste(SortiePath, date, "BarPlotCompo_BV_Roots_on_Non_AgaSalTrees--TUL.pdf",  sep=""), width = 10.2, height = 8, plot = last_plot())


#### 5. Spatial analyses ###################################################
#### 5.1. Fig. S10-A - Mantel test - ##########################################

SpatialDist <- readRDS("Path_to/24_05_06_SampleCoordinates-Mantel.rds")
Table <- TableA

Table <- Table[,which(colnames(Table)%in%Metadata_Bark$Ech_ID)]
Table <- Table[,which(!colnames(Table)%in%paste("E", seq(602,616), sep =""))] # samples without coordinates
Table <- Table[, which(apply(Table, 2, sum)>0)]
Table.dist.euc <- data.frame("SampleID" = colnames(Table),
                             "TreeID" = unlist(lapply(colnames(Table), function(x){Metadata_Bark$Individual_ID[which(Metadata_Bark$Ech_ID==x)]})))
Table.dist.euc$CoordX <- unlist(lapply(Table.dist.euc$TreeID, function(x){SpatialDist$xcoord[which(SpatialDist$Individual_ID==x)]}))
Table.dist.euc$CoordY <- unlist(lapply(Table.dist.euc$TreeID, function(x){SpatialDist$ycoord[which(SpatialDist$Individual_ID==x)]}))
rownames(Table.dist.euc) <- Table.dist.euc$SampleID
Table.dist.euc <- subset(Table.dist.euc, select = c(CoordX, CoordY))

Table.dist.euc <- as.matrix(dist(Table.dist.euc, method = "euclidean"))
Table.dist.euc <- Table.dist.euc[which(colnames(Table.dist.euc)%in%colnames(Table)),
                                 which(colnames(Table.dist.euc)%in%colnames(Table))]

mantel <- vegan::mantel 
mantel.test <- mantel(xdis = Table.dist.euc, ydis = Table.dist.ecol, permutations=999,  method = "spearman")
mantel.correl <- mantel.correlog(D.eco = Table.dist.ecol, D.geo=Table.dist.euc, 
                                 cutoff=TRUE, r.type="spearman", nperm=999, mult="bonferroni", 
                                 progressive=TRUE)
Mantel.Res <- as.data.frame(mantel.correl$mantel.res)
Mantel.Res <- Mantel.Res[which(is.na(Mantel.Res$Mantel.cor)==FALSE),]
Mantel.Res <- Mantel.Res[which(is.na(Mantel.Res$`Pr(corrected)`)==FALSE),]
Mantel.Res$Signif <- unlist(lapply(Mantel.Res$`Pr(corrected)`, function(x){ifelse(x<=0.05, 21,1)}))

ggplot(data = Mantel.Res, aes(x = class.index, y = Mantel.cor), color = "#218380") +
  geom_line(col = "#218380") +
  geom_point(aes(fill = factor(Signif)), shape = 21, colour = "black", size = 5) +
  scale_fill_manual(values = c("1"="white", "21"="#218380"))+
  theme_classic() +
  scale_x_continuous(breaks = seq(0, max(Mantel.Res$class.index)+50, 50))+
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  xlab("Class distance index (m)") +
  ylab("Mantel correlation index") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=28),
        legend.background = element_rect(size = .25, linetype = "solid", colour = "black"),
        legend.text = element_text(size = 20),
        legend.spacing.x = unit(.5, "cm"),
        legend.title = element_text(size = 0),
        legend.position = "none")
#ggsave(paste(SortiePath, date, "_Correl_TableA.pdf",  sep=""), plot = last_plot())

#### 5.2. Fig. S10-B VarPart analysis ###############################################

SpatialDist <- readRDS("Path_to/24_05_06_SampleCoordinates-Mantel.rds")
Table <- TableA

Table <- Table[,c(9:ncol(Table))]
Table <- Table[,which(colnames(Table)%in%Metadata_Bark$Ech_ID)]
Table <- Table[,which(!colnames(Table)%in%paste("E", seq(602,616), sep =""))] # samples with no coordinate
Table <- Table[which(apply(Table, 1, sum)>0),]
Table <- Table[,which(apply(Table, 2, sum)>0)]

Table.env <- subset(Metadata_Bark, select = c(Ech_ID, Tree_species, Plot, Zone)) # Individual redundant with some species


Table.env <- subset(Table.env, Ech_ID%in%colnames(Table))
rownames(Table.env) <- Table.env$Ech_ID
Table.env <- subset(Table.env, select = -c(Ech_ID))
Table <- as.data.frame(t(Table))

rda.var1 <- rda(Table~1, Table.env)
rda.var2 <- rda(Table~., Table.env)
step.forward <- ordistep(rda.var1, scope=formula(rda.var2),
                         direction = "forward", perm.max=200, pstep=999)
rda.env.var <- rda(Table ~ Tree_species + Zone + Plot , data=Table.env)

vif.cca(rda.env.var)
alias(rda.env.var) # no aliased terms)


## ## ## PCNM ## ## ##
Table.dist.euc <- data.frame("SampleID" = rownames(Table),
                             "TreeID" = unlist(lapply(rownames(Table), function(x){Metadata_Bark$Individual_ID[which(Metadata_Bark$Ech_ID==x)]})))
Table.dist.euc$CoordX <- unlist(lapply(Table.dist.euc$TreeID, function(x){SpatialDist$xcoord[which(SpatialDist$Individual_ID==x)]}))
Table.dist.euc$CoordY <- unlist(lapply(Table.dist.euc$TreeID, function(x){SpatialDist$ycoord[which(SpatialDist$Individual_ID==x)]}))
rownames(Table.dist.euc) <- Table.dist.euc$SampleID
Table.dist.euc <- subset(Table.dist.euc, select = c(CoordX, CoordY))

Table.dist.euc <- as.matrix(dist(Table.dist.euc, method = "euclidean"))
Table.dist.euc <- Table.dist.euc[which(colnames(Table.dist.euc)%in%rownames(Table)),
                                 which(colnames(Table.dist.euc)%in%rownames(Table))]

Table.dist.euc <- as.dist(Table.dist.euc)

Table.dist.euc <- as.matrix(Table.dist.euc)

Table.pcnm <- pcnm(Table.dist.euc)
Table.pcnm$values
barplot(Table.pcnm$values)
Table.pcnm.vec <- Table.pcnm$vectors 

# Modified from Blanchet, Legendre, & Borcard, 2008 (Selection of PCNM eigenvectors)
# The adespatial package is needed here

MEM_model <- "Positive"
MEM.select <- NULL
R2adj <- NULL
fsel <- NULL

if (MEM_model != "all") { # We consider only positively or negatively autocorrelated MEM
  if (anova.cca(rda(Table, Table.pcnm.vec), permutations = 9999)$Pr[1] <= 0.05) { # Global adjusted R-squared of the model
    R2adj <- RsquareAdj(rda(Table, Table.pcnm.vec))$adj.r.squared
    # FWD with two stopping criteria
    fsel <- forward.sel(Table, Table.pcnm.vec, adjR2thresh = R2adj, nperm = 999) # We order the selected MEM by decreasing eigenvalue 
    sorted_sel <- sort(fsel$order)
    # Object containing the selected MEM
    MEM.select <- as.data.frame(Table.pcnm.vec)[, c(sorted_sel)]
  } else print("No significant spatial autocorrelation was detected in the response")
  }else{ # We consider both positively and negatively autocorrelated predictors 
  # List to save the positively and negatively autocorrelated MEM separately 
  mem.sign <- vector("list", 2)
  signif <- c("FALSE", "FALSE")
  # We select the positive and negative MEM separately after testing the global
  # significance of both models at a corrected threshold value of null hypothesis # rejection (Sidak correction)
  for (i in 1:2) {
    if(i== 1){ # Positive MEM
      mem <- Table.pcnm.vec[, which(Table.pcnm$values > 0)]
    } else { # Negative MEM
      mem <- Table.pcnm.vec[, which(Table.pcnm$values < 0)]
      #mem <- MEM[, which(attributes(MEM)$values < 0)]
    }
    # Global test of significance with the Sidak correction for multiple tests
    if (anova.cca(rda(Table, mem), permutations = 9999)$Pr[1] <= (1-(1-0.05)^0.5)) { # Global adjusted R-squared of the model
      R2adj <- RsquareAdj(rda(Table, mem))$adj.r.squared
      # FWD with two stopping criteria
      fsel <- forward.sel(Table, mem, adjR2thresh = R2adj, nperm = 999) # We order the selected MEM by decreasing eigenvalue 
      sorted_sel <- sort(fsel$order)
      # We save the selection of MEM
      mem.sign[[i]] <- as.data.frame(mem)[, c(sorted_sel)]
      signif[i] <- "TRUE"
    }
  }
  # MEM.select will contain both positive and negative MEM, only positive or only # negative MEM, depending on the significance of the global tests.
  if (length(which(signif == "FALSE")) != 2) {
    if (length(which(signif == "TRUE")) == 2) {
      MEM.select <- cbind(mem.sign[[1]], mem.sign[[2]]) } else if (signif[1] == "TRUE") {
        MEM.select <- mem.sign[[1]]
      } else MEM.select <- mem.sign[[2]]
  } else print("No significant spatial autocorrelation was detected in the response") }

#Table.pcnm.vec # old
#MEM.select # new (after selection) 

# -- -- -- 
spe.part <- varpart(Table, as.factor(Table.env[,1]),as.factor(Table.env[,2]),as.factor(Table.env[,3]), MEM.select) # Check the choosen column ! <== == == == == == == 
spe.part
plot(spe.part, digits=2) 
spe.part$part$indfract$Adj.R.square[1] # 
rda.pcnm <- rda (Table, MEM.select) #[b+c]
rda.species <- rda (Table, as.factor(Table.env[,1])) # [a+b]
rda.pcnm.species <- rda (X = Table, Y = MEM.select, Z = as.factor(Table.env[,1])) # [b]
rda.species.pcnm <- rda (X = Table, Y = as.factor(Table.env[,1]), Z = MEM.select) # [a]
anova(rda.pcnm)
anova(rda.species)
anova(rda.pcnm.species)
anova(rda.species.pcnm)
## ## ## ## ## ##


#### 6. Rarefaction curves ################################################

Table <- TableB # (or TableC for Fig. S4-C)

# Need to change some samples' name
colnames(Table)[which(colnames(Table)=="RBV53bp1")] <- "RBV53b(P)_1"
colnames(Table)[which(colnames(Table)=="RBV53bp2")] <- "RBV53b(P)_2"
colnames(Table)[which(colnames(Table)=="RBV53bp3")] <- "RBV53b(P)_3"
colnames(Table)[which(colnames(Table)=="RBV53bp4")] <- "RBV53b(P)_4"
colnames(Table)[which(colnames(Table)=="RBV53bp5")] <- "RBV53b(P)_5"
colnames(Table)[which(colnames(Table)=="RBV94_1.1")] <- "RBV94_1"

# Provide a data.frame with samples in row and species in column.
Table <- Table[,c(9:ncol(Table))]

rarefy_for_sample_sd_option <- function(df, n, g){ 
  
  r_curve <- data.frame(it = NA, val = NA, type = NA, cat = NA)
  for (k in 1:nrow(df)){ 
    print(paste("Take", k, "/", nrow(df), "sample(s) into consideration"))
    ml <- as.numeric(n)
    for (i in 1:n){
      r_df <- as.data.frame(df[sample(1:nrow(df), k, replace = FALSE),]) 
      r_df <- r_df[,which(apply(r_df, MARGIN = 2, FUN = sum) > 0)] 
      if (is.data.frame(r_df)==FALSE){ml[i] <- 0}
      else (ml[i] <- ncol(r_df))
    }
    r_curve <- rbind(r_curve, data.frame(it = k, val = c(mean(ml),max(ml), min(ml)), type = c("m", "max", "min"), cat = g))
  }
  return(r_curve[-1,])
}

for (g in unique(Metadata_Roots$Site)){
  Ech_Grille <- Metadata_Roots$Sample_ID[which(Metadata_Roots$Site==g)]
  rTable <- Table[,which(colnames(Table)%in%Ech_Grille)]
  assign(paste("result_rarefy_for_sample",g, sep ="_"), rarefy_for_sample_sd_option(as.data.frame(t(rTable)), n = 10, g=g))
}

for (g in unique(Metadata_Bark$Tree_species)){
  if (!g %in%c("LabCal", "NuxVer", "AgaSal", "SyzCim", "SyzBor", "AntBor", "SidBor")){
    print(paste0("Did not take into consideration : ", g))
    next}
  Ech_Grille <- Metadata_Bark$Ech_ID[which(Metadata_Bark$Tree_species==g)]
  rTable <- Table[,which(colnames(Table)%in%Ech_Grille)]
  assign(paste("result_rarefy_for_sample",g, sep ="_"), rarefy_for_sample_sd_option(as.data.frame(t(rTable)), n = 10, g=g))
}

df <- rbind(result_rarefy_for_sample_BB, result_rarefy_for_sample_GE,
            result_rarefy_for_sample_BV, result_rarefy_for_sample_ML,
            result_rarefy_for_sample_Ed, result_rarefy_for_sample_Tr,
            result_rarefy_for_sample_AgaSal, result_rarefy_for_sample_AntBor,
            result_rarefy_for_sample_LabCal, result_rarefy_for_sample_NuxVer,
            result_rarefy_for_sample_SidBor, result_rarefy_for_sample_SyzBor)

Col <- wes_palette(n=5, name="Darjeeling1")
names(Col) <- c("ML", "BB", "Ed", "BV", "GE")

Col <- c(rep("#ec8931", 5), "#8f2d56", rep("#04a896", 5))
names(Col) <- c("ML", "BB", "Ed", "BV", "GE", "AgaSal", "LabCal", "NuxVer", "AntBor", "SidBor", "SyzBor")

AddLabels <- data.frame(Cat =c(), x=c(), y=c())
for (i in unique(df$cat)){
  dfi <- subset(df, cat==i & type=="m")
  yi <- max(dfi$val)
  xi <- max(dfi$it)
  AddLabels <- rbind(AddLabels, data.frame(Cat=i, x=xi, y=yi))
}
AddLabels <- subset(AddLabels, Cat != "Tr")

ggplot(subset(df, cat != "Tr"), aes(x=it, y = val, linetype = type, col = cat, size = type))+
  geom_point(alpha = 0) + theme_bw()+
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)))+
  scale_linetype_manual(values = c("m" = "solid", "max" = "dotted", "min"="dotted")) +
  scale_color_manual(values = Col)+
  scale_size_manual(values = c("m" = 1.8, "max" = .7, "min"=.7))+
  scale_y_continuous(n.breaks = 6)+
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.position = "none") +xlab("Number of samples")+ylab("Number of OTUs")+
  annotate(geom="text", x=AddLabels$x+10, y=AddLabels$y+10, label=AddLabels$Cat,
           col = "black", size = 4.8)

#ggsave(paste(SortiePath,date, "Raref_FULLDATA_Total-Species.pdf",  sep=""), width = 7, height = 6, plot = last_plot())


# One curve for all samples

Table1 <- Table[,colnames(Table)%in%Metadata_Bark$Ech_ID] # 373 samples
assign("result_rarefy_TotalBark", rarefy_for_sample_sd_option(as.data.frame(t(Table1)), n = 10, g="none"))
result_rarefy_TotalBark$cat <- "Bark"
Table2 <- Table[,colnames(Table)%in%Metadata_Roots$Sample_ID] # 368 samples
assign("result_rarefy_TotalRoots", rarefy_for_sample_sd_option(as.data.frame(t(Table2)), n = 10, g="none"))
result_rarefy_TotalRoots$cat <- "Roots"


ggplot(rbind(result_rarefy_TotalBark, result_rarefy_TotalRoots),
       aes(x=it, y = val, linetype = type, size = type, col = cat), black)+
  geom_point(alpha = 0) + theme_bw()+
  geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)))+
  scale_linetype_manual(values = c("m" = "solid", "max" = "dotted", "min"="dotted")) +
  scale_color_manual(values = c("Bark" = "#bc3250", "Roots" = "#ec8931"))+
scale_size_manual(values = c("m" = 1.8, "max" = .7, "min"=.7))+
  scale_y_continuous(n.breaks = 6)+
  theme(axis.title = element_text(size = 0),
        axis.text = element_text(size = 23),
        legend.position = "none") +xlab("Number of samples")+ylab("Number of OTUs")+
  annotate(geom="text", x = 100, y = 110, label = "Bark", size = 15, col = "#bc3250")+
  annotate(geom="text", x = 280, y = 65, label = "Roots",  size = 15, col = "#ec8931")
#ggsave(paste(SortiePath,date, "Raref_FULLDATA_Total2-Species.pdf",  sep=""), width = 7, height = 6, plot = last_plot())


#### 7. Fungal traits ######################################################

# Põlme, S., Abarenkov, K., Henrik Nilsson, R. et al. 
# FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles. Fungal Diversity 105, 1–16 (2020). 
# https://doi.org/10.1007/s13225-020-00466-2

FunTraitDB <- readRDS("Path_to/24_05_07_FunTraitDB.rds")

NSample_per_sp <- c()

for (TreeSpecies in c("All", "AgaSal", "AntBor", "LabCal", "NuxVer")){
  
  Table <- TableB
  BarkSamples <- colnames(Table)[colnames(Table)%in%Metadata_Bark$Ech_ID] # Bark only
  BarkSamples <- BarkSamples[!BarkSamples%in%paste0("E", 602:616)] # some out of sampling samples 
  
  
  if(TreeSpecies!="All"){
    BarkSamples <- BarkSamples[BarkSamples%in%Metadata_Bark$Ech_ID[Metadata_Bark$Tree_species==TreeSpecies]] # Select tree species
  }
  
  BarkSamples <- c("OTU_ID", "taxonomy", BarkSamples)
  Table <- Table[,colnames(Table)%in%BarkSamples] 
  print(paste(TreeSpecies, " ", ncol(Table)-2))
  NSample_per_sp <- c(NSample_per_sp, ncol(Table)-2)
  
  Table <- Table[apply(Table[,3:ncol(Table)], 1, sum)>0,]
  
  Table$sum <- apply(Table[,3:ncol(Table)], 1, sum)
  Table <- Table %>% relocate(sum, .after = "taxonomy")
  
  # Look for traits in data
  head(Table$taxonomy)
  ConCatLines <- function(x){
    l <- unique(x)
    if (length(l)>=1){l <- l[!is.na(l)]}
    paste(l, collapse = '; ')
  }
  
  
  
  ResFungalTraits <- data.frame()
  db <- FunTraitDB
  # These lines return FungTrait data for OTUs assigned at least at the order level 
  for (i in 1:nrow(Table)){
    TaxName <- getTaxoForAllOTUsFromGenus2(Table$taxonomy[i])
    df <- NULL
    
    # Genus
    GoodLine <- which(db$GENUS%in%TaxName)
    if (length(GoodLine)!=0){
      df <- db[GoodLine, 8:25]
      df <- cbind(OTU_ID = Table$OTU_ID[i], Identification = "Genus", Abundance = as.numeric(Table$sum[i]),
                  Taxonomy = TaxName, df)
      df <- df %>% group_by(OTU_ID) %>% mutate_at(vars(-group_cols()), .funs = ConCatLines) %>% distinct(.keep_all = TRUE)
      ResFungalTraits <- rbind(ResFungalTraits, df)
      next
    }
    
    # Family
    GoodLine <- which(db$Family%in%TaxName)
    if (length(GoodLine)!=0){
      df <- db[GoodLine, 8:25]
      df <- cbind(OTU_ID = Table$OTU_ID[i], Identification = "Family", Abundance = as.numeric(Table$sum[i]),
                  Taxonomy = TaxName, df)
      df <- df %>% group_by(OTU_ID) %>% mutate_at(vars(-group_cols()), .funs = ConCatLines) %>% distinct(.keep_all = TRUE)
      ResFungalTraits <- rbind(ResFungalTraits, df)
      next
    }
    
    # Order
    GoodLine <- which(db$Order%in%TaxName)
    if (length(GoodLine)!=0){
      df <- db[GoodLine, 8:25]
      df <- cbind(OTU_ID = Table$OTU_ID[i], Identification = "Order", Abundance = as.numeric(Table$sum[i]),
                  Taxonomy = TaxName, df)
      df <- df %>% group_by(OTU_ID) %>% mutate_at(vars(-group_cols()), .funs = ConCatLines) %>% distinct(.keep_all = TRUE)
      ResFungalTraits <- rbind(ResFungalTraits, df)
      next
    }
    
    # Class
    GoodLine <- which(db$Class%in%TaxName)
    if (length(GoodLine)!=0){
      df <- db[GoodLine, 8:25]
      df <- cbind(OTU_ID = Table$OTU_ID[i], Identification = "Class", Abundance = as.numeric(Table$sum[i]),
                  Taxonomy = TaxName, df)
      df <- df %>% group_by(OTU_ID) %>% mutate_at(vars(-group_cols()), .funs = ConCatLines) %>% distinct(.keep_all = TRUE)
      ResFungalTraits <- rbind(ResFungalTraits, df)
      next
    }
    
    # Phylum
    GoodLine <- which(db$Phylum%in%TaxName)
    if (length(GoodLine)!=0){
      df <- db[GoodLine, 8:25]
      df <- cbind(OTU_ID = Table$OTU_ID[i], Identification = "Phylum", Abundance = as.numeric(Table$sum[i]),
                  Taxonomy = TaxName, df)
      df <- df %>% group_by(OTU_ID) %>% mutate_at(vars(-group_cols()), .funs = ConCatLines) %>% distinct(.keep_all = TRUE)
      ResFungalTraits <- rbind(ResFungalTraits, df)
      next
    }
  }
  
  ResFungalTraits <- ResFungalTraits %>% mutate_if(is.character, list(~na_if(.,""))) # replace empty values by NA
  
  
  head(ResFungalTraits)

  ResFungalTraits <- ResFungalTraits[!is.na(ResFungalTraits$primary_lifestyle),] # 3/2
  # Retain only assignment with less than (or equal to) 3 values
  
  df.ft <- as.data.frame(table(ResFungalTraits$primary_lifestyle[ str_count(ResFungalTraits$primary_lifestyle, ";") <= 2 ]))
  df.ft <- rbind(df.ft,
                data.frame(Var1 = "NPE", Freq = nrow(ResFungalTraits)-sum(df.ft$Freq))) # NPE = not precise enough
  
  
  # Change Var1
  df.ft <- df.ft %>% mutate(Var1=recode_factor(Var1, 
                                               dung_saprotroph="saprotroph", 
                                               litter_saprotroph="saprotroph",
                                               wood_saprotroph="saprotroph",
                                               soil_saprotroph="saprotroph",
                                               `unspecified_saprotroph; soil_saprotroph`="saprotroph",
                                               `litter_saprotroph; unspecified`="saprotroph",
                                               `wood_saprotroph; soil_saprotroph; litter_saprotroph`="saprotroph",
                                               unspecified_saprotroph="saprotroph",
                                               `wood_saprotroph; litter_saprotroph`="saprotroph")) %>%
    group_by(Var1) %>% mutate(Freq = sum(Freq)) %>% distinct()
  
  df.ft <- df.ft %>% arrange(-Freq)
  
  
  TraitList <- c("NPE", "saprotroph", "plant_pathogen", "animal_parasite", "mycoparasite", "lichenized",
                 "lichen_parasite", "ectomycorrhizal", "root_endophyte")
  df.ft <- df.ft %>% mutate(Show = ifelse(Var1%in%TraitList, "Yes", "No")) %>% group_by(Show) %>%
    mutate(Freq = case_when(Show == "Yes" ~ Freq, Show == "No" ~ sum(Freq))) %>% 
    mutate(Var1 = case_when(Show == "Yes" ~ Var1, Show == "No" ~ "Other")) %>% distinct()
  df.ft$Freq
  
  df.ft$Freq <- df.ft$Freq/sum(df.ft$Freq)
  df.ft$TreeSpecies <- c(TreeSpecies)
  
  assign(paste0("df.ft_", TreeSpecies), df.ft)
}


df.ft.all <- rbind(df.ft_All, df.ft_AgaSal, df.ft_AntBor, df.ft_LabCal, df.ft_NuxVer)


Col.FungTrait <- c("NPE" = "#ffffff", 
                   "ectomycorrhizal" = "#bc3250", "saprotroph" = "#e05972", 
                   "lichenized" = "#fbb13c", "lichen_parasite" = "#ec8931",
                   "plant_pathogen" = "#218380", "animal_parasite"="#02c39a",
                   "mycoparasite" ="#06d6a0", 
                   "root_endophyte" = "#9bc9e3",
                   "Other" = "#d6d6d6")

df.ft.all <- subset(df.ft.all, Var1 != "NPE")

df.ft.all$TreeSpecies <- factor(df.ft.all$TreeSpecies, levels = c("AgaSal", "AntBor", "LabCal", "NuxVer", "All"))
df.ft.all$Var1 <- factor(df.ft.all$Var1, levels = rev(c("ectomycorrhizal", "saprotroph", "lichen_parasite", "lichenized",
                                                        "plant_pathogen", "animal_parasite", "mycoparasite", "root_endophyte",
                                                        "Other")))

ggplot(data = df.ft.all, aes(y = Freq, fill=Var1, x = TreeSpecies)) + 
  geom_bar(width = .7, color="black", stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = Col.FungTrait, breaks = names(Col.FungTrait), name = "Fungal guild")+
  ylab("Cumulative frequency")+ 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 24),
        axis.text.y = element_text(size = 20),
        axis.text.x=element_text(size = 20, angle = 90, hjust = .5))
#ggsave(paste(SortiePath, date, "FungTraits-Comp_PresenceAbsenceLeg.pdf",  sep=""), width = 9, height = 7, plot = last_plot())

