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

## Packages ####

library(dplyr)
library(reshape2)
library(stringr)

library(vegan) 

library(egg)
library(randomForest)


## Load raw tables and initial calculations ####

df <- readRDS("Path_to/OrchidSurvey_1.rds")
dataMLCount <- readRDS("Path_to/OrchidSurvey_2.rds")
dataMLMeta <- readRDS("Path_to/OrchidSurvey_3.rds")

## 1. Plot raw data ####
## 2. Plot orchid abundance ####

CountOrchids <- data.frame(Species = c(), Count = c())
for (i in unique(df$Species)){
  df.i <- subset(df, Species == i)
  s.i <- table(unlist(df.i[,6:ncol(df.i)]))["TRUE"]
  CountOrchids <- rbind(CountOrchids, data.frame(Species = i, Count = s.i))
  rm(df.i, s.i)
}
CountOrchids <- CountOrchids[order(CountOrchids$Count, decreasing = TRUE),]
CountOrchids$Species <- factor(CountOrchids$Species, levels = CountOrchids$Species)

ggplot(CountOrchids, aes(x = Species, y = Count))+
  geom_bar(stat = "identity", fill = "#04a896") + theme_classic() +
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 20)) + ylab("Abundance") +
  geom_text(aes(label = Count), vjust = -1, size = 6) +ylim(0, 800)
#ggsave(paste(SortiePath, date, "OrchidIndividualNumber.pdf",  sep=""), width = 8, height = 6, plot = last_plot())

## 3. Plot tree DBH ####

df_stats <- dataMLMeta[,c(1,3,2)]

# Get estimated surface of each tree, S=DBH*H
# DBH = diameter at breast height
df_stats$H <- df_stats$DBH*df_stats$H
colnames(df_stats)[3] <- "S"
df_stats$count <- c(1)

df_stats <- df_stats %>% group_by(Host) %>%
  mutate(M = mean(DBH),
         s = mean(S),
         count = sum(count)) %>%
  arrange(-s)

df_stats$Host <- factor(df_stats$Host, levels = unique(df_stats$Host)) 
df_stats$Host_plot <- paste0(df_stats$Host, " (", df_stats$count, ")")
df_stats$Host_plot <- factor(df_stats$Host_plot, levels = unique(df_stats$Host_plot)) 

# Keep only trees with >10 individuals
df_stats <- df_stats[df_stats$count >= 10,] 

ggplot(df_stats, aes(x=Host_plot, y=S)) +
  geom_boxplot(fill = "#04a896") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, size = 17, vjust = .5, hjust = 1),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 17),
        axis.title.x = element_text(size = 0)) + ylab("DBH * H")
#ggsave(paste(SortiePath, date, "DBHxH_ML_Sup10.pdf",  sep=""), width = 10, height = 6, plot = last_plot())

## 4. Plot network ####

df_count <- dataMLCount
df_count$Host <- dataMLMeta$Host

# Select tree species with > 10 individuals
KeepTreeSpecies <- names(table(df_count$Host))[table(df_count$Host) >= 10]

df_count <- df_count %>% filter(Host %in% KeepTreeSpecies) %>% 
  relocate(Host, 1)
# Store raw presence/absence data for randomization procedure
df_count_pa_raw <- df_count
df_count$count <- c(1)

df_count <- df_count %>% group_by(Host) %>%
  mutate_at(vars(-group_cols()), .funs = sum) %>%
  distinct(.keep_all = TRUE)
df_count$Host <- paste0(df_count$Host, " (", df_count$count, ")")
df_count$count <- NULL

for (k in 1:ncol(df_count)){
  if (colnames(df_count)[k] %in% c("Host")){next}
  colnames(df_count)[k] <- paste0(colnames(df_count)[k], " (", sum(df_count[,k]), ")")
}

df_count <- cbind(df_count[1], decostand(df_count[2:ncol(df_count)], MARGIN = 2, method = "total"))
df_count <- as.data.frame(df_count)

# Melt for use in ggplot
df_countm <- melt(df_count)

# Get Significance
# The significance of each tree species-orchid species interaction will be assessed here through a permutational analysis
# First the raw presence/absence data are randomized 10,000 times
# Then for each randomization of raw data, a new network of tree species-orchid species interactions is built
# Finally, for each tree species-orchid species interaction, the observed value of interaction is compared
# to the 10,000 values obtained after randomization. 
# This comparison includes tests for observed value > or < randomized values

NPerm <- 10000

# permatfull function of vegan allows randomizations of matrices (here raw presence/absence data)
# permatfull keeps connectance, here for rows and columns meaning for "orchids host range" and for "from poor-to-good-host"
p_df_count <- permatfull(df_count_pa_raw[,2:ncol(df_count_pa_raw)],
                         fixedmar = "both", mtype = "prab", time = NPerm) 


# Just check how many matrices we get here
mt <- c()
for (x in 1:NPerm){
  m1 <- paste(unlist(p_df_count$perm[x]), collapse = "")
  mt <- c(m1, mt)
}
length(unique(mt))
length(mt)
# ==> we get 10,000 different matrices out of 10,000 (100%)
rm(mt)

RandomMatrices <- list()
for (r in 1:NPerm){
  df_count_r <- as.data.frame(p_df_count$perm[r])
  df_count_r$Host <- df_count_pa_raw$Host
  df_count_r <- df_count_r %>% group_by(Host) %>%
    mutate_at(vars(-group_cols()), .funs = sum) %>%
    distinct(.keep_all = TRUE)
  df_count_r$Host <- NULL
  
  df_count_r <- decostand(df_count_r, MARGIN = 2, method = "total")

  RandomMatrices[[r]] <- as.data.frame(df_count_r)
  
}
rm(df_count_r)

# Results of comparisons for observed values > or < randomized values are stored in different matrices
MatrixSignifPos <- matrix(NA, nrow = nrow(df_count), ncol = ncol(df_count)-1)
MatrixSignifNeg <- matrix(NA, nrow = nrow(df_count), ncol = ncol(df_count)-1)
MatrixSignif <- matrix(NA, nrow = nrow(df_count), ncol = ncol(df_count)-1)

# Fill the Positive and Negative matrices
for (k in 1:(ncol(df_count)-1)){
  for (i in 1:nrow(df_count)){
    randomVal <- c() 
    
    for (n in 1:NPerm){randomVal <- c(randomVal, RandomMatrices[[n]][i,k])}
    MatrixSignifPos[i, k] <- sum(df_count[i, k+1] > randomVal)/length(randomVal)
    MatrixSignifNeg[i, k] <- sum(df_count[i, k+1] < randomVal)/length(randomVal)
  }
}

# Then take only the significative values (either positive or negative)
for (k in 1:ncol(MatrixSignif)){
  for (i in 1:nrow(MatrixSignif)){
    val1 <- MatrixSignifPos[i, k]
    val1 <- ifelse(val1>.95, val1, 0)
    val2 <- MatrixSignifNeg[i, k]
    val2 <- ifelse(val2>.95, val2, 0)
    
    if (val1 > 0 & val2 > 0){stop("error")}
    
    MatrixSignif[i, k] <- ifelse(val1!=0, val1, -val2)
  }
}

MatrixSignif_m <- melt(MatrixSignif)
df_countm$signif <- MatrixSignif_m$value

df_countm$value2 <- round(df_countm$value*100)
colVal <- unlist(lapply(df_countm$value2, function(x){ifelse(x>=1 & x < 80, "#212f3c", "white")}))

# Plot observed data
p1 <- ggplot(df_countm, aes(x=Host, y=variable, fill=value)) +
  geom_tile(color = "white", linetype = 1, lwd = 1.5) + theme_classic()+
  scale_fill_gradient(low = "white", high = "#04a896", name = "Frequency")+
  geom_text(aes(Host, variable, label = round(value*100)), color = colVal, size = 6)+
  theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.4,  size = 15),
        axis.title = element_text(size = 0),
        legend.position = "none",
        axis.text.y = element_text(size = 15))

BarPlot1 <- data.frame(Host = df_count$Host,
                       Sum = apply(df_count[,2:(ncol(df_count))], 1, function(x){sum(x>0)})) 
p2 <- ggplot(BarPlot1, aes(x = Host, y=Sum)) +
  geom_bar(stat="identity", fill = "#04a896") +theme_classic() + ylab("N orchid species") +
  theme(axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 0))

BarPlot2 <- data.frame(Orchid = colnames(df_count)[-1],
                       Sum = apply(df_count[,2:(ncol(df_count))], 2, function(x){sum(x>0)}))
p3 <- ggplot(BarPlot2, aes(x = Orchid, y=Sum)) +
  geom_bar(stat="identity", fill = "#04a896")+
  theme_classic() + ylab("N host tree species") + coord_flip() +
  theme(axis.title.y = element_text(size = 0),
        axis.title.x = element_text(size = 18),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 0))

p4 <- ggplot(NULL) + geom_blank() + theme(panel.background = element_blank())

p0 <- ggarrange(p2, p4, p1, p3, nrow = 2, ncol = 2, widths = c(2.5,1), heights = c(1,2.5))
#ggsave(paste(SortiePath, date, "DataDescription.pdf",  sep=""), width = 15, height = 12, plot = p0)
#ggsave(paste(SortiePath, date, "DataDescription_shorten.pdf",  sep=""), width = 10, height = 8, plot = p0)

# Plot for significance 
# New values for color scaling
df_countm$signif2 <- unlist(lapply(df_countm$signif, function(x){
  if (x==0){0}
  else (ifelse(x>0, -log10(1-x), log10(1-abs(x) ))) }))
df_countm$signif2[(df_countm$signif2=="Inf")] <- log10(NPerm)
df_countm$signif2[(df_countm$signif2=="-Inf")] <- -log10(NPerm)

ggplot(df_countm, aes(x=Host, y=variable, fill=signif2)) +
  geom_tile(color = "white", linetype = 1, lwd = 1.5) + theme_classic()+
  scale_fill_gradientn(colours = c("#bc3250", "white", "#04a896"))+
  geom_text(aes(Host, variable, label = round(value*100)), color = colVal, size = 6)+
  theme(axis.text.x = element_text(angle =90, hjust = 1, vjust = 0.4,  size = 15),
        axis.title = element_text(size = 0),
        axis.text.y = element_text(size = 15),
        legend.position = "none")
#ggsave(paste(SortiePath, date, "DataDescription-Significance.pdf",  sep=""), width = 9, height = 8, plot = last_plot())

## 5. Plot tree colonization ####

CountByHost <- dataMLCount
CountByHost$Host <- dataMLMeta$Host
CountByHost <- CountByHost[which(apply(CountByHost, 1, function(x){sum(is.na(x))})==0),]
CountByHost <- CountByHost %>% relocate(Host, 1)
CountByHost$count <- c(1)

CountByHost <- CountByHost %>% group_by(Host) %>%
  mutate_at(vars(-group_cols()), .funs = sum)
CountByHost <- as.data.frame(CountByHost[-which(duplicated(CountByHost)),])
# Main tree species : 
KeepTrees <- c("AgaSal", "CyaBor", "NuxVer", "LabCal", "SyzSp", "AntBor", "DorApe", "HomPan")
CountByHost <- subset(CountByHost, Host %in% KeepTrees)
CountByHost$Host <- paste0(CountByHost$Host, " (", CountByHost$count, ")")
for (i in 1:nrow(CountByHost)){
  CountByHost[i, 2:(ncol(CountByHost)-1)] <- CountByHost[i, 2:(ncol(CountByHost)-1)]/CountByHost[i, "count"]
}

CountByHost$count <- NULL
m_CountByHost <- melt(CountByHost)
m_CountByHost$value <- 100*m_CountByHost$value
ggplot(m_CountByHost, aes(x = Host, y = value, group = variable)) +
  geom_bar(aes(fill = variable), stat = "identity", position = position_dodge()) + 
  theme_classic() + ylab("% of colonized tree individuals")+ ylim(0,100)+
  #geom_hline(yintercept = c(25, 50, 75), linetype = "dashed", size = .3, col = "darkgrey")+
  geom_vline(xintercept = seq(1.5, 7.5), size = .3, col = "grey")+
  theme(axis.text.x = element_text(angle = 90, size = 15, vjust = .5, hjust = 1),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size=18),
        legend.title = element_text(colour = "white")) +
  scale_fill_manual(values = c(c("#07658b", "#00808f", "#04a896", "#02c39a", "#f1f4bd", "#9bc9e3", "#4b9cb9", "#0f2f45", "#f5b940", "#ec8931")))
#ggsave(paste(SortiePath, date, "ColonizationOnTreeIndividuals.pdf",  sep=""), width = 7, height = 5, plot = last_plot())

## 6. Random Forest - all analyses ####

df <- df %>% arrange(Species)
# The dataset is subsampled 100 times, each time selecting 50 random orchid individuals per species when possible,
# This prevents the model to overpredicts abundant species (B. prismaticum/A. mauritianum) compared to less abundant ones.  
rf_conf_total <- data.frame(Species = c(), Confusion = c(), NegConfusion = c(), text = c(), SubSampling = c())
importance.totale_total <- data.frame(MeanDecreaseGini = c(), Facteurs = c(), SubSampling = c())
importance.species_total <- data.frame(Var1 = c(), Species = c(), value = c(), SubSampling = c())
ppTot_Total <- data.frame(x = c(), y = c(), sp = c())
set.seed(19)

for (SubSampling in 1:100){
  
  df_sub <- df
  print(paste0("Subsampling n = ", SubSampling, " is being processed"))

  SelectAlea_F <- c()
  for (spa in unique(df_sub$Species)){
    SelectAlea <- NULL
    beg <- which(df_sub$Species==spa)[1]
    end <- rev(which(df_sub$Species==spa))[1]
    if (sum(df_sub$Species==spa) >= 50){SelectAlea <- sample(beg:end)[1:50]} # randomly here
    else {SelectAlea <- which(df_sub$Species==spa)}
    SelectAlea_F <- c(SelectAlea_F, SelectAlea)
  }
  df_sub <- df_sub[SelectAlea_F,]
  table(df_sub$Species)
  
  rf <- randomForest(df_sub[,c(2:3, 5:ncol(df_sub))], df_sub$Species, ntree = 2000)
  
  ## 2.1 Confusion ####
  rf_conf <- data.frame(Species = names(rf$confusion[,11]),
                        Confusion = rf$confusion[,11],
                        NegConfusion = 1-rf$confusion[,11])
  rf_conf$text <- unlist(lapply(rf_conf$NegConfusion, function(x){ifelse(x>0, round(x*100, 1), "")}))
  rf_conf$SubSampling <- SubSampling
  
  # Store results : 
  rf_conf_total <- rbind(rf_conf_total, rf_conf)
  
  ## 2.2 Total importance ####
  importance.totale <- as.data.frame(rf$importance)
  rownames(importance.totale)[rownames(importance.totale)=="co"] <- "Plot"
  importance.totale$Facteurs <- rownames(importance.totale)
  importance.totale <- importance.totale[order(importance.totale$MeanDecreaseGini, decreasing = TRUE),]
  importance.totale$Facteurs <- factor(importance.totale$Facteurs, levels = importance.totale$Facteurs)
  importance.totale$SubSampling <- SubSampling
  
  # Store results : 
  importance.totale_total <- rbind(importance.totale_total, importance.totale)
  
  
  ## 2.3 Importance for each species versus other species ####
  KeepTrees <- c("AgaSal", "CyaBor", "NuxVer", "LabCal", "SyzSp", "AntBor", "DorApe", "HomPan")
  
  importance.species <- data.frame(Var1=c(), Species = c(), value = c())
  for (sp in unique(df_sub$Species)){
    
    df.sp <- df_sub
    df.sp <- df.sp[,c("Species", "DBHxH", "Host", "co", KeepTrees)]
    df.sp$Species <- ifelse(df.sp$Species==sp, as.character(sp), "Other")
    df.sp$Species <- factor(df.sp$Species)
    rf.sp <- randomForest(df.sp[,c(2, 4:ncol(df.sp))], df.sp$Species, ntree = 2000)

    if (sp == "BulVar"){
      rf.sp.BulVar <- rf.sp
    } # For partial plot below
    
    #print(rf.sp) 
    importance.species.sp <- melt(rf.sp$importance)
    colnames(importance.species.sp)[2] <- "Species"
    importance.species.sp$Species <- c(sp)
    importance.species.sp$value <- decostand(importance.species.sp$value, method = "total", MARGIN = 2)
    
    importance.species <- rbind(importance.species, importance.species.sp)
  }
  
  importance.species <- importance.species %>% mutate(Var1=recode_factor(Var1, co="Plot"))
  
  importance.species$Var1 <- factor(importance.species$Var1, levels = c("DBH", "DBHxH", "Plot", "AgaSal", "CyaBor", "DorApe",
                                                                        "SyzSp", "AntBor", "HomPan", "LabCal", "NuxVer"))
  importance.species$SubSampling <- SubSampling
  rf_conf$NegConfusion
  importance.species_total <- rbind(importance.species_total, importance.species)
  
  ## 2.4 Partial plots ####

  pp1 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "AgaSal", plot = FALSE)
  pp2 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "CyaBor", plot = FALSE)
  pp3 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "DorApe", plot = FALSE)
  pp4 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "SyzSp", plot = FALSE)
  pp5 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "AntBor", plot = FALSE)
  pp6 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "HomPan", plot = FALSE)
  pp7 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "LabCal", plot = FALSE)
  pp8 <- partialPlot(rf.sp.BulVar, pred.data = df.sp, x.var = "NuxVer", plot = FALSE)
  
  ppTot <- data.frame(x=c(pp1$x, pp2$x, pp3$x, pp4$x, pp5$x, pp6$x, pp7$x, pp8$x),
                      y=c(pp1$y, pp2$y, pp3$y, pp4$y, pp5$y, pp6$y, pp7$y, pp8$y),
                      sp=c(rep("AgaSal", 2), rep("CyaBor", 2), rep("DorApe", 2), rep("SyzSp", 2),
                           rep("AntBor", 2), rep("HomPan", 2), rep("LabCal", 2), rep("NuxVer", 2)))
  ppTot$SubSampling <- SubSampling
  
  ppTot_Total <- rbind(ppTot_Total, ppTot)
  
}

## Plot section : 

# 1. Plot confusion
rf_conf_total 
Plot_conf <- rf_conf_total %>% group_by(Species) %>% mutate(M = mean(NegConfusion),
                                                           SD = sd(NegConfusion)) %>%
  select(Species, M, SD) %>% distinct()
Plot_conf$Mc <- round(Plot_conf$M*100, 1)


ggp1 <- ggplot(Plot_conf, aes(x=Species, y=M))+
  geom_bar(fill = "#04a896", col = "white", stat = "identity", alpha = .7) + ylab("Successful prediction") + theme_classic() + xlab("") +
  geom_errorbar( aes(ymin=M-SD, ymax=M+SD), col = "#323031", width = .3, size = .3)+
  theme(axis.text.x = element_text(size = 15, vjust = .5, angle = 90, hjust = 1),
        axis.text.y = element_text(size = 13, angle = 45, hjust = .5),
        axis.title.y = element_text(size=20),
        axis.title.x = element_text(size = 20),
        legend.position = "none")+
  scale_y_continuous(breaks = c(0, .5, 1))+
  geom_hline(yintercept = seq(0, 1, .25), linetype = "dashed", col = "darkgrey", size = .2)+
  geom_vline(xintercept = seq(1.5, 10.5, 1), linetype = "dashed", col = "darkgrey", size = .4)+
  annotate(geom = "text", x=1:10, y=Plot_conf$M+0.2, label = Plot_conf$Mc, size = 5, angle = 90)
  
ggp1

# 2. Plot total importance 
importance.totale_total

Plot_itt <- importance.totale_total %>% group_by(Facteurs) %>% mutate(SD = sd(MeanDecreaseGini),
                                                                      MeanDecreaseGini = mean(MeanDecreaseGini)) %>%
  select(-SubSampling) %>% distinct()
Plot_itt$MeanDecreaseGini <- Plot_itt$MeanDecreaseGini[order(Plot_itt$MeanDecreaseGini, decreasing = TRUE)]
Plot_itt <- Plot_itt %>% mutate(Facteurs=recode_factor(Facteurs, DBHxH="DBH*H"))

ggplot(Plot_itt, aes(x=Facteurs, y=MeanDecreaseGini))+
  geom_bar(stat = "identity", fill = "#04a896", col = "white") + 
  geom_errorbar( aes(ymin=MeanDecreaseGini-SD, ymax=MeanDecreaseGini+SD), col = "#323031", width = .3)+
  ylab("Importance\n(mean decrease Gini)") + theme_classic() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 15),
        axis.title.x = element_text(size = 0),
        axis.title.y = element_text(size = 20),
        axis.text.y = element_text(size = 15))
#ggsave(paste(SortiePath, date, "ImportanceTotale.pdf",  sep=""), width = 8, height = 6, plot = last_plot())


# 3. Plot per species importance
importance.species_total
mean(subset(importance.species_total, Species=="BulVar" & Var1=="AgaSal")$value)
sd(subset(importance.species_total, Species=="BulVar" & Var1=="AgaSal")$value)

importance.species_total$VarSpe <- paste(importance.species_total$Species, importance.species_total$Var1, sep = "_")
Plot_itt <- importance.species_total %>% group_by(VarSpe) %>% mutate(SD = sd(value), value = mean(value)) %>%
  select(Species, Var1, value, SD) %>% distinct()

ggp2 <- ggplot(Plot_itt, aes(x=Species, y=value, fill=as.factor(Var1)))+
  geom_bar(stat="identity", col = "white", position=position_dodge(width=0.9)) +
  geom_errorbar( aes(ymin=value, ymax=value+SD), col = "#323031", width = .5, position=position_dodge(width=0.9))+
  scale_fill_manual(values = c("DBH" = "#07658b", "DBHxH" = "#00808f","Plot" = "#04a896",
                               "AgaSal" ="#bc3250", "CyaBor"="#0f2f45", "DorApe"="#165577",
                               "SyzSp" ="#49a2c3", "AntBor"="#9bc9e3", "HomPan"="#cde5f2", "LabCal"="#ec8931", "NuxVer"="#fbb13c"))+
  theme_classic() +
  geom_vline(xintercept = seq(1.5, 14.5), linetype = "dashed", col = "darkgrey") +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(size = 15, angle = 90, hjust = .5),
        legend.title = element_text(size = 0),
        legend.text = element_text(size = 15),
        legend.spacing.x = unit(.5, "cm"),
        legend.background = element_rect(linewidth = .25, linetype = "solid", colour = "black")) + xlab("") +
  ylab("Relative importance")
ggp2

gg0 <- ggarrange(ggp2,ggp1, nrow = 2, heights = c(1,.2))
gg0
#ggsave(paste(SortiePath, date, "RelImportance_RF_DualPlot-2.pdf",  sep=""), width = 12, height = 6.5, plot = gg0)


# 4. Plot partial effects for BulVar
ppTot_Total

ppTot_Total$xsp <- paste(ppTot_Total$x, ppTot_Total$sp, sep = "_")
Plot_ME <- ppTot_Total %>% group_by(xsp) %>% mutate(SD = sd(y), y = mean(y)) %>%
  select(x, sp, y, SD, xsp) %>% distinct()

ggplot(Plot_ME, aes(x=sp, y=y, group = x))+
  geom_bar(aes(fill = x), stat="identity", col = "white", position=position_dodge(width=0.9)) +
  geom_errorbar( aes(ymin=y-SD, ymax=y+SD), col = "#323031", width = .5, position=position_dodge(width=0.9))+
  theme_classic() + #coord_flip()+
  scale_fill_manual(values = c("TRUE" = "#e05972", "FALSE" = "#218380")) +
  scale_x_discrete(position = "top")+
  geom_hline(yintercept = 0,  size = .2, col = "black") +
  geom_vline(xintercept = seq(1.5, 7.5, 1), size = .4, col = "grey") +
  ylab("Marginal effect logit(p)") +
  theme(axis.title.x = element_text(size = 0),
        axis.text.x = element_text(size = 15, angle = 90),
        axis.text.y = element_text(size = 15),
        axis.title.y =element_text(size = 17),
        legend.title = element_text(size = 0))
#ggsave(paste(SortiePath, date, "MarginalEffect-BulVar-AllTrees.pdf",  sep=""), width = 6, height = 4, plot = last_plot())

## 7. Chi square analysis - Goodness of fit ####

KeepTrees <- c("AgaSal", "CyaBor", "NuxVer", "LabCal", "SyzSp", "AntBor", "DorApe", "HomPan")

df_stats <- dataMLMeta[,c(1,3,2)]
df_stats <- df_stats %>% mutate(H = DBH*H) %>% rename("S"="H") %>% mutate(count = c(1))
df_stats$Host <- factor(df_stats$Host, levels = unique(df_stats$Host)) 
df_stats <- subset(df_stats, Host %in% KeepTrees)

# 
df_count <- dataMLCount
df_count$Host <- dataMLMeta$Host
df_count <- df_count %>% filter(Host %in% KeepTrees) %>% relocate(Host, 1) %>% mutate(count = c(1))

df_count <- df_count %>% group_by(Host) %>%
  mutate_at(vars(-group_cols()), .funs = sum) %>%
  distinct(.keep_all = TRUE)
df_count$Host <- paste0(df_count$Host, " (", df_count$count, ")")
df_count$count <- NULL
df_count_abs <- df_count 

# Three models considered : 
# Model 1 : equal probabilities for each species - no of tree individuals NOT considered
# Model 2 : equal probabilities for each species - no of tree individuals considered
# Model 3 : probabilities equal to total surface i.e. no of tree individuals considered

## Required : choose one model above and then run the code below. Save the plot and re-run the code from start. 

# 1. Model 1
df_stats <- df_stats[,c("Host", "count")] %>% group_by(Host) %>% mutate(count = 1) %>% distinct(.keep_all = TRUE)
# 2. Model 2
df_stats <- df_stats[,c("Host", "count")] %>% group_by(Host) %>% mutate(count = sum(count)) %>% distinct(.keep_all = TRUE)
# 3. Model 3
df_stats <- df_stats[,c("Host", "S")] %>% group_by(Host) %>% mutate(S = sum(S)) %>% distinct(.keep_all = TRUE)

Ref_factor <- sapply(df_count_abs$Host, function(x){
  x <- unlist(str_split(x, "\\ \\("))[1]
  df_stats[str_detect(df_stats$Host, x), 2] # column 2 i.e. count, S or M
})

if (length(Ref_factor)!=nrow(df_count_abs)){stop("problem")}

data.chi <- data.frame(Species = c(), val1 = c(), val2 = c(), val3 = c())
tot_bestChiNames <- data.frame(Species = c(), val1 = c(), val2 = c(), val3 = c())
ValFish <- c()
rm( list = paste0("res.", colnames(dataMLCount)) ) # just a security: a error is return at first use

for (OrchSpec in colnames(dataMLCount)){
  n <- which(str_detect(colnames(df_count_abs), OrchSpec))
  observed <- df_count_abs[,n]
  expected <- c(decostand(unlist(Ref_factor), 2, method = "total"))
  
  if (sum(expected) != 1){stop("Error")}
  test_chi <- chisq.test(x = observed, p = expected)
  # Test for significance
  test_Fish <- fisher.test(data.frame(x = observed, y = round(expected*sum(observed))), simulate.p.value=TRUE)
  # simulate.p.value : a logical indicating whether to compute p-values by Monte Carlo simulation, in larger than 2×22×2 tables.
  assign(paste0("res.", OrchSpec), test_chi)
  
  ValFish <- c(ValFish, test_Fish$p.value)
  
  test_chi.res <- as.data.frame(test_chi$residuals)
  colnames(test_chi.res) <- "V1"
  test_chi.res$Host <- seq(1, nrow(test_chi.res)) # Order matters
  test_chi.res <- test_chi.res[order(test_chi.res$V1, decreasing = TRUE),]
  
  bestChi <- test_chi.res$V1[1:3]
  bestChiNames <- c(as.numeric(rownames(test_chi.res)[1:3]))
  bestChiNames <- df_count_abs$Host[bestChiNames]
  bestChiNames <- sapply(bestChiNames, function(x){unlist(str_split(x, "\\ \\("))[1]})
  
  data.chi <- rbind(data.chi,
                    data.frame(Species = c(OrchSpec), val1 = c(bestChi[1]), val2 = c(bestChi[2]), val3 = c(bestChi[3])))
  tot_bestChiNames <- rbind(tot_bestChiNames,
                            data.frame(Species = c(OrchSpec), val1 = c(bestChiNames[1]), val2 = c(bestChiNames[2]), val3 = c(bestChiNames[3])))
}
#data.chi
m.data.chi <- reshape2::melt(data.chi)

#tot_bestChiNames
m.tot_bestChiNames <- reshape2::melt(tot_bestChiNames, id = "Species")
m.data.chi$Names <- m.tot_bestChiNames$value

FillSpeciesTree <- c("AgaSal" ="#bc3250", "CyaBor"="#0f2f45", "DorApe"="#165577", 
                     "SyzSp" ="#49a2c3", "AntBor"="#9bc9e3", "HomPan"="#cde5f2", "LabCal"="#ec8931", "NuxVer"="#fbb13c")

ValFish <- sapply(ValFish, function(x){
  if (x > 0.05){"n.s"}
  else if (x <= 0.05 & x > 0.01){"*"}
  else if (x <= 0.01 & x > 0.005){"**"}
  else if (x <= 0.005){"***"}
})

gg0 <- ggplot(m.data.chi, aes(x=Species, y=value, group = variable))+
  geom_bar(aes(fill = Names), stat = "identity", position = position_dodge(), width = .8) + 
  theme_classic() +
  scale_fill_manual(values = FillSpeciesTree, breaks = names(FillSpeciesTree)) +
  geom_vline(xintercept = seq(1.5, 10.5, 1), linetype = "dashed", col = "darkgrey")+
  geom_hline(yintercept = 5, linetype = "dotted", col = "darkgrey")+ # ligne with residuals > 5 == very signi
  ylab("Chisq-test's residuals") 

# Model 1
gg1 <- gg0 + ylim(-2, 20) + 
  annotate(geom = "text", x = c(1:10), y=rep(20, 10), label = ValFish, size = 7)+
  theme(axis.text.x = element_text(size=0),
        axis.title.x = element_text(size = 0),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20))

# Model 2
gg2 <- gg0 + ylim(-3, 30) + 
  annotate(geom = "text", x = c(1:10), y=rep(28, 10), label = ValFish, size = 7)+
  theme(axis.text.x = element_text(size=0),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 0))

# Model 3
gg3 <- gg0 + ylim(-3, 20) + 
  annotate(geom = "text", x = c(1:10), y=rep(20, 10), label = ValFish, size = 7)+
  theme(axis.text.x = element_text(size=18, angle = 90, vjust = .5, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 0))


p0 <- ggarrange(gg1, gg2, gg3, nrow = 3, ncol = 1, widths = c(1), heights = c(1,1,1))

#ggsave(paste(SortiePath, date, "GoodnessOfFit-S.pdf",  sep=""), width = 12, height = 12, plot = p0)
