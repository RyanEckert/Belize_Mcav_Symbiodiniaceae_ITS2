################################################################################
### Statistical analysis of Symbiodiniaceae ITS2 reads     version 09/29/19  ###
### Written by Ryan Eckert (ryan.j.eckert@gmail.com)                         ###
### adapted from pipeline developed by Carly Kenkel (Kenkel and Bay, 2018)   ###
### (http://doi.org/10.5281/zenodo.1208684)                                  ###
###--------------------------------------------------------------------------###
###                                                                          ###
### This is my adaptation of analyzing ITS2 amplicon sequence variants       ###
### This can definitely be done more elegantly, but this script works and    ###
### can be adapted to work with your data.                                   ###
###                                                                          ###
### BEFORE STARTING, replace:                                                ###
### - path/to/working/directory with the literal path to where you stored    ###
###   your processed sequencing reads                                        ###
### The lines beginning with hash marks (#) are explanations and additional  ###
### instructions, or commented out code that you may wish to un-comment and  ###
### run. Make sure to read the comments before running the code.             ###
###                                                                          ###
################################################################################
##----------SET WORKING DIRECTORY-----------------------------------------------
# First set the R working directory to the directory where your files are saved
# change the path below to your working directory path
# setwd("~/Desktop/ITS2_pipeline/test_ITS2_analysis")

setwd("~/path/to/working/directory")

##----------INSTALL AND LOAD PACKAGES-------------------------------------------
# un-comment the below lines to install packages for the first time

# install.packages("BiocManager")
# install.packages("devtools")
# library(BiocManager)
# library(devtools)
 
# install.packages("car")
# install.packages("cluster")
# install.packages("cowplot")
# install.packages("dunn.test")
# BiocManager::install("edgeR")
# install.packages("forcats")
# install.packages("ggfortify")
# install.packages("ggplot2")
# install.packages("ggpubr")
# install.packages("gridExtra")
# install.packages("labdsv")
# install.packages("MASS")
# install.packages("MCMC.OTU")
# install.packages("multcompView")
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# install.packages("pgirmess")
# install.packages("vegan")

library(car)
library(cluster)
library(cowplot)
library(dunn.test)
library(edgeR)
library(forcats)
library(ggfortify)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(labdsv)
library(MASS)
library(MCMC.OTU)
library(multcompView)
library(pgirmess)
library(pairwiseAdonis)
library(vegan)

#setting seed allows randomized processes to be repeated later
set.seed(3254)

##----------LOAD DATA INTO R----------------------------------------------------
# Load in the data from dada2 and LULU

curated_result_97 = readRDS("curated_result_97.rds")
ITS2data <- read.csv("ASVtable.csv")
head(ITS2data)

##----------CREATE A COLOR PALETTE----------------------------------------------
# We need 10 colors, I'm using ggPlot colors, but in a better order for 
# visualizing ASV bar chart

its2ColPal=c("#0081EF", "#00B0F6", "#00BFC4", "#39B600", "#9590FF", "#E76BF3", 
             "#FF62BC", "#F8766D", "#D89000", "#A3A500")

##----------DATA SETUP----------------------------------------------------------
# First, we need to set up the data we will use, and purge outlying OTUS and
# low quality samples.Pull out the curated ASV list, re-transpose and re-combine
# with additional sample data

ITS2data = cbind(ITS2data[, 1:3], data.frame(t(curated_result_97$curated_table)))

# now ordering sites from north to south so later data can be plotted shallow to
# deep; north to south
ITS2data$Depth = factor(ITS2data$Depth, levels = c("10", "16", "25", "35"))
levels(ITS2data$Site)
ITS2data$Site = factor(ITS2data$Site, levels(ITS2data$Site)[c(4, 2, 3, 1)])
ITS2data = ITS2data[order(ITS2data$Site, ITS2data$Depth), ]
ITS2data$newOrder = 1:nrow(ITS2data)
ITS2data = cbind(ITS2data[, length(ITS2data), drop = FALSE], 
                 ITS2data[, c(1:length(ITS2data) - 1)])
row.names(ITS2data) = ITS2data$Sample
head(ITS2data)

goods = purgeOutliers(ITS2data, count.columns = 5:length(ITS2data),
                      otu.cut = 0.001)
goods$newOrder = 1:nrow(goods)
row.names(goods) = goods$sample
head(goods)
# write.csv(goods, "goodsASV.csv", row.names = FALSE)
# un-comment the above line to save your good ASVs to the disk

goods = read.csv("goodsASV.csv")
row.names(goods) = goods$sample

# what is the proportion of samples with data for these ASVs?
withData = apply(goods[, 5:length(goods[1, ])], 2, 
                 function(x) 
                   {
	                 sum(x > 0)/length(x)
                   }
                   )
hist(withData, breaks = 50)

props = apply(goods[, 5:length(goods[1, ])], 2, 
              function(x)
                {
	              sum(x)/sum(goods[, 5:length(goods[1, ])])
                }
                )
barplot(props, xaxt = "n", log = "y")
props

# normalize otu counts with weighted trimmed mean of M-values (TMM; Robinson and
# Oshlack 2010)
itsGoodsTransposed = t(goods[, 5:length(goods[1, ])])
itsGoodsList = DGEList(counts = itsGoodsTransposed)
head(itsGoodsList$samples)
its2Norm =  calcNormFactors(itsGoodsList, method = "TMM")
head(its2Norm$samples)
its2TMM = t(cpm(its2Norm, normalized.lib.sizes = TRUE))
its2Norm = cbind(goods[,c(2:4)], its2TMM)
head(its2Norm)
#write.csv(its2Norm, "its2_NormalizedGoods.csv", row.names = FALSE)
its2Norm = read.csv("its2_NormalizedGoods.csv")
# renaming sequence coulmns to include Clade types identified with Blast-n
colnames(its2Norm)[4:ncol(its2Norm)] = c("sq01_C3", "sq10_C3", "sq11_C3g",
	       "sq14_C3g", "sq17_C3g", "sq18_C3g", "sq25_C3g", "sq05_C3z", "sq07_C3e", 
	       "sq08_C3g")
its2Norm$Site = factor(its2Norm$Site, levels(its2Norm$Site)[c(4, 2, 3, 1)])
its2Norm$Depth = as.factor(its2Norm$Depth)
head(its2Norm)
its2Norm$Site
its2Norm$Depth
# Everything looks good now!

##----------CALCULATE ALPHA DIVERSITY MEASURES---------------------------------
# Species richness per sample

its2Alpha = specnumber(its2Norm[,4:ncol(its2Norm)])
its2Shannon = diversity(its2Norm[,4:ncol(its2Norm)],index = "shannon")
its2Simpson = diversity(its2Norm[,4:ncol(its2Norm)],index = "simpson")
its2Div = cbind(its2Norm[,1:3], its2Alpha, its2Shannon, its2Simpson)
head(its2Div)
colnames(its2Div)[4:6] = c("alpha", "shannon", "simpson")
its2Div$effectiveShannon = exp(its2Div$shannon)
its2Div$effectiveSimpson = (1/(1-its2Div$simpson))
head(its2Div)

##----------TESTING DATA FOR NORMALITY AND HOMOSCHEDASTICITY--------------------


tapply(its2Div$alpha, its2Div$Site, shapiro.test)
tapply(its2Div$shannon, its2Div$Site, shapiro.test)
tapply(its2Div$simpson, its2Div$Site, shapiro.test)
tapply(its2Div$alpha, its2Div$Site, shapiro.test)
tapply(its2Div$shannon, its2Div$Site, shapiro.test)
tapply(its2Div$simpson, its2Div$Site, shapiro.test)

leveneTest(its2Div$alpha, its2Div$Depth)
leveneTest(its2Div$shannon, its2Div$Depth)
leveneTest(its2Div$simpson, its2Div$Depth)
leveneTest(its2Div$alpha, its2Div$Site)
leveneTest(its2Div$shannon, its2Div$Site)    
leveneTest(its2Div$simpson, its2Div$Site)
# Data are non-normally distributed and also heteroschedastic
# Using non-parametric tests moving forward

##----------KRUSKAL-WALLIS TEST ON SHANNON INDICES------------------------------


kruskal.test(shannon ~ Site, data = its2Div)
# 	Kruskal-Wallis rank sum test
# 
# data:  shannon by Site
# Kruskal-Wallis chi-squared = 2.7884, df = 3, p-value = 0.4254

kruskal.test(shannon ~ Depth, data = its2Div)
# 	 Kruskal-Wallis rank sum test
# 
# data:  shannon by Depth
# Kruskal-Wallis chi-squared = 3.0743, df = 3, p-value = 0.3803

# NOTHING SIGNIFICANT

##----------KRUSKAL-WALLIS TEST ON SIMPSON INDICES------------------------------


kruskal.test(simpson ~ Site, data = its2Div)
# 	Kruskal-Wallis rank sum test
# 
# data:  simpson by Site
# Kruskal-Wallis chi-squared = 3.4937, df = 3, p-value = 0.3216

kruskal.test(simpson ~ Depth, data = its2Div)
# 	 Kruskal-Wallis rank sum test
# 
# data:  simpson by Depth
# Kruskal-Wallis chi-squared = 3.8446, df = 3, p-value = 0.2787

# NOTHING SIGNIFICANT

##----------KRUSKAL-WALLIS TEST ON SPECIES RICHNESS-----------------------------


kruskal.test(alpha ~ Site, data = its2Div)
#     Kruskal-Wallis rank sum test
# 
# data:  alpha by Site
# Kruskal-Wallis chi-squared = 2.7628, df = 3, p-value = 0.4297

kruskal.test(alpha ~ Depth, data = its2Div)
# 	 Kruskal-Wallis rank sum test
# 
# data:  alpha by Depth
# Kruskal-Wallis chi-squared = 38.993, df = 3, p-value = 1.741e-08

##----------DUNN'S TEST ON SPECIES RICHNESS ACROSS DEPTH------------------------
# Using Dunn test to do multiple comparisons with bonferroni correction
dunn.test(its2Div$alpha, g = its2Div$Depth, method = "bonferroni", table = FALSE, 
          list = TRUE)

#  Kruskal-Wallis rank sum test
# 
# data: x and group
# Kruskal-Wallis chi-squared = 38.9935, df = 3, p-value = 0
# 
# 
#                          Comparison of x by group                            
#                                 (Bonferroni)                                  
# 
# List of pairwise comparisons: Z statistic (adjusted p-value)
#
# 10 - 16 :  1.449108 (0.4419)
# 10 - 25 :  4.989355 (0.0000)*
# 16 - 25 :  3.504135 (0.0014)*
# 10 - 35 :  5.057011 (0.0000)*
# 16 - 35 :  3.591117 (0.0010)*
# 25 - 35 :  0.133160 (1.0000)
# 
# alpha = 0.05
# Reject Ho if p <= alpha/2

##----------PLOTTING ALPHA METRICS IN A SINGLE FIGURE GRID----------------------
# Plot a-diversity by site
its2AlphaPlotS = ggplot(data = its2Div, aes(x = Site, y = alpha)) +
                 geom_point(aes(fill = Depth), shape = 21, color = "gray20", 
                            size = 2, position = position_jitterdodge()) +
                 geom_boxplot(alpha = 0, outlier.shape = NA, color = "gray30") +
                 scale_fill_manual(values = its2ColPal[c(1,4,6,2)], 
                                   labels = c("10 m", "16 m", "25 m", "35 m")) +  
                 xlab("Reef Site") + 
                 ylab("Richness") + 
                 stat_compare_means(size = 5, geom = "label", label.x = 1, 
                                    label.y = 10.7) + 
                 coord_cartesian( ylim = c(1, 11)) +
                 scale_y_continuous(breaks = seq(2, 10, 2))+
                 guides(fill=guide_legend(ncol=2))+
                 theme_bw()

its2AlphaPlotSite = its2AlphaPlotS +
                    theme(axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.title.y = element_text(color = "black", size = 14),
                          axis.text.y = element_text(color = "black", size = 12),
                          legend.title = element_text(color = "black", size = 14),
                          legend.text = element_text(color = "black", size = 12),
                          panel.border = element_rect(size = 1.1, color = "black"),
                          panel.background = element_rect(fill = "white"),
                          plot.background = element_rect(fill = "white"),
                          legend.position = "bottom",
                          legend.background = element_rect(color = "black"),
                          plot.margin = unit(c(0.5, 0.1, 0.1, 0.61), "cm"))

# its2AlphaPlotSite

# takes the legend and saves it as a separate graphical object (grob), this way we
# can have 1 legend per group of similar figures
get_legend <- function(its2AlphaPlotSite) 
                {
                tmp <- ggplot_gtable(ggplot_build(its2AlphaPlotSite))
                leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                legend <- tmp$grobs[[leg]]
                return(legend)
                }
legend.site = get_legend(its2AlphaPlotSite)

# re-plot without legend
its2AlphaPlotSite = its2AlphaPlotSite + theme(legend.position = "none")

# K-W test found significant differences, so we will prepare pairwise comparisons 
# to plot letters to denote simmilarity between groups on the plot
alphaDkmc = kruskalmc(its2Div$alpha ~ its2Div$Depth) # multiple-comparison test
print(alphaDkmc)

# Multiple comparison test after Kruskal-Wallis 
# p.value: 0.05 
# Comparisons
#         obs.dif critical.dif difference
# 10-16 14.040805     32.74642      FALSE
# 10-25 48.134463     32.60500       TRUE
# 10-35 49.441667     33.04241       TRUE
# 16-25 34.093659     32.88254       TRUE
# 16-35 35.400862     33.31630       TRUE
# 25-35  1.307203     33.17732      FALSE

alphaDkmcDiff = alphaDkmc$dif.com$difference # select logical vector
names(alphaDkmcDiff) = row.names(alphaDkmc$dif.com)# add comparison names

# create a list with "homogenous groups" coded by letter
alphaDsigLetters = multcompLetters(alphaDkmcDiff, compare="<", threshold=0.05,
  Letters=c(letters, LETTERS, "."), reversed = FALSE)
alphaLetters = as.data.frame(alphaDsigLetters$Letters)
alphaLetters$depth = row.names(alphaLetters)
colnames(alphaLetters)[1] = "id"
alphaLetters$y = c(1, 1, 1, 1)
head(alphaLetters)

its2AlphaPlotD = ggplot(data = its2Div, aes(x = Depth, y = alpha)) +
                 geom_point(aes(fill = Site), shape = 21, color = "gray 20", 
                            size = 2, position = position_jitterdodge()) +
                 geom_boxplot(alpha = 0, outlier.shape = NA, color = "gray30") +
                 scale_fill_manual(values = its2ColPal[c(8,5,3,9)], 
                                   labels = c("Tobacco Reef", "Raph's Wall", 
                                              "South Reef", "Gover's Reef")) +
                 stat_compare_means(size = 5, geom = "label", label.x = 1.07, 
                                    label.y = 10.7) +
                 geom_text(data = alphaLetters, aes(x = depth, y = y, label = id), 
                           size = 6) + 
                 coord_cartesian( ylim = c(1, 11)) +
                 scale_y_continuous(breaks = seq(2, 10, 2))+
                 xlab("Depth (m)") +
                 ylab("Richness") +
                 guides(fill=guide_legend(ncol=2))+
                 theme_bw()

its2AlphaPlotDepth = its2AlphaPlotD +
                     theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank(),
                           axis.title.y = element_blank(),
                           axis.text.y = element_blank(),
                           legend.title = element_text(color = "black", size = 14),
                           legend.text = element_text(color = "black", size = 12),
                           panel.border = element_rect(size = 1.1, color = "black"),
                           panel.background = element_rect(fill = "white"),
                           plot.background = element_rect(fill = "white"),
                           legend.position = "bottom",
                           legend.background = element_rect(color = "black"),
                           plot.margin = unit(c(0.5, 0.5, 0.1, 0.85), "cm"))

# its2AlphaPlotDepth

#takes the legend and saves it as a separate object (grob)
get_legend = function(its2AlphaPlotDepth) 
                {
                tmp <- ggplot_gtable(ggplot_build(its2AlphaPlotDepth))
                leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
                legend <- tmp$grobs[[leg]]
                return(legend)
                }
legend.depth = get_legend(its2AlphaPlotDepth)

# re-plot without legend
its2AlphaPlotDepth = its2AlphaPlotDepth + theme(legend.position = "none")

# Plot Shannon's index by site
its2ShannonPlotS = ggplot(data = its2Div, aes(x = Site, y = shannon)) +
                   geom_point(aes(fill = Depth), shape = 21, color = "gray20", 
                              size = 2, position = position_jitterdodge()) +
                   geom_boxplot(alpha = 0, outlier.shape = NA, color = "gray30") +	
                   scale_fill_manual(values = its2ColPal[c(1,4,6,2)], 
                                     labels = c("10 m", "16 m", "25 m", "35 m")) +
                   stat_compare_means(size = 5, geom = "label", label.x = 1, 
                                      label.y = 1.7) + 	
                   expand_limits(y = c(0, 1.75)) +
                   xlab("Reef Site") +
                   ylab("Shannon's index") +
                   theme_bw()

its2ShannonPlotSite = its2ShannonPlotS + 
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.title.y = element_text(color = "black", size = 14),
                            axis.text.y = element_text(color = "black", size = 12),
                            legend.title = element_text(color = "black", size = 14),
                            legend.text = element_text(color = "black", size = 12),
                            panel.border = element_rect(size = 1.1, color = "black"),
                            panel.background = element_rect(fill = "white"),
                            plot.background = element_rect(fill = "white"),
                            legend.position = "none",
                            legend.background = element_rect(color = "black"),
                            plot.margin = unit(c(0.1, 0.1, 0.1, 0.5), "cm"))

# its2ShannonPlotSite

# Plot Shannon's index by depth
its2ShannonPlotD = ggplot(data = its2Div, aes(x = Depth, y = shannon))+
                   geom_point(aes(fill = Site), shape = 21, color = "gray20", 
                              size = 2, position = position_jitterdodge()) +
                   geom_boxplot(alpha = 0, outlier.shape = NA, color = "gray30") +	
                   scale_fill_manual(values = its2ColPal[c(8,5,3,9)], 
                                     labels = c("Tobacco Reef", "Raph's Wall", 
                                     "South Reef", "Gover's Reef")) +
                   stat_compare_means(size = 5, geom = "label", label.x = 1, 
                                      label.y = 1.7) +
                   expand_limits(y = c(0, 1.75)) +
                   xlab("Depth (m)") +
                   ylab("Shannon's index") +
                   theme_bw()

its2ShannonPlotDepth = its2ShannonPlotD + 
                       theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.title.y = element_blank(),
                             axis.text.y = element_blank(),
                             legend.title = element_text(color = "black", size = 14),
                             legend.text = element_text(color = "black", size = 12),
                             panel.border = element_rect(size = 1.1, color = "black"),
                             panel.background = element_rect(fill = "white"),
                             plot.background = element_rect(fill = "white"),
                             legend.position = "none",
                             legend.background = element_rect(color = "black"),
                             plot.margin = unit(c(0.1, 0.5, 0.1, 0.85), "cm"))

# its2ShannonPlotDepth


# Plot Simpson's index by site
its2SimpsonPlotS = ggplot(data = its2Div, aes(x = Site, y = simpson)) +
                   geom_point(aes(fill = Depth), shape = 21, color = "gray20", 
                              size = 2, position = position_jitterdodge()) +
                   geom_boxplot(alpha = 0, outlier.shape = NA) +	
                   scale_fill_manual(values = its2ColPal[c(1,4,6,2)], 
                                     labels = c("10 m", "16 m", "25 m", "35 m")) +
                   stat_compare_means(size = 5, geom = "label", label.x = 1, 
                                      label.y = 0.77) + 	
                   expand_limits(y = c(0, 0.8)) +
                   xlab("Reef Site") +
                   ylab("Simpson's index") +
                   theme_bw()

its2SimpsonPlotSite = its2SimpsonPlotS + 
                      theme(axis.title.x = element_text(color = "black", size = 14),
                            axis.text.x = element_text(color = "black", size = 12),
                            axis.title.y = element_text(color = "black", size = 14),
                            axis.text.y = element_text(color = "black", size = 12),
                            legend.title = element_text(color = "black", size = 14),
                            legend.text = element_text(color = "black", size = 12),
                            panel.border = element_rect(size = 1.1, color = "black"),
                            panel.background = element_rect(fill = "white"),
                            plot.background = element_rect(fill = "white"),
                            legend.position = "none",
                            legend.background = element_rect(color = "black"),
                            plot.margin = unit(c(0.11, 0.1, 0, 0.5), "cm"))

# its2SimpsonPlotSite

# plot simpson's index by depth
its2SimpsonPlotD = ggplot(data = its2Div, aes(x = Depth, y = simpson)) +
                   geom_point(aes(fill = Site), shape = 21, color = "gray20",
                              size = 2, position = position_jitterdodge()) +
                   geom_boxplot(alpha = 0, outlier.shape = NA) +	
                   scale_fill_manual(values = its2ColPal[c(8,5,3,9)], 
                                     labels = c("Tobacco Reef", "Raph's Wall", 
                                     "South Reef", "Gover's Reef")) +
                   stat_compare_means(size = 5, geom = "label", 
                                      label.x = 1, label.y = 0.77) + 	
                   expand_limits(y = c(0, 0.8)) +
                   xlab("Depth (m)") +
                   ylab("Simpson's index") +
                   theme_bw()

its2SimpsonPlotDepth = its2SimpsonPlotD + 
  theme(axis.title.x = element_text(color = "black", size = 14),
        axis.text.x = element_text(color = "black", size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 12),
        panel.border = element_rect(size = 1.1, color = "black"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "none",
        legend.background = element_rect(color = "black"),
        plot.margin = unit(c(0.1, 0.5, 0, 0.85), "cm"))

# its2SimpsonPlotDepth

# Throw all those badboys into a single figure grid
divPlots = ggarrange(its2AlphaPlotSite, its2AlphaPlotDepth, its2ShannonPlotSite, 
                     its2ShannonPlotDepth, its2SimpsonPlotSite, 
                     its2SimpsonPlotDepth,legend.site, legend.depth, ncol = 2, 
                     nrow = 4, widths = c(3.5,3.4), heights = c(3,3,3.4,1),
                     labels = c("a","b", "c", "d", "e", "f","",""),
                     font.label = list(size = 18, color = "black", face =
                                       "bold")
                     )

divPlots

ggsave("its2_diversityPlots.eps", plot = divPlots, width = 8.25, height = 10, 
       unit = "in", dpi = 600)

##----------PCoA ON ALL SIGNIFICANT OTUS----------------------------------------
# first we need to set up data for PCoA analysis in R
# create a distance matrix with Bray-Curtis similarity using vegan

its2Dist = vegdist(its2Norm[, c(4:ncol(its2Norm))], method = "bray")
# write.csv(its2Dist, "its2Dist.csv") #saving the output matrix 

# this is the actual PCoA step, we want to calculate the eigenvalues 
# so later we can figure out % variation shown on each PCo
its2Mds = cmdscale(its2Dist, eig = TRUE, x.ret = TRUE)

# calculate how much variation each axis accounts for using Eigenvalues
its2Var = round(its2Mds$eig/sum(its2Mds$eig)*100, 1)

its2Var

# formatting the data so we can make it look nice with ggpPlot
its2Values = its2Mds$points

its2Values =as.data.frame(its2Values, sample = rownames(its2Values))
its2Pcoa = data.frame(sample = rownames(its2Values), 
	                     site = as.factor(its2Norm$Site),
                       depth = as.factor(its2Norm$Depth),
                       PCo1 = its2Values[,1],
                       PCo2 = its2Values[,2]
                       )

head(its2Pcoa)

# plotting PCoA with ggPlot
its2PcoaA = ggplot(its2Pcoa, aes(x = PCo1, y = PCo2)) +
        	  geom_hline(yintercept = 0, color = "gray90", size = 0.5)+
            geom_vline(xintercept = 0, color = "gray90", size = 0.5)+ 
            geom_point(aes(shape = factor(depth), size = depth, fill = site), 
                       color = "black") +
            scale_fill_manual(values = its2ColPal[c(8,1,3,9)], name = "Reef Site", 
            		              labels = c("Tobacco Reef", "Raph's Wall", 
            		              "South Reef", "Glover's Reef")) +
          	scale_shape_manual(values = c(24,23,22,21), name = "Depth",
          	                   labels = c("10 m", "16 m", "25 m", "35 m")) +
          	scale_size_manual(values = c(4, 4.75, 4.75, 4.75)) + 
          	guides(shape = guide_legend( override.aes = list(size = c(4, 4.75, 
          	       4.75, 4.75), fill = "white")), fill = 
          	       guide_legend(override.aes = list(shape = 22, size = 5.75, 
          	       color = "white")),	size = FALSE, color = FALSE)+
          	xlab(paste ("PCo 1 (", its2Var[1],"%)", sep = "")) +
          	ylab(paste ("PCo 2 (", its2Var[2],"%)", sep = "")) +
          	theme_bw()

its2Pcoa = its2PcoaA + 
            theme(axis.title.x = element_text(color = "black", size = 12),
                	axis.text.x = element_blank(),
                	axis.ticks.x = element_blank(),
                	axis.line.x = element_blank(),
                	axis.title.y = element_text(color = "black", size = 12),
                	axis.text.y = element_blank(),
                	axis.ticks.y = element_blank(),
                	axis.line.y = element_blank(),
                	legend.position = "right",
                	legend.title = element_text(color = "black", size = 12),
                	legend.text = element_text(color = "black", size = 12),
                	legend.key = element_blank(),
                	panel.border = element_rect(color = "black", size = 1.2),
                	panel.background = element_rect(fill = "white"),
                	plot.background = element_rect(fill = "white"),
                	panel.grid.major = element_blank(),
                	panel.grid.minor = element_blank()
                  )

its2PcoaA

# Save the generated PCoA plot
ggsave("ITS2_OTU_PCoA.eps", plot = its2PCoA, width = 10, height = 6.25, 
       dpi = 600, device="eps")	

##----------NMDS ON ITS2 OTU DATA-----------------------------------------------


its2Nmds = metaMDS(its2Norm[4:ncol(its2Norm)], try = 50)
its2Nmds

# Prepare the above data to plot with ggplot rather than base R
its2Scores = as.data.frame(scores(its2Nmds))
its2Scores$site = factor(its2Norm$Site)
its2Scores$depth = as.factor(its2Norm$Depth)
its2Scores$sample = row.names(its2Scores)
head(its2Scores)

its2Clades = as.data.frame(scores(its2Nmds,"species"))
its2Clades$seq = row.names(its2Clades)
its2Clades
its2Clades$seq

# Plot the data
its2NmdsPlotA = ggplot() + 
  	            geom_point(data = its2Scores, aes(x = NMDS1, y = NMDS2, 
  	                       shape = depth, size = depth, fill = site), 
  	                       color = "black") + # add the site points	
                scale_fill_manual(values = its2ColPal[c(8,5,3,9)], 
                                  name = "Reef Site", labels = c("Tobacco Reef",
                                  "Raph's Wall", "South Reef", "Glover's Reef")) +
              	scale_shape_manual(values = c(24, 23, 22, 21), name = "Depth", 
              	                   labels = c("10 m", "16 m", "25 m", "35 m")) +
              	scale_size_manual(values = c(3, 3.75, 3.75, 3.75)) + 
              	guides(shape = guide_legend( override.aes = list(size = 
              	       c(3, 3.75, 3.75, 3.75), fill = "white")), fill =
              		     guide_legend(override.aes = list(shape = 22, size = 3.75, 
              		     color = NA)), size = FALSE)+
              	geom_text(data = its2Clades, aes(x = NMDS1, y = NMDS2, label = seq), 
              		        color = "#000000", size = 4, fontface = "italic") +  
                          # add seq labels	
              	annotate("label", x = 1.25, y = 0.725, label = paste 
              		        ("Distance = Bray-Curtis\nStress = ", 
              		        round(its2Nmds$stress, 4), sep = ""), size = 4) +
                labs(x = "nMDS1", y = "nMDS2") +
              	coord_equal() +
              	theme_bw()
	
its2NmdsPlot = its2NmdsPlotA + 
               theme(axis.title.x = element_text(color = "black", size = 12),
              	     axis.text.x = element_blank(),
              	     axis.ticks.x = element_blank(),
              	     axis.title.y = element_text(color = "black", size = 12),
              	     axis.text.y = element_blank(),
              	     axis.ticks.y = element_blank(),
              	     legend.position = "right",
              	     legend.title = element_text(color = "black", size = 12),
              	     legend.text = element_text(color = "black", size = 12),
              	     legend.key = element_blank(),
              	     legend.background = element_blank(),
              	     panel.border = element_rect(color = "black", size = 1.2),
              	     panel.background = element_rect(fill = "white"),
              	     plot.background = element_blank()
                     )

its2NmdsPlot

# Saving nMDS plot
ggsave("nMDS_ITS2_OTU.eps", plot = its2NmdsPlot, width = 8, height = 5, dpi = 600, 
       device = "eps")

##----------REMOVE OUTLYING SAMPLES---------------------------------------------
## Samples with the C3g types skewing the data:
# 42, 80, 104, 98, 20, 57, 149, 152, 58, 163
# let's look at these samples specifically

itsOuts = its2Norm[its2Norm$sample %in% c("42", "80", "104", "98", "20", "57", 
                                          "149", "152", "58", "163"), ]
itsOuts$sums = apply(itsOuts[,4:13], 1, function(x) sum(x))
itsNoOut = its2Norm 
itsNoOut <- itsNoOut[!itsNoOut$sample %in% c("42", "80", "104", "98", "20", "57",
                                             "149", "152", "58", "163"), ]
itsNoOut$sums <- apply(itsNoOut[,4:13], 1, function(x) sum(x))

itsOuts
itsNoOut
max(itsOuts$sum)
max(itsNoOut$sum)
mean(itsNoOut$sum)
summary(itsOuts$sum)
summary(itsNoOut$sum)
# nothing unusual, neither the highest nor the lowest seq. depth of all samples

##----------NMDS WITHOUT OUTLYING SAMPLES---------------------------------------
# same idea as above, just with the new data lacking 10 C3g dominated samples

its2Nmds2 = metaMDS(itsNoOut[4:ncol(itsNoOut)], try = 50)
its2Nmds2

its2Scores2 = as.data.frame(scores(its2Nmds2))
its2Scores2$site = factor(itsNoOut$Site)
its2Scores2$depth = as.factor(itsNoOut$Depth)
its2Scores2$sample = row.names(its2Scores2)
head(its2Scores2)

its2Clades2 = as.data.frame(scores(its2Nmds2,"species"))
its2Clades2$seq = row.names(its2Clades2)
its2Clades2
its2Clades2$seq

# Plot the data
its2NmdsPlotB = ggplot() + 
                geom_point(data = its2Scores2, aes(x = NMDS1, y = NMDS2, 
                           shape = depth, size = depth, fill = site), 
                           color = "black") + # add the site points	
                scale_fill_manual(values = its2ColPal[c(8,1,3,9)], name = 
                                  "Reef Site",labels = c("Tobacco Reef", 
                                  "Raph's Wall", "South Reef", "Glover's Reef")) +
                scale_shape_manual(values = c(24, 23, 22, 21), name = "Depth", 
                                   labels = c("10 m", "16 m", "25 m", "35 m")) +
                scale_size_manual(values = c(3, 3.75, 3.75, 3.75)) + 
                guides(shape = guide_legend( override.aes = list(size = 
                       c(3, 3.75, 3.75, 3.75), fill = "white")), fill = 
                       guide_legend(override.aes = list(shape = 22, size = 3.75,
                       color = NA)), size = FALSE)+
                geom_text(data = its2Clades2, aes(x = NMDS1, y = NMDS2, 
                          label = seq), color = "#000000", size = 4, 
                          fontface = "italic") + 
                          # add seq labels	
                annotate("label", x = 0.95, y = 1.17, label = paste 
                         ("Distance = Bray-Curtis\nStress = ", 
                         round(its2Nmds2$stress, 4), sep = ""), size = 4) +
                labs(x = "nMDS1", y = "nMDS2") +
                coord_equal() +
                theme_bw()

its2NmdsPlot2 = its2NmdsPlotB + 
                theme(axis.title.x = element_text(color = "black", size = 12),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.y = element_text(color = "black", size = 12),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "right",
                      legend.title = element_text(color = "black", size = 12),
                      legend.text = element_text(color = "black", size = 12),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      panel.border = element_rect(color = "black", size = 1.2),
                      panel.background = element_rect(fill = "white"),
                      plot.background = element_blank()
                      )

its2NmdsPlot2

# saving new nMDS plot
ggsave("nMDS_ASVs_noOut.eps", plot = its2NmdsPlot2, width = 9.2, height = 5.75,
       dpi = 600, device = "eps")

##----------PLOTTING OTU RELATIVE ABUNDANCES------------------------------------
## PLOTTING OTU RELATIVE ABUNDANCES
# Calculate realtive abundances of each unique ITS2 

head(its2Norm)

its2NormPerc = its2Norm
its2NormPerc$sum = apply(its2NormPerc[, c(4:length(its2NormPerc[1, ]))], 1, 
                         function(x) 
                           {
	                         sum(x, na.rm = T)
                           }
                         )	

its2NormPerc = cbind(its2NormPerc[, c(1:3)], (its2NormPerc[, 
                     c(4:(ncol(its2NormPerc)-1))] / its2NormPerc$sum))
head(its2NormPerc)

# test that all are now 100% = 1
apply(its2NormPerc[, c(4:(ncol(its2NormPerc)))], 1, 
      function(x) 
        {
      	sum(x, na.rm = T)
        }
      )
  
#write.csv(its2NormPerc, "its2RelAbund.csv", row.names = FALSE)

# I added an additional column (in Excel) to sort better for the stacked barplot
# This was just a work around to get the facet_grid() to play nice with our data
# I added a coulumn "barPlotOrder" and for each population I filled in 1-16 for 
# each sample, so now there's no large blank expanses on the plot. The resulting
# file was named 'its2RelAbund2.csv' so it doesn't get overwritten if we run the 
# above steps again. I'll just load it in below and work on that file moving
# forward.

its2ra = read.csv("its2RelAbund2.csv")
rownames(its2ra) = its2ra$sample
colnames(its2ra)[5:ncol(its2ra)] = c("sq01_C3", "sq10_C3", "sq11_C3g",
                                     "sq14_C3g", "sq17_C3g", "sq18_C3g", 
                                     "sq25_C3g", "sq05_C3z", "sq07_C3e", 
                                     "sq08_C3g")
head(its2ra)

gss = otuStack(its2ra, count.columns = c(5:length(its2ra[1, ])),
               condition.columns = c(1:4))[1:2330,] # remove summ rows

levels(gss$otu)
gss$otu = factor(gss$otu, levels(gss$otu)[c(1, 2:8, 11, 9:10)])
gss$otu

levels(gss$Depth)
levels(gss$Depth) = c("10 m", "16 m", "25 m", "35 m")
levels(gss$Site)
levels(gss$Site) = c("Glover's Reef", "Raph's Wall", "South Reef", "Tobacco Reef")
gss$Site = factor(gss$Site, levels(gss$Site)[c(4, 2, 3, 1)])

OTUplotA = ggplot(gss, aes(x = barPlotOrder, y = count, fill = factor(otu))) +
          geom_bar(position = "stack", stat = "identity", color = "black",
                   size = 0.25) + 
        	ylab("Normalized proportion") +
        	scale_fill_manual(values=its2ColPal)+ 
        	labs(fill = "OTU_Clade type") +
	        facet_grid(Depth ~ Site, scales = "free_x") + 
	        theme_bw()
		  
OTUplot = OTUplotA +
           theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_text(color = "black", size = 12),
                 axis.text.y = element_text(color = "black", size = 12),
                 legend.position = "right",
                 legend.title = element_text(color = "black", size = 12),
                 legend.text = element_text(color = "black", size = 12),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_rect(fill = "white"),
                 plot.background = element_blank(),
                 strip.text.x = element_text(size = 12),
                 strip.text.y = element_text(size = 12),
                 strip.background = element_rect(fill = "white", size = 0.9)
                 )

OTUplot

# save the OTU barplot
ggsave("OTUbarPlot2.eps", plot = OTUplot, width = 8, height = 6, unit = "in", 
       dpi = 600)

##----------PERMANOVA ON OTU DATA-----------------------------------------------
# testing for dispersion

anova(betadisper(its2Dist, its2Norm$Depth))

#Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq  Mean Sq F value   Pr(>F)   
# Groups      3 0.5246 0.174860  5.1007 0.001971 **
#   Residuals 229 7.8505 0.034282                    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# dispersion is heteroschedastic, but PERMANOVA is robust to deviations in 
# homgeneity of variance

# Now let's see how different communities are from each other with PERMANOVA
its2Adonis = adonis2(its2Norm[, c(4:ncol(its2Norm))] ~ Depth*Site, 
                     data = its2Norm, permutations = 9999, method = "bray")
its2Adonis

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = its2Norm[, c(4:ncol(its2Norm))] ~ Depth * Site, 
# data = its2Norm,  permutations = 9999, method = "bray")
#             Df SumOfSqs      R2      F Pr(>F)   
# Depth        3   0.3946 0.04443 3.4741 0.0077 **
# Site         3   0.0349 0.00393 0.3073 0.9085   
# Depth:Site   9   0.2364 0.02662 0.6938 0.7970   
# Residual   217   8.2153 0.92502                 
# Total      232   8.8812 1.00000                 
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# see what is happening with depth using pairwise.adonis()
# these multiple comparisons are false dicovery rate (FDR) corrected 
its2PWAdonis = pairwise.adonis(its2Norm[, c(4:ncol(its2Norm))],
                  factors = its2Norm$Depth, sim.method = "bray", 
                  p.adjust.m = "BH", perm = 9999)
its2PWAdonis

#      pairs Df    SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1 10 vs 16  1 3.791992e-02 0.5212379 0.004473330  0.5145    0.61740    
# 2 10 vs 25  1 1.215016e-01 4.0807976 0.033703094  0.0002    0.00120   *
# 3 10 vs 35  1 1.209040e-01 3.9840789 0.033767936  0.0026    0.00680   *
# 4 16 vs 25  1 2.553940e-01 5.8424234 0.048347453  0.0034    0.00680   *
# 5 16 vs 35  1 2.551471e-01 5.7118059 0.048523645  0.0129    0.01935   .
# 6 25 vs 35  1 7.599081e-05 0.1802355 0.001592465  0.8289    0.82890  

##----------PERMANOVA WITHOUT 35 m SAMPLES--------------------------------------
goods2 = read.csv("goodsASVnoDeep.csv")
row.names(goods2) = goods2$sample

# what is the proportion of samples with data for these ASVs?
withData = apply(goods2[, 5:length(goods2[1, ])], 2, function(x) {
                 sum(x > 0)/length(x)
})
hist(withData, breaks = 50)

props = apply(goods2[, 5:length(goods2[1, ])], 2, function(x) 
              {
              sum(x)/sum(goods2[, 5:length(goods2[1, ])])
              }
              )

barplot(props, xaxt = "n", log = "y")
props

# Norm2alize otu counts with weighted trimmed mean of M-values (TMM; 
# Robinson and Oshlack 2010)
itsgoods2Transposed = t(goods2[, 5:length(goods2[1, ])])
itsgoods2List = DGEList(counts = itsgoods2Transposed)
head(itsgoods2List$samples)
its2Norm2 =  calcNormFactors(itsgoods2List, method = "TMM")
head(its2Norm2$samples)
its2TMM = t(cpm(its2Norm2, Normalized.lib.sizes = TRUE))
its2Norm2 = cbind(goods2[,c(2:4)], its2TMM)
head(its2Norm2)

#write.csv(its2Norm2, "its2_Norm2alizedgoods2.csv", row.names = FALSE)
# its2Norm2 = read.csv("its2_Normalizedgoods2.csv")
# renaming sequence coulmns to include Clade types identified with Blast-n
colnames(its2Norm2)[4:ncol(its2Norm2)] = c("sq01_C3", "sq10_C3", "sq11_C3g",
                                           "sq14_C3g", "sq17_C3g", "sq18_C3g",
                                           "sq25_C3g", "sq05_C3z", "sq07_C3e", 
                                           "sq08_C3g")
its2Norm2$Site = factor(its2Norm2$Site, levels(its2Norm2$Site)[c(4, 2, 3, 1)])
its2Norm2$Depth = as.factor(its2Norm2$Depth)
head(its2Norm2)
its2Norm2$Site
its2Norm2$Depth

its2Dist2 = vegdist(its2Norm2[, c(4:ncol(its2Norm2))], method = "bray")
anova(betadisper(its2Dist2, its2Norm2$Depth))
# Analysis of Variance Table
# 
# Response: Distances
# Df Sum Sq Mean Sq F value  Pr(>F)  
# Groups      2 0.3448 0.17240  3.9396 0.02122 *
#   Residuals 174 7.6143 0.04376                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

its2Adonis2 <- adonis2(its2Norm2[, c(4:ncol(its2Norm2))] ~ Depth*Site, 
                       data = its2Norm2, permutations = 9999, method = "bray")
its2Adonis2
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = its2Norm2[, c(4:ncol(its2Norm2))] ~ Depth * Site, 
#          data = its2Norm2, permutations = 9999, method = "bray")
#             Df SumOfSqs      R2      F Pr(>F)  
# Depth        2   0.2828 0.03326 2.9269 0.0388 *
# Site         3   0.0367 0.00432 0.2532 0.9348  
# Depth:Site   6   0.2111 0.02483 0.7284 0.677  
# Residual   165   7.9711 0.93759                
# Total      176   8.5017 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

its2PWAdonis2 = pairwise.adonis(its2Dist2, factors = its2Norm2$Depth,
                                perm = 9999, p.adjust.m = "BH")
its2PWAdonis2
#      pairs Df  SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1 10 vs 16  1 0.03164863 0.4486379 0.003852667  0.6011     0.6011    
# 2 10 vs 25  1 0.13184809 4.5463670 0.037404384  0.0001     0.0003  **
# 3 16 vs 25  1 0.26179196 6.1925677 0.051096926  0.0022     0.0033   *

##----------GENERALIZED LINEAR MIXED MODEL OF SIGNIFCANTLY CHANGING OTUS--------
# stacking the data table

its2Glm = its2Norm
colnames(its2Glm)[4:ncol(its2Glm)] = c("sq01_C3", "sq10_C3", "sq11_C3g", 
                                       "sq14_C3g", "sq17_C3g", "sq18_C3g", 
                                       "sq25_C3g", "sq05_C3z", "sq07_C3e", 
                                       "sq08_C3g")
glmStack=otuStack(data=its2Glm, count.columns = c(4:ncol(its2Glm)), 
                  condition.columns=c(1:3)) 

levels(glmStack$otu)
glmStack$otu = factor(glmStack$otu, levels(glmStack$otu)[c(1, 2:8, 11, 9:10)])
glmStack$otu
glmStack$count = round(glmStack$count, 0)
head(glmStack)

# fitting the model, using Depth as a fixed effect,
its2MmD = mcmc.otu(fixed = "Depth", data = glmStack,
                   nitt=100000, thin=25, burnin=5000)

# selecting the OTUs that were modeled reliably
acpassD = otuByAutocorr(its2MmD, glmStack)

# calculating effect sizes and p-values:
ssD = OTUsummary(its2MmD, glmStack, summ.plot = FALSE)

# correcting for mutliple comparisons (FDR) 
ssD = padjustOTU(ssD)

# getting significatly changing OTUs (FDR<0.05)
sigsD = signifOTU(ssD)

# plotting them
ss2D = OTUsummary(its2MmD, glmStack, otus = sigsD, whiskers = "sd", 
                  ptype = "mcmc")

its2MCMCglmD = ss2D$summary
levels(its2MCMCglmD$otu)
# write.csv(its2MCMCglmD, "its2MCMCglm_depth.csv", row.names = FALSE)
# its2MCMCglmD = read.csv("its2MCMCglm_depth.csv")
its2MCMCglmD$Depth = as.factor(its2MCMCglmD$Depth)
its2MCMCglmD$otu = factor(its2MCMCglmD$otu, levels = c("sq01_C3", "sq10_C3",
                          "sq14_C3g", "sq17_C3g", "sq18_C3g", "sq25_C3g", 
                          "sq08_C3g", "sq07_C3e"))
# displaying effect sizes and p-values for significant OTUs
ssD$otuWise[sigsD]

# plotting with custom colors in ggplot
pd = position_dodge(0.3)

its2glmPlotA = ggplot(its2MCMCglmD, aes(x = Depth, y = mean, group = otu, 
                      colour =otu))+
                      geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), 
                                    lwd = 0.6, width = 0.5, 	position = pd) +
                      geom_line(aes(group = otu), lwd = 0.6, position = pd)+
                      geom_point(aes(group = otu), position = pd, size = 2.5)+
                      scale_color_manual(values = its2ColPal[c(1,2,4:8,10)], 
                                         name = "OTU_Clade type") +
                      xlab("Depth")+
                      ylab("log10(proportion)") +
                      theme_bw()

its2glmPlotD = its2glmPlotA +	
               theme(axis.title.x = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12),
                     axis.title.y = element_text(color = "black", size = 12),
                     axis.text.y = element_text(color = "black", size = 12),
                     legend.position = "right",
                     legend.title = element_text(color = "black", size = 12),
                     legend.text = element_text(color = "black", size = 12),
                     legend.key = element_blank(),
                     legend.background = element_blank(),
                     panel.border = element_rect(color = "black", size = 1.2),
                     panel.background = element_rect(fill = "white"),
                     plot.background = element_blank()
                     )

its2glmPlotD

ggsave("its2glmDepth.eps", plot=its2glmPlotD, width = 8, height = 6, dpi = 600)
