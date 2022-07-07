#### Visualizing PCA in Metabolome and Proteome

library(VIM)
library(tidyverse)

## Data Cleaning ##

## Reading data
met_dat <- read.delim("Metabolites.txt", sep = " ")
pro_dat <- read.delim("Proteins.txt", sep = " ")

## Imputing values with knn, removing booleans, and log2 transform

knn_log2 <- function(dat){
  dat %>%
    kNN(k = 10) %>%
    select(1:64) %>%
    log2()
}

new_met <- knn_log2(met_dat); new_pro <- knn_log2(pro_dat)

## PCA Visualization ##

first2pc <- function(dat, omic_field, scaled){
  plot(prcomp(dat)$x[,1:2], 
       main = paste0("PC-1 vs PC-2 in ", omic_field, "\nScaled = ", scaled), 
       scale. = scaled)
}

par(mfrow = c(2,2))

first2pc(new_met, "Metabolome", T); first2pc(new_met, "Metabolome", F)
first2pc(new_pro, "Proteome", T); first2pc(new_pro, "Proteome", F)

## Trying to combine proteome and metabolome

## pulling out first two PCs
pro_pc12 <- prcomp(new_pro, scale. = T)$x[,1:2]
met_pc12 <- prcomp(new_met, scale. = T)$x[,1:2]

## labeling whether protein or metabolite
pro_pc12 <- cbind(pro_pc12, rep("protein", length(pro_pc12[,1])))
met_pc12 <- cbind(met_pc12, rep("metabolite", length(met_pc12[,1])))

## combining data
pc12 <- as.data.frame(rbind(pro_pc12, met_pc12))
pc12$PC1 <- as.numeric(pc12$PC1)
pc12$PC2 <- as.numeric(pc12$PC2)

ggplot(pc12, mapping = aes(x = PC1,
                           y = PC2, 
                           col = pc12[, 3])) + geom_point()