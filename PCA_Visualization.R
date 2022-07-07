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

new_met <- knn_log2(met_dat)
new_pro <- knn_log2(pro_dat)

## PCA Visualization ##

## Metabolome first two principal components
#  No scaling
pr.met <- prcomp(new_met)
plot(pr.met$x[,1:2], main = "PC-1 vs PC-2 in Metabolome")

## Proteome first two principal components
#No scaling
pr.pro <- prcomp(new_pro)
plot(pr.pro$x[,1:2])

## Metabolome first two principal components
#Scaled
pr.met <- prcomp(new_met, scale. = T)
plot(pr.met$x[,1:2])

## Proteome first two principal components
#Scaled
pr.pro <- prcomp(new_pro, scale. = T)
plot(pr.pro$x[,1:2])

## Trying to combine proteome and metabolome

## pulling out first two PCs
pro_pc12 <- pr.pro$x[,1:2]
met_pc12 <- pr.met$x[,1:2]

## labeling whether protein or metabolite
pro_pc12 <- cbind(pro_pc12, rep("protein", length(pro_pc12[,1])))
met_pc12 <- cbind(met_pc12, rep("metabolite", length(met_pc12[,1])))

## combining data
pc12 <- rbind(pro_pc12, met_pc12)
ggplot(data.frame(pc12), mapping = aes(x = pc12[, 1],
                                  y = pc12[, 2], col = pc12[, 3])) + 
  geom_point()

##WTF??????