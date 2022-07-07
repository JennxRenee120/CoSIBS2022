library(VIM)
library(tidyverse)

## Data Cleaning ##

#Reading data
met_dat <- read.delim("../Metabolites.txt", sep = " ")
pro_dat <- read.delim("../Proteins.txt", sep = " ")

#knn imputation
met_dat_imputed <- kNN(data = met_dat, k = 10)
pro_dat_imputed <- kNN(data = pro_dat, k = 10)

#log base 2 normalization
met_dat_imp_norm <- log2(met_dat_imputed)
pro_dat_imp_norm <- log2(pro_dat_imputed)

#excluding boolean columns
new_met <- met_dat_imp_norm[, 1:64]
new_pro <- pro_dat_imp_norm[, 1:64]

## PCA Visualization ##

#Metabolome first two principal components
#No scaling
pr.met <- prcomp(new_met)
plot(pr.met$x[,1:2])

#Proteome first two principal components
#No scaling
pr.pro <- prcomp(new_pro)
plot(pr.pro$x[,1:2])

#Metabolome first two principal components
#Scaled
pr.met <- prcomp(new_met, scale. = T)
plot(pr.met$x[,1:2])

#Proteome first two principal components
#Scaled
pr.pro <- prcomp(new_pro, scale. = T)
plot(pr.pro$x[,1:2])

#Trying to combine proteome and metabolome

#pulling out first two PCs
pro_pc12 <- pr.pro$x[,1:2]
met_pc12 <- pr.met$x[,1:2]

#labeling whether protein or metabolite
pro_pc12 <- cbind(pro_pc12, rep("protein", length(pro_pc12[,1])))
met_pc12 <- cbind(met_pc12, rep("metabolite", length(met_pc12[,1])))

#combining data
pc12 <- rbind(pro_pc12, met_pc12)
ggplot(data.frame(pc12), mapping = aes(x = pc12[, 1],
                                  y = pc12[, 2], col = pc12[, 3])) + 
  geom_point()

##WTF??????