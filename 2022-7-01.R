library(VIM)
# knn imputation
# log base 2 normalization
met_dat <- read.delim("Desktop/COSIBS2022/Project/Metabolites.txt", sep = " ")
pro_dat <- read.delim("Desktop/COSIBS2022/Project/Proteins.txt", sep = " ")

met_dat_imputed <- kNN(data = met_dat, k = 10)
pro_dat_imputed <- kNN(data = pro_dat, k = 10)

met_dat_imp_norm <- log2(met_dat_imputed)
pro_dat_imp_norm <- log2(pro_dat_imputed)

## -- Spectral Clustering Prototype -- ##

new_met <- met_dat_imp_norm[, 1:64]
new_pro <- pro_dat_imp_norm[, 1:64]

library(kernlab)
sc_met <- specc(new_met, centers=2)
sc_pro <- specc(new_pro, centers=2)

table(sc_met, sc_pro)
