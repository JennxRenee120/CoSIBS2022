
library(tidyverse)
library(VIM)

source("../CharlieCode/CoSIBS.R")

metabolites <- read.table("../Data/Metabolites.txt")
protiens <- read.table("../Data/Proteins.txt")



# check missingness by column
# missingness_metabolites <- data.frame(apply(metabolites, 2, FUN = function(x){sum(is.na(x))/length(x)}))
# missingness_protiens <- data.frame(apply(protiens, 2, FUN = function(x){sum(is.na(x))/length(x)}))



# imputation
met_imp <- kNN(data = metabolites, k = 10)
pro_imp <- kNN(data = protiens, k = 10)


# log base 2 normalization
pro <- log(pro_imp, base = 2)
met <- log(met_imp, base = 2)

# checking if we need quantile normalization
# dim(metabolites_normalized)
# 
# (
#   metabolites_normalized
#   %>% rowid_to_column("ID")
#   %>% mutate(ID = as.character(ID))
#   %>% select(-ends_with("_imp"))
#   %>% gather(-ID, key = 'metabolites', value = 'value')
# ) -> metabolites_gathered
# 
# 
# (
#   ggplot(data = metabolites_gathered)
#   + geom_boxplot(aes(x = ID, y = value))
#   + theme(axis.title.x=element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# )
# 
# 
# dim(protiens_normalized)
# 
# (
#   protiens_normalized
#   %>% rowid_to_column("ID")
#   %>% mutate(ID = as.character(ID))
#   %>% select(-ends_with("_imp"))
#   %>% gather(-ID, key = 'protiens', value = 'value')
# ) -> protiens_gathered
# 
# 
# (
#   ggplot(data = protiens_gathered)
#   + geom_boxplot(aes(x = ID, y = value))
#   + theme(axis.title.x=element_blank(),
#          axis.text.x=element_blank(),
#          axis.ticks.x=element_blank())
# )

