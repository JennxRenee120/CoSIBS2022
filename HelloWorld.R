
library(tidyverse)
library(VIM)
library(ggthemes)

theme_set(theme_fivethirtyeight())

metabolites <- read.delim("C:/Users/wz8878wh/Desktop/CoSIBS/Project/Data/Metabolites.txt", sep = " ")
protiens <- read.delim("C:/Users/wz8878wh/Desktop/CoSIBS/Project/Data/Proteins.txt", sep = " ")



# check missingness by column
missingness_metabolites <- data.frame(apply(metabolites, 2, FUN = function(x){sum(is.na(x))/length(x)}))
missingness_protiens <- data.frame(apply(protiens, 2, FUN = function(x){sum(is.na(x))/length(x)}))



# imputation
metabolites_imputed <- kNN(data = metabolites, k = 10)
protiens_imputed <- kNN(data = protiens, k = 10)


# log base 2 normalization
protiens_normalized <- log(protiens_imputed, base = 2)
metabolites_normalized <- log(metabolites_imputed, base = 2)


dim(metabolites_normalized)
(
  metabolites_normalized
  %>% rowid_to_column("ID")
  %>% mutate(ID = as.character(ID))
  %>% select(1:65)
  %>% gather(2:65, key = 'metabolites', value = 'value')
) -> metabolites_gathered


(ggplot(data = metabolites_gathered)
  + geom_boxplot(aes(x = ID, y = value))
  )
