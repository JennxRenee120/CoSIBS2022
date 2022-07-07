#Project 6

#Libraries
library(VIM)
library(FCPS)

#Reading Datasets
MData = read.delim("Metabolites.txt", sep = " ")
PData = read.delim("Proteins.txt", sep = " ")
#Finding Amount of Data Missing for Metabolites Data
MissingM = numeric(length(MData))

for(i in 1:length(MData)){
  MissingM[i] = sum(is.na(MData[,i]))/107
}
summary(MissingM)

#Finding Amount of Data Missing for Protein Data
MissingP = numeric(length(PData))
for(i in 1:length(PData)){
  MissingP[i] = sum(is.na(PData[,i]))/107
}
summary(MissingP)

#KNN Imputation for Missingness
MDataFull = kNN(MData, k = 10)
PDataFull = kNN(PData, k = 10)
str(MDataFull)

#Normalization of Data
MDataNorm = log2(MDataFull[,1:64])
PDataNorm = log2(PDataFull[,1:90])

#Quantile Normalization


#Boxplots
par(mfrow = c(1,2))
boxplot(MDataNorm)
boxplot(PDataNorm)

#Spectral Clustering 
SpectralClustering(MDataNorm, ClusterNo = 3, PlotIt = TRUE)
SpectralClustering(PDataNorm, ClusterNo = 3, PlotIt = TRUE)
SpectralClustering(c(MDataNorm, PDataNorm), ClusterNo = 3, Plotit= TRUE)


