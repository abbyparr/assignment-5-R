library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

#A load the data
read_data1 <- read_excel("Gene_Expression_Data.xlsx")
read_data2 <- read_csv("Gene_Information.csv")
read_data3 <- read.table("Sample_Information.tsv")

#B change the sample names 
change <- merge(read_data1, read_data3, by = "SampleID")
change$SampleID <- change$Phenotype
print(change)

#C split the merged data 
split_merge <- split(read_data1, read_data1$Phenotype)
tumor <- split_merge[['tumor']]
normal <- split_merge[['normal']]

#D compute the average expression
avg_tumor <- rowMeans(tumor[, -1])
avg_normal <- rowMeans(normal[, -1])
print(avg_tumor)
print(avg_normal)

#E determine fold change 
fold_change <- log2((avg_tumor - avg_normal)/ avg_normal)
print(fold_change)

#F identify fold change magnitude 
read_data2 <- read.csv("Gene_Information.csv")
read_data1<- merge(read_data1, fold_change)
identify <- merge[abs(merge$fold_change) >5]

#G add a column 
exp <- ifelse(avg_normal > avg_tumor, "normal", "tumor")
newdata$exp <- exp

#Part 2
#A perform EDA 
head(read_data1)
head(read_data2)
head(read_data3)
summary(read_data1)
summary(read_data2)
summary(read_data3)
dim(read_data1)
dim(read_data2)
dim(read_data3)
cor(read_data1)
cor(read_data2)
cor(read_data3)

#B create a histogram 
hist(read_data1, main = 'DEGs by chromosome', xlab = 'DEGs', ylab = 'chromosome')

#C make another histogram 
hist(read_data3, main = 'DEGs by chromosome segregated by sample type', xlab = 'DEGs by chromosome', ylab = 'sample type(normal or tumor)')

#D create a bar chart 
barplot(read_data3, main = 'Percentage of upregulated and downregulated DEGs', xlab ='DEGs', ylab = 'Upregulated or downregulated')

#E create a heatmap 
#gene expression by sample
heatmap(read_data1, scale = 'row', main="gene expression by sample")

#F create a clustermap 
#visualizing gene expression by sample
clustermap(read_data1, main= "visualizing gene expression by sample")

#E analysis: 
#While analyzing the code, I noticed that DEGs by chromosome have a wide distribution. I also noticed that upregulation and downregulation in tumor samples also has a wide distribution and effect. 
