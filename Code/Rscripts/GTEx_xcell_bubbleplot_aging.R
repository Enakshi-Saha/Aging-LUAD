# Must set seed
rm(list = ls())
set.seed(123)

# Load GTEx Male & Female Indegree

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("variancePartition")
#BiocManager::install("doParallel")

library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(recount)

# Load immune cell composition
GTEx_xcell <- read.csv("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/xcell/GTEx_xcell.txt", sep="")
rownames(GTEx_xcell) = GTEx_xcell$cell_type
GTEx_xcell = GTEx_xcell[,-1]

GTEx_IDs = colnames(GTEx_xcell)
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

# Get phenotypes: GTEx
GTEx_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_phenotypes.txt", sep=""))
GTEx_phenotypes = GTEx_phenotypes[match(GTEx_IDs, GTEx_phenotypes$SUBJID),]
GTEx_phenotypes = GTEx_phenotypes[,c("gtex.sex", "AGE", "RACE", "MHSMKSTS", "TRISCH", "gtex.smrin", "gtex.smnabtcht")]
colnames(GTEx_phenotypes) = c("gender", "age", "race", "smoking", "ischemic_time", "rna_degrad", "batch")

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("MALE", "FEMALE"))

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"

smoking_status = GTEx_phenotypes$smoking
smoking_status[-which(smoking_status %in% c("No", "Yes"))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})


# Linear model coefficient of each cell type
GTEx_lm = apply(GTEx_xcell, MARGIN = 1, function(x){summary(lm(x ~ age + gender + race + ischemic_timeH + smoking_status*gender))$coefficients[2,c(1,4)]})
rownames(GTEx_lm) = c("age_coefficient", "p-value")
GTEx_lm = t(GTEx_lm)
GTEx_lm

# Bubbleplot plot by group
# importing the ggplot2 library
library(ggplot2)
# creating the dataframe from the above columns
fulldata = data.frame(GTEx_lm)
fulldata = fulldata[which(fulldata$p.value < 0.05),]
fulldata$cell_type = factor(rownames(fulldata))
fulldata$trend = rep("decreasing", nrow(fulldata))
fulldata$trend[which(fulldata$age_coefficient >= 0)] = "increasing"
fulldata$trend = factor(fulldata$trend, levels = c("increasing", "decreasing"))
head(fulldata)

# Plot bubble plot
ggplot(fulldata, aes(x=p.value, y=cell_type, size=age_coefficient, color = trend))+
  geom_point(alpha=0.8) + labs(title="GTEx: Change in Cell Composition with Age",y="cell type", x = "pvalue")

# Save immune deconv table
immune_table = GTEx_lm
immune_table = format(immune_table, digits=4)
colnames(immune_table)[1] = "age coefficient"
write.csv(immune_table, "/home/esaha/aging_plots/aging_plots2023/GTEx_immuneInfiltration_aging.csv", row.names = T, quote = F)
