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

# Load indegree matrices
indegree_GTEx_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_male.txt"))
indegree_GTEx_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_female.txt"))

# Load indegree matrices
indegree_TCGA_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_TCGA_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

###### softplus indegrees
# indegree_GTEx_male = data.frame(fread("/home/esaha/Lung_lioness/softplus_indegree/softmax_indegree_tpm_GTEx_male.txt"))
# indegree_GTEx_female = data.frame(fread("/home/esaha/Lung_lioness/softplus_indegree/softmax_indegree_tpm_GTEx_female.txt"))
# indegree_TCGA_male = data.frame(fread("/home/esaha/Lung_lioness/softplus_indegree/softmax_indegree_tpm_TCGA_male.txt"))
# indegree_TCGA_female = data.frame(fread("/home/esaha/Lung_lioness/softplus_indegree/softmax_indegree_tpm_TCGA_female.txt"))
# indegree_TCGA_female[4753, 127] = 97580

head(indegree_GTEx_male[,1:4])
head(indegree_GTEx_female[,1:4])
head(indegree_TCGA_male[,1:4])
head(indegree_TCGA_female[,1:4])


# Check if male and female edges are in the same order
sum(indegree_GTEx_male$V1 != indegree_GTEx_female$V1)
sum(indegree_GTEx_male$V1 != indegree_TCGA_male$V1)
sum(indegree_GTEx_female$V1 != indegree_TCGA_female$V1)

genes = indegree_GTEx_male$V1

sample_male.GTEx = colnames(indegree_GTEx_male)[-1]
sample_female.GTEx = colnames(indegree_GTEx_female)[-1]
sample_male.TCGA = colnames(indegree_TCGA_male)[-1]
sample_female.TCGA = colnames(indegree_TCGA_female)[-1]

c(length(sample_male.GTEx), length(sample_male.TCGA), 
  length(sample_female.GTEx), length(sample_female.TCGA))

indegree_male.GTEx = indegree_GTEx_male[,-1]
indegree_male.TCGA = indegree_TCGA_male[,-1]
indegree_female.GTEx = indegree_GTEx_female[,-1]
indegree_female.TCGA = indegree_TCGA_female[,-1]

indegree.GTEx = cbind(indegree_male.GTEx, indegree_female.GTEx)
indegree.TCGA = cbind(indegree_male.TCGA, indegree_female.TCGA)

sample.GTEx = c(sample_male.GTEx, sample_female.GTEx)
sample.TCGA = c(sample_male.TCGA, sample_female.TCGA)

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
names(gene_info.TCGA)
head(gene_info.TCGA)
gene_name = rep(0, length(genes))
chr_loc = rep(0, length(genes))
for (i in 1:length(genes)){
  index = which(gene_info.TCGA$gene_id == genes[i])[1]
  gene_name[i]=gene_info.TCGA$gene_name[index]
  chr_loc[i]=as.character(gene_info.TCGA$seqnames[index])
}
gene_name[1:100]
chr_loc[1:100]

# Don't remove Y genes in aging analysis
rownames(indegree.GTEx) = genes
rownames(indegree.TCGA) = genes

# Format male & female IDs to match with phenotypic data: GTEx
GTEx_IDs = sample.GTEx
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

# Get phenotypes: GTEx
# Get phenotypes: GTEx
GTEx_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_phenotypes.txt", sep=""))
GTEx_phenotypes = GTEx_phenotypes[match(GTEx_IDs, GTEx_phenotypes$SUBJID),]
GTEx_phenotypes = GTEx_phenotypes[,c("gtex.sex", "AGE", "RACE", "MHSMKSTS", "TRISCH", "gtex.smrin", "gtex.smnabtcht")]
colnames(GTEx_phenotypes) = c("gender", "age", "race", "smoking", "ischemic_time", "rna_degrad", "batch")
head(GTEx_phenotypes)

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
GTEx_gender = factor(gender, levels = c("MALE", "FEMALE"))

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)
GTEx_age = age

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"
GTEx_race <- as.factor(race)

smoking_status = GTEx_phenotypes$smoking
smoking_status[-which(smoking_status %in% c("Yes", "No"))] = "Unknown"
GTEx_smoking_status = smoking_status

# remove X from the starting of TCGA IDs
sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")] = substring(sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")], 2, length(sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(sample.TCGA, split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))

# Get phenotypes: TCGA
TCGA_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", sep=""))
TCGA_phenotypes = TCGA_phenotypes[match(TCGA_IDs, rownames(TCGA_phenotypes)),]
TCGA_phenotypes = TCGA_phenotypes[,c("tcga.gdc_cases.demographic.gender", "tcga.gdc_cases.demographic.race", 
                                     "tcga.xml_age_at_initial_pathologic_diagnosis", "tcga.gdc_cases.samples.sample_type", 
                                     "tcga.gdc_cases.diagnoses.tumor_stage", "tcga.xml_tobacco_smoking_history")]
colnames(TCGA_phenotypes) = c("gender", "race", "age", "sample_type", "tumor_stage", "smoking_status")
head(TCGA_phenotypes)
dim(TCGA_phenotypes)

# Define the covariates: gender, age, race, smoking, tumor_stage
TCGA_phenotypes$gender[which(TCGA_phenotypes$gender == "male")] = "MALE" 
TCGA_phenotypes$gender[which(TCGA_phenotypes$gender == "female")] = "FEMALE" 
TCGA_phenotypes$gender = factor(TCGA_phenotypes$gender, levels = c("MALE", "FEMALE"))

TCGA_phenotypes$race[which(TCGA_phenotypes$race != "black or african american" & TCGA_phenotypes$race != "white")] = "others"
TCGA_phenotypes$race = factor(TCGA_phenotypes$race)

TCGA_phenotypes$age <- as.numeric(as.character(TCGA_phenotypes$age))
TCGA_phenotypes$age[which(is.na(TCGA_phenotypes$age))] = mean(TCGA_phenotypes$age,na.rm=TRUE)

TCGA_phenotypes$smoking_status[which(is.na(TCGA_phenotypes$smoking_status))] = "Unknown"
TCGA_phenotypes$smoking_status[which(TCGA_phenotypes$smoking_status == 1)] = "No"
TCGA_phenotypes$smoking_status[which(TCGA_phenotypes$smoking_status %in% 2:5)] = "Yes"

###### Combined Phenotypes ######
gender = c(GTEx_gender, TCGA_phenotypes$gender)
gender = factor(gender, levels = c("MALE", "FEMALE"))
age = c(GTEx_age, TCGA_phenotypes$age)
race = c(GTEx_race, TCGA_phenotypes$race)
smoking = c(GTEx_smoking_status, TCGA_phenotypes$smoking_status)
smoking = factor(smoking, levels = c("No", "Yes", "Unknown"))
sample_type = factor(c(rep("healthy", length(GTEx_gender)), rep("tumor", length(TCGA_phenotypes$gender))))
sample_type = factor(sample_type, levels = c("healthy", "tumor"))

age_binary = age
age_binary[which(age<=55)] = "below55"
age_binary[-which(age<=55)] = "above55"
age_binary = factor(age_binary, levels = c("below55", "above55"))

###### Combined Indegree ######
indegree = cbind(indegree.GTEx, indegree.TCGA)
head(indegree[,1:4])

#### Design Matrix
design = model.matrix(~ age + gender + race + sample_type*age_binary + smoking)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Load aging genes
tb <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_aging.txt", sep="")
increasing_gene = tb$gene_name[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = tb$gene_name[which(tb$P.Value<0.05 & tb$t<0)]

# Save table for sex difference: sample_type tumor
tb = topTable(fit,coef="sample_typetumor",number=Inf)
tb$chr = c()
tb$gene_name = c()
for (i in 1:nrow(tb)){
  index = match(rownames(tb)[i], gene_info.TCGA$gene_id)
  tb$chr[i]=gene_info.TCGA$seqnames[index]
  tb$gene_name[i]=gene_info.TCGA$gene_name[index]
}
head(tb)

interaction_increasingGenes = tb$t[which(tb$gene_name %in% increasing_gene)]
interaction_decreasingGenes = tb$t[which(tb$gene_name %in% decreasing_gene)]

boxplot(interaction_increasingGenes, interaction_decreasingGenes, main = "below 55", names = c("increasing", "decreasing"))
