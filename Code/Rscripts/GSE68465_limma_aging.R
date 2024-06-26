#########################################
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
library(plyr)

# Load indegree
indegree_male = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_male.txt")
indegree_female = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_female.txt")
indegree = cbind(indegree_male, indegree_female[,-1])
head(indegree[,1:4])
genes = indegree$V1

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(genes, gene_info.TCGA$gene_name)]
chr_loc = gene_info.TCGA$seqnames[match(genes, gene_info.TCGA$gene_name)]

# Remove Y genes
indegree = indegree[-which(chr_loc == "chrY"),]
genes = genes[-which(chr_loc == "chrY")]
indegree = data.frame(indegree[,-1])
rownames(indegree) = genes

# Get clinical data
clinical = read.csv("~/validation_data/validation_newdata2023/GSE68465_phenotypes.txt", sep="")
clinical = clinical[match(colnames(indegree),clinical$geo_accession),]
head(clinical[,1:5])

gender = clinical$Sex.ch1
gender = factor(gender, levels = c("Male", "Female"))

race = clinical$race.ch1
race[which(race == "Not Reported")] = "Unknown"
race[-which(race %in% c("White", "Black or African American", "Unknown"))] = "others"
race = factor(race)

age = clinical$age.ch1

smoking = clinical$smoking_history.ch1
smoking[-which(smoking %in% c("Never smoked", "Currently smoking", "Smoked in the past"))] = "Unknown"
smoking = factor(smoking, levels = c("Never smoked", "Currently smoking", "Smoked in the past", "Unknown"))
table(smoking)

stage = clinical$disease_stage.ch1
tumor_stage = substr(stage, nchar(stage), nchar(stage))
tumor_stage = factor(tumor_stage, levels = c("1", "2", "3", "4", "p"))

#### Design Matrix
design = model.matrix(~ age + gender + tumor_stage + race + smoking)
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderfemale
tb = topTable(fit,coef="age",number=Inf)

############# GSEA ###################
# Rank genes in limma table
indegree_rank <- setNames(object=tb[,"t"], rownames(tb))
head(indegree_rank)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

fgseaRes[which(fgseaRes$pathway %in% c("KEGG_PATHWAYS_IN_CANCER", "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_P53_SIGNALING_PATHWAY", "KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_SMALL_CELL_LUNG_CANCER", "KEGG_NON_SMALL_CELL_LUNG_CANCER")),]
fgseaRes[c(grep("kine", fgseaRes$pathway, ignore.case = T),grep("immu", fgseaRes$pathway, ignore.case = T)),]

write.table(tb, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GSE68465_aging.txt", col.names = T, row.names = T)
save(fgseaRes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_GSE68465_aging.RData")
