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
genes = indegree_GTEx_male$V1

indegree_GTEx = cbind(indegree_GTEx_male[,-1], indegree_GTEx_female[,-1])
rownames(indegree_GTEx) = genes

GTEx_IDs = colnames(indegree_GTEx)
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree_GTEx), gene_info.TCGA$gene_id)]
chr_loc = gene_info.TCGA$seqnames[match(rownames(indegree_GTEx), gene_info.TCGA$gene_id)]

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

batch = factor(GTEx_phenotypes$batch)

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

rin = as.numeric(as.character(GTEx_phenotypes$rna_degrad))

###### Combined Indegree ######
indegree = indegree_GTEx
head(indegree[,1:4])

#### Design Matrix
design = model.matrix(~ age + gender + race + smoking_status + ischemic_timeH )
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderFEMALE
tb = topTable(fit,coef="age",number=Inf)
tb$chr = c()
tb$gene_name = c()
for (i in 1:nrow(tb)){
  index = match(rownames(tb)[i], gene_info.TCGA$gene_id)
  tb$chr[i]=gene_info.TCGA$seqnames[index]
  tb$gene_name[i]=gene_info.TCGA$gene_name[index]
}
head(tb)

############# GSEA ###################
# Rank genes in limma table
indegree_rank <- setNames(object=tb[,"t"], tb$gene_name)
head(indegree_rank)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

# Cancer Pathways
fgseaRes[which(fgseaRes$pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_PATHWAYS_IN_CANCER", "KEGG_SMALL_CELL_LUNG_CANCER",
                                       "KEGG_NON_SMALL_CELL_LUNG_CANCER", "KEGG_WNT_SIGNALING_PATHWAY", "KEGG_ERBB_SIGNALING_PATHWAY",
                                       "KEGG_P53_SIGNALING_PATHWAY"))]


# Immune Pathways

rbind(fgseaRes[grep("immun", fgseaRes$pathway, ignore.case = T),], fgseaRes[grep("kine", fgseaRes$pathway, ignore.case = T),])

sig_pathways = fgseaRes$pathway[which(fgseaRes$padj < 0.05)]

write.table(tb, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_aging.txt", col.names = T, row.names = T)
save(fgseaRes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_GTEx_aging.RData")
