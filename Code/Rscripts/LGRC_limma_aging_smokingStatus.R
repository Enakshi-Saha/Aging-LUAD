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
indegree = fread("/home/esaha/Lung_lioness/validation_networks/GSE47460/GSE47460_inDegree.txt")
head(indegree[,1:4])
genes = indegree$V1

# Get chromosome information
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
names(gene_info.TCGA)
head(gene_info.TCGA)
gene_name = gene_info.TCGA$gene_name
chr_loc = rep(0, length(genes))
for (i in 1:length(genes)){
  index = which(gene_name == genes[i])[1]
  chr_loc[i]=as.character(gene_info.TCGA$seqnames[index])
}
chr_loc[1:100]

# Remove Y genes if any
length(which(chr_loc == "chrY"))
genes = indegree$V1[-which(chr_loc == "chrY")]
indegree = indegree[-which(chr_loc == "chrY"),-1]
rownames(indegree) = genes

# Load phenotypes
phenotypes = read.csv("/home/esaha/validation_data/GSE47460/GSE47460_phenotypes.txt", sep="")
phenotypes = phenotypes[match(colnames(indegree), rownames(phenotypes)),]

age = phenotypes$age
gender = factor(phenotypes$sex, levels = c("Male", "Female"))
smoking = phenotypes$smoking_status
smoking[is.na(smoking)] = "NA"
smoking = factor(smoking, levels = c("Never", "Ever", "NA"))

design = model.matrix(~ age*smoking + gender + smoking)

fit <- lmFit(indegree, design)
fit <- eBayes(fit)
rownames(fit$coefficients) = rownames(indegree)

# Save table for sex difference: genderFEMALE
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

sig_pathways = fgseaRes$pathway[which(fgseaRes$padj < 0.05)]

write.table(tb, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_LGRC_nonsmoker_aging.txt", col.names = T, row.names = T)
save(fgseaRes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_LGRC_nonsmoker_aging.RData")

#######################################
colnames(design) = make.names(colnames(design))
contrast.matrix <- makeContrasts(age_smoker = age + age.smokingEver, levels=design)
colnames(fit$coefficients) = rownames(contrast.matrix)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
tb = topTable(fit2,coef="age_smoker",number=Inf)

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

sig_pathways = fgseaRes$pathway[which(fgseaRes$padj < 0.05)]

write.table(tb, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_LGRC_smoker_aging.txt", col.names = T, row.names = T)
save(fgseaRes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_LGRC_smoker_aging.RData")

