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
TCGA_xcell <- read.csv("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/xcell/TCGA_xcell.txt", sep="")
rownames(TCGA_xcell) = TCGA_xcell$cell_type
TCGA_xcell = TCGA_xcell[,-1]

# Load pathway scores
pathway_score_pca <- read.csv("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/pathway_score_pca.txt", sep="")
TCGA_xcell = TCGA_xcell[,match(colnames(pathway_score_pca), colnames(TCGA_xcell))]

sample.TCGA = colnames(TCGA_xcell)
# remove X from the starting of TCGA IDs
sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")] = substring(sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")], 2, length(sample.TCGA[which(substring(sample.TCGA, 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(sample.TCGA, split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))

TCGA_immune_tumor = TCGA_xcell

head(TCGA_immune_tumor[,1:4])

pathway_score_tumor = pathway_score_pca

# Correlation of each cell type with each pathway score
TCGA_cor_pathway = list()
TCGA_cor_pathway$cor = list()
TCGA_cor_pathway$cells = list()
for (i in 1:nrow(pathway_score_tumor)){
  pathway_name = rownames(pathway_score_tumor)[i]
  score = as.numeric(pathway_score_tumor[i,])
  cor_pathway = apply(TCGA_immune_tumor, MARGIN = 1, function(x){rbind(cor.test(x,score)$estimate,cor.test(x,score)$p.value)})
  rownames(cor_pathway) = c("correlation", "p-value")
  cor_pathway = t(cor_pathway)
  TCGA_cor_pathway$cells[[pathway_name]] = cor_pathway[which(cor_pathway[,2]<0.05),]
  TCGA_cor_pathway$cor[[pathway_name]] = cor_pathway
}

# Correlation of immune score with each pathway score
TCGA_cor_pathway = matrix(0, ncol = 2, nrow = nrow(pathway_score_tumor))
colnames(TCGA_cor_pathway) = c("correlation", "p-value")
rownames(TCGA_cor_pathway) = rownames(pathway_score_tumor)
for (i in 1:nrow(pathway_score_tumor)){
  pathway_name = rownames(pathway_score_tumor)[i]
  score = as.numeric(pathway_score_tumor[i,])
  immune_score = as.numeric(TCGA_immune_tumor["immune score",])
  cor_test = cor.test(immune_score,score)
  TCGA_cor_pathway[i,] = c(cor_test$estimate,cor_test$p.value)
}
TCGA_cor_pathway = data.frame(TCGA_cor_pathway)
save(TCGA_cor_pathway, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/TCGA_xcell_pathway_correlation.RData")

TCGA_fgsea = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_TCGA_aging.RData"))
sig_paths_TCGA = TCGA_fgsea$pathway[which(TCGA_fgsea$padj < 0.05)]

sig_paths = sig_paths_TCGA

tab = TCGA_cor_pathway[which(TCGA_cor_pathway$p.value < 0.05),]
sigtab = tab[which(rownames(tab) %in% sig_paths),]
tab_immune = sigtab
rownames(tab_immune) = stringr::str_to_sentence(gsub("KEGG_", "", rownames(tab_immune)))
tab_immune

write.csv(tab_immune, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/TCGA_xcell_immuneScore_pathway_correlation.csv", row.names = T, quote = F)
