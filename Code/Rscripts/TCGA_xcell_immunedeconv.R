# remotes::install_github("grst/immunedeconv")

library(recount)
library(recount3)
library(immunedeconv)

# Gene expression should be TPM, not log transformed, row names should be HGNC symbols
TCGA_tpm = read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_TPM.txt", sep="")
# Transform data to non log scale
TCGA_tpm = 2^(TCGA_tpm)-1

info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info = data.frame(info.TCGA@rowRanges)
gene_names = gene_info$gene_name[match(rownames(TCGA_tpm), gene_info$gene_id)]

# For rows with same gene name, take the mean as expression
TCGA_TPM = aggregate(TCGA_tpm, by=list(Category=gene_names), FUN=mean)
rownames(TCGA_TPM) = TCGA_TPM$Category
TCGA_TPM = TCGA_TPM[,-1]

xcell_tcga <- data.frame(immunedeconv::deconvolute(TCGA_TPM, "xcell"))

write.table(xcell_tcga, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/xcell/TCGA_xcell.txt", row.names = T, col.names = T)
