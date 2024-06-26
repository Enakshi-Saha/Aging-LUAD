library(ggplot2)

tb <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_aging.txt", sep="")

all_cosmic = read.csv("/home/ubuntu/cancer_genes_cosmic_2023.txt")
oncogene = all_cosmic$Gene.Symbol[which(all_cosmic$Role.in.Cancer == "oncogene")]
TSG = all_cosmic$Gene.Symbol[which(all_cosmic$Role.in.Cancer == "TSG")]

maintitle = "Cosmic Gene Targeting: t-statistic (age coefficient in GTEx)"
score_tsg = tb$t[which(tb$gene_name %in% TSG)]
score_onco = tb$t[which(tb$gene_name %in% oncogene)]
score_others = tb$t[-which(tb$gene_name %in% all_cosmic$Gene.Symbol)]
score = c(score_tsg, score_onco, score_others)
stype = c(rep("TSG", length(score_tsg)), rep("Oncogene", length(score_onco)), rep("Non-cancer gene", length(score_others)))
fulldata = data.frame(score)
fulldata$sample_type = stype
ggplot(fulldata, aes(x = sample_type, y = score)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed", 
                                                                                                                                                                     