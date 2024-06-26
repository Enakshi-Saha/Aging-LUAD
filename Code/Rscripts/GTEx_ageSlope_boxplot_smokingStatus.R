library(ggplot2)
library(fgsea)

tb <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_aging.txt", sep="")

tb_smoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_smoker_aging.txt", sep="")
tb_nonsmoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_nonsmoker_aging.txt", sep="")

increasing_gene = tb$gene_name[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = tb$gene_name[which(tb$P.Value<0.05 & tb$t<0)]

maintitle = "Gene score in GTEx (genes with increasing targeting with age)"
score_smoker = tb_smoker$t[which(tb_smoker$gene_name %in% increasing_gene)]
score_nonsmoker = tb_nonsmoker$t[which(tb_nonsmoker$gene_name %in% increasing_gene)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("Ever", length(score_smoker)), rep("Never", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$smoking_status = factor(stype, levels = c("Ever", "Never"))
ggplot(fulldata, aes(x = smoking_status, y = score, fill = smoking_status)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

wilcox.test(score_smoker, score_nonsmoker, alternative = "greater")  

maintitle = "Gene score in GTEx (genes with decreasing targeting with age)"
score_smoker = tb_smoker$t[which(tb_smoker$gene_name %in% decreasing_gene)]
score_nonsmoker = tb_nonsmoker$t[which(tb_nonsmoker$gene_name %in% decreasing_gene)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("Ever", length(score_smoker)), rep("Never", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$smoking_status = factor(stype, levels = c("Ever", "Never"))
ggplot(fulldata, aes(x = smoking_status, y = score, fill = smoking_status)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

wilcox.test(score_smoker, score_nonsmoker, alternative = "less")  


