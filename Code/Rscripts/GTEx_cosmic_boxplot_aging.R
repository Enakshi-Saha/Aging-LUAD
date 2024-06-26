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
                                                                                                                                                                     color = "red") 
#################################
wilcox.test(score_tsg, alternative = "greater")
wilcox.test(score_onco, alternative = "greater")
wilcox.test(score_others, alternative = "greater")

wilcox.test(score_tsg, score_others, alternative = "greater")
wilcox.test(score_onco, score_others, alternative = "greater")

###########################################################
##### COSMIC Gene Boxplot by Smoking Status ###############

library(ggplot2)

tb_smoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_smoker_aging.txt", sep="")
tb_nonsmoker <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_nonsmoker_aging.txt", sep="")

maintitle = "Oncogene Targeting in GTEx: t-statistic"
score_smoker = tb_smoker$t[which(tb_smoker$gene_name %in% oncogene)]
score_nonsmoker = tb_nonsmoker$t[which(tb_nonsmoker$gene_name %in% oncogene)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("smoker", length(score_smoker)), rep("nonsmoker", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$sample_type = factor(stype, levels = c("smoker", "nonsmoker"))
ggplot(fulldata, aes(x = sample_type, y = score, fill = sample_type)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

pval_oncogene = wilcox.test(score_smoker, score_nonsmoker, alternative = "greater")                                     
                                                                                                                                
maintitle = "TSG Targeting in GTEx: t-statistic"
score_smoker = tb_smoker$t[which(tb_smoker$gene_name %in% TSG)]
score_nonsmoker = tb_nonsmoker$t[which(tb_nonsmoker$gene_name %in% TSG)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("smoker", length(score_smoker)), rep("nonsmoker", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$sample_type = factor(stype, levels = c("smoker", "nonsmoker"))
ggplot(fulldata, aes(x = sample_type, y = score, fill = sample_type)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

pval_tsg = wilcox.test(score_smoker, score_nonsmoker, alternative = "greater")  

maintitle = "Non-cancer gene Targeting in GTEx: t-statistic"
score_smoker = tb_smoker$t[-which(tb_smoker$gene_name %in% all_cosmic$Gene.Symbol)]
score_nonsmoker = tb_nonsmoker$t[-which(tb_nonsmoker$gene_name %in% all_cosmic$Gene.Symbol)]
score = c(score_smoker, score_nonsmoker)
stype = c(rep("smoker", length(score_smoker)), rep("nonsmoker", length(score_nonsmoker)))
fulldata = data.frame(score)
fulldata$sample_type = factor(stype, levels = c("smoker", "nonsmoker"))
ggplot(fulldata, aes(x = sample_type, y = score, fill = sample_type)) + geom_boxplot() + theme_bw() + ggtitle(maintitle) + xlab("") + ylab("t-statistic (age coefficient)") + geom_hline(yintercept=0, linetype="dashed",color = "red")  

pval_noncancer = wilcox.test(score_smoker, score_nonsmoker, alternative = "greater")  
