#### Heatmap 1: GTEx vs LGRC ####

GTEx = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_GTEx_aging.RData"))
LGRC = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_LGRC_aging.RData"))
LGRC = LGRC[match(GTEx$pathway, LGRC$pathway),]

tumor = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_GTEx_TCGA_aging.RData"))
tumor = tumor[match(GTEx$pathway, tumor$pathway),]

tab = cbind(GTEx$NES, LGRC$NES, tumor$NES)

rownames(tab) = GTEx$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((GTEx$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("GTEx", "LGRC", "TCGA vs GTEx")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/aging_plots/aging_plots2023/heatmaps/GTEx_LGRC_tumor_aging_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 45, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

###################################################
#### Heatmap 1: TCGA vs GSE68465 ####

TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_TCGA_aging.RData"))
GSE68465 = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_GSE68465_aging.RData"))
GSE68465 = GSE68465[match(TCGA$pathway, GSE68465$pathway),]

tab = cbind(TCGA$NES, GSE68465$NES)

rownames(tab) = TCGA$pathway
head(tab)

# Significant Pathways
sig_pathways = rownames(tab)[which((TCGA$padj<0.05))]
tab_subset = tab[which(rownames(tab) %in% sig_pathways),]

colnames(tab_subset) = c("TCGA", "GSE68465")
rownames(tab_subset) = stringr::str_to_title(lapply(strsplit(rownames(tab_subset), split = "_"), function(x){paste(x[-1], collapse = " ")}))

# Make heatmap comparing all groups

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/aging_plots/aging_plots2023/heatmaps/TCGA_GSE68465_tumor_aging_heatmap.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()


# Heatmap of pathways validated in GSE68465
tab_subset0 = tab_subset[which(sign(tab_subset[,1]*tab_subset[,2]) == 1),]

library(ggplot2)
library(gplots)
library(RColorBrewer)

mycol2 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(50)
pdf("/home/esaha/aging_plots/aging_plots2023/heatmaps/TCGA_GSE68465_tumor_aging_heatmap_validatedPathsOnly.pdf",h=25,w=16)
heatmap_tcga_stages = heatmap.2(as.matrix(tab_subset0),density.info="none",trace="none",col=mycol2,symbreaks=T,symkey=T, 
                                cexRow=1.5, cexCol=1.5, srtCol = 0, keysize=0.5, mar=c(20,50), key.title=NULL, key.xlab="NES", Colv = F)
dev.off()

validated_paths = rownames(tab_subset0)
write.table(validated_paths, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/TCGA_sigPaths_validatedPaths.txt", col.names = F, row.names = F)

