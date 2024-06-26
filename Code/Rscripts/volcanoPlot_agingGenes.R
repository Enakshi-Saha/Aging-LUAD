library(ggplot2)
library(ggrepel)

# Load limma tables
tb_GTEx = read.csv("~/aging_plots/limma_tables/GTEx_age_limma.txt", sep="")
tb_TCGA = read.csv("~/aging_plots/limma_tables/TCGA_age_limma.txt", sep="")

tb = tb_GTEx
# cutoff for logFC
FClim = 10
# add a column of NAs
tb$direction <- "Not significant"
# if log2Foldchange > 25 and -log10(pvalue) > 5, set as "UP" 
tb$direction[tb$logFC > FClim & tb$P.Value < 0.05] <- "increasing"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
tb$direction[tb$logFC < -FClim & tb$P.Value < 0.05] <- "decreasing"

mycolors <- c("red", "blue", "grey")
names(mycolors) <- c("increasing", "decreasing", "Not significant")

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
tb$delabel <- NA
tb$delabel[tb$direction != "Not significant"] <- tb$gene_name[tb$direction != "Not significant"]

# Remove gene names that copntain "-" or "." in name
tb$delabel[grep("-",tb$delabel)] <- NA
tb$delabel[grep(".",tb$delabel, fixed = T)] <- NA

# Re-plot but this time color the points with "diffexpressed"
ggplot(data=tb, aes(x=logFC, y=-log10(P.Value), col=direction, label = delabel)) +
  geom_point() + theme_minimal() +
  geom_vline(xintercept=c(-FClim, FClim), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") +
  scale_colour_manual(values = mycolors) +
  ggtitle("GTEx: Volcanoplot of Differentially Targeted Genes with Age") +
  geom_text_repel(size = 3)

