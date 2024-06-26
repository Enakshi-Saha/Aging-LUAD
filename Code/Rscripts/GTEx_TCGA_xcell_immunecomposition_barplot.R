library(ggplot2)
library(data.table)

GTEx_xcell <- read.csv("~/aging_plots/xcell_immuneInfiltration/GTEx_immuneInfiltration_aging.csv")
TCGA_xcell <- read.csv("~/aging_plots/xcell_immuneInfiltration/TCGA_immuneInfiltration_aging.csv")

sig_cells = unique(c(GTEx_xcell$X[which(GTEx_xcell$p.value<0.05)],GTEx_xcell$X[which(GTEx_xcell$p.value<0.05)]))

age_coefficient = c(GTEx_xcell$age.coefficient, TCGA_xcell$age.coefficient)
sample_type = c(rep("GTEx healthy", nrow(GTEx_xcell)), rep("TCGA tumor", nrow(TCGA_xcell)))
sample_type = factor(sample_type, levels = c("TCGA tumor", "GTEx healthy"))
cell_type = c(GTEx_xcell$X, TCGA_xcell$X)
fulldata = data.frame(cell_type, sample_type, age_coefficient)
# Select cells that significantly change with age in either GTEx or TCGA
# fulldata = fulldata[which(fulldata$cell_type %in% sig_cells),]
head(fulldata)

ggplot(fulldata, aes(fill=sample_type, x=age_coefficient, y=cell_type)) + 
  geom_bar(position="dodge", stat="identity") + 
  labs(title="Change in Cell Composition with Age",y="cell type", x = "t-statistic of age coefficient") +
  geom_vline(xintercept=qnorm(0.025), linetype="dashed",color = "red") +
  geom_vline(xintercept=qnorm(0.975), linetype="dashed",color = "red")
  
  
