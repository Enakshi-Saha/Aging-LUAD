# Load aging genes
tb <- read.csv("/home/esaha/aging_plots/limma_tables/GTEx_age_limma.txt", sep="")
increasing_gene = tb$gene_name[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = tb$gene_name[which(tb$P.Value<0.05 & tb$t<0)]


write.table(increasing_gene, file = "/home/esaha/aging_plots/limma_tables/GTEx_increasing_genes_0.05.txt", row.names = F, col.names = F, quote = F)
write.table(decreasing_gene, file = "/home/esaha/aging_plots/limma_tables/GTEx_decreasing_genes_0.05.txt", row.names = F, col.names = F, quote = F)
