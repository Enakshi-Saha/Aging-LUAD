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

# Load indegree matrices
indegree_TCGA_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_TCGA_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

genes = indegree_TCGA_male$V1

indegree_TCGA = cbind(indegree_TCGA_male[,-1], indegree_TCGA_female[,-1])
rownames(indegree_TCGA) = genes

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree_TCGA), gene_info.TCGA$gene_id)]
chr_loc = gene_info.TCGA$seqnames[match(rownames(indegree_TCGA), gene_info.TCGA$gene_id)]

# remove X from the starting of TCGA IDs
colnames(indegree_TCGA)[which(substring(colnames(indegree_TCGA), 1,1) == "X")] = substring(colnames(indegree_TCGA)[which(substring(colnames(indegree_TCGA), 1,1) == "X")], 2, length(colnames(indegree_TCGA)[which(substring(colnames(indegree_TCGA), 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(colnames(indegree_TCGA), split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))

# Get phenotypes: TCGA
TCGA_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", sep=""))
TCGA_phenotypes = TCGA_phenotypes[match(TCGA_IDs, rownames(TCGA_phenotypes)),]
TCGA_phenotypes = TCGA_phenotypes[,c("tcga.gdc_cases.demographic.gender", "tcga.gdc_cases.demographic.race", 
                                     "tcga.xml_age_at_initial_pathologic_diagnosis", "tcga.gdc_cases.samples.sample_type", 
                                     "tcga.gdc_cases.diagnoses.tumor_stage", "tcga.xml_tobacco_smoking_history")]
colnames(TCGA_phenotypes) = c("gender", "race", "age", "sample_type", "tumor_stage", "smoking_status")
head(TCGA_phenotypes)
dim(TCGA_phenotypes)

# Define the covariates: gender, age, race, smoking, tumor_stage
gender = TCGA_phenotypes$gender
gender[which(gender == "male")] = "MALE" 
gender[which(gender == "female")] = "FEMALE" 
gender = factor(gender, levels = c("MALE", "FEMALE"))

race = TCGA_phenotypes$race
race[which(race != "black or african american" & race != "white")] = "others"
race = factor(race)

age <- as.numeric(as.character(TCGA_phenotypes$age))
age[which(is.na(age))] = mean(age,na.rm=TRUE)

# Get 4 levels of smoking
smoking_4levels = TCGA_phenotypes$smoking_status
smoking_4levels[which(is.na(smoking_4levels))] = "Unknown"
smoking_4levels[which(smoking_4levels == 1)] = "nonsmoker"
smoking_4levels[which(smoking_4levels == 2)] = "current smoker"
smoking_4levels[which(smoking_4levels == 3)] = "past smoker (>15)"
smoking_4levels[which(smoking_4levels == 4)] = "past smoker (<=15)"
smoking_4levels[which(smoking_4levels == 5)] = "past smoker (unspecificed)"
smoking_4levels = factor(smoking_4levels, levels = c("nonsmoker", "past smoker (>15)", "past smoker (<=15)", "current smoker", "past smoker (unspecificed)", "Unknown"))
smoking_status = smoking_4levels

#smoking_status = TCGA_phenotypes$smoking_status
#smoking_status[which(is.na(smoking_status))] = "Unknown"
#smoking_status[which(smoking_status == 1)] = "No"
#smoking_status[which(smoking_status %in% 2:5)] = "Yes"
#smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

tumor_stage = TCGA_phenotypes$tumor_stage
tumor_stage[which(tumor_stage == "stage i" | tumor_stage == "stage ia" | tumor_stage == "stage ib")] = "stageI"
tumor_stage[which(tumor_stage == "stage ii" | tumor_stage == "stage iia" | tumor_stage == "stage iib")] = "stageII"
tumor_stage[which(tumor_stage == "Stage iii" | tumor_stage == "stage iiia" | tumor_stage == "stage iiib")] = "stageIII"
tumor_stage[which(tumor_stage == "stage iv")] = "stageIV"
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "stageIV", "not reported"))

design = model.matrix(~ age + gender + race + smoking_status + tumor_stage)
indegree = indegree_TCGA
fit <- lmFit(indegree, design)
fit <- eBayes(fit)

# Save table for sex difference: genderFEMALE
tb = topTable(fit,coef="age",number=Inf)
tb$chr = c()
tb$gene_name = c()
for (i in 1:nrow(tb)){
  index = match(rownames(tb)[i], gene_info.TCGA$gene_id)
  tb$chr[i]=gene_info.TCGA$seqnames[index]
  tb$gene_name[i]=gene_info.TCGA$gene_name[index]
}
head(tb)

############# GSEA ###################
# Rank genes in limma table
indegree_rank <- setNames(object=tb[,"t"], tb$gene_name)
head(indegree_rank)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways, indegree_rank, minSize=15, maxSize=500)
head(fgseaRes)

# Cancer Pathways
fgseaRes[which(fgseaRes$pathway %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_PATHWAYS_IN_CANCER", "KEGG_SMALL_CELL_LUNG_CANCER",
                                       "KEGG_NON_SMALL_CELL_LUNG_CANCER", "KEGG_WNT_SIGNALING_PATHWAY",
                                       "KEGG_P53_SIGNALING_PATHWAY"))]


# Immune Pathways

rbind(fgseaRes[grep("immun", fgseaRes$pathway, ignore.case = T),], fgseaRes[grep("kine", fgseaRes$pathway, ignore.case = T),])

sig_pathways = fgseaRes[which(fgseaRes$padj < 0.05),]

write.table(tb, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_TCGA_aging.txt", col.names = T, row.names = T)
save(fgseaRes, file = "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/gsea_TCGA_aging.RData")
