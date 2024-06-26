# Survival Analysis with Pathway Score
set.seed(0)

library(fgsea)
library(limma)
library(Biobase)
library(ggplot2)
library(igraph)
library(data.table)
library(recount3)
library(recount)

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")
sigPaths = read.table("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/TCGA_sigPaths_validatedPaths.txt", quote="\"", comment.char="")$V1
# Load indegree matrices
indegree_TCGA_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_male.txt"))
indegree_TCGA_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/TCGA_lung_lioness_indegree_female.txt"))

genes = indegree_TCGA_male$V1

indegree = cbind(indegree_TCGA_male[,-1], indegree_TCGA_female[,-1])
rownames(indegree) = genes

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
gene_name = gene_info.TCGA$gene_name[match(rownames(indegree), gene_info.TCGA$gene_id)]

# remove X from the starting of TCGA IDs
colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")] = substring(colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")], 2, length(colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(colnames(indegree), split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))

# Get phenotypes: TCGA
TCGA_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/TCGA_phenotypes.txt", sep=""))
TCGA_phenotypes = TCGA_phenotypes[match(TCGA_IDs, rownames(TCGA_phenotypes)),]
TCGA_phenotypes = as.data.frame(cbind(TCGA_phenotypes$tcga.gdc_cases.demographic.gender, 
                                      TCGA_phenotypes$tcga.gdc_cases.demographic.race,
                                      TCGA_phenotypes$tcga.xml_age_at_initial_pathologic_diagnosis,
                                      TCGA_phenotypes$tcga.gdc_cases.samples.sample_type, 
                                      TCGA_phenotypes$tcga.gdc_cases.diagnoses.tumor_stage,
                                      TCGA_phenotypes$tcga.tcga_barcode,
                                      TCGA_phenotypes$tcga.xml_tobacco_smoking_history,
                                      TCGA_phenotypes$tcga.xml_days_to_death,
                                      TCGA_phenotypes$tcga.xml_days_to_last_followup,
                                      TCGA_phenotypes$tcga.xml_vital_status,
                                      TCGA_phenotypes$tcga.cgc_drug_therapy_pharmaceutical_therapy_type))
colnames(TCGA_phenotypes) = c("gender", "race", "age", "sample_type", 
                              "tumor_stage", "barcode", "smoking_status",
                              "death_days", "last_contact_days", "vital_status", "therapy")
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

tumor_stage = TCGA_phenotypes$tumor_stage
tumor_stage[which(tumor_stage == "stage i" | tumor_stage == "stage ia" | tumor_stage == "stage ib")] = "stageI"
tumor_stage[which(tumor_stage == "stage ii" | tumor_stage == "stage iia" | tumor_stage == "stage iib")] = "stageII"
tumor_stage[which(tumor_stage == "Stage iii" | tumor_stage == "stage iiia" | tumor_stage == "stage iiib")] = "stageIII"
tumor_stage[which(tumor_stage == "stage iv")] = "stageIV"
tumor_stage = factor(tumor_stage, levels = c("stageI", "stageII", "stageIII", "stageIV", "not reported"))

therapy = TCGA_phenotypes$therapy
therapy[which(is.na(therapy))] = "NA"
therapy[-which(therapy %in% c("Chemotherapy", "NA"))] = "other"

death_days <- as.numeric(as.character(TCGA_phenotypes$death_days))
last_contact_days <- as.numeric(as.character(TCGA_phenotypes$last_contact_days))
vital_status <- TCGA_phenotypes$vital_status
vital_status = factor(vital_status)

# Survival Analysis: Cox Regression

time <-  ifelse(is.na(death_days), yes=last_contact_days, no=death_days)
time <- time/365
status <- vital_status
status <- sub("Alive", 0, status)
status <- sub("Dead", 1, status)
status <- sub("[Discrepancy]", NA, status)
status <- as.numeric(status)

# Pathways significantly changing with age and validated on GSE68465

## Load survival package
library(survival)
SurvObj = Surv(time, status)

##### Get pathway score
gene_pca = function(indegree_matrix){
  gene.pca = prcomp(t(indegree_matrix))
  gene.pca$x[,1]
}
pathway_score_median = c()
pathway_score_pca = c()
for (i in 1:length(pathways)){
  indegrees = indegree[which(gene_name %in% pathways[i][[1]]),]
  pathway_score_median = rbind(pathway_score_median,unlist(apply(indegrees, MARGIN = 2, FUN = median)))
  pathway_score_pca = rbind(pathway_score_pca, gene_pca(indegrees))
}
rownames(pathway_score_median) = names(pathways)
rownames(pathway_score_pca) = names(pathways)

pathway_score = pathway_score_pca

write.table(pathway_score_pca, "/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/pathway_score_pca.txt", row.names = T, col.names = T)

# Fit Cox regression for every pathway
score_coeff = rep(0, length(pathways))
pval_coeff = score_coeff
for (i in 1:length(pathways)){
  #score = unlist(pathway_score[[sigPaths[i]]])
  score = pathway_score[i,]
  cox_tcga <- coxph(SurvObj ~ age + gender + race + tumor_stage + smoking_status + therapy 
                    + score)
  coeffs = summary(cox_tcga)$coefficients
  score_coeff[i] = coeffs["score",4]
  pval_coeff[i] = coeffs["score",5]
}

score_coeff = data.frame(score_coeff)
rownames(score_coeff) = stringr::str_to_title(lapply(strsplit(names(pathways), split = "_"), function(x){paste(x[-1], collapse = " ")}))
score_coeff$pvalue = pval_coeff
head(score_coeff)

TCGA_survival_tab = score_coeff[which(score_coeff$pvalue < 0.1),]
TCGA_survival_tab = TCGA_survival_tab[which(rownames(TCGA_survival_tab) %in% sigPaths),]

write.table(TCGA_survival_tab, "/home/esaha/aging_plots/aging_plots2023/survival_aging/TCGA_survival_pathways.csv", quote = F, row.names = T, col.names = T)

# pathways interacting with chemotherapy
therapy = factor(therapy, levels = c("NA", "Chemotherapy", "other"))
score_coeff = matrix(0, nrow = length(pathways), ncol = 2)
pval_coeff = score_coeff
for (i in 1:length(pathways)){
  score = pathway_score[i,]
  cox_tcga <- coxph(SurvObj ~ age + gender + race + tumor_stage + smoking_status + therapy 
                    + score*therapy)
  coeffs = summary(cox_tcga)$coefficients
  score_coeff[i,] = coeffs[c("score", "therapyChemotherapy:score"),4]
  pval_coeff[i,] = coeffs[c("score", "therapyChemotherapy:score"),5]
}
score_coeff = data.frame(score_coeff)
rownames(score_coeff) = stringr::str_to_title(lapply(strsplit(names(pathways), split = "_"), function(x){paste(x[-1], collapse = " ")}))
colnames(score_coeff) = c("no_therapy", "chemotherapy")
score_coeff$pval_interaction = pval_coeff[,2]
score_coeff[,2] = score_coeff[,1]+score_coeff[,2]
head(score_coeff)

TCGA_survival_therapy_tab = score_coeff[which(score_coeff$pval_interaction < 0.05),]
TCGA_survival_therapy_tab = TCGA_survival_therapy_tab[which(rownames(TCGA_survival_therapy_tab) %in% sigPaths),]

write.csv(TCGA_survival_therapy_tab, "/home/esaha/aging_plots/aging_plots2023/survival_aging/TCGA_survival_therapy_pathways.csv", quote = F, row.names = T)

