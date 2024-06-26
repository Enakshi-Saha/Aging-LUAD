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
indegree_male = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_male.txt")
indegree_female = fread("/home/esaha/validation_data/validation_newdata2023/validation_indegree/GSE68465_lioness_indegree_female.txt")
indegree = cbind(indegree_male, indegree_female[,-1])
genes = indegree$V1
indegree = data.frame(indegree[,-1])
rownames(indegree) = genes
head(indegree[,1:4])

# Get clinical data
clinical = read.csv("~/validation_data/validation_newdata2023/GSE68465_phenotypes.txt", sep="")
clinical = clinical[match(colnames(indegree),clinical$geo_accession),]
head(clinical[,1:5])

gender = clinical$Sex.ch1
gender = factor(gender, levels = c("Male", "Female"))

race = clinical$race.ch1
race[which(race == "Not Reported")] = "Unknown"
race[-which(race %in% c("White", "Black or African American", "Unknown"))] = "others"
race = factor(race)

age = clinical$age.ch1

smoking = clinical$smoking_history.ch1
smoking[-which(smoking %in% c("Never smoked", "Currently smoking", "Smoked in the past"))] = "Unknown"
smoking = factor(smoking, levels = c("Never smoked", "Currently smoking", "Smoked in the past", "Unknown"))
table(smoking)

stage = clinical$disease_stage.ch1
tumor_stage = substr(stage, nchar(stage), nchar(stage))
tumor_stage = factor(tumor_stage, levels = c("1", "2", "3", "4", "p"))

chemo = factor(clinical$clinical_treatment_adjuvant_chemo.ch1)
radiation = factor(clinical$clinical_treatment_adjuvant_rt.ch1)

# Compute pathway scores
gene_pca = function(indegree_matrix){
  gene.pca = prcomp(t(indegree_matrix))
  gene.pca$x[,1]
}
pathway_score_median = c()
pathway_score_pca = c()
for (i in 1:length(pathways)){
  indegrees = indegree[which(genes %in% pathways[i][[1]]),]
  pathway_score_median = rbind(pathway_score_median,unlist(apply(indegrees, MARGIN = 2, FUN = median)))
  pathway_score_pca = rbind(pathway_score_pca, gene_pca(indegrees))
}
rownames(pathway_score_median) = names(pathways)
rownames(pathway_score_pca) = names(pathways)

pathway_score = pathway_score_pca
# Survival Analysis: Cox Regression

time <-  as.numeric(as.character((clinical$months_to_last_contact_or_death.ch1)))
time[is.na(time)] = 0
vital_status <- clinical$vital_status.ch1
vital_status = factor(vital_status)

status <- vital_status
status <- sub("Alive", 0, status)
status <- sub("Dead", 1, status)
status <- as.numeric(status)

# Pathways significantly changing with age and validated on GSE68465

## Load survival package
library(survival)
SurvObj = Surv(time, status)

# pathways interacting with chemotherapy
therapy = rep(0, length(chemo))
therapy[which(chemo=="No" & radiation =="No")] = "No"
therapy[which(chemo=="Yes" & radiation =="No")] = "chemo"
therapy[which(chemo=="No" & radiation =="Yes")] = "radiation"
therapy[which(chemo=="Yes" & radiation =="Yes")] = "both"
therapy = factor(therapy, levels = c("No", "chemo", "radiation", "both", "Unknown"))

chemo_score_coeff = matrix(0, nrow = length(pathways), ncol = 2)
chemo_pval_coeff = chemo_score_coeff
radio_score_coeff = matrix(0, nrow = length(pathways), ncol = 2)
radio_pval_coeff = radio_score_coeff
for (i in 1:length(pathways)){
  score = pathway_score[i,]
  cox_tcga <- coxph(SurvObj ~ age + gender + race + tumor_stage + smoking + therapy +
                    + score*therapy)
  coeffs = summary(cox_tcga)$coefficients
  chemo_score_coeff[i,] = coeffs[c("score", "therapychemo:score"),4]
  chemo_pval_coeff[i,] = coeffs[c("score", "therapychemo:score"),5]
  radio_score_coeff[i,] = coeffs[c("score", "therapyradiation:score"),4]
  radio_pval_coeff[i,] = coeffs[c("score", "therapyradiation:score"),5]
}
rownames(chemo_score_coeff) = stringr::str_to_title(lapply(strsplit(names(pathways), split = "_"), function(x){paste(x[-1], collapse = " ")}))
colnames(chemo_score_coeff) = c("no_therapy", "chemotherapy")
chemo_score_coeff = data.frame(chemo_score_coeff)
chemo_score_coeff$pval_interaction = chemo_pval_coeff[,2]
chemo_score_coeff[,2] = chemo_score_coeff[,1]+chemo_score_coeff[,2]

head(chemo_score_coeff)
survival_chemo_tab = chemo_score_coeff[which(chemo_score_coeff$pval_interaction < 0.05),]
survival_chemo_tab = survival_chemo_tab[which(rownames(survival_chemo_tab) %in% sigPaths),]

rownames(radio_score_coeff) = stringr::str_to_title(lapply(strsplit(names(pathways), split = "_"), function(x){paste(x[-1], collapse = " ")}))
colnames(radio_score_coeff) = c("no_therapy", "radiotherapy")
radio_score_coeff = data.frame(radio_score_coeff)
radio_score_coeff$pval_interaction = radio_pval_coeff[,2]
radio_score_coeff[,2] = radio_score_coeff[,1]+radio_score_coeff[,2]

head(radio_score_coeff)
survival_radio_tab = radio_score_coeff[which(radio_score_coeff$pval_interaction < 0.05),]
survival_radio_tab = survival_radio_tab[which(rownames(survival_radio_tab) %in% sigPaths),]
