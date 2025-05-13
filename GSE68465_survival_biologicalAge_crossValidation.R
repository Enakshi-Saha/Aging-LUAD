# Survival Analysis with Pathway Score
# First run TCGA_survival_biologicalAge_crossValidation.R
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

therapy = clinical$clinical_treatment_adjuvant_chemo.ch1

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

pathway_score_selected = scale(pathway_score[colnames(sig_pathways_score0_clean_test),])

####################
biological_age <- predict(fit, newx = t(pathway_score_selected), type = "link")
names(biological_age) = rownames(data_test)

biological_age_group = biological_age
biological_age_group[which(biological_age <= median(biological_age))] = "lower aging"
biological_age_group[-which(biological_age <= median(biological_age))] = "higher aging"
biological_age_group = factor(biological_age_group)
fit_test = survfit(SurvObj ~ biological_age_group)
diff <- survdiff(SurvObj ~ biological_age_group)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)

library(survival)
library(ggplot2)
library(GGally)
library(stringr)

ggsurv(fit_test, CI=FALSE,  plot.cens = F, cens.col = "gray", surv.col=c("red", "orange"),
       xlab = "Time (years)", ylab = "Overall survival probability",
       main="Survival probability by Aging Signature (GSE68465)") +
  coord_cartesian(ylim=c(0,1)) +
  theme(text=element_text(size=12), axis.text=element_text(size=10), plot.title=element_text(size=16),
        legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.9,0.8)) +
  geom_text(x=5, y=0.05, label=paste0("p=",pval), size=5, family="Helvetica")

# Check if biological age or chronological age is significant in Cox Model on Test Data
cox_all = coxph(SurvObj ~ age + gender + race + tumor_stage + smoking +
                + biological_age)
summary(cox_all)

# Effect on Chemotherapy-response
cox_all_therapy = coxph(SurvObj ~ biological_age*therapy)
summary(cox_all_therapy)


