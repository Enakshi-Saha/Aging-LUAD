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
# Load reactome pathways
# pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.reactome.v2022.1.Hs.symbols.gmt")
# Load GO pathways
# pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c5.go.bp.v2022.1.Hs.symbols.gmt")

# Load indegree matrices
indegree_male = data.frame(fread("/home/esaha/Lung_lioness/TPM/TCGA/TCGA_indegree_tpm_male.txt"))
indegree_female = data.frame(fread("/home/esaha/Lung_lioness/TPM/TCGA/TCGA_indegree_tpm_female.txt"))

head(indegree_male[,1:4])
head(indegree_female[,1:4])

# Check if male and female edges are in the same order
sum(indegree_male$V1 != indegree_female$V1)

genes = indegree_male$V1

# Get chromosome and gene names
info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_RSE.RData"))
gene_info.TCGA = rowRanges(info.TCGA)
gene_info.TCGA = data.frame(gene_info.TCGA)
gene_info.TCGA$gene_id = gsub("\\..*","",gene_info.TCGA$gene_id)
names(gene_info.TCGA)
head(gene_info.TCGA)
gene_name = rep(0, length(genes))
chr_loc = rep(0, length(genes))
for (i in 1:length(genes)){
  index = which(gene_info.TCGA$gene_id == genes[i])[1]
  gene_name[i]=gene_info.TCGA$gene_name[index]
  chr_loc[i]=as.character(gene_info.TCGA$seqnames[index])
}
gene_name[1:100]
chr_loc[1:100]

rownames(indegree_male) = genes
rownames(indegree_female) = genes

indegree = cbind(indegree_male[,-1], indegree_female[,-1])

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

# write.table(pathway_score_pca, "/home/esaha/cibersort/pathway_score_pca.txt", col.names = T, row.names = T)

# remove X from the starting of TCGA IDs
colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")] = substring(colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")], 2, length(colnames(indegree)[which(substring(colnames(indegree), 1,1) == "X")]))
# Format male & female IDs to match with phenotypic data: TCGA
TCGA_IDs = unlist(lapply(strsplit(colnames(indegree), split=".", fixed = T),
                         function(x){paste(x, collapse ="-")}))


# Get phenotypes: TCGA
pheno_info.TCGA = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/TCGA_eset.RData"))
pheno_info.TCGA = pData(pheno_info.TCGA)

TCGA_phenotypes = as.data.frame(cbind(pheno_info.TCGA$external_id,
                                      pheno_info.TCGA$tcga.gdc_cases.demographic.gender, 
                                      pheno_info.TCGA$tcga.gdc_cases.demographic.race,
                                      pheno_info.TCGA$tcga.xml_age_at_initial_pathologic_diagnosis,
                                      pheno_info.TCGA$tcga.gdc_cases.samples.sample_type, 
                                      pheno_info.TCGA$tcga.gdc_cases.diagnoses.tumor_stage,
                                      pheno_info.TCGA$tcga.tcga_barcode,
                                      pheno_info.TCGA$tcga.xml_tobacco_smoking_history,
                                      pheno_info.TCGA$tcga.xml_days_to_death,
                                      pheno_info.TCGA$tcga.xml_days_to_last_followup,
                                      pheno_info.TCGA$tcga.xml_vital_status,
                                      pheno_info.TCGA$tcga.cgc_drug_therapy_pharmaceutical_therapy_type))
TCGA_phenotypes =  TCGA_phenotypes[match(TCGA_IDs, rownames(pheno_info.TCGA)),]
colnames(TCGA_phenotypes) = c("id","gender", "race", "age", "sample_type", 
                              "tumor_stage", "barcode", "smoking_status",
                              "death_days", "last_contact_days", "vital_status", "therapy")
head(TCGA_phenotypes)
dim(TCGA_phenotypes)

# Define the covariates: gender, age, race, smoking, tumor_stage
TCGA_phenotypes$gender[which(TCGA_phenotypes$gender == "male")] = "MALE" 
TCGA_phenotypes$gender[which(TCGA_phenotypes$gender == "female")] = "FEMALE" 
TCGA_phenotypes$gender = factor(TCGA_phenotypes$gender, levels = c("MALE", "FEMALE"))

TCGA_phenotypes$race[which(TCGA_phenotypes$race != "black or african american" & TCGA_phenotypes$race != "white")] = "others"
TCGA_phenotypes$race = factor(TCGA_phenotypes$race)

TCGA_phenotypes$age <- as.numeric(as.character(TCGA_phenotypes$age))
TCGA_phenotypes$age[which(is.na(TCGA_phenotypes$age))] = mean(TCGA_phenotypes$age,na.rm=TRUE)

TCGA_phenotypes$smoking_status[which(is.na(TCGA_phenotypes$smoking_status))] = "Unknown"
TCGA_phenotypes$smoking_status[which(TCGA_phenotypes$smoking_status == 1)] = "No"
TCGA_phenotypes$smoking_status[which(TCGA_phenotypes$smoking_status %in% 2:5)] = "Yes"

tumor_stage = TCGA_phenotypes$tumor_stage
tumor_stage[which(tumor_stage == "stage i" | tumor_stage == "stage ia" | tumor_stage == "stage ib")] = "Stage I"
tumor_stage[which(tumor_stage == "stage ii" | tumor_stage == "stage iia" | tumor_stage == "stage iib")] = "Stage II"
tumor_stage[which(tumor_stage == "Stage iii" | tumor_stage == "stage iiia" | tumor_stage == "stage iiib")] = "Stage III"
tumor_stage[which(tumor_stage == "stage iv")] = "Stage IV"
tumor_stage = factor(tumor_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "not reported"))
TCGA_phenotypes$tumor_stage = factor(tumor_stage)

therapy = TCGA_phenotypes$therapy
therapy[which(is.na(therapy))] = "NA"
therapy[-which(therapy %in% c("Chemotherapy", "NA"))] = "other"
TCGA_phenotypes$therapy = factor(therapy)

# Separate tumor and adjacent normal samples
TCGA_indegree_normal = indegree[, which(TCGA_phenotypes$sample_type == "Solid Tissue Normal")]
TCGA_indegree_tumor = indegree[, which(TCGA_phenotypes$sample_type == "Primary Tumor")]

TCGA_phenotypes_normal = TCGA_phenotypes[which(TCGA_phenotypes$sample_type == "Solid Tissue Normal"),]
TCGA_phenotypes_tumor = TCGA_phenotypes[which(TCGA_phenotypes$sample_type == "Primary Tumor"),]

pathway_score_median_normal = pathway_score_median[, which(TCGA_phenotypes$sample_type == "Solid Tissue Normal")]
pathway_score_median_tumor = pathway_score_median[, which(TCGA_phenotypes$sample_type == "Primary Tumor")]

pathway_score_pca_normal = pathway_score_pca[, which(TCGA_phenotypes$sample_type == "Solid Tissue Normal")]
pathway_score_pca_tumor = pathway_score_pca[, which(TCGA_phenotypes$sample_type == "Primary Tumor")]

gender = TCGA_phenotypes_tumor$gender
race = TCGA_phenotypes_tumor$race
age = TCGA_phenotypes_tumor$age
smoking_status = TCGA_phenotypes_tumor$smoking_status
tumor_stage = TCGA_phenotypes_tumor$tumor_stage
therapy = TCGA_phenotypes_tumor$therapy

death_days <- as.numeric(as.character(TCGA_phenotypes_tumor$death_days))
#death_days[which(is.na(death_days))] <- mean(death_days,na.rm=TRUE)

last_contact_days <- as.numeric(as.character(TCGA_phenotypes_tumor$last_contact_days))
#last_contact_days[which(is.na(last_contact_days))] <- mean(last_contact_days,na.rm=TRUE)

vital_status <- TCGA_phenotypes_tumor$vital_status
vital_status = factor(vital_status)

# Survival Analysis: Cox Regression

time <-  ifelse(is.na(death_days), yes=last_contact_days, no=death_days)
time <- time/365
status <- vital_status
status <- sub("Alive", 0, status)
status <- sub("Dead", 1, status)
status <- sub("[Discrepancy]", NA, status)
status <- as.numeric(status)

sig_paths_GTEx = read.csv("/home/esaha/aging_plots/KEGG/sigpathsList_NES/GTEx_significantPathways_05.txt", sep="")$x
TCGA_fgsea = get(load("/home/esaha/aging_plots/KEGG/GSEA_tables/TCGA_age_fgseaRes.RData"))
sig_paths_TCGA = TCGA_fgsea$pathway[which(TCGA_fgsea$padj < 0.05)]

sig_paths = sig_paths_TCGA

# Get pathways significantly different in male and female
pathway_score_tumor = pathway_score_pca_tumor
pathway_score_tumor = data.frame(pathway_score_tumor)
sig_pathways_score =  t(pathway_score_tumor[match(sig_paths, rownames(pathway_score_tumor)),])

## Load survival package
library(survival)
SurvObj = Surv(time, status)
age_group = age
age_group[which(age <= median(age))] = "younger"
age_group[-which(age <= median(age))] = "older"
age_group = factor(age_group)
diff = survdiff(Surv(time,status)~age_group)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)
plot(survfit(Surv(time,status)~age_group), main = "Plot of Survival Curves by Age", xlab = "Length of Survival",ylab="Proportion of Individuals who have Survived",col=c("blue","red"))
legend("topright", legend=levels(age_group),fill=c("blue","red"),bty="n")

# Fit Cox regression for every pathway
#age_group = rep("old", length(age))
#age_group[which(age < median(age))] = "young"
# age_group = factor(age_group, levels = c("old", "young"))
score_coeff = rep(0, length(sig_paths))
pval_coeff = score_coeff
for (i in 1:length(sig_paths)){
  score = sig_pathways_score[,i]
  cox_tcga <- coxph(SurvObj ~ age + gender + race + tumor_stage + smoking_status + therapy 
                    + score)
  coeffs = summary(cox_tcga)$coefficients
  score_coeff[i] = coeffs["score",4]
  pval_coeff[i] = coeffs["score",5]
}
names(score_coeff) = colnames(sig_pathways_score)
score_coeff = data.frame(score_coeff)
score_coeff$pvalue = pval_coeff
head(score_coeff)

######## Biological age score
sig_paths_cox = rownames(score_coeff[which(score_coeff$pvalue < 0.19),])

sig_pathways_score_subset = sig_pathways_score[,which(colnames(sig_pathways_score) %in% sig_paths_cox)]
sig_pathways_score0 = sig_pathways_score[,-which(colnames(sig_pathways_score) %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_RETINOL_METABOLISM"))]
# cox_tcga <- coxph(SurvObj ~ sig_pathways_score0)

# save(cox_tcga, file="/home/esaha/aging_plots/kaplan_meier/cox_tcga_biologicalAge.RData")

# cox_results = summary(cox_tcga)$coefficients
# cox_results[which(cox_results[,5] < 0.05),]

# biological_age = predict(cox_tcga, data.frame(sig_pathways_score))
# biological_age0 = biological_age
# names(biological_age0) = TCGA_phenotypes_tumor$id

###########################################################
############## Regularize Cox Regression ###################
set.seed(1000)
library(glmnet)
library(survival)
SurvObj = Surv(time, status)

# standardize pathway scores
sig_pathways_score0 = scale(sig_pathways_score0)

# Combine into a data frame to remove rows with any NAs
data_combined <- cbind(SurvObj, sig_pathways_score0)
data_combined <- data_combined[complete.cases(data_combined), ]

training_samples = sample(1:nrow(data_combined), size = floor(0.5*nrow(data_combined)), replace = F)
data_training = data_combined[training_samples,]
data_test = data_combined[-training_samples,]

# Re-extract response and predictor
SurvObj_clean <- data_training[, 1:2] # Assuming SurvObj has 2 columns: time and status
sig_pathways_score0_clean <- as.matrix(data_training[, -(1:2)])
time_shifted <- SurvObj_clean[,1] + 0.01
SurvObj_clean2 <- Surv(time = time_shifted, event = SurvObj_clean[,2])

cvfit <- cv.glmnet(sig_pathways_score0_clean, SurvObj_clean2, family = "cox", type.measure = "deviance", 
                   lambda = seq(from = 0, to = 0.01, by = 0.001))
cvfit$lambda.1se

fit <- glmnet(sig_pathways_score0_clean, SurvObj_clean2, family = "cox", lambda = cvfit$lambda.1se)


# Check selected coefficients
coefs <- coef(fit)
# Convert to a numeric vector
coefs_vec <- as.matrix(coefs)
# See non-zero coefficients (i.e., selected variables)
selected_vars <- rownames(coefs_vec)[which(coefs_vec != 0)]

########################################################
# Predict linear predictor (risk score) on training data
SurvObj_clean_test <- data_test[, 1:2] # Assuming SurvObj has 2 columns: time and status
sig_pathways_score0_clean_test <- as.matrix(data_test[, -(1:2)])
time_shifted_test <- SurvObj_clean_test[,1] + 0.01
SurvObj_clean2_test <- Surv(time = time_shifted_test, event = SurvObj_clean_test[,2])
biological_age <- predict(fit, newx = sig_pathways_score0_clean_test, type = "link")
names(biological_age) = rownames(data_test)

# test_samples = rownames(SurvObj_clean_test)
# biological_age = TCGA_phenotypes_tumor$age[match(test_samples, rownames(sig_pathways_score0))]
biological_age_group = biological_age
biological_age_group[which(biological_age <= median(biological_age))] = "lower aging"
biological_age_group[-which(biological_age <= median(biological_age))] = "higher aging"
biological_age_group = factor(biological_age_group)
fit_test = survfit(SurvObj_clean2_test ~ biological_age_group)
diff <- survdiff(SurvObj_clean2_test ~ biological_age_group)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)

library(survival)
library(ggplot2)
library(GGally)
library(stringr)

ggsurv(fit_test, CI=FALSE,  plot.cens = F, cens.col = "gray", surv.col=c("red", "orange"),
       xlab = "Time (years)", ylab = "Overall survival probability",
       main="Survival probability by Aging Signature (TCGA)") +
  coord_cartesian(ylim=c(0,1)) +
  theme(text=element_text(size=12), axis.text=element_text(size=10), plot.title=element_text(size=16),
        legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.9,0.8)) +
  geom_text(x=1.2, y=0.05, label=paste0("p=",pval), size=5, family="Helvetica")

##########################
# Check if biological age or chronological age is significant in Cox Model on Test Data
# Effect on Survival
biological_age_all <- predict(fit, newx = sig_pathways_score0, type = "link")
cox_all = coxph(SurvObj ~ age + gender + race + tumor_stage + smoking_status
      + biological_age_all)
summary(cox_all)
# Effect on Chemotherapy-response
cox_all_therapy = coxph(SurvObj ~ biological_age_all*therapy)
summary(cox_all_therapy)

##########################
c_age = TCGA_phenotypes_tumor$age[match(test_samples, rownames(sig_pathways_score0))]
cor.test(c_age, biological_age)

tumor_stage_test = tumor_stage[match(test_samples, rownames(sig_pathways_score0))]
anova(lm(biological_age ~ tumor_stage_test))
