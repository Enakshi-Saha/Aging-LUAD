#########################################
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
library(plyr)

# Load indegree
indegree_male = fread("/home/esaha/Lung_lioness/validation_networks/GSE68465/indegree_GSE68465_male.txt")
indegree_female = fread("/home/esaha/Lung_lioness/validation_networks/GSE68465/indegree_GSE68465_female.txt")
indegree = cbind(indegree_male, indegree_female[,-1])
head(indegree[,1:4])
genes = indegree$V1

# Load KEGG pathways
pathways <- gmtPathways("/home/ubuntu/GSEA_pathways/c2.cp.kegg.v2022.1.Hs.symbols.gmt")

genes = indegree$V1
indegree = indegree[,-1]
rownames(indegree) = genes

# Get clinical data
clinical = fread("/home/esaha/validation_data/GSE68465/GSE68465_phenotype.txt")
clinical = clinical[match(colnames(indegree),clinical$geo_accession),]
head(clinical[,1:5])

gender = clinical$`Sex:ch1`
gender = factor(gender, levels = c("Male", "Female"))
race = factor(clinical$`race:ch1`)
age = clinical$`age:ch1`

smoking = clinical$`smoking_history:ch1`
smoking[which(smoking == "--")] = "Unknown"
smoking = factor(smoking)

stage = clinical$`disease_stage:ch1`
tumor_stage = substr(stage, nchar(stage), nchar(stage))

batch = clinical$source_name_ch1

chemo = clinical$`clinical_treatment_adjuvant_chemo:ch1`
radiation = clinical$`clinical_treatment_adjuvant_rt:ch1`

therapy = chemo
therapy[which(chemo == "Unknown" & radiation == "Unknown")] = "NA"
therapy[which(chemo == "Yes" & radiation == "No")] = "chemo"
therapy[which(chemo == "No" & radiation == "Yes")] = "radiation"
therapy[which(chemo == "Yes" & radiation == "Yes")] = "both"
therapy[which(chemo == "No" & radiation == "No")] = "noTherapy"

therapy = factor(therapy, levels = c("noTherapy", "NA", "chemo", "radiation", "both"))


# Compute pathway scores
pathway_score = c()
for (i in 1:length(pathways)){
  indegrees = indegree[which(genes %in% pathways[i][[1]]),]
  pathway_score = rbind(pathway_score,unlist(apply(indegrees, MARGIN = 2, FUN = mean)))
}
rownames(pathway_score) = names(pathways)

# Survival Analysis: Cox Regression

time <-  as.numeric(as.character((clinical$`months_to_last_contact_or_death:ch1`)))
time[is.na(time)] = 0
vital_status <- clinical$`vital_status:ch1`
vital_status = factor(vital_status)

status <- vital_status
status <- sub("Alive", 0, status)
status <- sub("Dead", 1, status)
status <- as.numeric(status)

TCGA_fgsea = get(load("/home/esaha/aging_plots/KEGG/GSEA_tables/TCGA_age_fgseaRes.RData"))
sig_paths_TCGA = TCGA_fgsea$pathway[which(TCGA_fgsea$padj < 0.05)]
GSE68465_fgsea = get(load("~/aging_plots/KEGG/GSEA_tables/GSE68465_age_fgseaRes.RData"))
sig_path_GSE68465 = GSE68465_fgsea$pathway[which(GSE68465_fgsea$padj < 0.05)]

sig_paths = sig_path_GSE68465

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
  indegrees = indegree[which(genes %in% pathways[i][[1]]),]
  pathway_score_median = rbind(pathway_score_median,unlist(apply(indegrees, MARGIN = 2, FUN = median)))
  pathway_score_pca = rbind(pathway_score_pca, gene_pca(indegrees))
}
rownames(pathway_score_median) = names(pathways)
rownames(pathway_score_pca) = names(pathways)

# Get pathways significantly different in male and female
sig_pathways_score =  pathway_score_pca[which(rownames(pathway_score_pca) %in% sig_paths_TCGA),]
therapy = rep("Unknown", length(radiation))
therapy[which(chemo == "Yes" & radiation == "Yes")] = "both"
therapy[which(chemo == "Yes" & radiation == "No")] = "chemo_only"
therapy[which(chemo == "No" & radiation == "Yes")] = "radio_only"
therapy = factor(therapy, levels = c("Unknown", "chemo_only", "radio_only", "both"))

sig_pathways_score0 = sig_pathways_score[-which(rownames(sig_pathways_score) %in% c("KEGG_DRUG_METABOLISM_CYTOCHROME_P450", "KEGG_RETINOL_METABOLISM")),]

cox_tcga = get(load("/home/esaha/aging_plots/kaplan_meier/cox_tcga_biologicalAge.RData"))
biological_age = cox_tcga$coefficients %*% sig_pathways_score0

library(ggplot2)
fulldata = data.frame(age)
fulldata$aging_score = as.numeric(biological_age)
tumor_stage[which(tumor_stage == "p")] = "NA"
fulldata$tumor_stage = tumor_stage
ggplot(fulldata, aes(x = age , y = aging_score, col = tumor_stage)) +
  geom_point(size=1.5) + labs(title="Network-informed Aging score in Tumor (GSE68465)",
                              y="aging score", x = "age")

library(survival)
SurvObj = Surv(time, status)
biological_age_group = biological_age
biological_age_group[which(biological_age <= median(biological_age))] = "lower aging"
biological_age_group[-which(biological_age <= median(biological_age))] = "higher aging"
biological_age_group = factor(biological_age_group)
fit = survfit(SurvObj ~ biological_age_group)
diff <- survdiff(SurvObj ~ biological_age_group)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)

library(survival)
library(ggplot2)
library(GGally)
library(stringr)

ggsurv(fit, CI=FALSE,  plot.cens = F, cens.col = "gray", surv.col=c("red", "orange"),
       xlab = "Time (months)", ylab = "Overall survival probability",
       main="Survival probability by Aging Signature (GSE68465)") +
  coord_cartesian(ylim=c(0,1)) +
  theme(text=element_text(size=12), axis.text=element_text(size=10), plot.title=element_text(size=16),
        legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.9,0.8)) +
  geom_text(x=6, y=0.05, label=paste0("p=",pval), size=5, family="Helvetica")

####### Kaplan Meier Plot by Clinical Variables ###########
SurvObj0 = SurvObj[which(gender == "Female"),]
biological_age_group0 = biological_age_group[which(gender == "Female")]
fit = survfit(SurvObj0 ~ biological_age_group0)
diff <- survdiff(SurvObj0 ~ biological_age_group0)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)
ggsurv(fit, CI=FALSE,  plot.cens = F, cens.col = "gray", surv.col=c("red", "orange"),
       xlab = "Time (months)", ylab = "Overall survival probability",
       main="Survival by Aging Signature (GSE68465 Female)") +
  coord_cartesian(ylim=c(0,1)) +
  theme(text=element_text(size=12), axis.text=element_text(size=10), plot.title=element_text(size=16),
        legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.9,0.8)) +
  geom_text(x=4, y=0.05, label=paste0("p=",pval), size=5, family="Helvetica")

############################
library(survival)
library(ggplot2)
library(GGally)
library(stringr)

age_group = age
age_group[which(age <= median(age))] = "younger"
age_group[-which(age <= median(age))] = "older"
age_group = factor(age_group)

fit = survfit(SurvObj ~ age_group)
diff <- survdiff(SurvObj ~ age_group)
pval <- round(1 - pchisq(diff$chisq, length(diff$n) - 1),3)

ggsurv(fit, CI=FALSE,  plot.cens = F, cens.col = "gray", surv.col=c("red", "orange"),
       xlab = "Time (months)", ylab = "Overall survival probability",
       main="Survival probability by Chronological Age (GSE68465)") +
  coord_cartesian(ylim=c(0,1)) +
  theme(text=element_text(size=12), axis.text=element_text(size=10), plot.title=element_text(size=16),
        legend.title=element_blank(), legend.justification=c(1,0), legend.position=c(0.9,0.8)) +
  geom_text(x=6, y=0.05, label=paste0("p=",pval), size=5, family="Helvetica")

############### Survival under Therapy #################
biological_age = as.numeric(biological_age)
# cox_tcga <- coxph(SurvObj ~ biological_age*chemo + tumor_stage + gender + smoking + race)
cox_tcga <- coxph(SurvObj ~ biological_age*chemo)
cox_results = summary(cox_tcga)$coefficients

# cox_tcga <- coxph(SurvObj ~ age*chemo + tumor_stage + gender + smoking + race)
cox_tcga <- coxph(SurvObj ~ age*chemo)
cox_results = summary(cox_tcga)$coefficients

