library(data.table)
library(ggplot2)
library(dplyr)

# Load indegree
indegree = data.frame(fread("/home/esaha/Lung_lioness/validation_networks/GSE47460/GSE47460_inDegree.txt"))
head(indegree[,1:4])
genes = indegree$V1

rownames(indegree) = genes
indegree = indegree[,-1]

# Load phenotypes
phenotypes = read.csv("/home/esaha/validation_data/GSE47460/GSE47460_phenotypes.txt", sep="")
phenotypes = phenotypes[match(colnames(indegree), rownames(phenotypes)),]

age = phenotypes$age
gender = factor(phenotypes$sex, levels = c("Male", "Female"))
smoking = phenotypes$smoking_status
smoking[is.na(smoking)] = "NA"
smoking = factor(smoking, levels = c("Never", "Ever", "NA"))

#########################################
# GTEx age trajectory

tb <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_LGRC_aging.txt", sep="")

increasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t<0)]

indegree_increasing = apply(indegree[which(rownames(indegree) %in% increasing_gene),], MARGIN=2, mean)
indegree_decreasing = apply(indegree[which(rownames(indegree) %in% decreasing_gene),], MARGIN=2, mean)

fulldata = data.frame(age)
fulldata$smoking_status = as.character(smoking)
fulldata$smoking_status[which(fulldata$smoking_status == "Ever")] = "smoker"
fulldata$smoking_status[which(fulldata$smoking_status == "Never")] = "nonsmoker"

fulldata$gender = gender
fulldata$indegree_increasing = indegree_increasing
fulldata$indegree_decreasing = indegree_decreasing

age_cuts = quantile(age, probs = seq(0, 1, 0.05))
age_cuts[1] = age_cuts[1]-1
age_cuts[length(age_cuts)] = age_cuts[length(age_cuts)] + 1
mean_interval = function(x){
  y = rep(0, (length(x)-1))
  for (i in 2:length(x)){
    y[i-1] = (x[i-1] + x[i])/2
  }
  return(y)
}
fulldata$age_group = cut(fulldata$age,
                         breaks=age_cuts,
                         labels=mean_interval(age_cuts))
fulldata$age_group = as.numeric(as.character(fulldata$age_group))
fulldata = fulldata[-which(fulldata$smoking_status == "NA"),]
head(fulldata)

fulldata = fulldata[-which(fulldata$age > 80 | fulldata$age < 44),]

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, smoking_status)  %>%
  summarise(mean_pathScore = mean(indegree_increasing))
head(fulldata_aggregated)

fulldata_aggregated$smoking_status = factor(fulldata_aggregated$smoking_status, levels = c("smoker", "nonsmoker"))

# Age trajectory by smoking status
ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = smoking_status)) +
  geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title="LGRC: Genes increasingly targeted with age ",
                                                             y="gene score", x = "age")

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, smoking_status)  %>%
  summarise(mean_pathScore = mean(indegree_decreasing))
head(fulldata_aggregated)

fulldata_aggregated$smoking_status = factor(fulldata_aggregated$smoking_status, levels = c("smoker", "nonsmoker"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = smoking_status)) +
  geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title="LGRC: Genes decreasingly targeted with age ",
                                                             y="gene score", x = "age")

