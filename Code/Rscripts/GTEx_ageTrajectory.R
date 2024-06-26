library(data.table)
library(ggplot2)
library(dplyr)

# Load indegree matrices
indegree_GTEx_male = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_male.txt"))
indegree_GTEx_female = data.frame(fread("/home/esaha/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/indegree/GTEx_lung_lioness_indegree_female.txt"))
genes = indegree_GTEx_male$V1

indegree = cbind(indegree_GTEx_male[,-1], indegree_GTEx_female[,-1])
rownames(indegree) = genes

GTEx_IDs = colnames(indegree)
# Format GTEx IDs to match with phenotypic data: GTEx
GTEx_IDs = unlist(lapply(strsplit(GTEx_IDs, split=".", fixed = T),
                         function(x){paste(x[1:2], collapse= "-")}))

# Get phenotypes: GTEx
GTEx_phenotypes <- data.frame(read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023/GTEx_phenotypes.txt", sep=""))
GTEx_phenotypes = GTEx_phenotypes[match(GTEx_IDs, GTEx_phenotypes$SUBJID),]
GTEx_phenotypes = GTEx_phenotypes[,c("gtex.sex", "AGE", "RACE", "MHSMKSTS", "TRISCH", "gtex.smrin", "gtex.smnabtcht")]
colnames(GTEx_phenotypes) = c("gender", "age", "race", "smoking", "ischemic_time", "rna_degrad", "batch")

# Define the covariates: gender, age, race etc
gender <- GTEx_phenotypes$gender
gender[which(gender == 1)] = "MALE"
gender[which(gender == 2)] = "FEMALE"
gender = factor(gender, levels = c("MALE", "FEMALE"))

age <- as.numeric(as.character(GTEx_phenotypes$age))
age[which(is.na(age))] <- mean(age,na.rm=TRUE)

# Categorize race 1, 4, 99 as single class "others", 2: "black", 3: "white"
race <- as.numeric(as.character(GTEx_phenotypes$race))
race[which(race != 2 & race != 3)] = "others"
race[which(race == 2)] = "black or african american"
race[which(race == 3)] = "white"

smoking_status = GTEx_phenotypes$smoking
smoking_status[-which(smoking_status %in% c("No", "Yes"))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

#########################################
# GTEx age trajectory

tb <- read.csv("~/Recount3_GTEx_TCGA_Data/TPM_expression/newdata2023_aging/limma_GSEA_aging/limma_GTEx_aging.txt", sep="")

increasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t>0)]
decreasing_gene = rownames(tb)[which(tb$P.Value<0.05 & tb$t<0)]

indegree_increasing = apply(indegree[which(rownames(indegree) %in% increasing_gene),], MARGIN=2, mean)
indegree_decreasing = apply(indegree[which(rownames(indegree) %in% decreasing_gene),], MARGIN=2, mean)

fulldata = data.frame(age)
fulldata$smoking_status = as.character(smoking_status)
#fulldata$smoking_status[which(fulldata$smoking_status == "Yes")] = "smoker"
#fulldata$smoking_status[which(fulldata$smoking_status == "No")] = "nonsmoker"

fulldata$gender = gender
fulldata$race = race
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
fulldata = fulldata[-which(fulldata$smoking_status == "Unknown"),]
head(fulldata)

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, smoking_status)  %>%
  summarise(mean_pathScore = mean(indegree_increasing))
head(fulldata_aggregated)

fulldata_aggregated$smoking_status = factor(fulldata_aggregated$smoking_status, levels = c("Yes", "No"))

# Age trajectory by smoking status
ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = smoking_status)) +
  geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title="GTEx: Genes increasingly targeted with age ",
                                                             y="gene score", x = "age")

# Summarize data for each age
fulldata_aggregated = fulldata %>% group_by(age_group, smoking_status)  %>%
  summarise(mean_pathScore = mean(indegree_decreasing))
head(fulldata_aggregated)

fulldata_aggregated$smoking_status = factor(fulldata_aggregated$smoking_status, levels = c("Yes", "No"))

ggplot(fulldata_aggregated, aes(x = age_group , y = mean_pathScore, color = smoking_status)) +
  geom_point(size = 0.5) + geom_smooth(method = "lm") + labs(title="GTEx: Genes decreasingly targeted with age ",
                                                             y="gene score", x = "age")

