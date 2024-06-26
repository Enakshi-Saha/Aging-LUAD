# Must set seed
# rm(list = ls())
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

##########################
gene_name = "MET"

GTEx_male = get(load(paste("/home/esaha/Lung_lioness/TPM/GTEx/", gene_name, "_subnet_GTEx_male.RData", sep = "")))
GTEx_female = get(load(paste("/home/esaha/Lung_lioness/TPM/GTEx/", gene_name, "_subnet_GTEx_female.RData", sep = "")))

########## limma analysis to find significant sex-difference in TFs targeting EGFR ####
edges = GTEx_male[,1]
gene = rep(gene_name, length(edges))
TF = unlist(lapply(strsplit(edges, split = "_"), function(x){x[1]}))

sample_male.GTEx = colnames(GTEx_male)[-1]
sample_female.GTEx = colnames(GTEx_female)[-1]

sample.GTEx = c(sample_male.GTEx, sample_female.GTEx)

GTEx_IDs = unlist(lapply(strsplit(sample.GTEx, split=".", fixed = T),
                         function(x){paste(x[-length(x)], collapse ="-")}))

# Get phenotypes: GTEx
pheno_info.GTEx = get(load("/home/esaha/Recount3_GTEx_TCGA_Data/GTEx_eset.RData"))
pheno_info.GTEx = pData(pheno_info.GTEx)

GTEx_phenotypes = as.data.frame(cbind(pheno_info.GTEx$SEX, pheno_info.GTEx$AGE, pheno_info.GTEx$RACE, pheno_info.GTEx$MHSMKSTS,
                                      pheno_info.GTEx$SMNABTCHT, pheno_info.GTEx$TRISCH,
                                      pheno_info.GTEx$SMGEBTCHD, pheno_info.GTEx$SMRIN,
                                      pheno_info.GTEx$MHSMKYRS))
GTEx_phenotypes =  GTEx_phenotypes[match(GTEx_IDs, rownames(pheno_info.GTEx)),]
colnames(GTEx_phenotypes) = c("gender", "age", "race", "smoking", "batch", "ischemic_time", "batch_date", "rna_degrad", "smoke_years")
GTEx_phenotypes$batch_date = unlist(lapply(strsplit(GTEx_phenotypes$batch_date, split = "/"), function(x){x[3]}))

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
smoking_status[which(is.na(smoking_status))] = "Unknown"
smoking_status = factor(smoking_status, levels = c("No", "Yes", "Unknown"))

batch = GTEx_phenotypes$batch

ischemic_time = GTEx_phenotypes$ischemic_time
hour_minute = strsplit(ischemic_time, split=" ")
ischemic_timeH = sapply(hour_minute, function(x){as.numeric(as.character(x[1]))+as.numeric(as.character(x[1]))/60})

rin = as.numeric(as.character(GTEx_phenotypes$rna_degrad))

design = model.matrix(~ age + gender + race + smoking_status*gender + batch + rin +ischemic_timeH)
GTEx_edges = cbind(GTEx_male[,-1], GTEx_female[,-1])
rownames(GTEx_edges) = TF
fit <- lmFit(GTEx_edges, design)
fit <- eBayes(fit)

# Save table for age coefficients
tb = topTable(fit,coef="age",number=Inf)
head(tb)

#### Draw Graphs #########
library(visNetwork)

subnet = data.frame(cbind(tb$t, tb$P.Value))
subnet$TF = TF
subnet$gene = gene
colnames(subnet) = c("GTEx_t", "GTEx_pval", "TF", "gene")
head(subnet)

# Network Plot
edges = subnet[which(tb$P.Value < 0.05),c("TF", "gene", "GTEx_t")]
edges$arrows = "to" 
colnames(edges) <- c("from","to","force","arrows")

edges0=edges
edges = edges[1:50,]

nodes <- data.frame(id = unique(as.vector(as.matrix(edges[,c(1,2)]))) , 
                    label=unique(as.vector(as.matrix(edges[,c(1,2)]))))
nodes$group <- ifelse(nodes$id %in% edges$from, "TF", "gene")

# Color male edges blue and female edges red
edge_col = rep("red", nrow(edges))
edge_col[which(sign(edges[,3]) == -1)] = "blue"

edges$color = edge_col

net <- visNetwork(nodes, edges, width = "100%", 
                  main = paste("GTEx: Change in ", gene_name, " Network with Age", sep = ""))
net <- visGroups(net, groupname = "TF", shape = "dot",
                 color = list(background = "teal", border="black"))
net <- visGroups(net, groupname = "gene", shape = "square",       
                 color = list(background = "gold", border="black"))
visLegend(net, main="Legend", position="right", ncol=1) 

# Save plot as PNG
html_name <- tempfile(fileext = ".html")
visSave(net, html_name)

#library(webshot); #webshot::install_phantomjs() #in case phantomjs was not installed 

#webshot(html_name, zoom = 2, file = "/home/esaha/paper_plots/network_plots/TCGA_EGFR_limma_sigEdges.png")
