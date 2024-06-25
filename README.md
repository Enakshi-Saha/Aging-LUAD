# Aging-associated Alterations in the Gene Regulatory Network Landscape associate with Risk, Prognosis and Response to Therapy in Lung Adenocarcinoma

This repository contains code for replicating the analysis presented in the paper titled "Aging-associated Alterations in the Gene Regulatory Network Landscape associate with Risk, Prognosis and Response to Therapy in Lung Adenocarcinoma": doi: TBD

## Discovery Data
We downloaded uniformly processed RNA-Seq data from the Recount3 database for the following datasets on May 26, 2022: (i) healthy lung tissue samples from the Genotype Tissue Expression (GTEx) Project (version 8) and (ii) lung adenocarcinoma (LUAD) samples from The Cancer Genome Atlas (TCGA). Phenotypic data for GTEx including age, sex, race, smoking status etc. were downloaded from the dbGap website (https://dbgap.ncbi.nlm.nih.gov/) under study accession phs000424.v8.p2. Extensive phenotypic data for TCGA and clinical information for both the discovery datasets GTEx and TCGA were obtained from Recount3.

## Validation Data
We used the following datasets for validating our findings from the discovery datasets mentioned above. Two independent studies from the Gene Expression Omnibus (GEO) were used as validation datasets: GSE47460 (or, LGRC) and GSE68465. We downloaded the LGRC datset on Feb 12, 2023. From this dataset, we used only 108 samples (59 female and 49 male) that had been annotated as “control” samples.

## Constructing Individual-specific Gene Regulatory Networks
We implemented the PANDA and LIONESS algorithms using Python package netzooPy (version 0.9.10) to construct individual sample-specific gene regulatory networks from all the discovery and validation datasets. Along with the gene expression data, two other data sources were integrated for constructing the networks: Transcription factor/target gene regulatory prior (derived by mapping Transcription factor motifs from the Catalog of Inferred Sequence Binding Preferences (CIS-BP) to the promoter of target genes) and protein-protein interaction (using the interaction scores from StringDb v11.5 between all Transcription factor in the regulatory prior).

The networks are stored in an Amazon Web Services s3 bucket and will be made available upon reasonable request.

## Code
R code for replicating the analysis are documented in README.txt
