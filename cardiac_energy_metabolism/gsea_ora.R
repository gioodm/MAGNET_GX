#=============================================================================#
# GSEA and ORA                                                                #
#=============================================================================#
library(BiocManager)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(wordcloud)
library(ggnewscale)

setwd("~/Desktop/Uni/MSB/Period 3")

#GSEA
de_results <- read.table("FvsM_PathwayAnalysis_filtered.txt", header = TRUE, sep = "\t")

male_geneList <- de_results$DCMvsControl_male_logFC
names(male_geneList) <- de_results$Ensembl_GeneID

female_geneList <- de_results$DCMvsControl_female_logFC
names(female_geneList) <- de_results$Ensembl_GeneID

#sort the lists in decreasing order (needed for GSEA)
male_geneList = sort(male_geneList, decreasing = TRUE)
female_geneList = sort(female_geneList, decreasing = TRUE)

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#perform GSEA

gse_male_filtered <- gseGO(geneList=male_geneList, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

dotplot(gse_male, showCategory=36, split=".sign") + facet_grid(.~.sign)

gse_male2 <- pairwise_termsim(gse_male_filtered)
emapplot(gse_male2)

gse_female_filtered <- gseGO(geneList=female_geneList, 
                  ont ="BP", 
                  keyType = "ENSEMBL", 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = 0.05, 
                  verbose = TRUE, 
                  OrgDb = organism, 
                  pAdjustMethod = "BH")

dotplot(gse_female, showCategory=10, split=".sign") + facet_grid(.~.sign)

gse_female2 <- pairwise_termsim(gse_female_filtered)
emapplot(gse_female2)

# retrieving the pathways that are enriched in females and males
pathways_female <- gse_female_filtered@result$Description
pathways_male <- gse_male_filtered@result$Description

combined_pathways <- cbind(pathways_female, pathways_male)

#find pathways that are only differentially expressed in females/ males
only_females <- setdiff(combined_pathways[,1], combined_pathways[, 2])
only_males <- setdiff(combined_pathways[,2], combined_pathways[, 1])

#create network plots for pathways only enriched in females/ males
emapplot(gse_male2, showCategory = only_males)
emapplot(gse_female2, showCategory = only_females)

#ORA
male_de_genes <- read.table("ORA_M_significant.txt", header = FALSE, sep = "\t")



male_geneList <- data.frame(de_results$Ensembl_GeneID ,de_results$DCMvsControl_male_logFC, de_results$DCMvsControl_male_adj.P.Val)
male_geneList = male_geneList[order(male_geneList$de_results.DCMvsControl_male_logFC, decreasing = TRUE), ]

#Extract significant genes
male_sig_genes_df = subset(male_geneList, male_geneList$de_results.DCMvsControl_male_adj.P.Val< 0.05)

male_genes <- male_sig_genes_df$de_results.DCMvsControl_male_logFC

names(male_genes) <- male_sig_genes_df$de_results.Ensembl_GeneID

go_enrich <- enrichGO(gene = male_de_genes,
                      universe = male_geneList$de_results.Ensembl_GeneID,
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH")

library(enrichplot)
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
