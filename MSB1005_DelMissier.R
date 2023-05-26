#=============================================================================#
# MSB1005_DelMissier.R                                                        #
#										                                      									  # 															                                #
# Date: Dec 17, 2021											                                    #
# Author: Giorgia Del Missier, i6292403, Maastricht University                #
#														                                                  # 
#=============================================================================#


#-----------------------------------------------------------------------------#
# Assignment 1: Data import
#-----------------------------------------------------------------------------#


# 1.a: Import the data.

# Current working directory
setwd("~/Desktop/Systems Biology/First year/2.Experimental Design & Data Management/Skills/MAGNET_GX")

# Load data
gxData <- read.table(file = "MAGNET_GeneExpressionData_CPM_19112020.txt", 
                     row.names = 1, header = TRUE)
sampleInfo <- read.table(file = "MAGNET_SampleData_19112020.txt", 
                         header = TRUE)
geneTotExonLengths <- read.table(file = "MAGNET_exonLengths.txt", 
                                 row.names = 1, header = TRUE)


# 1.b: Export a summary table and/or figure(s) of participant characteristics.

library(arsenal)
library(grid)
library(gridExtra)
library(ggplot2)
library(wesanderson)
library(hrbrthemes)
library(tidyr)

# General info about the sampleInfo dataset
summary(sampleInfo)
str(sampleInfo)


# Create summary table
table <- tableby(Disease ~ Sex + Age + Ethnicity, 
                 data = sampleInfo, 
                 test=FALSE, 
                 numeric.stats=c("mean","median", "range"))

summary_table <- as.data.frame(summary(table, text = TRUE))

# Export table as pdf file
pdf(file = "1b - Summary table.pdf", width = 10)

title <- "Summary table - Participants' characteristics"
grid::grid.text(title, x = (0.5), y = (0.8),
                gp = gpar(fontsize = 18, fontface = "bold"))
my_theme <- ttheme_default(core = list(fg_params = list(hjust=0, x=0.1)))
grid.table(summary_table, theme = my_theme, rows = NULL)

dev.off()


# Export png figure of age distribution in the samples
png("1b - Age distribution.png")

ggplot(data = sampleInfo, aes(x=Age, color=Disease, fill=Disease)) +
  geom_histogram(alpha=0.8, binwidth = 6) +
  ggtitle("Age distribution by disease") +
  xlab("Age (years)") +
  ylab("Frequency") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("GrandBudapest2")) +
  scale_color_manual(values = wes_palette("GrandBudapest2")) +
  facet_wrap(~Disease)

dev.off()


# 1.c: Transform the data to FPKM values.

# Function to convert CPM values to FPKM values
cpm2fpkm <- function(x) {
  geneTotExonLengths_kb <- geneTotExonLengths[, 1] / 1E3
  .t <- 2^(x) / geneTotExonLengths_kb
}


# Check whether the order of genes is the same in the two objects
all(rownames(geneTotExonLengths) == rownames(gxData)) # TRUE

# Conversion of dataset CPM values to FPKM values
gxData_fpkm <- cpm2fpkm(gxData)

# Density plot of FPKM values
plotData_fpkm <- gather(gxData_fpkm, key = "SampleID", value = "FPKM")

ggplot(plotData_fpkm, aes(x = FPKM, color = SampleID)) +
  geom_density(alpha = 0.2) + 
  theme(legend.position="none") # the right-skewed distribution indicates 
                                # that most of the values lie near zero

# log2 transformation allows a better visualization of the entire distribution
gxData_log2fpkm <- log2(gxData_fpkm)

# Density plot of log2 transformed FPKM values
plotData_log2fpkm <- gather(gxData_log2fpkm, key = "SampleID", value = "log2_FPKM")

ggplot(plotData_log2fpkm, aes(x = log2_FPKM, color = SampleID)) + 
  geom_density(alpha = 0.2) + 
  theme(legend.position="none") # now the distribution is closer to a normal one


#-----------------------------------------------------------------------------#
# Assignment 2: Diagnostic plots
#-----------------------------------------------------------------------------#


# 2.a: Create at least one figure containing readable boxplots of the gene expression values 
#      and export it.

library(dplyr)

# Check whether the order of samples is the same in the two objects
all(sampleInfo[, 1] == colnames(gxData)) # TRUE


# Create plot data
plotData1 <- gather(gxData, key = "SampleID", value = "CPM")
plotData2 <- sampleInfo %>% slice(rep(1:n(), each = nrow(gxData)))

# Check again whether the order of samples is the same in the two objects and merge them
all(plotData1[, 1] == plotData2[, 1]) # TRUE
plotData_box <- cbind(plotData1, plotData2[, 2:ncol(plotData2)])


# Create boxplot of gene expression values
box_disease_sex <- ggplot(plotData_box, aes(x = Sex, y = CPM, fill = Disease)) + 
  geom_boxplot() +
  ggtitle("Gene expression values by sex for each disease group") +
  xlab("Sex") +
  ylab("CPM value") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("Royal1")) +
  scale_color_manual(values = wes_palette("Royal1"))

box_disease_eth <- ggplot(plotData_box, aes(x = Ethnicity, y = CPM, fill = Disease)) + 
  geom_boxplot() +
  ggtitle("Gene expression values by ethnicity for each disease group") +
  xlab("Ethnicity") +
  ylab("CPM value") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("Royal1")) +
  scale_color_manual(values = wes_palette("Royal1"))

# Export them as png figure
ggsave("2a - Gene expression boxplots.png", 
       arrangeGrob(box_disease_sex, box_disease_eth), 
       height = 10)


# 2.b: Create at least one PCA plot showing the sample clustering colored by all available 
#      co-variates and export it.

library(pcaMethods)

# Perform PCA on gene expression values
pcaRes <- pca(t(gxData), nPcs = 10)
pcaRes@R2
plot(pcaRes, main = "Cumulative explained variance")

# Check again that the order of samples is the same and create plot data
all(rownames(pcaRes@scores) == sampleInfo[, 1]) # TRUE
plotData_pca <- cbind(data.frame(pcaRes@scores), sampleInfo)


# Create PCA plots of gene expression values for different PCs and export them
pca_disease_sex12 <- ggplot(plotData_pca, aes(x = PC1, y = PC2, color = Disease, shape = Sex)) + 
  geom_point() + 
  ggtitle("PC1 vs PC2, by sex and disease group") +
  xlab("PC1 (27.26% of variance explained)") +
  ylab("PC2 (7.35% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  scale_color_manual(values = wes_palette("BottleRocket2"))

pca_disease_eth12 <- ggplot(plotData_pca, aes(x = PC1, y = PC2, color = Disease, shape = Ethnicity)) + 
  geom_point() + 
  ggtitle("PC1 vs PC2, by ethnicity and disease group") +
  xlab("PC1 (27.26% of variance explained)") +
  ylab("PC2 (7.35% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  scale_color_manual(values = wes_palette("BottleRocket2"))

pca_disease_age12 <- ggplot(plotData_pca, aes(x = PC1, y = PC2, color = Age, shape=Disease)) + 
  geom_point() + 
  ggtitle("PC1 vs PC2, by age and disease group") +
  xlab("PC1 (27.26% of variance explained)") +
  ylab("PC2 (7.35% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_color_gradientn(colours = wes_palette("BottleRocket2", 100, type = "continuous"))

ggsave("2b - PCA plots by all covariates (PC1 vs PC2).png", 
       arrangeGrob(pca_disease_sex12, pca_disease_eth12, pca_disease_age12), height = 15)


pca_disease_sex13 <- ggplot(plotData_pca, aes(x = PC1, y = PC3, color = Disease, shape = Sex)) + 
  geom_point() + 
  ggtitle("PC1 vs PC3, by sex and disease group") +
  xlab("PC1 (27.26% of variance explained)") +
  ylab("PC3 (5.63% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  scale_color_manual(values = wes_palette("BottleRocket2"))

pca_disease_eth13 <- ggplot(plotData_pca, aes(x = PC1, y = PC3, color = Disease, shape = Ethnicity)) + 
  geom_point() + 
  ggtitle("PC1 vs PC3, by ethnicity and disease group") +
  xlab("PC1 (27.26% of variance explained)") +
  ylab("PC3 (5.63% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  scale_color_manual(values = wes_palette("BottleRocket2"))

pca_disease_age13 <- ggplot(plotData_pca, aes(x = PC1, y = PC3, color = Age, shape=Disease)) + 
  geom_point() + 
  ggtitle("PC1 vs PC3, by age and disease group") +
  xlab("PC1 (27.26% of variance explained)") +
  ylab("PC3 (5.63% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_color_gradientn(colours = wes_palette("BottleRocket2", 100, type = "continuous"))

ggsave("2b - PCA plots by all covariates (PC1 vs PC3).png", 
       arrangeGrob(pca_disease_sex13, pca_disease_eth13, pca_disease_age13), height = 15)


pca_disease_sex23 <- ggplot(plotData_pca, aes(x = PC2, y = PC3, color = Disease, shape = Sex)) + 
  geom_point() + 
  ggtitle("PC2 vs PC3, by sex and disease group") +
  xlab("PC2 (7.35% of variance explained)") +
  ylab("PC3 (5.63% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  scale_color_manual(values = wes_palette("BottleRocket2"))

pca_disease_eth23 <- ggplot(plotData_pca, aes(x = PC2, y = PC3, color = Disease, shape = Ethnicity)) + 
  geom_point() + 
  ggtitle("PC2 vs PC3, by ethnicity and disease group") +
  xlab("PC2 (7.35% of variance explained)") +
  ylab("PC3 (5.63% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_fill_manual(values = wes_palette("BottleRocket2")) +
  scale_color_manual(values = wes_palette("BottleRocket2"))

pca_disease_age23 <- ggplot(plotData_pca, aes(x = PC2, y = PC3, color = Age, shape=Disease)) + 
  geom_point() + 
  ggtitle("PC2 vs PC3, by age and disease group") +
  xlab("PC2 (7.35% of variance explained)") +
  ylab("PC3 (5.63% of variance explained)") +
  theme_ipsum() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)) +
  scale_color_gradientn(colours = wes_palette("BottleRocket2", 100, type = "continuous"))

ggsave("2b - PCA plots by all covariates (PC2 vs PC3).png", 
       arrangeGrob(pca_disease_sex23, pca_disease_eth23, pca_disease_age23), height = 15)


# None of the co-variates on the first three PCs seem to separate clearly the samples,
# except for disease state.


#-----------------------------------------------------------------------------#
# Assignment 3: Additional gene annotation.
#-----------------------------------------------------------------------------#


# 3.a: Add at least gene symbols, gene names and gene descriptions to the data based on the 
#      provided Ensembl gene identifiers.

library(biomaRt)

# Connect to Ensembl database, version of February 2014
# needed in order to retrieve the annotation of genes whose ID is not in use anymore
ensembl_v75 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version=75)

# Connect to Ensembl database, most recent version (December 2021)
ensembl_v105 <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# Gene names
gene_id <- rownames(gxData)

# Retrieve gene annotation from version 75
gene_info_v75 <- getBM(filters= "ensembl_gene_id", 
                       attributes= c("ensembl_gene_id", "chromosome_name", "band", "description"),
                       values= gene_id,
                       mart= ensembl_v75)

# Check for duplicates
gene_info_v75[duplicated(gene_info_v75$ensembl_gene_id), ] # none

# Final check
length(gene_id[ !gene_id %in% gene_info_v75$ensembl_gene_id]) # no missing genes


# Retrieve gene annotation from version 105
gene_info_v105 <- getBM(filters= "ensembl_gene_id", 
                        attributes= c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id", "external_gene_name"),
                        values= gene_id,
                        mart= ensembl_v105)

# Missing genes from dataset
missing <- gene_id[ !gene_id %in% gene_info_v105$ensembl_gene_id]
length(missing) # 739 genes missing

# Check for duplicates and remove them
gene_info_v105[duplicated(gene_info_v105$ensembl_gene_id), ]
gene_info_v105 <- gene_info_v105[!duplicated(gene_info_v105$ensembl_gene_id), ]
gene_info_v105[duplicated(gene_info_v105$ensembl_gene_id), ] # no more duplicates


# Merge the two datasets, using NAs for missing values
gene_info <- full_join(gene_info_v105, gene_info_v75, by="ensembl_gene_id")
all(gene_info[,1] == rownames(gxData)) # FALSE : full_join() reorders rows, putting at the end the ones 
                                       #         containing missing values (NAs)

# Reorder of rows as in the original gene expression dataset and final check
gene_info <- gene_info %>% slice(match(gene_id, ensembl_gene_id))
all(gene_info[,1] == rownames(gxData)) # TRUE


#-----------------------------------------------------------------------------#
# Assignment 4: Statistical analysis
#-----------------------------------------------------------------------------#


# 4.a: Perform a differential gene expression analysis comparing DCM patients and healthy donors.

library(Biobase)
library(limma)

# Create ExpressionSet object
expr_mat <- as.matrix(gxData)
colnames(expr_mat) <- rownames(sampleInfo)
rownames(gene_info) <- gene_info[, 1]

eset <- ExpressionSet(assayData = expr_mat,
                      phenoData = AnnotatedDataFrame(sampleInfo),
                      featureData = AnnotatedDataFrame(gene_info))

# View the number of genes (rows) and samples (columns)
dim(eset)

# Create design matrix with no intercept
design_mat <- model.matrix(~0 + Disease, data = pData(eset))

# Count the number of samples modeled by each coefficient
colSums(design_mat)

# Create a contrasts matrix
contrasts_mat <- makeContrasts(DCMvDonor = DiseaseDCM - DiseaseDonor,
                               HCMvDonor = DiseaseHCM - DiseaseDonor,
                               DCMvHCM = DiseaseDCM - DiseaseHCM,
                               levels = design_mat)

# View the contrasts matrix
contrasts_mat

# Fit the model
model_fit <- lmFit(eset, design_mat)

# Fit the contrasts
contrasts_fit <- contrasts.fit(model_fit, contrasts = contrasts_mat)

# Calculate the t-statistics for the contrasts
stats <- eBayes(contrasts_fit)

# Summarize results
results <- decideTests(stats)
summary(results)

vennDiagram(results)


# Obtain the summary statistics for the contrast DCM vs Donor
stats_DCMvDonor <- topTable(stats, coef = "DCMvDonor", number = nrow(stats),
                         sort.by = "none")

# Obtain the summary statistics for the contrast HCM vs Donor
stats_HCMvDonor  <- topTable(stats, coef = "HCMvDonor", number = nrow(stats),
                            sort.by = "none")

# Obtain the summary statistics for the contrast DCM vs HCM
stats_DCMvHCM  <- topTable(stats, coef = "DCMvHCM", number = nrow(stats),
                              sort.by = "none")

# Create histograms of the p-values for each contrast
hist(stats_DCMvDonor[, "P.Value"]) 
hist(stats_HCMvDonor[, "P.Value"])
hist(stats_DCMvHCM[, "P.Value"]) 

# A large density of low p-values indicates many differentially expressed genes, 
# while a uniformly distributed histogram indicates there are few. 
# As expected, many genes are differentially expressed in the diseased categories (here DCM and HCM) 
# when compared to the healthy controls; few genes are differentially expressed when comparing two 
# diseased groups (here DCM vs HCM).


# 4.b: Correct for relevant co-variates and add comments to the scripts explaining your decision.

# Create design matrix with no intercept, including all co-variates
design_mat_covar <- model.matrix(~0 + Disease + Age + Sex + Ethnicity, data = pData(eset))

# Create a contrasts matrix
contrasts_mat_covar <- makeContrasts(DCMvDonor = DiseaseDCM - DiseaseDonor,
                                     HCMvDonor = DiseaseHCM - DiseaseDonor,
                                     DCMvHCM = DiseaseDCM - DiseaseHCM,
                                     levels = design_mat_covar)

# Fit the model
model_fit_covar <- lmFit(eset, design_mat_covar)

# Fit the contrasts
contrasts_fit_covar <- contrasts.fit(model_fit_covar, contrasts = contrasts_mat_covar)

# Calculate the t-statistics for the contrasts
stats_covar <- eBayes(contrasts_fit_covar)

# Summarize results
results_covar <- decideTests(stats_covar)
summary(results_covar)

vennDiagram(results_covar)


# Obtain the summary statistics for the contrast DCM vs Donor
stats_DCMvDonor_covar <- topTable(stats_covar, coef = "DCMvDonor", number = nrow(stats_covar),
                            sort.by = "none")

# Obtain the summary statistics for the contrast HCM vs Donor
stats_HCMvDonor_covar <- topTable(stats_covar, coef = "HCMvDonor", number = nrow(stats_covar),
                             sort.by = "none")

# Obtain the summary statistics for the contrast DCM vs HCM
stats_DCMvHCM_covar <- topTable(stats_covar, coef = "DCMvHCM", number = nrow(stats_covar),
                           sort.by = "none")

# Create histograms of the p-values for each contrast
hist(stats_DCMvDonor_covar[, "P.Value"]) 
hist(stats_HCMvDonor_covar[, "P.Value"])
hist(stats_DCMvHCM_covar[, "P.Value"]) 

# The number of genes differentially expressed is slightly different from the one previously
# obtained, when only Disease was used to create the model. However, from the PCA analysis,
# none of the covariates indicated a clear sample separation and hence the need to correct 
# for them. Thus, these results agree with the previous hypothesis and show that including 
# the covariates in the model does not produce significant differences.


# Additional analyses with model only fitted on Disease

# Extract the HGNC gene IDs
gene_hgnc <- stats$genes[, "hgnc_symbol"]

# Create a volcano plot for the contrast DCM vs Donor
volcanoplot(stats, coef = "DCMvDonor", highlight = 5, names = gene_hgnc)

# Create a volcano plot for the contrast HCM vs Donor
volcanoplot(stats, coef = "HCMvDonor", highlight = 5, names = gene_hgnc)

# Create a volcano plot for the contrast DCM vs HCM
volcanoplot(stats, coef = "DCMvHCM", highlight = 5, names = gene_hgnc)


# Extract the Entrez gene IDs
gene_entrez <- stats$genes[, "entrezgene_id"]

# Test for enriched KEGG Pathways for contrast DCM vs Donor
enrich_DCMvDonor <- kegga(stats, coef = "DCMvDonor", geneid = gene_entrez, species = "Hs")
# View the top 5 enriched KEGG pathways
topKEGG(enrich_DCMvDonor, number = 5)

# Test for enriched KEGG Pathways for contrast HCM vs Donor
enrich_HCMvDonor <- kegga(stats, coef = "HCMvDonor", geneid = gene_entrez, species = "Hs")
# View the top 5 enriched KEGG pathways
topKEGG(enrich_HCMvDonor, number = 5)

# Test for enriched KEGG Pathways for contrast DCM vs HCM
enrich_DCMvHCM <- kegga(stats, coef = "DCMvHCM", geneid = gene_entrez, species = "Hs")
# View the top 5 enriched KEGG pathways
topKEGG(enrich_DCMvHCM, number = 5)


#-----------------------------------------------------------------------------#
# Assignment 5: Relative expression levels
#-----------------------------------------------------------------------------#


# 5.a: Assess for each gene in the dataset whether it is expressed above background (noise) level.

# Compute background (noise) level as the average expression of Y chromosome genes in female subjects
expr_bg <- mean(as.matrix(gxData_log2fpkm[gene_info$chromosome_name == "Y", sampleInfo$Sex == "Female"]))

# Compute average expression of each gene
gene_means <- rowMeans(gxData_log2fpkm)

# Determine relative expression value and store it
relative_expr <- c()

for (gene_mean in gene_means) {
  if (gene_mean > expr_bg) {
    relative_expr <- c(relative_expr, "Yes")} 
  else {
    relative_expr <- c(relative_expr, "No")}
}

# Summary
table(relative_expr)


#-----------------------------------------------------------------------------#
# Assignment 6: Export the results
#-----------------------------------------------------------------------------#


# Compute average gene expression value (in CPM) for each disease group category (rounded)
avg_DCM <- round(rowMeans(gxData[ ,sampleInfo$Disease == "DCM"]), 3)
avg_Donors <- round(rowMeans(gxData[ ,sampleInfo$Disease == "Donor"]), 3)
avg_HCM <- round(rowMeans(gxData[ ,sampleInfo$Disease == "HCM"]), 3)
avg_PPCM <- round(rowMeans(gxData[ ,sampleInfo$Disease == "PPCM"]), 3)

# Retrieve Log2 Fold Change and P-value of each differential gene expression performed (rounded)
logFC_DCMvDonor <- round(stats_DCMvDonor$logFC, 3)
pval_DCMvDonor <- round(stats_DCMvDonor$P.Value, 3)

logFC_HCMvDonor <- round(stats_HCMvDonor$logFC, 3)
pval_HCMvDonor <- round(stats_HCMvDonor$P.Value, 3)

logFC_DCMvHCM <- round(stats_DCMvHCM$logFC, 3)
pval_DCMvHCM <- round(stats_DCMvHCM$P.Value, 3)


# Create final table and export it as tab-separated txt file
table_final <- cbind(gene_info, 
                     avg_DCM, avg_Donors, avg_HCM, avg_PPCM,
                     logFC_DCMvDonor, pval_DCMvDonor, 
                     logFC_HCMvDonor, pval_HCMvDonor, 
                     logFC_DCMvHCM, pval_DCMvHCM, 
                     relative_expr)

colnames(table_final) <- c("Gene ID, Ensembl", "Gene ID, HGNC", "Gene ID, Entrez", "External gene name",
                           "Chr", "Band", "Gene description",
                           "DCM mean expression (CPM)", "Donors mean expression (CPM)", 
                           "HCM mean expression (CPM)", "PPCM mean expression (CPM)",
                           "Log2 Fold Change, DCM vs Donor", "p-value, DCM vs Donor",
                           "Log2 Fold Change, HCM vs Donor", "p-value, HCM vs Donor",
                           "Log2 Fold Change, DCM vs HCM", "p-value, DCM vs HCM",
                           "Relative expression level (Yes = expressed above noise level)")
                        
write.table(table_final, file = "6 - Final results.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = FALSE, col.names = TRUE)


#=============================================================================#