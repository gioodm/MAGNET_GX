# =============================================================================#
                   #
# =============================================================================#



#-----------------------------------------------------------------------------#
# Setup                                                                       #
#-----------------------------------------------------------------------------#
# Clear workspace
rm(list = ls())

# !# Replace below (in quotes) the working directory (wd) with the location
# of MAGNET_GX folder on your computer
wd <- "/Users/michalskawinski/Desktop/Nauka/Studies/Systems Biology/Period 3 - research project 1/"
setwd(wd)

# Install and load required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("Biobase")) BiocManager::install("Biobase")
if (!require("pcaMethods")) BiocManager::install("pcaMethods")
if (!require("limma")) BiocManager::install("limma")
if (!require("biomaRt")) BiocManager::install("biomaRt")
if (!require("edgeR")) BiocManager::install("edgeR")





#-----------------------------------------------------------------------------#
# 1. Data import                                                              #
#-----------------------------------------------------------------------------#
# a. Import the data

# Import gene expression data as data frame
gxData <- read.table("MAGNET_GX/MAGNET_GeneExpressionData_CPM_19112020.txt",
                     sep = "\t",
                     header = TRUE, row.names = 1
)

# Import sample information as data frame
sampleInfo <- read.table("MAGNET_GX/MAGNET_SampleData_19112020.txt", sep = "\t", header = TRUE)

# Import gene total gene exon lengths
geneTotExonLengths <- read.delim("MAGNET_GX/MAGNET_exonLengths.txt",
                                 as.is = T,
                                 row.names = 1
)

# b. Export a summary table and/or figure(s) of participant characteristics

sampleInfo_DCM <- sampleInfo[sampleInfo$Disease == "DCM",]
sampleInfo_Donor <- sampleInfo[sampleInfo$Disease == "Donor",]

#Male
mean(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Male"])
mean(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Male"])

hist(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Male"])
hist(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Male"])

sd(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Male"])
sd(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Male"])

#Female
mean(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Female"])
mean(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Female"])

hist(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Female"])
hist(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Female"])

sd(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Female"])
sd(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Female"])

#Significant Males Age DCM vs Donor?   -> not significant
shapiro.test(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Male"]) # not
shapiro.test(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Male"]) # not
wilcox.test(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Male"], sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Male"]) #not significant 

#Significant Females Age DCM vs Donor?  -> significant
shapiro.test(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Female"]) # not
shapiro.test(sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Female"]) # not
wilcox.test(sampleInfo_DCM$Age[sampleInfo_DCM$Sex == "Female"], sampleInfo_Donor$Age[sampleInfo_Donor$Sex == "Female"]) #significant 

# Create a summary table of participants
participant_characteristics <- sampleInfo %>%
  group_by(Disease) %>%
  summarize(
    Total_n = n(),
    Female_n = sum(Sex == "Female"),
    Male_n = sum(Sex == "Male"),
    African.American_n = sum(Ethnicity == "African.American"),
    Caucasian_n = sum(Ethnicity == "Caucasian"),
    #Age_mean_f = mean(),
    #Age_mean_m = mean(sampleInfo$Age[sampleInfo$Sex == "Male"]),
    Age_mean = mean(Age),
  )

# View table
#participant_characteristics

# Export table as csv file
#write.table(participant_characteristics, file = "1b_participant_characteristics.csv")

# Create a summary plot of participants in DCM and Donors
sampleInfo_DCM_Donor <- sampleInfo[sampleInfo$Disease %in% c("DCM","Donor"),]
participant_characteristics_plot <- ggplot(sampleInfo_DCM_Donor, aes(x = Age, y = SampleID, color = Sex)) +
  geom_point(position = "stack", alpha = 0.7) +
  facet_wrap(~Disease, scales = "free") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  ylab("Samples")

# View plot
participant_characteristics_plot
#For Females: Donors are older


# Export plot as png file
#png(file = "1b_participant_characteristics_plot.png")
#participant_characteristics_plot
#dev.off()

# c. Transform the data to FPKM values (use the function from the skills trainings)
all(rownames(geneTotExonLengths) == rownames(gxData)) # TRUE (just a check)
cpm2fpkm <- function(x) {
  .t <- 2^(x) * 1E3 / geneTotExonLengths[, 1]
}
gxData_fpkm <- cpm2fpkm(gxData)



#-----------------------------------------------------------------------------#
# 2A. Statistical analysis - NOT FILTERED                                         #
#-----------------------------------------------------------------------------#
# Add rownames to sampleInfo, to be able to later match by samples IDs
rownames(sampleInfo) <- sampleInfo$SampleID

# Create ExpressionSet object with CPM and sampleInfo data
eset <- ExpressionSet(
  assayData = data.matrix(gxData),
  phenoData = AnnotatedDataFrame(sampleInfo)
)

# Delete rownames from sampleInfo
rownames(sampleInfo) <- NULL

# DIFFERENTIAL EXPRESSION ANALYSIS
#Create 2x2 groups
group <- with(pData(eset),
              paste(Disease, Sex, sep = "."))
group <- factor(group)

# Create design matrix for all groups
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Make interesting contrasts 
contr.matrix <- makeContrasts(Female_Disease = DCM.Female - Donor.Female,
                              Male_Disease = DCM.Male - Donor.Male,
                              levels = design)

# Fit coefficients
fit <- lmFit(eset, design)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrasts = contr.matrix)

# Calculate t-statistics
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
summary(results)

#-----------------------------------------------------------------------------#
# 2B. Statistical analysis -  FILTERED                                         #
#-----------------------------------------------------------------------------#
rownames(sampleInfo) <- sampleInfo$SampleID

eset <- ExpressionSet(
  assayData = data.matrix(gxData),
  phenoData = AnnotatedDataFrame(sampleInfo)
)

rownames(sampleInfo) <- NULL

group <- with(pData(eset),
              paste(Disease, Sex, sep = "."))
group <- factor(group)

# Create design matrix for all groups
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

#FILTERED
gxData_DCM_Donor <- gxData[,sampleInfo$Disease %in% c("DCM","Donor")]
dge <- DGEList(counts = gxData_DCM_Donor)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
#create filtered eset
gxData_DCM_Donor_filtered <- gxData_DCM_Donor %>%
  mutate(geneID = rownames(gxData_DCM_Donor)) %>%
  filter(geneID %in% rownames(dge$counts)) %>%
  dplyr::select(-geneID)

sampleInfo_DCM_Donor <- sampleInfo %>%
  filter(Disease %in% c("DCM","Donor"))
rownames(sampleInfo_DCM_Donor) <- sampleInfo_DCM_Donor$SampleID
sampleInfo_DCM_Donor <- sampleInfo_DCM_Donor %>%
  dplyr::select(-SampleID)

eset_filtered <- ExpressionSet(
  assayData = data.matrix(gxData_DCM_Donor_filtered),
  phenoData = AnnotatedDataFrame(sampleInfo_DCM_Donor)
)

group_filtered <- with(pData(eset_filtered),
              paste(Disease, Sex, sep = "."))
group_filtered <- factor(group_filtered)

design_filtered <- model.matrix(~ 0 + group_filtered)
colnames(design_filtered) <- levels(group_filtered)

# Make interesting contrasts 
contr.matrix <- makeContrasts(Female_Disease = DCM.Female - Donor.Female,
                              Male_Disease = DCM.Male - Donor.Male,
                              levels = design_filtered)

# Fit coefficients
fit <- lmFit(eset_filtered, design_filtered)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrasts = contr.matrix)

# Calculate t-statistics
fit2 <- eBayes(fit2)

# Summarize results
results <- decideTests(fit2)
summary(results)



###################################################
# 3. CREATE OUTPUT FILES
###################################################
# Get results
results_df <- as.data.frame(results@.Data)
results_df <- results_df %>%
  mutate(geneID = rownames(results@.Data)) 

###Create output files
#FOR ORA
top_table_F <- topTable(fit2, coef = c(1), n=Inf) 
top_table_F <- top_table_F %>%
  mutate(geneID = rownames(top_table_F))
top_table_F_significant <- top_table_F %>%
  filter(adj.P.Val < 0.05 )
top_table_F_not_significant <- top_table_F %>%
  filter(adj.P.Val > 0.05 )

top_table_M <- topTable(fit2, coef = c(2), n=Inf) 
top_table_M <- top_table_M %>%
  mutate(geneID = rownames(top_table_M))
top_table_M_significant <- top_table_M %>%
  filter(adj.P.Val < 0.05 )
top_table_M_not_significant <- top_table_M %>%
  filter(adj.P.Val > 0.05 )

top_table_significant <- top_table_F_significant %>%
  full_join(top_table_M_significant, by=c("geneID"), suffix=c(".F", ".M"))


#Sig_F, nSig_M
top_table_F_significant <- top_table_F_significant %>% 
  inner_join(top_table_M_not_significant, by = c("geneID"), suffix = c(".F", ".M")) 
cat("", file="ORA_F_significant.txt", append=FALSE)
for(i in 1:length(top_table_F_significant$geneID)) {
  cat(paste(top_table_F_significant$geneID[[i]]),sep="\n", file="ORA_F_significant.txt", append=TRUE)
}

#Sig_M, nSig_F
top_table_M_significant <- top_table_M_significant %>% 
  inner_join(top_table_F_not_significant, by = c("geneID"), suffix = c(".F", ".M")) 

cat("", file="ORA_M_significant.txt", append=FALSE)
for(i in 1:length(top_table_M_significant$geneID)) {
  cat(paste(top_table_M_significant$geneID[[i]]),sep="\n", file="ORA_M_significant.txt", append=TRUE)
}

#Sig_F, all_M
top_table_F_significant <- top_table_F %>%
  filter(adj.P.Val < 0.05 )
top_table_F_significant <- top_table_F_significant %>% 
  left_join(top_table_M_not_significant, by = c("geneID"), suffix = c(".F", ".M"))
cat("", file="ORA_F_significant_M_all.txt", append=FALSE)
for(i in 1:length(top_table_F_significant$geneID)) {
  cat(paste(top_table_F_significant$geneID[[i]]),sep="\n", file="ORA_F_significant_M_all.txt", append=TRUE)
}

#GSEA: 

ensembl <- useEnsembl(
  dataset = "hsapiens_gene_ensembl",
  biomart = "genes"
)
gene_ids <- top_table_significant$geneID
atribbuteNames <- c("ensembl_gene_id", "external_gene_name")
genes <- getBM(
  attributes = atribbuteNames,
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)
if( (nrow(genes) == nrow(top_table_significant)) == FALSE) {
  full_genes <- as.data.frame(gene_ids) %>%
    full_join(genes, by = c("gene_ids" = "ensembl_gene_id"))
}

top_table_significant <- top_table_significant %>% 
  mutate(GeneSymbol = full_genes$external_gene_name) %>%
  dplyr::select(
    Ensembl_GeneID = geneID,
    GeneSymbol,
    DCMvsControl_male_logFC	= logFC.M,
    DCMvsControl_male_P.Value	= P.Value.M,
    DCMvsControl_male_adj.P.Val	= adj.P.Val.M,
    DCMvsControl_female_logFC	= logFC.F,
    DCMvsControl_female_P.Value	= P.Value.F,
    DCMvsControl_female_adj.P.Val = adj.P.Val.F
  )
write.table(top_table_significant, file = "FvsM_PathwayAnalysis_significant.txt", quote = FALSE, sep = "\t", row.names = FALSE)





#SINGLE FILES
#Only Female genes
diff_genes_in_F <- data.frame(Female_Disease = numeric(),
                                 Male_Disease = numeric(),
                                 geneID = character())
diff_genes_in_F_Up <- data.frame(Female_Disease = numeric(),
                              Male_Disease = numeric(),
                              geneID = character())
diff_genes_in_F_Down <- data.frame(Female_Disease = numeric(),
                                 Male_Disease = numeric(),
                                 geneID = character())

#Only Males genes
diff_genes_in_M <- data.frame(Female_Disease = numeric(),
                                 Male_Disease = numeric(),
                                 geneID = character())
diff_genes_in_M_Up <- data.frame(Female_Disease = numeric(),
                              Male_Disease = numeric(),
                              geneID = character())
diff_genes_in_M_Down <- data.frame(Female_Disease = numeric(),
                                 Male_Disease = numeric(),
                                 geneID = character())

#Change (-1f/1m)
diff_genes_change <- data.frame(Female_Disease = numeric(),
                              Male_Disease = numeric(),
                              geneID = character())
diff_genes_change_Mhigher <- data.frame(Female_Disease = numeric(),
                                Male_Disease = numeric(),
                                geneID = character())
diff_genes_change_Fhigher <- data.frame(Female_Disease = numeric(),
                                Male_Disease = numeric(),
                                geneID = character())



row_n1 <- 1
row_n2 <- 1
row_n3 <- 1
row_n4 <- 1
row_n5 <- 1
row_n6 <- 1
row_n7 <- 1
row_n8 <- 1
row_n9 <- 1

for(i in 1:nrow(results_df)) {
  if(results_df$Female_Disease[i] == 1 && results_df$Male_Disease[i] == 0) {
    diff_genes_in_F[row_n1,] <- results_df[i,]
    diff_genes_in_F_Up[row_n2,] <- results_df[i,]
    row_n1 <- row_n1 + 1
    row_n2 <- row_n2 + 1
  } 
  else if (results_df$Female_Disease[i] == -1 && results_df$Male_Disease[i] == 0) {
    diff_genes_in_F[row_n1,] <- results_df[i,]
    diff_genes_in_F_Down[row_n3,] <- results_df[i,]
    row_n1 <- row_n1 + 1
    row_n3 <- row_n3 + 1
  }
  else if (results_df$Female_Disease[i] == 0 && results_df$Male_Disease[i] == 1) {
    diff_genes_in_M[row_n4,] <- results_df[i,]
    diff_genes_in_M_Up[row_n5,] <- results_df[i,]
    row_n4 <- row_n4 +1
    row_n5 <- row_n5 +1
  }
  else if (results_df$Female_Disease[i] == 0 && results_df$Male_Disease[i] == -1) {
    diff_genes_in_M[row_n4,] <- results_df[i,]
    diff_genes_in_M_Down[row_n6,] <- results_df[i,]
    row_n4 <- row_n4 +1
    row_n6 <- row_n6 +1
  }
  else if(results_df$Female_Disease[i] == 1 && results_df$Male_Disease[i] == -1) {
    diff_genes_change[row_n7,] <- results_df[i,]
    diff_genes_change_Fhigher[row_n8,] <- results_df[i,]
    row_n7 <- row_n7 +1
    row_n8 <- row_n8 +1
  } 
  else if (results_df$Female_Disease[i] == -1 && results_df$Male_Disease[i] == 1) {
    diff_genes_change[row_n7,] <- results_df[i,]
    diff_genes_change_Mhigher[row_n9,] <- results_df[i,]
    row_n7 <- row_n7 +1
    row_n9 <- row_n9 +1
  }
}

# Plot differential expression
#plotMD(fit2, column = 3, status = results[, 3])

#Save all the files 
cat("", file="diff_genes_in_F.txt", append=FALSE)
cat("", file="diff_genes_in_F_Up.txt", append=FALSE)
cat("", file="diff_genes_in_F_Down.txt", append=FALSE)
cat("", file="diff_genes_in_M.txt", append=FALSE)
cat("", file="diff_genes_in_M_Up.txt", append=FALSE)
cat("", file="diff_genes_in_M_Down.txt", append=FALSE)
cat("", file="diff_genes_change.txt", append=FALSE)
cat("", file="diff_genes_change_Fhigher.txt", append=FALSE)
cat("", file="diff_genes_change_Mhigher.txt", append=FALSE)

for(i in 1:length(diff_genes_in_F$geneID)) {
  cat(paste(diff_genes_in_F$geneID[[i]]),sep="\n", file="diff_genes_in_F.txt", append=TRUE)
}
for(i in 1:length(diff_genes_in_F_Up$geneID)) {
  cat(paste(diff_genes_in_F_Up$geneID[[i]]),sep="\n", file="diff_genes_in_F_Up.txt", append=TRUE)
}
for(i in 1:length(diff_genes_in_F_Down$geneID)) {
  cat(paste(diff_genes_in_F_Down$geneID[[i]]),sep="\n", file="diff_genes_in_F_Down.txt", append=TRUE)
}
for(i in 1:length(diff_genes_in_M$geneID)) {
  cat(paste(diff_genes_in_M$geneID[[i]]),sep="\n", file="diff_genes_in_M.txt", append=TRUE)
}
for(i in 1:length(diff_genes_in_M_Up$geneID)) {
  cat(paste(diff_genes_in_M_Up$geneID[[i]]),sep="\n", file="diff_genes_in_M_Up.txt", append=TRUE)
}
for(i in 1:length(diff_genes_in_M_Down$geneID)) {
  cat(paste(diff_genes_in_M_Down$geneID[[i]]),sep="\n", file="diff_genes_in_M_Down.txt", append=TRUE)
}
for(i in 1:length(diff_genes_change$geneID)) {
  cat(paste(diff_genes_change$geneID[[i]]),sep="\n", file="diff_genes_change.txt", append=TRUE)
}
for(i in 1:length(diff_genes_change_Fhigher$geneID)) {
  cat(paste(diff_genes_change_Fhigher$geneID[[i]]),sep="\n", file="diff_genes_change_Fhigher.txt", append=TRUE)
}
for(i in 1:length(diff_genes_change_Mhigher$geneID)) {
  cat(paste(diff_genes_change_Mhigher$geneID[[i]]),sep="\n", file="diff_genes_change_Mhigher.txt", append=TRUE)
}


#-----------------------------------------------------------------------------#
# 4. Additional gene annotation                                               #
#-----------------------------------------------------------------------------#
# a.Add at least gene symbols, gene names and gene descriptions to the data
# based on the provided Ensembl gene identifiers.

# Create biomart object with genes database from Hsapiens database with human genes
# from GRCh38.p13 version of human genome
ensembl <- useEnsembl(
  dataset = "hsapiens_gene_ensembl",
  biomart = "genes"
)
# Retrive genes IDs
gene_ids1 <- diff_genes_in_DCM[,3]



# Choose attributes
atribbuteNames <- c("ensembl_gene_id", "external_gene_name", "chromosome_name", "description")

# Retrieve gene annotations from biomart
genes <- getBM(
  attributes = atribbuteNames,
  filters = "ensembl_gene_id",
  values = gene_ids,
  mart = ensembl
)

# Check: do all genes have annotations?
nrow(genes) == nrow(diff_genes_in_DCM)

# If check is FALSE, we need to impute missing annotations for genes with NA values
if( (nrow(genes) == nrow(diff_genes_in_DCM)) == FALSE) {
  full_genes <- as.data.frame(gene_ids) %>%
    full_join(genes, by = c("gene_ids" = "ensembl_gene_id"))
}

# Glimpse into genes annotations
glimpse(full_genes)


#Print all genes ID
for(i in 1:length(diff_genes_in_DCM$geneID)) {
  cat(paste(diff_genes_in_DCM$geneID[[i]]),sep="\n")
}




#-----------------------------------------------------------------------------#
# 6. Export the results                                                       #
#-----------------------------------------------------------------------------#
# a. The exported file with results should be a tab-delimited text file.
# b. The file should contain all the additionally generated data with clear column names.
# c. Do not include the original log2-transformed CPM values (or FPKM values) for each sample,
# to keep the results file manageable. However, for visualization purposes (e.g. pathway analysis)
# it is useful to have an average expression value for DCM and an average expression value
# for controls for each gene.

# Add additional columns with average expression of each gene for DCM and control groups
DCMvsDonor_1 <- DCMvsDonor %>%
  mutate(
    AveExpr_DCM = rowMeans(exprs(eset)[, sampleInfo$Disease == "DCM"]),
    AveExpr_Donor = rowMeans(exprs(eset)[, sampleInfo$Disease == "Donor"])
  )

# Add a columns with differential expression outcome (Up=1/Down=-1/NotSign=0)
DCMvsDonor_2 <- DCMvsDonor_1 %>%
  mutate(DCMvsDonor = as.numeric(results))

# Add a column with genes IDs to DCMvsDonor_2 to be able to access them
DCMvsDonor_2$gene_ids <- rownames(DCMvsDonor_2)

# Remove row names from DCMvsDonor_2
rownames(DCMvsDonor_2) <- NULL

# Add annotations for each gene
results_output <- full_genes_1 %>%
  # join with annotations
  inner_join(DCMvsDonor_2, by = "gene_ids") %>%
  # and select some columns (and change their names)
  dplyr::select(
    EnsemblId = gene_ids,
    GeneName = external_gene_name,
    Chromosome = chromosome_name,
    Description = description,
    DCMvsDonor, logFC, P.Value, adj.P.Val,
    AveExpr, ExprAboveBackground, AveExpr_DCM, AveExpr_Donor
  )

# Glimpse into created object
glimpse(results_output)

# Export object
#write.table(results_output, file = "6_results.txt", sep = "\t")
