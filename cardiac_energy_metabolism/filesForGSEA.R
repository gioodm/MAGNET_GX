

#Load files
gxData <- read.table("MAGNET_GX/MAGNET_GeneExpressionData_CPM_19112020.txt",
                     sep = "\t",
                     header = TRUE, row.names = 1
)

sampleInfo <- read.table("MAGNET_GX/MAGNET_SampleData_19112020.txt", sep = "\t", header = TRUE)
top_table <- read.table("FvsM_PathwayAnalysis.txt", sep="\t", header = TRUE)


#Create databases (.txt)
gxData_DCM_Donor_F <- gxData[,sampleInfo$Sex == "Female" & sampleInfo$Disease %in% c("DCM","Donor")]
gxData_DCM_Donor_M <- gxData[,sampleInfo$Sex == "Male" & sampleInfo$Disease %in% c("DCM","Donor")]


gxData_DCM_Donor_F <- gxData_DCM_Donor_F %>%
  mutate(NAME = rownames(gxData_DCM_Donor_F),
         DESCRIPTION = NA) %>%
  dplyr::select(NAME, DESCRIPTION, 1:length(colnames(gxData_DCM_Donor_F)))

gxData_DCM_Donor_M <- gxData_DCM_Donor_M %>%
  mutate(NAME = rownames(gxData_DCM_Donor_M),
         DESCRIPTION = NA) %>%
  dplyr::select(NAME, DESCRIPTION, 1:length(colnames(gxData_DCM_Donor_M)))

write.table(gxData_DCM_Donor_F,
           file = "GSEA/gxData_DCM_Donor_F.txt",
           sep = "\t",
           row.names = FALSE)
write.table(gxData_DCM_Donor_M,
            file = "GSEA/gxData_DCM_Donor_M.txt",
            sep = "\t",
            row.names = FALSE)

#Print Sample Names to  create phenotype data (.cls)
sampleInfo_DCM_Donor_F <- sampleInfo[sampleInfo$SampleID %in% colnames(gxData_DCM_Donor_F),]
for(i in 1:nrow(sampleInfo_DCM_Donor_F)) {
  cat(paste(sampleInfo_DCM_Donor_F$Disease[i], " "))
}

sampleInfo_DCM_F <- sampleInfo_DCM_Donor_F[sampleInfo_DCM_Donor_F$Disease == "DCM",]
for(i in 1:nrow(sampleInfo_DCM_F)) {
  cat(paste("\"",sampleInfo_DCM_F$SampleID[i], "\"\n",sep=""))
}
sampleInfo_Donor_F <- sampleInfo_DCM_Donor_F[sampleInfo_DCM_Donor_F$Disease == "Donor",]
for(i in 1:nrow(sampleInfo_Donor_F)) {
  cat(paste("\"",sampleInfo_Donor_F$SampleID[i], "\"\n",sep=""))
}



sampleInfo_DCM_Donor_M <- sampleInfo[sampleInfo$SampleID %in% colnames(gxData_DCM_Donor_M),]
for(i in 1:nrow(sampleInfo_DCM_Donor_M)) {
  cat(paste(sampleInfo_DCM_Donor_M$Disease[i], " "))
}

sampleInfo_DCM_M <- sampleInfo_DCM_Donor_M[sampleInfo_DCM_Donor_M$Disease == "DCM",]
for(i in 1:nrow(sampleInfo_DCM_M)) {
  cat(paste("\"",sampleInfo_DCM_M$SampleID[i], "\"\n",sep=""))
}
sampleInfo_Donor_M <- sampleInfo_DCM_Donor_M[sampleInfo_DCM_Donor_M$Disease == "Donor",]
for(i in 1:nrow(sampleInfo_Donor_M)) {
  cat(paste("\"",sampleInfo_Donor_M$SampleID[i], "\"\n",sep=""))
}



#Create lists for WebGestalt
webGestalt_F <- data.frame(Ensembl_GeneID = top_table$Ensembl_GeneID, 
                           DCMvsControl_female_logFC = top_table$DCMvsControl_female_logFC) %>%
  arrange(desc(DCMvsControl_female_logFC))

write.table(webGestalt_F, 
            file="webGestalt_F.txt",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

