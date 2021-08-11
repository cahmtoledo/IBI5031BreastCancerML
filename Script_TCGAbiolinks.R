library(TCGAbiolinks)
library(DT)
library(SummarizedExperiment)
library(tidyverse)

query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts",
                  experimental.strategy = "RNA-Seq",
                  legacy = FALSE)
GDCdownload(query, method = "api", files.per.chunk = 10)
data <- GDCprepare(query)

count_data <- as.data.frame(assay(data))
metadata <- as.data.frame(colData(data))

write.csv(count_data, "IBI5031/TCGA-BRCA_1222_samples_raw_counts.csv")
write.csv(metadata, "IBI5031/TCGA-BRCA_1222_samples_metadata.csv")

saveRDS(count_data, "IBI5031/TCGA-BRCA_1222_samples_raw_counts.RDS")
saveRDS(metadata, "IBI5031/TCGA-BRCA_1222_samples_metadata.RDS")

write.csv(count_data, "IBI5031/TCGA-BRCA_1222_samples_raw_counts.csv")
write.csv(metadata, "IBI5031/TCGA-BRCA_1222_samples_metadata.csv")

primary_tumor_metadata <- subset(metadata, definition == "Primary solid Tumor")
primary_tumor_count_data <- count_data %>% select(intersect(rownames(primary_tumor_metadata), colnames(count_data)))
primary_tumor_subtype <- primary_tumor_metadata %>% select(patient, sample, barcode, paper_BRCA_Subtype_PAM50)
  
saveRDS(primary_tumor_count_data, "IBI5031/TCGA-BRCA_1102_primary_tumors_raw_counts.RDS")
saveRDS(primary_tumor_metadata, "IBI5031/TCGA-BRCA_1102_primary_tumors_metadata.RDS")
saveRDS(primary_tumor_subtype, "IBI5031/TCGA-BRCA_1102_primary_tumors_subtypes.RDS")

write.csv(primary_tumor_count_data, "IBI5031/TCGA-BRCA_1102_primary_tumors_raw_counts.csv")
write.csv(primary_tumor_metadata, "IBI5031/TCGA-BRCA_1102_primary_tumors_metadata.csv")
write.csv(primary_tumor_subtype, "IBI5031/TCGA-BRCA_1102_primary_tumors_subtypes.csv")