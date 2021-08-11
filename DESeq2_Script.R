counts <- readRDS("D:/pan_cancer_analysis/TCGA_BRCA/data/TCGA-BRCA_1222_samples_raw_counts.RDS")
metadata <- readRDS("D:/pan_cancer_analysis/TCGA_BRCA/data/TCGA-BRCA_1222_samples_metadata.RDS")

train_samples <- read.csv("D:/pan_cancer_analysis/TCGA_BRCA/data/TCGA-BRCA_Primary_Train_Tumors_vs_Solid_Normal_Tissues.csv",
                          header = TRUE,
                          row.names = 1,
                          sep = ";",
                          check.names = FALSE)

library(dplyr)

metadata <- metadata %>% select(sample_type)

metadata <- subset(metadata, rownames(metadata) %in% rownames(train_samples))

counts <- counts %>% select(intersect(colnames(counts), rownames(metadata)))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ sample_type)

dds <- DESeq(dds)
saveRDS(dds, "dds.RDS")

