raw_data <- read.csv("D:/IBI5031_project/X_raw_counts.csv",
                     header = TRUE,
                     row.names = 1,
                     sep = ",")

only_train_raw_data <- read.csv("D:/IBI5031_project/X_train_raw_counts.csv",
                               header = TRUE,
                               row.names = 1,
                               sep = ",")
  
library(edgeR)

dge <- DGEList(raw_data)

dge_TMM <- calcNormFactors(dge, method = "TMM")
normalized_data <- cpm(dge_TMM)

#write.csv(normalized_data, "D:/IBI5031_project/Train+Test_TMM-normalized_counts.csv")
saveRDS(normalized_data, "D:/IBI5031_project/Train+Test_TMM-normalized_counts.RDS")

dge_train <- DGEList(only_train_raw_data)

dge_train_TMM <- calcNormFactors(dge_train, method = "TMM")
normalized_train_data <- cpm(dge_train_TMM)

#write.csv(normalized_train_data, "D:/IBI5031_project/Train_TMM-normalized_counts.csv")
saveRDS(normalized_train_data, "D:/IBI5031_project/Train_TMM-normalized_counts.RDS")

means <- rowMeans(normalized_train_data)
df_means <- as.data.frame(means)
rownames(df_means) <- rownames(normalized_train_data)

quantile(df_means$means, 0)    
quantile(df_means$means, 0.25) 
quantile(df_means$means, 0.5)  
quantile(df_means$means, 0.75) # 3.126453
quantile(df_means$means, 1)    

sub_means <- subset(df_means, means > 3.126453)

library(matrixStats)
library(ggplot2)

vars <- rowVars(normalized_train_data)     
df_vars <- as.data.frame(vars)
rownames(df_vars) <- rownames(normalized_train_data)

quantile(df_vars$vars, 0)    
quantile(df_vars$vars, 0.25) 
quantile(df_vars$vars, 0.5)  
quantile(df_vars$vars, 0.75) # 18.4733 
quantile(df_vars$vars, 1)    

sub_vars <- subset(df_vars, vars > 18.4733)

length(intersect(rownames(sub_means), rownames(sub_vars)))

higher_mean_and_var <- intersect(rownames(sub_means), rownames(sub_vars))

sub_normalized_train_data <- subset(normalized_train_data, rownames(normalized_train_data) %in% higher_mean_and_var)

write.csv(as.data.frame(sub_normalized_train_data), "D:/IBI5031_project/Train_TMM-normalized_count_data_13212_most_expressed_and_variable_genes.csv")

sub_normalized_data <- subset(normalized_data, rownames(normalized_data) %in% intersect(rownames(normalized_data), rownames(sub_normalized_train_data)))

write.csv(as.data.frame(sub_normalized_data), "D:/IBI5031_project/Train+Test_TMM-normalized_count_data_13212_most_expressed_and_variable_genes.csv")
