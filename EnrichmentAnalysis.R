library(biomaRt)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)

ensembl <- useEnsembl(biomart='ensembl')
ensembl <- useDataset('hsapiens_gene_ensembl', mart=ensembl)
probe_mutinfo <- read.table('mutinfogenes.txt')

# Retriving genes anotations
annotation_mutinfo <- getBM(attributes = c('ensembl_gene_id',
                                   'external_gene_name','entrezgene_id', 'gene_biotype'), 
                    filters = 'ensembl_gene_id', #choose filter based in the plataform used
                    values = probe_mutinfo,
                    mart = ensembl,
                    useCache = T)

# Testing for enhanced pathways in Network of Cancer Gene (NCG)
# NCGenhanced_mutinfo <- enrichNCG(annotation_mutinfo$entrezgene_id)
# GOenhanced_mutinfo <- enrichGO(annotation_mutinfo$external_gene_name, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
KEGGenhanced_mutinfo <- enrichKEGG(annotation_mutinfo$entrezgene_id)

probe_deg <- read.table('deggenes.txt')

# Retriving genes anotations
annotation_deg <- getBM(attributes = c('ensembl_gene_id',
                                   'external_gene_name','entrezgene_id', 'gene_biotype'), 
                    filters = 'ensembl_gene_id', #choose filter based in the plataform used
                    values = probe_deg,
                    mart = ensembl,
                    useCache = T)

# NCGenhanced <- enrichNCG(annotation_deg$entrezgene_id)
# GOenhanced <- enrichGO(annotation_deg$external_gene_name, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")
KEGGenhanced_deg <- enrichKEGG(annotation_deg$entrezgene_id)

barplot_mutinfo <- barplot(KEGGenhanced_mutinfo) + xlab("Count") + xlim(0,10) +
        ggtitle("Mutual Information")
barplot_deg <- barplot(KEGGenhanced_deg) + xlab("count") + xlim(0,10) + 
        ggtitle("DEGs")
 
cowplot::plot_grid(barplot_mutinfo,barplot_deg, ncol = 1)

