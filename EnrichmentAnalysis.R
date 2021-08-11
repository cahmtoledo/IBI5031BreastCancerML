library(biomaRt)
library(clusterProfiler)

ensembl <- useEnsembl(biomart='ensembl')
ensembl <- useDataset('hsapiens_gene_ensembl', mart=ensembl)
probe_ids <- read.table('~/carmensmind/IBI5031/IBI5031 - Projeto/mutinfogenes.txt')

annotation <- getBM(attributes = c('ensembl_gene_id',
                                   'external_gene_name', 'gene_biotype'), 
                    filters = 'ensembl_gene_id', #choose filter based in the plataform used
                    values = probe_ids,
                    mart = ensembl,
                    useCache = T)

enrichNGC(annotation$external_gene_name)
enri