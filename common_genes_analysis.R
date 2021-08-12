common_genes <- 0:100

genes_in_common <- function(n){
        choose(100,n)*choose(13112,100-n)/choose(13212, 100)
        
}

probabilities <- sapply(common_genes, genes_in_common)

choose(13212, 100)

sum(probabilities[15:101]) # p_value of 14

sum(probabilities[10:101]) # p_value of 9
