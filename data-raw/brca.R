## code to prepare `brca` dataset goes here

options(timeout = 240) # my internet is sloooowww apparently

brca0 <- TCGA2STAT::getTCGA(disease = "BRCA", data.type = "RNASeq")
brca1 <- as.data.frame(t(brca0$dat))
gene_median <- apply(brca1, 2, median)
gene_median_order <- order(gene_median, decreasing = TRUE)

example_brca     <- brca1[, head(gene_median_order, 10000)]
example_patients <- substring(row.names(example_brca), first = 1, last = 12)
example_genes    <- colnames(example_brca)

usethis::use_data(example_brca, overwrite = TRUE)
usethis::use_data(example_patients, overwrite = TRUE)
usethis::use_data(example_genes, overwrite = TRUE)
