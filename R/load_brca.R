#' @export
load_brca <- function(d = NULL) {
  brca0 <- TCGA2STAT::getTCGA(disease = "BRCA", data.type = "RNASeq")
  brca1 <- as.data.frame(t(brca0$dat))
  gene_median <- apply(brca1, 2, median)
  gene_median_order <- order(gene_median, decreasing = TRUE)

  if (is.null(d))
    d <- length(gene_median)

  brca  <- brca1[, head(gene_median_order, d)]
  patients <- substring(row.names(brca), first = 1, last = 12)

  list(brca = brca, patients = patients)
}
