library(bigsimr)
size <- 4
prob <- 3e-04
type <- "spearman"
n <- 1e3

margins <- alist(
        qnbinom(size = size, prob = prob),
        qnbinom(size = size, prob = prob)
)

tmp_bounds <- cor_bounds(margins, type = type)
