set.seed(2020-02-25)

box::use(
  dplyr[...],
  bigsimr[...],
  ./R/utils[mom_nbinom]
)

Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
bs <- bigsimr_setup(pkg_check = FALSE)
dist <- distributions_setup()

make_nbinom_margins <- function(sizes, probs) {
  margins <- lapply(1:length(sizes), function(i) {
    dist$NegativeBinomial(sizes[i], probs[i])
  })
  do.call(c, margins)
}

load("data/example_brca.rda")

brca1000 <- example_brca %>%
  select(all_of(1:1000)) %>%
  mutate(across(everything(), as.double))

R_S <- bs$cor(as.matrix(brca1000), bs$Spearman)
R_X <- bs$cor_convert(R_S, bs$Spearman, bs$Pearson)
R_X <- bs$cor_nearPD(R_X)

nbinom_fit <- apply(brca1000, 2, mom_nbinom)
sizes <- nbinom_fit["size",]
probs <- nbinom_fit["prob",]
nb_margins <- make_nbinom_margins(sizes, probs)
sim_nbinom <- bs$rvec(10000, R_X, nb_margins)
colnames(sim_nbinom) <- colnames(brca1000)

saveRDS(sim_nbinom, file = "results/brca_1000_sim.rds")
