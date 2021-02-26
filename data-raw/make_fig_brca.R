box::use(
  dplyr[...],
  ggplot2[...],
  patchwork[...],
  bigsimr[...],
  ./R/utils[mom_nbinom]
)

Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
bs <- bigsimr_setup(pkg_check = FALSE)
dist <- distributions_setup()

sim_nbinom <- readRDS("results/brca_1000_sim.rds")
load("data/example_brca.rda")

brca1000 <- example_brca %>%
  select(all_of(1:1000)) %>%
  mutate(across(everything(), as.double))
R_S <- bs$cor(as.matrix(brca1000), bs$Spearman)
nbinom_fit <- apply(brca1000, 2, mom_nbinom)
sizes <- nbinom_fit["size",]
probs <- nbinom_fit["prob",]

R_S_hat <- bs$cor(as.matrix(sim_nbinom), bs$Spearman)
sim_nbinom_fit <- apply(sim_nbinom, 2, mom_nbinom)
corr_compare <- tibble(
  `Target Correlation` = R_S[lower.tri(R_S, diag = FALSE)],
  `Estimated Correlation` = R_S_hat[lower.tri(R_S_hat, diag = FALSE)]
)
target_compare <- tibble(
  `Target Size` = sizes,
  `Estimated Size` = sim_nbinom_fit[1,]
)
logprob_compare <- tibble(
  `Target Log Probability` = log(probs),
  `Estimated Log Probability` = log(sim_nbinom_fit[2,])
)

p1 <- ggplot(corr_compare, aes(`Target Correlation`, `Estimated Correlation`)) +
  geom_point(alpha = 0.05) +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  theme_bw() + coord_equal()

p2 <- ggplot(target_compare, aes(`Target Size`, `Estimated Size`)) +
  geom_point(alpha = 0.85) +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  theme_bw() + coord_equal()

p3 <- ggplot(logprob_compare,
             aes(`Target Log Probability`, `Estimated Log Probability`)) +
  geom_point(alpha = 0.85) +
  geom_abline(slope = 1, col = "red", lty = "dashed") +
  theme_bw() + coord_equal()

p <- p1 + p2 + p3

ggsave(filename = "fig/ch050-figBRCA.png",
       plot = p,
       device = "png",
       width = 8,
       height = 4,
       units = "in")
