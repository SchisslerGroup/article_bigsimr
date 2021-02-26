if (!file.exists("data/example_brca.rda")) {
  e <- new.env()
  source("data-raw/brca.R", local = e)
  rm(e)
}
if (!file.exists("data/bivariate_normal_sims.rda")) {
  e <- new.env()
  source("data-raw/bivariate_normal_sims.R", local = e)
  rm(e)
}
if (!file.exists("data/bivariate_gamma_sims.rda")) {
  e <- new.env()
  source("data-raw/bivariate_gamma_sims.R", local = e)
  rm(e)
}
if (!file.exists("data/bivariate_nbinom_sims.rda")) {
  e <- new.env()
  source("data-raw/bivariate_nbinom_sims.R", local = e)
  rm(e)
}
if (!file.exists("data/benchmark_dependences.rda")) {
  e <- new.env()
  source("data-raw/benchmark_dependences.R", local = e)
  rm(e)
}
if (!file.exists("results/brca_1000_sim.rds")) {
  e <- new.env()
  source("data-raw/sim_nbinom.R", local = e)
  rm(e)
}
if (!file.exists("fig/ch050-figBRCA.png")) {
  e <- new.env()
  source("data-raw/make_fig_brca.R", local = e)
  rm(e)
}
if (!file.exists("results/brca_1000_frobenius.rds")) {
  e <- new.env()
  source("data-raw/sim_nbinom_corr.R", local = e)
  rm(e)
}
