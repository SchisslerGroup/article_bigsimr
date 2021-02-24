## code to prepare `benchmark_dependences` dataset goes here

library(tidyverse)
library(bigsimr)
Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
bs <- bigsimr_setup(pkg_check = FALSE)
dist <- distributions_setup()

# "preheat" the bigsimr functions
bs$cor(rnorm(10), rnorm(10), bs$Pearson)
bs$cor(rnorm(10), rnorm(10), bs$Spearman)
bs$cor(rnorm(10), rnorm(10), bs$Kendall)

bs$cor_fast(matrix(rnorm(100), 10, 10), bs$Pearson)
bs$cor_fast(matrix(rnorm(100), 10, 10), bs$Spearman)
bs$cor_fast(matrix(rnorm(100), 10, 10), bs$Kendall)

bs$cor_convert(matrix(rnorm(100), 10, 10), bs$Spearman, bs$Pearson)
bs$cor_convert(matrix(rnorm(100), 10, 10), bs$Kendall, bs$Pearson)

d <- dist$NegativeBinomial(20, 0.2)
bs$pearson_match(0.5, d, d)

d <- dist$NegativeBinomial(20, 0.002)
bs$pearson_bounds(d, d)
bs$pearson_match(0.5, d, d)

d <- dist$Gamma(2, 0.002)
bs$pearson_bounds(d, d)
bs$pearson_match(0.5, d, d)

s <- bs$cor_randPSD(100)
p <- bs$cor_convert(s, bs$Spearman, bs$Pearson)
r <- bs$cor_nearPD(p)
JuliaCall::julia_call("Bigsimr.iscorrelation", r)


# Prepare the data
load("data/example_brca.rda")
brca_mat <- do.call(cbind, lapply(example_brca, as.double))

mom_gamma <- function(x) {
  m <- mean(x)
  s <- sd(x)
  c(shape = m^2 / s^2, rate = m / s^2)
}

make_margins <- function(par1, par2) {
  margins <- lapply(1:length(par1), function(i) {
    dist$Gamma(par1[i], par2[i])
  })
  do.call(c, margins)
}

params <- apply(brca_mat, 2, mom_gamma)
margins <- make_margins(params[1, 1:3], params[2, 1:3])
bs$pearson_bounds(margins)

# steps
# 1. estimate pearson correlation: cor(data, cor_type)
# 2. calculate adjusted pearson correlation: match(corr, margins)
# 3. check admissibility:
#    3a) isvalidcorrelation(adjusted_corr)
#    3b) cor_nearPD(adjusted_corr)
# 4. simulate data: rvec(1000, adjusted_corr_pd, margins)

bench <- function(type, d, n=1000) {
  # step 0
  m <- brca_mat[,1:d]
  margins <- make_margins(params[1, 1:d], params[2, 1:d])

  # step 1
  if (type == "Kendall") {
    corr_time <- system.time(corr <- bs$cor_fast(m, bs[[type]]))
  } else {
    corr_time <- system.time(corr <- bs$cor(m, bs[[type]]))
  }

  # step 2
  if (type == "Pearson") {
    adj_time <- system.time(adj_corr <- bs$pearson_match(corr, margins))
  } else {
    adj_time <- system.time(adj_corr <- bs$cor_convert(corr, bs[[type]], bs$Pearson))
  }

  # step 3
  admiss_time <- system.time({
    adj_corr_pd <- adj_corr
    needed_near_pd <- FALSE
    if (!JuliaCall::julia_call("Bigsimr.iscorrelation", adj_corr)) {
      needed_near_pd <- TRUE
      adj_corr_pd <- bs$cor_nearPD(adj_corr)
    }
  })

  # step 4
  sim_time <- system.time(x <- bs$rvec(n, adj_corr_pd, margins))

  list(
    corr_type = type,
    dim = d,
    n_sim = n,
    corr_time = unname(corr_time["elapsed"]),
    adj_time = unname(adj_time["elapsed"]),
    admiss_time = unname(admiss_time["elapsed"]),
    sim_time = unname(sim_time["elapsed"]),
    needed_near_pd = needed_near_pd
  )
}

dims <- c(100, 250, 500, 1000)
pearson_bench <- lapply(dims, function(d) bench("Pearson", d)) %>%
  do.call(what = bind_rows)

dims <- c(100, 250, 500, 1000, 2500, 5000, 10000)
spearman_bench <- lapply(dims, function(d) bench("Spearman", d)) %>%
  do.call(what = bind_rows)

dims <- c(100, 250, 500, 1000, 2500, 5000, 10000)
kendall_bench <- lapply(dims, function(d) bench("Kendall", d)) %>%
  do.call(what = bind_rows)

benchmark_dependences <- bind_rows(pearson_bench, spearman_bench, kendall_bench) %>%
  mutate(corr_type = factor(corr_type, levels = c("Pearson", "Spearman", "Kendall")),
         total_time = corr_time + adj_time + admiss_time + sim_time)

benchmark_dependences

usethis::use_data(benchmark_dependences, overwrite = TRUE)
