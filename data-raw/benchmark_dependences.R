## code to prepare `benchmark_dependences` dataset goes here

box::use(
  bigsimr[...],
  dplyr[...],
  JuliaCall[julia_call],
  beepr[beep],
  ./R/utils[mom_gamma]
)

Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
bs <- bigsimr_setup(pkg_check = FALSE)
dist <- distributions_setup()


make_margins <- function(par1, par2) {
  margins <- lapply(1:length(par1), function(i) {
    dist$Gamma(par1[i], par2[i])
  })
  do.call(c, margins)
}


# "preheat" the bigsimr functions
bs$cor(rnorm(10), rnorm(10), bs$Pearson)
bs$cor(rnorm(10), rnorm(10), bs$Spearman)
bs$cor(rnorm(10), rnorm(10), bs$Kendall)

bs$cor(matrix(rnorm(100), 10, 10), bs$Pearson)
bs$cor(matrix(rnorm(100), 10, 10), bs$Spearman)
bs$cor(matrix(rnorm(100), 10, 10), bs$Kendall)

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
bs$cor_nearPD(p, 1e-6)
bs$cor_nearPD(p, 1e-6, tol=1e-6)
bs$cor_fastPD(p)

k <- bs$cor_randPSD(10000)
p <- bs$cor_convert(s, bs$Kendall, bs$Pearson)
julia_call("Bigsimr.iscorrelation", p)
system.time(r1 <- bs$cor_nearPD(p, 1e-6, tol=1e-2))
system.time(r2 <- bs$cor_nearPD(p, 1e-6, tol=1e-6))
system.time(r3 <- bs$cor_fastPD(p))
julia_call("Bigsimr.iscorrelation", r1)
julia_call("Bigsimr.iscorrelation", r2)
julia_call("Bigsimr.iscorrelation", r3)

d <- dist$Gamma(3, 0.01)
bs$rvec(10, bs$cor_randPD(2), c(d, d))

d <- 20
tmp_cor <- bs$cor_randPD(d)
shapes <- runif(d, 1, 10)
rates  <- rexp(d, 1/5)
tmp_margins <- make_margins(shapes, rates)
bs$pearson_match(tmp_cor, tmp_margins)

rm(p, r, r1, r2, r3, s, d, k, tmp_cor, tmp_margins, shapes, rates)
gc()


# Prepare the data
if (!file.exists("data/synthetic_data.rda")) {
  d <- 10000
  shapes <- runif(d, 1, 10)
  rates  <- rexp(d, 1/5)
  margins <- make_margins(shapes, rates)
  system.time(corr <- bs$cor_randPD(d))
  margin_bounds <- bs$pearson_bounds(margins)
  lb <- margin_bounds$lower
  ub <- margin_bounds$upper

  corr <- pmin(corr, ub)
  corr <- pmax(corr, lb)
  diag(corr) <- 1.0
  if (!julia_call("Bigsimr.iscorrelation", corr))
    corr <- bs$cor_nearPD(corr)

  data <- bs$rvec(1000, corr, margins)

  save(data,   file = "data/synthetic_data.rda")
  save(shapes, file = "data/synthetic_shapes.rda")
  save(rates,  file = "data/synthetic_rates.rda")
  rm(d, data, corr, shapes, rates, margins, margin_bounds, lb, ub)
  gc()
}

load("data/synthetic_data.rda")
gamma_params_hat <- apply(data, 2, mom_gamma)
shapes <- gamma_params_hat["shape",]
rates  <- gamma_params_hat["rate",]

# steps
# 1. estimate pearson correlation: cor(data, cor_type)
# 2. calculate adjusted pearson correlation: match(corr, margins)
# 3. check admissibility:
#    3a) isvalidcorrelation(adjusted_corr)
#    3b) cor_nearPD(adjusted_corr)
# 4. simulate data: rvec(1000, adjusted_corr_pd, margins)

bench <- function(type, d, n=1000) {
  # step 0
  margins <- make_margins(shapes[1:d], rates[1:d])
  m <- data[, 1:d]

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
    if (!julia_call("Bigsimr.iscorrelation", adj_corr)) {
      needed_near_pd <- TRUE
      adj_corr_pd <- bs$cor_nearPD(adj_corr, 1e-6, tol=0.005)
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


dims <- c(100, 250, 500, 1000, 2500, 5000, 10000)

pearson_bench <- local({
  x <- lapply(dims, function(d) bench("Pearson", d)) %>%
    do.call(what = bind_rows)
  beep(8)
  x
})
save(pearson_bench, file = "data/pearson_bench.rda")


spearman_bench <- local({
  x <- lapply(dims, function(d) bench("Spearman", d)) %>%
    do.call(what = bind_rows)
  beep(8)
  x
})
save(spearman_bench, file = "data/spearman_bench.rda")

kendall_bench <- local({
  x <- lapply(dims, function(d) bench("Kendall", d)) %>%
    do.call(what = bind_rows)
  beep(8)
  x
})
save(kendall_bench, file = "data/kendall_bench.rda")


benchmark_dependences <- bind_rows(pearson_bench, spearman_bench, kendall_bench) %>%
  mutate(corr_type = factor(corr_type, levels = c("Pearson", "Spearman", "Kendall")),
         total_time = corr_time + adj_time + admiss_time + sim_time)

save(benchmark_dependences, file = "data/benchmark_dependences.rda")
