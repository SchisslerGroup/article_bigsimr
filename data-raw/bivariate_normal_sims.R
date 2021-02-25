## code to prepare `bivariate_normal_sims` dataset goes here
box::use(
  bigsimr[...],
  usethis[use_data],
  ./R/method_of_moments[mom_norm]
)

Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
bs <- bigsimr_setup(pkg_check = FALSE)
dist <- distributions_setup()

mom_norm <- function(x) {
  m <- mean(x)
  s <- sd(x)
  list(mean = m, sd = s)
}

mu <- 0
sigma <- 1
margins <- c(dist$Normal(mu, sigma),
             dist$Normal(mu, sigma))
types <- c("Pearson", "Spearman", "Kendall")
n <- c(1e3, 1e4, 1e5)
eps <- 1e-2
grid_steps <- 100

sim_pars <- expand.grid(
  type = types,
  n = n,
  stringsAsFactors = FALSE
)

res <- data.frame()

for (i in 1:nrow(sim_pars)) {
  type <- sim_pars$type[i]
  n <- sim_pars$n[i]

  tmp_bounds <- bs$cor_bounds(margins, bs[[type]])

  cor_lo <- tmp_bounds$lower[1, 2] + eps
  cor_hi <- tmp_bounds$upper[1, 2] - eps
  cor_seq <- seq(cor_lo, cor_hi, length.out = grid_steps)

  for (rho in cor_seq) {
    Rho <- matrix(rho, 2, 2)
    diag(Rho) <- 1.0
    if (type == "Pearson") {
      Rho <- bs$pearson_match(Rho, margins)
    } else {
      Rho <- bs$cor_convert(Rho, bs[[type]], bs$Pearson)
    }

    time_data <- system.time(x <- bs$rvec(n, Rho, margins))

    id <- paste0(
      "d", 2,
      "-N", n,
      "-c", Sys.getenv("JULIA_NUM_THREADS"),
      "-r", rho,
      "-Cor", type,
      "-dev", "cpu",
      "-lib", "bigsimr"
    )

    ## Estimate statistics
    Rho_hat <- bs$cor(x, bs[[type]])
    rho_hat <- Rho_hat[1, 2]

    norm_args_hat <- mom_norm(x[, 1])
    mu_hat <- norm_args_hat$mean
    sigma_hat <- norm_args_hat$sd

    ## Save the results
    res <- rbind(res, data.frame(
      method = "bigsimr",
      device = "CPU",
      type = type,
      cores = Sys.getenv("JULIA_NUM_THREADS"),
      margins = "norm",
      d = 2,
      N = n,
      rho = rho,
      rho_hat = rho_hat,
      mean = mu,
      sd = sigma,
      mean_hat = mu_hat,
      sd_hat = sigma_hat,
      sim_time = unname(time_data["elapsed"])
    ))
  }
}
res$type <- factor(res$type, levels = c("Pearson", "Spearman", "Kendall"))

bivariate_normal_sims <- res

usethis::use_data(bivariate_normal_sims, overwrite = TRUE)
