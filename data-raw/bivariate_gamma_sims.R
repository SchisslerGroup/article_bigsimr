## code to prepare `bivariate_gamma_sims` dataset goes here

box::use(
  bigsimr[...],
  ./R/utils[mom_gamma]
)

Sys.setenv(JULIA_NUM_THREADS = parallel::detectCores())
bs <- bigsimr_setup(pkg_check = FALSE)
dist <- distributions_setup()

shape <- 10
rate <- 1
margins <- c(dist$Gamma(shape, rate),
             dist$Gamma(shape, rate))

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

    gamma_args_hat <- mom_gamma(x[,1])
    shape_hat <- gamma_args_hat["shape"]
    rate_hat  <- gamma_args_hat["rate"]

    ## Save the results
    res <- rbind(res, data.frame(
      method = "bigsimr",
      device = "CPU",
      type = type,
      cores = Sys.getenv("JULIA_NUM_THREADS"),
      margins = "gamma",
      d = 2,
      N = n,
      rho = rho,
      rho_hat = rho_hat,
      shape = shape,
      rate = rate,
      shape_hat = shape_hat,
      rate_hat = rate_hat,
      sim_time = unname(time_data["elapsed"])
    ))
  }
}

res$type <- factor(res$type, levels = c("Pearson", "Spearman", "Kendall"))
bivariate_gamma_sims <- res
save(bivariate_gamma_sims, file = "data/bivariate_gamma_sims.rda")
