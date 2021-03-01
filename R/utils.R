#' @export
mom_norm <- function(x) {
  box::use(stats[sd])
  c(mean = mean(x), sd = sd(x))
}

#' @export
mom_nbinom <- function(x) {
  box::use(stats[sd])
  m <- mean(x)
  s <- sd(x)
  c(size = m^2 / (s^2 - m), prob = m / s^2)
}

#' @export
mom_gamma <- function(x) {
  box::use(stats[sd])
  m <- mean(x)
  s <- sd(x)
  c(shape = m^2 / s^2, rate = m / s^2)
}
