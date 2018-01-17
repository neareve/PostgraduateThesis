indep_test <- function(u1, u2) {
  tau <- cor(u1, u2,method = "kendall")
  N <- length(u1)
  f <- sqrt((9 * N * (N - 1))/(2 * (2 * N + 5))) * abs(tau)
  return( p.value = 2 * (1 - pnorm(f)))
} 