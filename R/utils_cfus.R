find_poisson_cutoff <- function(counts, dil, N, V=1){
  mask <- counts<N
  r_p <- sum(counts[mask])/(sum(dil[mask])*V)
  r_p_std <- sqrt(r_p*r_p/sum(counts[mask]))
  tibble(r = r_p, stderr = r_p_std, type = "poisson_cutoff")
}

find_MLE <- function(counts, dil, N, V=1){
  f <- function(x) {
    sum(counts*dil/(N*(1-exp(-x*dil*V/N))))-sum(dil)
  }
  if (any(N<counts)) {
    tibble(r = Inf, stderr = Inf, type = "ML_estimator")
  } else {
    root <- uniroot(f, interval = c(1, 1e9))
    r_mle <- root$root
    p0 <- exp(-r_mle*dil*V/N)
    invVar <- sum(dil*dil*V*V*counts*p0/(N*N*(1-p0)**2))
    estVar <- 1/invVar
    r_mle_std <- sqrt(estVar)
    tibble(r = r_mle, stderr = r_mle_std, type = "ML_estimator")
  }
}

find_estimators <- function(counts, dil, N, NMLE) {
  bind_rows(find_poisson_cutoff(counts, dil, N),
            find_MLE(counts, dil, NMLE))
}