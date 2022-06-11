# Function to implement Shannon transfer entropy.
# Used internally by transfer_entropy; same arguments.
#
te_shannon <- function(x,
                       lx,
                       y,
                       ly,
                       shuffles,
                       type,
                       quantiles,
                       bins,
                       limits,
                       nboot,
                       burn,
                       quiet,
                       jidt = FALSE,
                       jidtCompat = FALSE,
                       jidtDebug = FALSE) {

  # Code time series
  x <- code_sample(x, type, quantiles, bins, limits)
  y <- code_sample(y, type, quantiles, bins, limits)

  # Lead = y
  if (!quiet) cat("  [calculate] X->Y transfer entropy\n")
  texy <- calc_te_shannon(y, lx = ly, x, ly = lx, jidt = jidt, jidtCompat = jidtCompat, jidtDebug = jidtDebug)
  consty <- shuffle_shannon(
    x = y,
    lx = ly,
    y = x,
    ly = lx,
    shuffles = shuffles,
    jidt = jidt,
    jidtCompat = jidtCompat,
    jidtDebug = jidtDebug
  )
  stexy <- texy - consty

  # Lead = x
  if (!quiet) cat("  [calculate] Y->X transfer entropy\n")
  teyx <- calc_te_shannon(x, lx = lx, y, ly = ly, jidt = jidt, jidtCompat = jidtCompat, jidtDebug = jidtDebug)
  constx <- shuffle_shannon(
    x = x,
    lx = lx,
    y = y,
    ly = ly,
    shuffles = shuffles,
    jidt = jidt,
    jidtCompat = jidtCompat,
    jidtDebug = jidtDebug
  )
  steyx <- teyx - constx

  # Bootstrap
  if (nboot > 0) {
    if (!quiet) {
      cat(sprintf(
        "  [bootstrap] %s time%s\n",
        nboot, mult_s(nboot)
      ))
    }

    boot <- future.apply::future_sapply(seq_len(nboot), function(i) {
      bootstrap_shannon(
        x = x,
        lx = lx,
        y = y,
        ly = ly,
        burn = burn
      )
    }, future.seed = TRUE)
  } else {
    boot <- NA
  }

  return(list(
    teyx = teyx,
    texy = texy,
    steyx = steyx,
    stexy = stexy,
    boot = boot
  ))
}
