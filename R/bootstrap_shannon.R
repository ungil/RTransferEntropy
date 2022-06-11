# Function for bootstrapping Shannon transfer entropy under H0 of independence
# between time series x and y. Used internally by transfer_entropy; same
# arguments.
#
bootstrap_shannon <- function(x,
                              lx,
                              y,
                              ly,
                              burn = 50,
                              jidt = FALSE,
                              jidtCompat = FALSE,
                              jidtDebug = FALSE) {
  bootx <- markov_boot_step(x, lx, burn)
  booty <- markov_boot_step(y, ly, burn)

  # Lead = x
  dteyx <- calc_te_shannon(
    x = bootx,
    lx = lx,
    y = y,
    ly = ly,
    jidt = jidt,
    jidtCompat = jidtCompat,
    jidtDebug = jidtDebug
  )

  # Lead = y
  dtexy <- calc_te_shannon(
    x = booty,
    lx = ly,
    y = x,
    ly = lx,
    jidt = jidt,
    jidtCompat = jidtCompat,
    jidtDebug = jidtDebug    
  )

  teboot <- c(dtexy, dteyx)
  names(teboot) <- c("dtexy", "dteyx")

  return(teboot)
}
