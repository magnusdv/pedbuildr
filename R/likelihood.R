
# Utility: Compute total loglikelihood
#' @importFrom pedprobr likelihood
loglikTotal = function(x) {

  # Catch potential fails (return NA)
  x = tryBreakLoops(x)
  if(is.null(x))
    return(NA_real_)

  #logliks = vapply(seq_len(nMark),
  #   function(i) likelihood(x, i, logbase = exp(1)), FUN.VALUE = 1)
  logliks = likelihood(x, logbase = exp(1))

  sum(logliks)
}


tryBreakLoops = function(x) {
  if(is.pedList(x))
    y = lapply(x, function(comp)
      breakLoops(comp, verbose = FALSE, errorIfFail = FALSE))
  else
    y = breakLoops(x, verbose = FALSE, errorIfFail = FALSE)

  # In unsuccessful: Return NULL
  if(hasUnbrokenLoops(y))
    return(NULL)

  y
}
