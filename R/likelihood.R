
#' @importFrom pedprobr likelihood
loglikTotal = function(x) {
  sum(likelihood(x, logbase = exp(1)))
}


# Main wrapper used in reconstruct()
loglikWrapper = function(ped, amatList, loci, founderInb = 0) {
  # Attach marker data
  x = setMarkersFAST(ped, amatList, loci)

  # Set founder inbreeding
  if(founderInb > 0)
    x = setFounderInbreeding(x, value = founderInb)

  # Return total loglik
  tryCatch(loglikTotal(x), error = function(e) NA_real_)
}

# Outdated code -------------------------------------------------------------------------------

loglikTotalOLD = function(x) {

  # Catch potential fails (return NA)
  x = tryBreakLoops(x)
  if(is.null(x))
    return(NA_real_)

  logliks = likelihood(x, logbase = exp(1))
  sum(logliks)
}

tryBreakLoops = function(x) {
  if(is.pedList(x))
    y = lapply(x, function(comp)
      pedprobr:::.breakLoops(comp, verbose = FALSE, errorIfFail = FALSE))
  else
    y = pedprobr:::.breakLoops(x, verbose = FALSE, errorIfFail = FALSE)

  # In unsuccessful: Return NULL
  if(hasUnbrokenLoops(y))
    return(NULL)

  y
}
