
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

