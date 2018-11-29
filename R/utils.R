stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

as_int = function(m) {
  structure(as.integer(m), dim = dim(m))
}

indent = function(depth) {
  strrep(" ", 2 * (depth - 1))
}

.mysetdiff = function(x, y)
  unique.default(x[match(x, y, 0L) == 0L])


# Utility: Compute total loglikelihood
loglikTotal = function(x) {
  if(is.pedList(x)) {
    nm = vapply(x, nMarkers, 0L)
    stopifnot(all(nm == nm[1]))
    nMark = nm[1]
  }
  else {
    nMark = nMarkers(x)
  }
  if(nMark == 0)
    return(0)

  loglik_per_marker = vapply(seq_len(nMark),
    function(i) pedprobr::likelihood(x, i, logbase = exp(1)), FUN.VALUE = 1)

  sum(loglik_per_marker)
}
