stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

stopifnot2 = function(...) {
  exprs = list(...)

  for (i in seq_along(exprs)) {
    expri = .subset2(exprs, i)
    if (length(expri) != 1L || is.na(expri) || !expri) {
      full_call = match.call()
      call = deparse(full_call[[i + 1]])
      stop(sQuote(call), " is not TRUE", call. = FALSE, domain = NA)
    }
  }
}

as_int = function(m) {
  structure(as.integer(m), dim = dim(m))
}

indent = function(depth) {
  strrep(" ", 2 * (depth - 1))
}

.mysetdiff = function(x, y)
  unique.default(x[match(x, y, 0L) == 0L])

tryBreakLoops = function(x) {
  if(is.pedList(x))
    y = lapply(x, function(comp)
      breakLoops(comp, verbose = F, errorIfFail = F))
  else
    y = breakLoops(x, verbose = F, errorIfFail = F)

  # In unsuccessful: Return NULL
  if(has_unbroken_loops(y))
    return(NULL)

  y
}

# Utility: Compute total loglikelihood
#' @importFrom pedprobr likelihood
loglikTotal = function(x) {
  if(is.pedList(x)) {
    nm = vapply(x, nMarkers, 0L)
    if(!all(nm == nm[1]))
      stop2("Components have unequal numbers of attached markers: ", nm)
    nMark = nm[1]
  }
  else
    nMark = nMarkers(x)

  if(nMark == 0)
    return(0)

  x = tryBreakLoops(x)
  if(is.null(x))
    return(NA_real_)

  loglik_per_marker = vapply(seq_len(nMark),
    function(i) likelihood(x, i, logbase = exp(1)), FUN.VALUE = 1)

  sum(loglik_per_marker)
}
