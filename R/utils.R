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

