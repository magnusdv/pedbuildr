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


`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

`%e%` = function(x, y) {
  if(x == "") y else x
}

as_int = function(m) {
  structure(as.integer(m), dim = dim(m))
}

indent = function(depth) {
  strrep(" ", 2 * (depth - 1))
}

# Equivalent to t.default(combn(x, 2)), but ~5 times faster.
.comb2 = function(x) {
  L = length(x)
  vec = L > 1
  n = if(vec) L else x

  if (n < 2)
    return(matrix(nrow = 0, ncol = 2))
  v1 = rep.int(seq_len(n - 1), (n - 1):1)
  v2 = NULL
  for (i in 2:n) v2 = c(v2, i:n)
  res = cbind(v1, v2, deparse.level = 0)

  if(vec)
    res[] = x[res]
  res
}

# Faster version of sort.int, especially for vectors of size 1 and 2
.mysortInt = function(v) {
  L = length(v)
  if(L == 1)
    return(v)
  if(L == 2) {
    if(v[1] > v[2])
      return(v[c(2, 1)])
    else return(v)
  }
  sort.int(v, method = "shell")
}

.mysetdiff = function(x, y)
  unique.default(x[match(x, y, 0L) == 0L])

