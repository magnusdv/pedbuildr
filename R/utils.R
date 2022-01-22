stop2 = function(...) {
  a = lapply(list(...), toString)
  a = append(a, list(call. = FALSE))
  do.call(stop, a)
}

isCount = function(x, minimum = 1, maximum = NA) {
    isTRUE(length(x) == 1 &&
             (is.integer(x) || (is.numeric(x) && x == as.integer(x))) &&
             x >= minimum && (is.na(maximum) || x <= maximum))
}

# Sample from Dirichlet distribution with mean p (where b controls variance)
#' @importFrom stats rgamma
rdirich = function(n, p, b = 1) {
  if(!is.numeric(b) || length(b) != 1 || b <= 0)
    stop("`b` must be a positive integer")

  alpha = p * b
  s = vapply(alpha, function(a) rgamma(n, a, 1), numeric(n))
  if (n == 1)
    dim(s) = c(1, length(alpha))

  s/rowSums(s)
}

powerset = function(x, force = NULL) {

  if(!is.null(force)) {
    y = .mysetdiff(x, force, makeUnique = FALSE)
    return(lapply(powerset(y), function(u) c(force, u)))
  }

  n = length(x)
  if(n == 0)
    return(list(x[0]))
  if(n == 1)
    return(list(x[0], x))
  if(n == 2)
    return(list(x[0], x[1], x[2], x))

  unlist(lapply(0:length(x), fast.combn, x = x), recursive = FALSE)
}

`%||%` = function(x, y) {
  if(is.null(x)) y else x
}

`%e%` = function(x, y) {
  if(x == "") y else x
}


indent = function(depth) {
  strrep(" ", 2 * (depth - 1))
}

ftime = function(st, digits = 3)
  format(Sys.time() - st, digits = digits)


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

.mysetdiff = function(x, y, makeUnique = TRUE) {
  if(is.null(y)) {
    if(makeUnique)
      unique.default(x)
    else
      x
  }
  else{
    if(makeUnique)
      unique.default(x[match(x, y, 0L) == 0L])
    else
      x[match(x, y, 0L) == 0L]
  }
}

# Fast intersection. NB: assumes no duplicates!
.myintersect = function (x, y) {
  y[match(x, y, 0L)]
}

# Stripped version of expand.grid
fast.grid = function(argslist, as.list = FALSE) {
  nargs = length(argslist)
  orep = nr = prod(lengths(argslist))
  if (nargs == 0L || nr == 0L)
    return(if(as.list) list() else matrix(ncol = 0, nrow = 0))

  rep.fac = 1L
  res = NULL
  for (x in argslist) {
    nx = length(x)
    orep = orep/nx
    res = c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)])  #this is res[, i]
    rep.fac = rep.fac * nx
  }

  dim(res) = c(nr, nargs)
  if (as.list)
    res = lapply(seq_len(nr), function(r) res[r, ])

  res
}


# Stripped version of utils::combn(x, m, FUN = NULL, simplify = FALSE)
# In particular: Never converts x to 1:x
fast.combn = function(x, m) {
  n = length(x)
  if (n < m)
    stop2("n < m")
  m <- as.integer(m)
  e <- 0
  h <- m
  a <- seq_len(m)

  out <- vector("list", choose(n, m))
  out[[1L]] <- x[a]

  if(m == 0)
    return(out)

  i <- 2L
  nmmp1 <- n - m + 1L
  while (a[1L] != nmmp1) {
    if (e < n - h) {
      h <- 1L
      e <- a[m]
      j <- 1L
    }
    else {
      e <- a[m - h]
      h <- h + 1L
      j <- 1L:h
    }
    a[m - h + j] <- e + j
    out[[i]] <- x[a]
    i <- i + 1L
  }
  out
}


