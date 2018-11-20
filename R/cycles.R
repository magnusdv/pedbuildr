# Determine if a directed 0-1 matrix has a cycle
hasCycle = function(m, maxlen = NA) {
  if(all(m == 0))
    return(FALSE)
  if(any(diag(m) == 1))
    return(TRUE)

  if(is.na(maxlen))
    maxlen = sum(m)

  mpow = m
  for(i in 1:maxlen) {
    mpow = mpow %*% m
    if(any(diag(mpow) == 1))
      return(TRUE)
  }

  FALSE
}

hasCycle2 = function(m) {
  if(all(m == 0))
    return(FALSE)
  if(any(diag(m) == 1))
    return(TRUE)

  for(i in 1:ncol(m)) {
    anci = pedmatrixAncestors(m, i)
    if(i %in% anci) {
      #message(i, " is own anc")
      return(TRUE)
    }
  }

  FALSE
}


pedmatrixAncestors = function(m, i) {
  anc = integer(0)
  newanc = i
  while(TRUE) {
    newanc = which(m[, newanc, drop=F] == 1, arr.ind = TRUE)[, 1]
    if(length(newanc) == 0 || all(newanc %in% anc))
      break
    anc = c(anc, newanc)

    if(i %in% newanc)
      break
  }
  anc
  #sort.default(unique.default(anc))
}

rmat = function(N = 5) {
  m = matrix(0, ncol=N, nrow=N)
  for(i in which(sample(c(T, F), N, replace = T)))
    m[sample.int(N, size=2), i] = 1
  diag(m) = 0
  m
}
