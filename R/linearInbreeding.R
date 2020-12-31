# Check for linear inbreeding (of distance at least minDist)
# Original version, very slow unless `descList` provided
linearInb = function(a, descList = NULL, minDist = 1) {
  n = dim(a)[1]
  if(n < 3 || minDist >= n)
    return(FALSE)

  nseq = seq_len(n)

  if(is.null(descList)) {
    descList = lapply(nseq, function(id)
      dagDescendants(a, i = id, minDist = minDist))
  }

  hasMultipleKids = nseq[.rowSums(a, n, n) > 1]
  for(i in hasMultipleKids) {
    kids = nseq[a[i, ] == 1]
    for(k in kids)
      if(any(kids %in% descList[[k]]))
        return(TRUE)
  }

  FALSE
}

# Faster version of the above
linearInb2 = function(a, minDist = 1) {
  n = dim(a)[1]
  if(n < 3 || minDist >= n)
    return(FALSE)

  nseq = seq_len(n)

  getParents = function(ids) {
    nids = length(ids)
    if(nids == 0)
      return(integer(0))
    if(nids == 1)
      return(if(ids == 0) integer(0) else nseq[a[, ids]])
    nseq[.rowSums(a[, ids], n, nids) > 0]
  }

  getAncs = function(id) {
    ANCS = vector("list", n - 1)
    newanc = getParents(id)
    j = 0
    while(j < n - 1 && length(newanc)) {
      j = j + 1
      ANCS[[j]] = newanc
      newanc = getParents(newanc)
    }
    if(j < minDist)
      return(integer(0))
    unlist(ANCS[minDist:j])
  }

  A = vector("list", n)
  for(id in nseq) {
    par = getParents(id)
    if(length(par) < 2)
      next
    p1 = par[1]; p2 = par[2]

    anc1 = A[[p1]]
    if(is.null(anc1))
      anc1 = A[[p1]] = getAncs(p1)
    if(p2 %in% anc1)
      return(TRUE)

    anc2 = A[[p2]]
    if(is.null(anc2))
      anc2 = A[[p2]] = getAncs(p2)
    if(p1 %in% anc2)
      return(TRUE)
  }

  FALSE
}
