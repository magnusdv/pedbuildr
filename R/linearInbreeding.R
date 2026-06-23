# Check for linear inbreeding (of distance at least minDist).
# Used in `addMissingParents()`.
# NB: Very slow unless `descList` provided.
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

# Faster version of the above, used in "extra" branch.
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


linealParents = function(par, anc) {
  length(par) == 2L &&
    (par[1] %in% anc[[par[2]]] || par[2] %in% anc[[par[1]]])
}

hasLinealMating = function(a, anc, ids) {
  # Check selected individuals for lineal mating between their parents.
  # Used after adding the new individual. It is enough to inspect `ids`,
  # since only these individuals received new ancestors.

  for(id in ids) {
    par = which(a[, id])

    if(length(par) == 2L &&
       (par[1] %in% anc[[par[2]]] || par[2] %in% anc[[par[1]]]))
      return(TRUE)
  }

  FALSE
}


# Pairs of positions in `ids` where one individual descends from the other.
# A partition is illegal if both positions are assigned to the same new parent.
linealPairs = function(ids, descList) {

  if(length(ids) < 2L)
    return(matrix(integer(0), ncol = 2L))

  out = lapply(seq_along(ids), function(i) {
    j = match(descList[[ids[i]]], ids, nomatch = 0L)
    j = j[j > 0L]
    if(length(j))
      cbind(i, j)
  })

  out = out[lengths(out) > 0L]

  if(length(out))
    do.call(rbind, out)
  else
    matrix(integer(0), ncol = 2L)
}

linealPartition = function(group, pairs) {
  if(!nrow(pairs))
    return(FALSE)

  group = as.integer(group)
  any(group[pairs[, 1]] == group[pairs[, 2]])
}
