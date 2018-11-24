# Determine if a directed 0-1 matrix has a cycle
hasCycle = function(m, maxlen = NA) {
  if(all(m == 0))
    return(FALSE)
  if(any(diag(m) == 1))
    return(TRUE)

  if(is.na(maxlen)) {
    rowRank = sum(rowSums(m) > 0)
    colRank = sum(colSums(m) > 0)
    maxlen = min(rowRank, colRank)
  }

  mpow = m
  for(i in 1:maxlen) {
    mpow = mpow %*% m
    if(any(diag(mpow) == 1))
      return(TRUE)
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

# All descendants of a node (or set of nodes) in a DAG
dagDescendants = function(adj, i, minDist = 1, maxDist = Inf) {
  desc = if(minDist == 0) i else integer(0)

  current = i # current set of individuals
  dist = 0 # generation number

  while(dist < maxDist && length(current) > 0) {
    dist = dist + 1

    # Next generation: All offspring of current
    current = which(adj[current, , drop = F], arr.ind = T)[, 2]

    # Add to storage if appropriate
    if(dist >= minDist)
      desc = c(desc, current)
  }

  # Return
  as.integer(desc)
}


isConnected = function(adj) {
  # If size 1: connected
  if(nrow(adj) == 1)
    return(TRUE)

  # Adj of underlying (undirected) graph
  a = adj | t.default(adj)

  # Count in+out edges for each node
  edgs = rowSums(a)
  if(any(edgs == 0))
    return(FALSE)

  # Start at node with few edges
  comp = integer()
  new = which.min(edgs)
  while(length(new) > 0) {
    comp = c(comp, new)
    neigh = which(a[new, , drop = F], arr.ind = T)[, 2]
    new = setdiff(neigh, comp)
  }

  # Connected iff all are included
  length(comp) == nrow(adj)
}

rmat = function(N = 5) {
  m = matrix(0, ncol=N, nrow=N)
  for(i in which(sample(c(T, F), N, replace = T)))
    m[sample.int(N, size=2), i] = 1
  diag(m) = 0
  m
}


