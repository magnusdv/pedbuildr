# Determine if a directed 0-1 matrix has a cycle
hasCycle = function(m, maxlen = NA) {
  if(all(m == 0))
    return(FALSE)
  if(any(diag(m) == 1))
    return(TRUE)

  n = dim(m)[1]
  if(is.na(maxlen)) {
    rowRank = sum(.rowSums(m, n, n) > 0)
    colRank = sum(.colSums(m, n, n) > 0)
    maxlen = min(rowRank, colRank)
  }

  if(maxlen < 2)
    return(FALSE)

  mpow = m
  for(i in 2:maxlen) {
    mpow = mpow %*% m
    if(any(diag(mpow) > 0))
      return(TRUE)
  }

  FALSE
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

# Check connectedness of a directed adjacency matrix
isConnected = function(adj) {
  n = dim(adj)[1]

  # If size 1: connected
  if(n == 1)
    return(TRUE)

  # Adj of undirected graph
  a = adj | t.default(adj)

  # Count in+out edges for each node
  edgs = .rowSums(a, n, n)
  if(any(edgs == 0))
    return(FALSE)

  nseq = seq_len(n)
  getNeigh = function(ids) {
    nids = length(ids)
    if(nids == 1)
      return(nseq[a[, ids]])
    nseq[.rowSums(a[, ids], n, nids) > 0]
  }

  # Start at node with few edges
  comp = integer()
  new = which.min(edgs)
  while(length(new) > 0) {
    comp = c(comp, new)
    neigh = getNeigh(new) # which(a[new, , drop = F], arr.ind = T)[, 2]
    new = .mysetdiff(neigh, comp, makeUnique = FALSE)
  }

  # Connected iff all are included
  length(comp) == n
}



