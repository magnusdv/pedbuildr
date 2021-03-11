# Convert adjacency matrix to ped format
adj2ped = function(adj, labs = NULL) {
  sex = as.integer(attr(adj, 'sex'))
  n = length(sex)
  nseq = seq_len(n)

  # Fix labels
  if(is.null(labs))
    labs = as.character(nseq)
  else {
    labs = as.character(labs)
    origN = length(labs)
    if(n > origN) {
      exSeq = seq_len(n - origN)
      labs[origN + exSeq] = paste0("e", exSeq)
    }
  }

  # Find fidx and midx
  fid = mid = integer(n)
  parents = nseq[.rowSums(adj, n, n) > 0]

  for(i in parents) {
    kids = nseq[adj[i, ]]
    if(sex[i] == 1)
      fid[kids] = i
    else
      mid[kids] = i
  }

  # If known to be connected, go straight to newPed()
  if(isTRUE(attr(adj, "connected")))
    return(newPed(labs, fid, mid, sex, ""))

  p = ped(id = nseq, fid = fid, mid = mid, sex = sex, reorder = FALSE, validate = FALSE)
  relabelFast(p, labs)
}


# Stripped-down version of pedtools::relabel()
relabelFast = function(x, newlabs) {

  if(is.pedList(x)) {
    y = lapply(x, function(comp) {
      comp$ID = newlabs[as.integer(comp$ID)]
      comp
    })
    class(y) = "pedList"
    return(y)
  }

  x$ID = newlabs

  x
}

# Not used
relabelAddedParents = function(x, origN) {
  if(is.pedList(x)) {
    y = lapply(x, relabelAddedParents, origN)
    class(y) = "pedList"
    return(y)
  }

  n = length(x$ID)
  if(n > origN) {
    exSeq = seq_len(n - origN)
    x$ID[origN + exSeq] = paste0("e", exSeq)
  }

  x
}


# Convert pedigree to adjacency matrix
ped2adj = function(ped) {
  if(is.pedList(ped)) {
    return(lapply(ped, ped2adj))
  }

  adj = matrix(0L, ncol = pedsize(ped), nrow = pedsize(ped),
               dimnames = list(labels(ped), labels(ped)))

  for(nf in nonfounders(ped))
    adj[parents(ped, nf), nf] = 1L

  adjMatrix(adj, sex = getSex(ped))
}
