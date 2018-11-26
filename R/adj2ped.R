# Convert adjacency matrix to ped format
adj2ped = function(adj, origSize = ncol(adj)) {
  sex = attr(adj, 'sex')
  N = length(sex)

  fid = mid = integer(N)
  parents = which(rowSums(adj) > 0)

  # TODO: make more efficient
  for(i in parents) {
    kids = which(adj[i, ] == 1)
    switch(sex[i], {fid[kids] = i}, {mid[kids] = i})
  }

  p = ped(id = 1:N, fid = fid, mid = mid, sex = sex, reorder = F, validate = F)

  # Relabel aded indivs to p1, p2, ...
  relabelAddedParents(p, origSize)
}

relabelAddedParents = function(x, origSize) {
  if(is.pedList(x)) {
    y = lapply(x, relabelAddedParents, origSize)
    class(y) = "pedList"
    return(y)
  }

  added = as.numeric(x$ID) > origSize
  if(any(added)) {
    x$ID[added] = paste0("p", seq_len(sum(added)))
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
