
# Build pedigrees between a set of individuals.
# Missing parents are added and allowed to interact in all possible ways.
# Note: knownPO, notPO, noChildren must be internal numerics
buildPedsParents = function(labs, sex, ageMat = NULL, knownPO = NULL, allKnown = FALSE,
                     notPO = NULL, noChildren = NULL, connected = TRUE,
                     maxLinearInb = Inf, sexSymmetry = TRUE, verbose = FALSE) {

  st = Sys.time()
  N = length(labs)

  if(verbose) cat("\nBuilding pedigree list:\n")

  # List possible sets of parent-offspring
  POsets = listPOsets(N, knownPO = knownPO, allKnown = allKnown, notPO = notPO)

  # Convert POlists to undirected adjacency matrices
  UA = lapply(POsets, function(po) po2adj(po, N))
  if(verbose) cat("  Undirected adjacency matrices:", length(UA), "\n")

  # For each undirAdj, build list of directed adjacency matrices
  DA = lapply(UA, function(ua)
    directedAdjs(ua, sex, ageMat = ageMat, noChildren = noChildren, verbose = FALSE))
  DA = unlist(DA, recursive = FALSE)
  if(verbose) cat("  Directed adjacency matrices:", length(DA), "\n")

  # Extend each matrix by adding parents in all possible ways
  DA_EXT = lapply(DA, function(da)
    addMissingParents(da, maxLinearInb = maxLinearInb, sexSymmetry = sexSymmetry))
  DA_EXT = unlist(DA_EXT, recursive = FALSE)
  if(verbose) cat("  After adding parents:", length(DA_EXT), "\n")

  # If `connected = TRUE`, remove disconnected matrices
  if(connected) {
    DA_EXT[!sapply(DA_EXT, isConnected)] = NULL
    if(verbose) cat("  Connected solutions:", length(DA_EXT), "\n")
  }

  # Convert to list of pedigrees
  peds = lapply(DA_EXT, function(a) adj2ped(a, labs))

  invisible(peds)
}

# Auxiliary function for listing all possible PO sets
listPOsets = function(N, knownPO = NULL, allKnown = FALSE, notPO = NULL) {
  if(allKnown)
    return(if(is.null(knownPO)) list() else list(knownPO))

  knownPO = lapply(knownPO, .mysortInt)
  notPO = lapply(notPO, .mysortInt)

  # Check for illegal overlaps
  knownPOchar = sapply(knownPO, paste, collapse="-")
  notPOchar = sapply(notPO, paste, collapse = "-")
  if(length(err <- intersect(knownPOchar, notPOchar)))
    stop2("`knownPO` and `notPO` must be disjoint: ", err)

  # Potential extra parent-offspring: All pairs except "knownPO" and "notPO"
  allPO = fast.combn(1:N, 2)
  allPOchar = sapply(allPO, paste, collapse = "-")
  potentialPO = allPO[!allPOchar %in% c(knownPOchar, notPOchar)]

  if(length(potentialPO) == 0)
    return(list(knownPO))

  if(length(potentialPO) > 22)
    stop2("Problem is too complex; number of potential PO relations = ", length(potentialPO))

  # Loop over all subsets of potentialPO and add to knownPO
  subs = fast.grid(rep(list(c(FALSE, TRUE)), length(potentialPO)))
  lapply(1:nrow(subs), function(i) c(knownPO, potentialPO[subs[i, ]]))
}


# Convert list of PO pairs to undirected (symmetric) adjacency matrix
po2adj = function(po, N) {
  unl = unlist(po, use.names = FALSE)
  NR = max(unl, N)
  m = logical(NR*NR)
  even = seq_along(po) * 2
  unlEven = unl[even]
  unlOdd = unl[even - 1]
  m[NR * (unlEven - 1) + unlOdd] = TRUE
  m[NR * (unlOdd - 1) + unlEven] = TRUE
  dim(m) = c(NR, NR)
  m
}


directedAdjs = function(undirAdj, sex, ageMat, noChildren, verbose = TRUE) {

  # Environment for holding the identified pedigrees
  storage = new.env()
  storage$dirAdjs = list()
  storage$target = sum(undirAdj)/2 # the target number of edges
  storage$ageData = lapply(unique.default(ageMat[, 1]), function(i)
    list(id = i, older = as.numeric(ageMat[ageMat[,1] == i, 2])))

  # Empty (directed) adjacency matrix
  adj = adjMatrix(sex = sex)

  # Apply age info
  undirAdj[ageMat] = FALSE

  # Apply age info
  undirAdj[noChildren, ] = FALSE

  # Start with an individual of max degree
  id = which.max(colSums(undirAdj))

  # Potential fathers and mothers
  rownum = seq_along(sex)
  PF = rownum[undirAdj[, id] & sex == 1]
  PM = rownum[undirAdj[, id] & sex == 2]

  # Add edges recursively
  for(f in c(0, PF)) for(m in c(0, PM)) {
    addEdge(adj, id, f, m, remaining = undirAdj, storage = storage, verbose = verbose)
  }

  storage$dirAdjs
}

addEdge = function(adj, id, father, mother, remaining, storage, verbose = TRUE, depth = 1) {
  if(verbose)
    cat(sprintf("%sDepth %d. ID = %d, father = %d, mother = %d, ",
                indent(depth), depth, id, father, mother))

  parents = c(father, mother)

  SEX = attr(adj, 'sex')
  rownum = seq_len(dim(remaining)[1])

  # parents of id are f and m
  adj[parents, id] = TRUE

  # id is not parent of f,m
  remaining[id, parents] =  FALSE

  # remaining adjacent are kids of id
  kids = rownum[remaining[id, ]]
  if(verbose)
    cat("kids =", toString(kids), "\n")

  if(length(kids)) {
    adj[id, kids] = TRUE

    # no-one else with same sex as id is parent of kids
    remaining[SEX == SEX[id], kids] = FALSE
    remaining[id, kids] = TRUE # restore these

    # no kids are parent of c(id,f,m) # TODO: utvid til alle ancs!!!
    remaining[kids, c(id, parents)] = FALSE
  }

  # Remove from "remaining" edges involving id
  remaining[id, ] = remaining[, id] = FALSE

  # No remaining edges?
  if(!any(remaining)) {

    if(sum(adj) == storage$target &&           # correct number of edges
       !hasCycle(adj) &&                       # no cycles
       !ageViolation(adj, storage$ageData)     # no older descendants
      ) {
      if(verbose) cat(indent(depth+1), "Adj matrix found\n")
      storage$dirAdjs = c(storage$dirAdjs, list(adj))
    }
    else {
      if(verbose) cat(indent(depth+1), "No solution here\n")
    }
    # Solution or not - go to next (in the outside loop)
    return()
  }

  # Remaining edges? If so - continue
  id = which.max(colSums(remaining))
  PF = rownum[remaining[, id] & SEX == 1]
  PM = rownum[remaining[, id] & SEX == 2]

  for(f in c(0, PF)) for(m in c(0, PM)) {
    addEdge(adj, id, f, m, remaining, storage, verbose = verbose, depth = depth + 1)
  }
}

# Check if any individual has descendants who are supposed to be older
ageViolation = function(adj, ageData) {
  # ageData: List of lists: list(id = 1, older = 2:4)
  for(dat in ageData) {
    if(any(dat$older %in% dagDescendants(adj, dat$id, minDist = 2)))
      return(TRUE)
  }
  return(FALSE)
}


