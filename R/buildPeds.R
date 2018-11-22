#' Build all possible pedigrees
#'
#' @param ids A vector of ID labels
#' @param sex A vector of the same length as `ids`, with entries 1 (male) or 2
#'   (female)
#' @param knownPO A list of pairs of ID labels: The known parent-offspring pairs
#' @param allKnown A logical. If TRUE, no other pairs than `knownPO` will be
#'   assigned as parent-offspring. If FALSE (default), all pairwise combinations
#'   of ids - except those in `notPO` - will be treated as potential
#'   parent-offspring pairs.
#' @param notPO A list of pairs of ID labels: Pairs known not to be parent-offspring.
#' @param maxLinearInbreeding A nonnegative integer, or `Inf` (default). If this
#'   is a finite number, it disallows mating between pedigree members X and Y if
#'   X is a linear descendant of Y separated by more than the given number. For
#'   example, setting `maxLinearInbreeding = 0` forbids mating between
#'   parent-child, grandparent-grandchild, a.s.o. If `maxLinearInbreeding = 1`
#'   then parent-child matings are allowed, but not grandparent-grandchild or
#'   higher.
#' @param verbose A logical
#'
#' @return A list
#'
#' @examples
#' p = buildPeds(1:3, sex = c(1,2,1), knownPO = list(c(1,3), c(2,3)))
#' stopifnot(length(p) == 26)
#' # plotPeds(p)
#'
#' # Remove pedigrees with linear inbreeding (e.g. parent-child)
#' p2 = buildPeds(1:3, sex = c(1,2,1), knownPO = list(c(1,3), c(2,3)),
#'               maxLinearInbreeding = 0)
#' stopifnot(length(p2) == 8)
#' # plotPeds
#'
#' @export
buildPeds = function(ids, sex, knownPO = NULL, allKnown = F, notPO = NULL,
                     maxLinearInbreeding = Inf, verbose = F) {
  N = length(ids)
  stopifnot(length(sex) == N, setequal(ids, 1:N))

  # List possible sets of parent-offspring
  POsets = listPOsets(knownPO = knownPO, allKnown = allKnown, notPO = notPO, N)

  # Convert POlists to undirected adjacency matrices
  UA = lapply(POsets, po2adj, n = N)
  if(verbose) cat("Undirected adjacency matrices:", length(UA), "\n")

  # For each undirAdj, build list of directed adjacency matrices
  DA = lapply(UA, function(ua) directedAdjs(ua, sex, verbose = F))
  DA = unlist(DA, recursive = F)
  if(verbose) cat("Directed adjacency matrices:", length(DA), "\n")

  # Extend each matrix by adding parents in all possible ways
  DA_EXT = lapply(DA, addMissingParents, maxLinearInbreeding = maxLinearInbreeding)
  DA_EXT = unlist(DA_EXT, recursive = F)
  if(verbose) cat("After adding parents:", length(DA_EXT), "\n")

  # Convert to list of pedigrees
  peds = lapply(DA_EXT, adj2ped, origSize = N)
  if(verbose) cat("Pedigrees:", length(peds), "\n")

  invisible(peds)
}

# Auxiliary function for listing all possible PO sets
#' @importFrom utils combn
listPOsets = function(knownPO = NULL, allKnown = F, notPO = NULL, n) {
  if(allKnown) {
    if(!is.null(notPO))
      stop2("When `allKnown`is TRUE, `notPO` must be NULL")
    return(if(is.null(knownPO)) list() else list(knownPO))
  }

  knownPO = lapply(knownPO, sort.int)
  notPO = lapply(notPO, sort.int)

  # Check for illegal overlaps
  knownPOchar = sapply(knownPO, paste, collapse="-")
  notPOchar = sapply(notPO, paste, collapse="-")
  if(length(err <- intersect(knownPOchar, notPOchar)))
    stop2("`knownPO` and `notPO` must be disjoint: ", err)

  # Potential extra parent-offspring: All pairs except "knownPO" and "notPO"
  allPO = combn(n, 2, simplify = F)
  allPOchar = sapply(allPO, paste, collapse="-")
  potentialPO = allPO[!allPOchar %in% c(knownPOchar, notPOchar)]

  if(length(potentialPO) == 0)
    return(list(knownPO))

  # Loop over all subsets of potentialPO and add to knownPO
  subs = expand.grid(rep(list(c(F,T)), length(potentialPO)))
  totalPO = apply(subs, 1, function(s) c(knownPO, potentialPO[s]))

  totalPO
}


# Convert list of PO pairs to undirected (symmetric) adjacency matrix
po2adj = function(po, n) {
  m = matrix(F, ncol = n, nrow = n)
  for(e in po)
    m[e[1], e[2]] = m[e[2], e[1]] = T
  m
}

directedAdjs = function(undirAdj, sex, verbose = T) {

  # Environment for holding the identified pedigrees
  storage = new.env()
  storage$dirAdjs = list()
  storage$target = sum(undirAdj)/2 # the target number of edges

  # Empty (directed) adjacency matrix
  adj = adjMatrix(sex = sex)

  # Start with an individual of max degree
  id = which.max(colSums(undirAdj))

  # Potential fathers and mothers
  PF = which(undirAdj[, id] & sex == 1)
  PM = which(undirAdj[, id] & sex == 2)

  # Add edges recursively
  for(f in c(0, PF)) for(m in c(0, PM)) {
    addEdge(adj, id, f, m, remaining = undirAdj, storage = storage, verbose = verbose)
  }

  storage$dirAdjs
}

addEdge = function(adj, id, father, mother, remaining, storage, verbose = T, depth = 1) {
  if(verbose)
    cat(sprintf("%sDepth %d. ID = %d, father = %d, mother = %d, ",
                indent(depth), depth, id, father, mother))

  SEX = attr(adj, 'sex')

  # parents of id are f and m
  adj[c(father, mother), id] = T

  # id is not parent of f,m
  remaining[id, c(father, mother)] =  F

  # remaining adjacent are kids of id
  kids = which(remaining[id, ])
  if(verbose)
    cat("kids =", toString(kids), "\n")

  if(length(kids)) {
    adj[id, kids] = T

    # no-one else with same sex as id is parent of kids
    remaining[SEX == SEX[id], kids] = F
    remaining[id, kids] = T # restore these

    # no kids are parent of c(id,f,m) # utvid til alle ancs!!!
    remaining[kids, c(id, father, mother)] = F
  }

  # Remove from "remaining" edges involving id
  remaining[id, ] = remaining[, id] = F

  # No remaining edges?
  if(!any(remaining)) {

    if(sum(adj) == storage$target && !hasCycle(adj)) { # Correct total number of edges, and no cycles
      if(verbose) cat(indent(depth+1), "Adj matrix found\n")
      storage$dirAdjs = c(storage$dirAdjs, list(adj))
    }
    else {
      if(verbose) cat(indent(depth+1), "No solution here\n")
    }
    # Solution of not - go to next (in the outside loop)
    return()
  }

  # Remaining edges? If so - continue
  id = which.max(colSums(remaining))
  PF = which(remaining[, id] & SEX == 1)
  PM = which(remaining[, id] & SEX == 2)

  for(f in c(0, PF)) for(m in c(0, PM)) {
    addEdge(adj, id, f, m, remaining, storage, verbose = verbose, depth = depth + 1)
  }
}
