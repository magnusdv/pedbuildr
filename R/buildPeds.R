#' Build all possible pedigrees
#'
#' @param ids A vector of ID labels.
#' @param sex A vector of the same length as `ids`, with entries 1 (male) or 2
#'   (female).
#' @param age A numeric vector of the same length as `ids`. If `age[i] <
#'   age[j]`, individual `i` will not be an ancestor of individual `j`. The
#'   numbers themselves are irrelevant, only the partial ordering is used.
#' @param knownPO A list of pairs of ID labels: The known parent-offspring
#'   pairs.
#' @param allKnown A logical. If TRUE, no other pairs than `knownPO` will be
#'   assigned as parent-offspring. If FALSE (default), all pairwise combinations
#'   of ids - except those in `notPO` - will be treated as potential
#'   parent-offspring pairs.
#' @param notPO A list of pairs of ID labels: Pairs known not to be
#'   parent-offspring.
#' @param connected A logical. If TRUE (default), only connected pedigrees are
#'   returned.
#' @param maxLinearInb A nonnegative integer, or `Inf` (default). If this
#'   is a finite number, it disallows mating between pedigree members X and Y if
#'   X is a linear descendant of Y separated by more than the given number. For
#'   example, setting `maxLinearInb = 0` forbids mating between
#'   parent-child, grandparent-grandchild, a.s.o. If `maxLinearInb = 1`
#'   then parent-child matings are allowed, but not grandparent-grandchild or
#'   higher.
#' @param genderSym A logical. If TRUE (default), pedigrees which are equal
#'   except for the gender distribution of the *added* parents, are regarded as
#'   equivalent, and only one of each equivalence class is returned. Example:
#'   paternal vs. maternal half sibs.
#'
#' @param verbose A logical.
#'
#' @return A list of pedigrees.
#'
#' @examples
#' p = buildPeds(1:3, sex = c(1,2,1), knownPO = list(c(1,3), c(2,3)))
#' stopifnot(length(p) == 25)
#' # plotPeds(p)
#'
#' # By default, one pedigree was removed above.
#' # Keep all gender combinations (in this case paternal/maternal half sibs):
#' p2 = buildPeds(1:3, sex = c(1,2,1), knownPO = list(c(1,3), c(2,3)),
#'               genderSym = FALSE)
#' stopifnot(length(p2) == 26)
#' # plotPeds(p2)
#'
#' # Remove pedigrees with linear inbreeding (e.g. parent-child)
#' p3 = buildPeds(1:3, sex = c(1,2,1), knownPO = list(c(1,3), c(2,3)),
#'               maxLinearInb = 0)
#' stopifnot(length(p3) == 7)
#' # plotPeds(p3)
#'
#' @export
buildPeds = function(ids, sex, age = NULL,
                     knownPO = NULL, allKnown = FALSE, notPO = NULL,
                     connected = TRUE, maxLinearInb = Inf,
                     genderSym = TRUE, verbose = FALSE) {
  N = length(ids)
  if(!setequal(ids, 1:N))
    stop2("`ids` must be an integer vector of the form `1:N`")
  if(length(sex) != N)
    stop2("`ids` and `sex` must have the same length\nids: ", ids, "\nsex: ", sex)

  # Sort by ids
  if(!identical(ids, 1:N)) {
    ids = ids[order(ids)]
    sex = sex[order(ids)]
  }

  if(allKnown && is.null(knownPO))
    stop2("`knownPO` cannot be NULL when `allKnown = TRUE`")

  # Convert age vector into matrix with all ordered pairs (works with NULL)
  if(!is.null(age) && (length(age) != N || !is.numeric(age)))
    stop2("`ids` and `age` must have the same length\nids: ", ids, "\nage: ", age)
  ageMat = which(outer(age, age, `<`), arr.ind = TRUE)

  if(verbose) {
    .knownPO = toString(vapply(knownPO, paste, collapse = "-", FUN.VALUE="")) %e% "-"
    .notPO = toString(vapply(notPO, paste, collapse = "-", FUN.VALUE="")) %e% "-"
    .age = toString(paste(ageMat[, 1], ageMat[, 2], sep = "<")) %e% "-"

    print(glue::glue("
      Pedigree parameters:
        ID labels: {toString(ids)}
        Sex: {toString(sex)}
        Age info: {.age}
        Known PO: {.knownPO}
        Known non-PO: {.notPO}
        Connected only: {connected}
        Symmetry filter: {genderSym}
        Linear inbreeding: {maxLinearInb}"
    ))
  }

  if(verbose) cat("\nBuilding pedigree list:\n")

  # List possible sets of parent-offspring
  POsets = listPOsets(knownPO = knownPO, allKnown = allKnown, notPO = notPO, N)

  # Convert POlists to undirected adjacency matrices
  UA = lapply(POsets, function(po) po2adj(po, N))
  if(verbose) cat("  Undirected adjacency matrices:", length(UA), "\n")

  # For each undirAdj, build list of directed adjacency matrices
  DA = lapply(UA, function(ua) directedAdjs(ua, sex, ageMat = ageMat, verbose = FALSE))
  DA = unlist(DA, recursive = FALSE)
  if(verbose) cat("  Directed adjacency matrices:", length(DA), "\n")

  # Extend each matrix by adding parents in all possible ways
  DA_EXT = lapply(DA, function(da)
    addMissingParents(da, maxLinearInb = maxLinearInb, genderSym = genderSym))
  DA_EXT = unlist(DA_EXT, recursive = FALSE)
  if(verbose) cat("  After adding parents:", length(DA_EXT), "\n")

  # If `connected = TRUE`, remove disconnected matrices
  if(connected) {
    DA_EXT = DA_EXT[sapply(DA_EXT, isConnected)]
    if(verbose) cat("  Connected solutions:", length(DA_EXT), "\n")
  }

  # Convert to list of pedigrees
  peds = lapply(DA_EXT, function(a) adj2ped(a, origSize = N))

  invisible(peds)
}

# Auxiliary function for listing all possible PO sets
#' @importFrom utils combn
listPOsets = function(knownPO = NULL, allKnown = FALSE, notPO = NULL, n) {
  if(allKnown) {
    if(!is.null(notPO))
      stop2("When `allKnown`is TRUE, `notPO` must be NULL")
    return(if(is.null(knownPO)) list() else list(knownPO))
  }

  knownPO = lapply(knownPO, .mysortInt)
  notPO = lapply(notPO, .mysortInt)

  # Check for illegal overlaps
  knownPOchar = sapply(knownPO, paste, collapse="-")
  notPOchar = sapply(notPO, paste, collapse = "-")
  if(length(err <- intersect(knownPOchar, notPOchar)))
    stop2("`knownPO` and `notPO` must be disjoint: ", err)

  # Potential extra parent-offspring: All pairs except "knownPO" and "notPO"
  allPO = combn(n, 2, simplify = FALSE)
  allPOchar = sapply(allPO, paste, collapse = "-")
  potentialPO = allPO[!allPOchar %in% c(knownPOchar, notPOchar)]

  if(length(potentialPO) == 0)
    return(list(knownPO))

  # Loop over all subsets of potentialPO and add to knownPO
  subs = expand.grid(rep(list(c(FALSE, TRUE)), length(potentialPO)))
  totalPO = apply(subs, 1, function(s) c(knownPO, potentialPO[s]))

  totalPO
}


# Convert list of PO pairs to undirected (symmetric) adjacency matrix
po2adj = function(po, n) {
  m = matrix(FALSE, ncol = n, nrow = n)
  for(e in po)
    m[e[1], e[2]] = m[e[2], e[1]] = TRUE
  m
}


directedAdjs = function(undirAdj, sex, ageMat, verbose = TRUE) {

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

  SEX = attr(adj, 'sex')
  rownum = seq_len(dim(remaining)[1])

  # parents of id are f and m
  adj[c(father, mother), id] = TRUE

  # id is not parent of f,m
  remaining[id, c(father, mother)] =  FALSE

  # remaining adjacent are kids of id
  kids = rownum[remaining[id, ]]
  if(verbose)
    cat("kids =", toString(kids), "\n")

  if(length(kids)) {
    adj[id, kids] = TRUE

    # no-one else with same sex as id is parent of kids
    remaining[SEX == SEX[id], kids] = FALSE
    remaining[id, kids] = TRUE # restore these

    # no kids are parent of c(id,f,m) # utvid til alle ancs!!!
    remaining[kids, c(id, father, mother)] = FALSE
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
    # Solution of not - go to next (in the outside loop)
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
