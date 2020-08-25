#' Add missing parents to a pedigree adjacency matrix
#'
#' @param a An adjMatrix object
#' @param maxLinearInbreeding A nonnegative integer, or `Inf` (default). If this
#'   is a finite number, it disallows mating between pedigree members X and Y if
#'   X is a linear descendant of Y separated by more than the given number. For
#'   example, setting `maxLinearInbreeding = 0` forbids mating between
#'   parent-child, grandparent-grandchild, a.s.o. If `maxLinearInbreeding = 1`
#'   then parent-child matings are allowed, but not grandparent-grandchild or
#'   higher.
#' @param genderSym A logical. If TRUE, pedigrees which are equal except for the
#'   gender distribution of the *added* parents, are regarded as equivalent, and
#'   only one of each equivalence class is returned. Example: paternal vs
#'   maternal half sibs.
#' @return A list of adjMatrix objects where all columns sum to either 0 or 2.
#'
#' @examples
#' a = adjMatrix(c(1,0,0,0), sex=c(1,1))
#' addMissingParents(a)
#'
#' b = adjMatrix(rbind(rep(1,4), 0,0,0), sex=c(1,1,1,1))
#' addMissingParents(b)
#'
#' @export
addMissingParents = function(a, maxLinearInbreeding = Inf, genderSym = FALSE) {
  sex = attr(a, "sex")
  n = ncol(a)
  idvec = seq_len(n)

  missingFa = idvec[colSums(a[sex == 1, , drop = F]) == 0]
  missingMo = idvec[colSums(a[sex == 2, , drop = F]) == 0]

  nMissFa = length(missingFa)
  nMissMo = length(missingMo)

  if(nMissFa > 7) stop2("More than 7 extra fathers needed: Too many possible combinations")
  if(nMissMo > 7) stop2("More than 7 extra mothers needed: Too many possible combinations")

  # Founders (needed in a later step)
  fou = idvec[colSums(a) == 0]

  # Setup for gender symmetry restriction
  if(genderSym) {
    observedInvariants = character()
    pows = 2^(0:(n - 1))
  }

  # List of descendants
  checkInb = maxLinearInbreeding < Inf
  if(checkInb)
    descList = lapply(1:n, function(id)
      dagDescendants(a, id, minDist = maxLinearInbreeding + 1))

  # All set partitions for the fathers and the mothers
  if(nMissFa > 0)
    pFa = partitions[[nMissFa]]
  else
    pFa = list(matrix(0L, ncol=1, nrow=1))

  if(nMissMo > 0)
    pMo = partitions[[nMissMo]]
  else
    pMo = list(matrix(0L, ncol=1, nrow=1))

  # Loop over all combinations
  res = vector(length(pFa) * length(pMo), mode = "list")
  for(i in seq_along(pFa)) for(j in seq_along(pMo)) {
    pf = pFa[[i]]
    pm = pMo[[j]]
    newfa = max(pf)
    newmo = max(pm)
    newpars = newfa + newmo

    # Matrix blocks to be added at the bottom of a
    FA = matrix(0L, ncol = n, nrow = newfa)
    FA[cbind(pf, missingFa)] = 1L

    MO = matrix(0L, ncol = n, nrow = newmo)
    MO[cbind(pm, missingMo)] = 1L

    bottom = rbind(FA, MO)

    # Gender restriction
    if(genderSym && newpars > 1) {
      inv_vec = .rowSums(bottom * rep(pows, each = newpars), m = newpars, n = n)
      inv = paste(.mysortInt(inv_vec), collapse = "-")
      if(inv %in% observedInvariants)
        next
      observedInvariants = c(observedInvariants, inv)
    }

    # Add blocks and create adjMatrix object
    adjExp = rbind(a, bottom)
    adjExp = cbind(adjExp, matrix(0L, ncol = newpars, nrow = nrow(adjExp)))

    # Check linear inbreeding
    if(checkInb && linearInbreeding(adjExp, descList))
      next

    # Create adjMatrix object
    sexExp = c(sex, rep(1L, newfa), rep(2L, newmo))
    A = adjMatrix(adjExp, sexExp, validate = F)

    # Remove superfluous added parents
    A = removeFounderParents(A, fou)

    res[[length(pMo) * (i-1) + j]] = A
  }

  res[!unlist(lapply(res, is.null))]
}

linearInbreeding = function(adj, descList = NULL, dist = 1) {
  n = dim(adj)[1]
  idvec = seq_len(n)

  if(is.null(descList)) {
    descList = lapply(idvec, function(id)
      dagDescendants(adj, i = id, minDist = dist))
  }

  hasMultipleKids = idvec[.rowSums(adj, n, n) > 1]
  for(i in hasMultipleKids) {
    kids = idvec[adj[i, ] == 1]
    for(k in kids)
      if(any(kids %in% descList[[k]]))
        return(TRUE)
  }

  FALSE
}


# Remove parents of original founders, unless these parents have other children
removeFounderParents = function(adj, fou) {
  if(length(fou) == 0)
    return(adj)

  idvec = seq_len(dim(adj)[1])

  remov = integer()
  for(id in fou) {
    pars = adj[, id]
    if(!any(adj[pars, -id])) # If parents have no other children
      remov = c(remov, idvec[pars])
  }

  if(length(remov) > 0) {
    sex = attr(adj, 'sex')
    adj = newAdjMatrix(adj[-remov, -remov], sex[-remov])
  }

  adj
}
