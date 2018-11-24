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
#'
#' @return A list of adjMatrix objects where all columns sum to either 0 or 2.
#'
#' @examples
#' a = adjMatrix(c(1,0,0,0), sex=c(1,1))
#' addMissingParents(a)
#'
#' b = adjMatrix(rbind(rep(1,4), 0,0,0), sex=c(1,1,1,1))
#' addMissingParents(b)
#'
#' @importFrom partitions setparts
#' @export
addMissingParents = function(a, maxLinearInbreeding = Inf, partitions = NULL) {
  sex = attr(a, "sex")

  missingFa = which(colSums(a[sex == 1, , drop = F]) == 0)
  missingMo = which(colSums(a[sex == 2, , drop = F]) == 0)

  nMissFa = length(missingFa)
  nMissMo = length(missingMo)

  if(nMissFa > 7) stop2("More than 7 extra fathers needed: Too many possible combinations")
  if(nMissMo > 7) stop2("More than 7 extra mothers needed: Too many possible combinations")

  # Founders (needed in a later step)
  fou = which(colSums(a) == 0)

  # List of descendants
  checkInb = maxLinearInbreeding < Inf
  if(checkInb)
    descList = lapply(1:ncol(a), function(id)
      dagDescendants(a, id, minDist = maxLinearInbreeding + 1))

  # All set partitions for the fathers and the mothers
  if(nMissFa > 0)
    pFa = if(is.null(partitions)) setPartitions(nMissFa) else partitions[[nMissFa]]
  else
    pFa = list(matrix(0L, ncol=1, nrow=1))

  if(nMissMo > 0)
    pMo = if(is.null(partitions)) setPartitions(nMissMo) else partitions[[nMissMo]]
  else
    pMo = list(matrix(0L, ncol=1, nrow=1))

  # Loop over all combinations
  res = apply(expand.grid(pFa, pMo), 1, function(r) {
    pf = r$Var1
    pm = r$Var2
    newfa = max(pf)
    newmo = max(pm)

    # Matrix blocks to be added at the bottom of a
    FA = matrix(0L, ncol = ncol(a), nrow = newfa)
    FA[cbind(pf, missingFa)] = 1L

    MO = matrix(0L, ncol = ncol(a), nrow = newmo)
    MO[cbind(pm, missingMo)] = 1L

    # Add blocks and create adjMatrix object
    adjExp = rbind(a, FA, MO)
    adjExp = cbind(adjExp, matrix(0L, ncol = newfa + newmo, nrow = nrow(adjExp)))

    # Check linear inbreeding
    if(checkInb && linearInbreeding(adjExp, descList))
      return(NULL)

    sexExp = c(sex, rep(1L, newfa), rep(2L, newmo))
    res = adjMatrix(adjExp, sexExp, validate = F)
    removeFounderParents(res, fou)
  })

  res[!sapply(res, is.null)]
}

linearInbreeding = function(adj, descList = NULL, dist = 1) {
  if(is.null(descList)) {
    descList = lapply(1:nrow(adj), function(id)
      dagDescendants(adj, i = id, minDist = dist))
  }

  hasMultipleKids = which(rowSums(adj) > 0)
  for(i in hasMultipleKids) {
    kids = which(adj[i, ] == 1)
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

  remov = integer()
  for(id in fou) {
    pars = adj[, id]
    if(!any(adj[pars, -id])) # If parents have no other children
      remov = c(remov, which(pars))
  }

  if(length(remov) > 0) {
    sex = attr(adj, 'sex')
    adj = newAdjMatrix(adj[-remov, -remov], sex[-remov])
  }

  adj
}
