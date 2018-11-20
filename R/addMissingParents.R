#' Add missing parents to an adjacency matrix
#'
#' @param a An adjMatrix object
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
addMissingParents = function(a) {
  sex = attr(a, "sex")

  missingFa = which(colSums(a[sex == 1, , drop = F]) == 0)
  missingMo = which(colSums(a[sex == 2, , drop = F]) == 0)

  nMissFa = length(missingFa)
  nMissMo = length(missingMo)

  if(nMissFa > 7) stop2("More than 7 extra fathers needed: Too many possible combinations")
  if(nMissMo > 7) stop2("More than 7 extra mothers needed: Too many possible combinations")

  # Singletons (needed in a later step)
  singls = which(rowSums(a) + colSums(a) == 0)

  # Compute all set partitions for the fathers and the mothers
  if(nMissFa > 0) {
    pFa_matr = partitions::setparts(nMissFa)
    pFa = lapply(seq_len(ncol(pFa_matr)), function(i) pFa_matr[, i])
  }
  else
    pFa = list(matrix(0L, ncol=1, nrow=1))

  if(nMissMo > 0) {
    pMo_matr = partitions::setparts(nMissMo)
    pMo = lapply(seq_len(ncol(pMo_matr)), function(i) pMo_matr[, i])
  }
  else
    pMo = list(matrix(0L, ncol=1, nrow=1))

  # Loop over all combinations
  apply(expand.grid(pFa, pMo), 1, function(r) {
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

    sexExp = c(sex, rep(1L, newfa), rep(2L, newmo))
    res = adjMatrix(adjExp, sexExp)
    removeSingletonParents(res, singls)
  })
}

removeSingletonParents = function(adj, singls) {
  if(length(singls) == 0) return(adj)
  remov = integer()
  for(id in singls) {
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
