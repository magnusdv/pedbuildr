#' Pedigree adjacency matrix
#'
#' Construct adjacency matrices corresponding to pedigrees.
#'
#' A pedigree can be thought of a a directed acyclic graph (DAG), with the
#' pedigree members as nodes, and an edge from A to B if and only if A is a
#' parent of B. Since no pedigree member have more than 2 parents, any node of
#' the corresponding graphs has at most 2 incoming edges.
#'
#' The adjacency matrix of a graph G is a square logical matrix with TRUE in
#' entry `(i,j)` if and only if there is an edge from node `i` to node `j`.
#'
#' Since sex is not encoded in the adjacency matrix itself, we add the sex
#' information as an attribute `sex`. This is an integer vector with elements 0,
#' 1 or 2.
#'
#' @param adj A square matrix, converted internally to be of "logical" type. If
#'   `adj` is a vector it is recycled and converted to matrix, using `sex` to
#'   deduce the dimensions. If missing, a matrix with all elements FALSE is
#'   returned.
#' @param sex A vector of length `ncol(adj)` indicating the sex of each
#'   individual. All elements should be either 1 (male) or 2 (female) or 0
#'   (unknown).
#' @param connected A logical. If the underlying graph is known to be
#'   (dis)connected, this information can be stored here.
#' @param validate A logical indicating if the validity of the created object
#'   should be checked.
#'
#' @return An object of class `adjMatrix`. This is a square logical matrix, with
#'   an attribute `sex` which is an integer vector with elements in {0,1,2}.
#'
#' @examples
#' # An adjacency matrix corresponding to a family trio
#' # Member order: Father - mother - son
#' pedbuildr:::adjMatrix(c(0,0,0,0,0,0,1,1,0), sex = c(1,2,1))
#'
#' #' # An empty adjacency matrix for 3 male individuals
#' pedbuildr:::adjMatrix(sex = c(1,1,1))
#'
#' @noRd
adjMatrix = function(adj, sex, connected = NA, validate = TRUE) {
  sex = as.integer(sex)

  if(missing(adj))
    adj = matrix(F, nrow = length(sex), ncol = length(sex))

  if(is.null(dim(adj)))
    adj = matrix(adj, nrow = length(sex))

  if(!is.matrix(adj))
    adj = as.matrix(adj)

  mode(adj) = "logical"

  obj = newAdjMatrix(adj, sex, connected)

  if(validate)
    obj = validateAdjMatrix(obj)
  obj
}


# Constructor for class "adjMatrix"
newAdjMatrix = function(adj, sex, connected = NA) {
  if(!all(is.logical(adj), is.integer(sex), is.logical(connected)))
    stop2("Type error in adjMatrix constructor")
  dm = dim(adj)
  if(!all(length(dm) == 2, dm[1] == dm[2], dm[1] == length(sex), length(connected) == 1))
    stop2("Wrong input to adjMatrix constructor")

  attr(adj, "sex") = sex
  attr(adj, "connected") = connected
  class(adj) = "adjMatrix"
  adj
}


# Validator for "adjMatrix" objects
validateAdjMatrix = function(adj) {
  sex = attr(adj, "sex")
  if(!all(sex %in% 0:2))
    stop2("Illegal elements found in `sex` attribute vector: ", setdiff(sex, 0:2))

  if(nrow(adj) != (ncol(adj)))
    stop2("Adjacency matrix must be square")

  if(nrow(adj) != length(sex))
    stop2("Length of `sex` attribute must equal the number of columns")

  if(any(diag(adj)))
    stop2("Individual is its own parent: ", which(diag(adj)))

  # Number of parents for each indiv
  nPar = colSums(adj)
  if(any(nPar > 2))
    stop2("More than two parents identified in column ", which(nPar > 2))

  # Number of fathers/mothers for each indiv
  nFa = colSums(adj[sex == 1, , drop = FALSE])
  if(any(nFa > 1))
    stop2("More than one father identified in column ", which(nFa > 1))
  nMo = colSums(adj[sex == 2, , drop = FALSE]) # if(all(sex > 0)) nPar - nFa
  if(any(nMo > 1))
    stop2("More than one mother identified in column ", which(nMo > 1))

  adj
}
