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
#' Since sex is not encoded in the adjacency matrix itself, we add the vector
#' of genders as an attribute `sex`.
#'
#'
#' @param adj A square matrix, converted internally to be of "logical" type. If
#'   `adj` is a vector it is recycled and converted to matrix, using `sex` to
#'   deduce the dimensions. If missing, a matrix with all elements FALSE is
#'   returned.
#' @param sex A gender vector, i.e. a vector where all elements are 1 (male) or
#'   2 (female).
#' @param validate A logical indicating if the validity of the created object
#'   should be checked.
#'
#' @return An object of class `adjMatrix`. This is a square logical matrix, with
#'   an attribute `sex` which is an integer vector of gender codes (1 = male; 2
#'   = female).
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
adjMatrix = function(adj, sex, validate = TRUE) {
  sex = as.integer(sex)

  if(missing(adj))
    adj = matrix(F, nrow = length(sex), ncol = length(sex))

  if(is.null(dim(adj)))
    adj = matrix(adj, nrow = length(sex))

  if(!is.matrix(adj))
    adj = as.matrix(adj)

  mode(adj) = "logical"

  obj = newAdjMatrix(adj, sex)
  if(validate)
    obj = validateAdjMatrix(obj)
  obj
}


# Constructor for class "adjMatrix"
newAdjMatrix = function(adj, sex) {
  if(!all(c(is.logical(adj), is.integer(sex), is.matrix(adj),
            nrow(adj) == ncol(adj), nrow(adj) == length(sex))))
    stop2("Wrong input to adjMatrix constructor")

  attr(adj, "sex") = sex
  class(adj) = "adjMatrix"
  adj
}


# Validator for "adjMatrix" objects
validateAdjMatrix = function(adj) {
  sex = attr(adj, "sex")
  if(!all(sex %in% 1:2))
    stop2("Illegal elements found in `sex` attribute vector: ", setdiff(sex, 1:2))

  if(nrow(adj) != (ncol(adj)))
    stop2("Adjacency matrix must be square")

  if(nrow(adj) != length(sex))
    stop2("Length of `sex` attribute must equal the number of columns")

  # Number of fathers for each indiv
  nFathers = colSums(adj[sex == 1, , drop = FALSE])
  nMothers = colSums(adj[sex == 2, , drop = FALSE])
  if(any(nFathers > 1))
    stop2("More than one father identified in column ", which(nFathers > 1)[1])
  if(any(nMothers > 1))
    stop2("More than one father identified in column ", which(nMothers > 1)[1])

  adj
}
