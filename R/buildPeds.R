#' Build a list of pedigrees
#'
#' Build all pedigrees between a set of individuals, subject to given
#' restrictions.
#'
#' The parameter `extra` controls which of two algorithms are used to create the
#' pedigree list.
#'
#' If `extra` is a nonnegative integer, it determines the number of extra
#' individuals allowed in the iterative pedigree construction. These extras
#' start off with undetermined sex, meaning that both males and females are
#' used. It should be noted that the final pedigrees may contain additional
#' extras, since missing parents are added at the end.
#'
#' If `extra` is the word "parents", the algorithm is not iterative. It first
#' generates all directed acyclic graphs between the original individuals. Then
#' their parents are added and merged in all possible ways. This option has the
#' advantage of not requiring an explicit/ad hoc number of "extras", but works
#' best in smaller cases.
#'
#' @param labs A character vector of ID labels.
#' @param sex A vector of the same length as `labs`, with entries 1 (male) or 2
#'   (female).
#' @param extra Either the word "parents", or a nonnegative integer. See
#'   details.
#' @param age A numeric or character vector. If numeric, and `age[i] < age[j]`,
#'   then individual `i` will not be an ancestor of individual `j`. The numbers
#'   themselves are irrelevant, only the partial ordering. Note that no
#'   interpretation is made about individuals of equal age.
#'
#'   Alternatively `age` may be a character vector of inequalities, e.g., `age =
#'   c("1>2", "1>3")`. This syntax allows finer control than the numeric
#'   version.
#' @param knownPO A list of vectors of length 2, containing the ID labels of
#'   pairs known to be parent-offspring. By default, both directions are
#'   considered; use `age` to force a specific direction.
#' @param allKnown A logical. If TRUE, no other pairs than `knownPO` will be
#'   assigned as parent-offspring. If FALSE (default), all pairs except those in
#'   `notPO` are treated as potential parent-offspring.
#' @param notPO A list of vectors of length 2, containing the ID labels of pairs
#'   known *not* to be parent-offspring.
#' @param noChildren A vector of ID labels, indicating individuals without
#'   children of their own.
#' @param connected A logical. If TRUE (default), only connected pedigrees are
#'   returned.
#' @param linearInb Either TRUE (allow any linear inbreeding), FALSE (disallow
#'   linear inbreeding) or a nonnegative integer indicating the maximum
#'   separation linearly related spouses. For example, `linearInb = 1` allows
#'   mating between parent and child, but not between grandparent and grandchild
#'   (or more distant).
#' @param maxLinearInb Deprecated; replaced by `linearInb`.
#' @param sexSymmetry A logical. If TRUE (default), pedigrees which are equal
#'   except for the gender distribution of the *added* parents, are regarded as
#'   equivalent, and only one of each equivalence class is returned. Example:
#'   paternal vs. maternal half sibs.
#' @param verbose A logical.
#'
#' @return A list of pedigrees. Each element is a `ped` object or a list of
#'   such.
#'
#' @examples
#'
#' # Showing off a few of the options
#' plist = buildPeds(1:3, sex = c(1,2,1), extra = 1, knownPO = list(1:2),
#'                   age = "1 > 2", linearInb = FALSE)
#' stopifnot(length(plist) == 12)
#' # plotPeds(plist)
#'
#'
#' # Slightly different output with `extra = "parents"`
#' plist2 = buildPeds(1:3, sex = c(1,2,1), extra = "parents", knownPO = list(1:2),
#'                    age = "1 > 2", linearInb = FALSE)
#' stopifnot(length(plist2) == 8)
#' # plotPeds(plist2)
#'
#'
#' @export
buildPeds = function(labs, sex, extra = "parents", age = NULL, knownPO = NULL, allKnown = FALSE,
                     notPO = NULL, noChildren = NULL, connected = TRUE, linearInb = TRUE,
                     maxLinearInb = NULL, sexSymmetry = TRUE, verbose = TRUE) {

  if(!is.null(maxLinearInb)) {
    warning("The argument `maxLinearInb` has been replaced with `linearInb`. Use this in new code.")
    linearInb = maxLinearInb
  }

  N = length(labs)
  ids = seq_along(labs)

  if(length(sex) != N)
    stop2("`labs` and `sex` must have the same length\nlabs: ", labs, "\nsex: ", sex)

  if(!(identical(extra, "parents") || isCount(extra, minimum = 0)))
     stop2('`extra` must be either the word "parents" or a nonnegative integer: ', extra)

  if(allKnown) {
    if(is.null(knownPO))
      stop2("`knownPO` cannot be NULL when `allKnown = TRUE`")
    if(!is.null(notPO))
      stop2("`notPO` must be NULL when `allKnown = TRUE`")
  }

  checkLabs = c(unlist(knownPO), unlist(notPO), noChildren)
  if(!all(checkLabs %in% labs))
    stop2("Unknown pedigree member: ", setdiff(labs, checkLabs))

  knownPO_int = lapply(knownPO, function(p) match(p, labs))
  notPO_int = lapply(notPO, function(p) match(p, labs))
  noChildren_int = match(noChildren, labs)

  # Convert age vector into matrix with all ordered pairs (works with NULL)
  if(is.numeric(age))
    age = convertNumAge(age, labs)
  ageMat = parseAge(age, labs)

  # Check noChildren, and convert to internal
  if(!all(noChildren %in% labs))
    stop2("Unknown label in `noChildren`: ", setdiff(noChildren, labs))

  # Convert `linearInb` to `maxLinearInb`
  maxLinearInb = if(isTRUE(linearInb)) Inf else if(isFALSE(linearInb)) 0 else linearInb

  if(verbose) {
    .knownPO = sapply(knownPO, paste, collapse = "-")
    .notPO = sapply(notPO, paste, collapse = "-")
    .age = paste(labs[ageMat[,2]], labs[ageMat[,1]], sep = ">")

    toStr = function(...) toString(...) %e% "-"

    print(glue::glue("
      Pedigree parameters:
        ID labels: {toString(labs)}
        Sex: {toString(sex)}
        Extra: {extra}
        Age info: {toStr(.age)}
        Known PO: {toStr(.knownPO)}
        Known non-PO: {toStr(.notPO)}
        No children: {toStr(noChildren)}
        Connected only: {connected}
        Symmetry filter: {sexSymmetry}
        Linear inbreeding: {linearInb}"
    ))
  }

  if(extra == "parents") {
    buildPedsParents(labs = labs, sex = sex, ageMat = ageMat, knownPO = knownPO_int,
                     allKnown = allKnown, notPO = notPO_int, noChildren = noChildren_int,
                     connected = connected, maxLinearInb = maxLinearInb,
                     sexSymmetry = sexSymmetry, verbose = verbose)
  }
  else {
    buildPedsExtra(labs = labs, sex = sex, extra = extra, ageMat = ageMat, knownPO = knownPO_int,
                   allKnown = allKnown, notPO = notPO_int, noChildren = noChildren_int,
                   connected = connected, maxLinearInb = maxLinearInb,
                   sexSymmetry = sexSymmetry, verbose = verbose)
  }
}



# Convert numeric age vector to string inequalities: "A>B,C,D", B>C"
convertNumAge = function(age, labs) {
  if(!is.numeric(age))
    stop2("`age` is not numeric")

  N = length(labs)

  # If age is named, convert to full vector
  if(!is.null(nms <- names(age))) {
    if(!all(nms %in% labs))
      stop2("Unknown name in `age`: ", setdiff(nms, labs))
    age.full = rep(NA_real_, N)
    names(age.full) = labs
    age.full[nms] = age
    age = age.full
  }

  if(length(age) != N)
    stop2("When `age` is an unnamed numeric, its length must match `labs`")

  s = lapply(seq_len(N), function(i) {
    younger = age < age[i]
    younger = younger & !is.na(younger)

    if(any(younger))
      sprintf("%s>%s", labs[i], labs[younger])
  })

  unlist(s)
}

# Convert vector to string inequalities: "A>B,C,D", B>C" to matrix
parseAge = function(a, labs, output = c("matrix", "list")) {
  if(is.null(a) || all(is.na(a)))
    return(NULL)

  if(!all(good <- grepl(">", a)))
    stop2("Character '>' missing in `age` entry: ", a[!good])

  lst = lapply(strsplit(a, ">"), function(par) {
    par1 = lapply(strsplit(par, ","), trimws)
    names(par1) = c("o", "y")
    par1
  })

  # Bind to a single matrix
  # Note column order: y - o. (This is the natural choice in `directedAdjs`)
  res = do.call(rbind, lapply(lst, function(l) {
    cbind(younger = l$y,
          older = rep(l$o, each = length(l$y)))
  }))

  if(!all(res %in% labs))
    stop2("Unknown ID label in `age`: ", sort(setdiff(res, labs)))

  # Convert to internal ordering
  cbind(younger = match(res[, "younger"], labs),
        older = match(res[, "older"], labs))
}

