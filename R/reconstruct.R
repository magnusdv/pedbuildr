#' Pedigree reconstruction
#'
#' Reconstruct the most likely pedigree from genotype data.
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
#' @param x A `pedtools::ped` object or a list of such.
#' @param ids A vector of ID labels from `x`. By default, the genotyped members
#'   of `x` are used.
#' @param alleleMatrix A matrix with two columns for each marker.
#' @param loci A list of marker attributes.
#' @param pedlist A list of pedigrees. If NULL, [buildPeds()] is used to
#'   generate a list.
#' @param inferPO A logical. If TRUE, an initial stage of pairwise IBD
#'   estimation is done, in order to infer certain parent-child pairs, as well
#'   as certain *non*-parent-child pairs. When this option is used, arguments to
#'   `knownPO` and `notPO` are ignored.
#' @param founderInb A number in the interval `[0,1]`, used as background
#'   inbreeding level in all founders.
#' @param sortResults A logical. If TRUE, the output is sorted so that the most
#'   likely pedigree comes first.
#' @param verbose A logical.
#' @inheritParams buildPeds
#'
#' @return An object of class `pedrec`, which is essentially list with the
#'   following entries:
#'
#'   * `pedlist`: A list of pedigrees, either built by [buildPeds()] or as
#'   supplied in the input argument `pedlist`. If `sortResults = TRUE`, the list
#'   is sorted so that the most likely pedigrees come first
#'
#'   * `logliks`: A numerical vector of pedigree log-likelihoods
#'
#'   * `kappa`: A data frame with pairwise estimates (if `inferPO = TRUE`)
#'
#'   * `alleleMatrix`: A matrix of marker alleles
#'
#'   * `loci`: A list of marker locus attributes
#'
#'   * `errPeds`: A list of pedigrees for which the likelihood calculation
#'   failed
#'
#'   * `errIdx`: The indices of pedigrees in `errPeds` as elements of `pedlist`
#'
#' @examples
#' #-----------------
#' # Example 1: Trio
#' #-----------------
#'
#' data(trioData)
#'
#' x = as.ped(trioData, locusAttributes = "snp-12")
#' summary(x)
#'
#' res = reconstruct(x, inferPO = TRUE, age = "1 > 2", linearInb = FALSE)
#'
#' # Plot most likely pedigrees
#' plot(res, top = 6)
#'
#'
#' ### Alternative workflow: Extract data manually...
#' als  = getAlleles(x)
#' loci = getLocusAttributes(x)
#' sex  = getSex(x)
#'
#' # ...and then reconstruct
#' res2 = reconstruct(alleleMatrix = als, loci = loci, sex = sex,
#'                    inferPO = TRUE, age = "1 > 2", linearInb = FALSE)
#'
#' stopifnot(identical(res, res2))
#'
#' #--------------------
#' # Example 2: Siblings
#' #--------------------
#'
#' ids = c("s1", "s2")
#' y = nuclearPed(children = ids)
#'
#' # Simulate data
#' y = forrel::markerSim(y, N = 50, ids = ids, seed = 123)
#'
#' # Reconstruct and plot
#' res3 = reconstruct(y, connected = FALSE)
#' res3
#' plot(res3)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom forrel showInTriangle
#' @export
reconstruct = function(x, ids, extra = "parents", alleleMatrix = NULL, loci = NULL,
                       pedlist = NULL, inferPO = FALSE, sex = NULL,
                       age = NULL, knownPO = NULL, allKnown = FALSE,
                       notPO = NULL, noChildren = NULL, connected = TRUE,
                       linearInb = TRUE, maxLinearInb = NULL, sexSymmetry = TRUE,
                       sortResults = TRUE, founderInb = 0,
                       verbose = TRUE) {


  st = Sys.time()

  if(!missing(x)) {
    if(!is.ped(x) && !is.pedList(x))
      stop2("Argument `x` must be a `ped` object, or a list of such")

    if(!is.null(sex))
      stop2("When a ped/pedlist `x` is given as input, `sex` must be NULL (it is extracted from `x`)")
    if(!is.null(alleleMatrix))
      stop2("When a ped/pedlist `x` is given as input, `alleleMatrix` must be NULL. (Data is extracted from `x`)")

    if(missing(ids))
      ids = typedMembers(x)
    ids = as.character(ids) # in case of pedlist

    sex = getSex(x, ids)
    alleleMatrix = getAlleles(x, ids)
    loci = getLocusAttributes(x, markers = loci)
  }
  else {
    ids = rownames(alleleMatrix)
    if(is.null(ids))
      ids = rownames(alleleMatrix) = as.character(1:nrow(alleleMatrix))
  }

  kappa = NULL

  ### Build pedigree list
  if(is.null(pedlist)) {

    # (Optional) pairwise analysis to establish parent-child relationships
    if(inferPO) {

      if(length(knownPO) > 0 | length(knownPO) > 0)
        stop2("`knownPO` and `notPO` must be NULL when `inferPO = TRUE`")

      POresult = inferPO(alleleMatrix, loci, list = TRUE, verbose = FALSE)
      knownPO = POresult$PO
      notPO = POresult$notPO
      kappa = POresult$kappa

      if(verbose) {
        po = toString(sapply(knownPO, paste, collapse = "-")) %e% "None identified"
        notpo = toString(sapply(notPO, paste, collapse = "-")) %e% "None identified"
        cat("Pairwise estimation:\n")
        cat("  PO:", po, "\n")
        cat("  non-PO:", notpo, "\n\n")
      }
    }

    pedlist = buildPeds(labs = ids, sex = sex, extra = extra, age = age,
                        knownPO = knownPO, allKnown = allKnown,
                        notPO = notPO, noChildren = noChildren,
                        connected = connected, linearInb = linearInb, maxLinearInb = maxLinearInb,
                        sexSymmetry = sexSymmetry, verbose = verbose)
  }

  npeds = length(pedlist)
  if(npeds == 0)
    stop2("Empty pedigree list")

  if(verbose)
    cat("\nComputing the likelihood of", npeds, "pedigrees.\n")

  # Progress bar
  if(progbar <- verbose && interactive())
    pb = txtProgressBar(min = 0, max = npeds, style = 3)

  # Compute likelihoods
  logliks = vapply(seq_len(npeds), function(i) {

    # Update progressbar
    if(progbar) setTxtProgressBar(pb, i)

    ped = pedlist[[i]]

    # Attach marker data
    x = setMarkers(ped, alleleMatrix = alleleMatrix, locusAttributes = loci)

    # Founder inbreeding
    if(founderInb > 0) {
      if(is.pedList(x))
        x = lapply(x, function(comp) `founderInbreeding<-`(comp, founders(comp), value = founderInb))
      else
        founderInbreeding(x, founders(x)) = founderInb
    }

    # Compute loglikelihood
    tryCatch(loglikTotal(x), error = function(e) NA_real_)
  },
  FUN.VALUE = 0)

  # Deal with failed likelihoods
  errs = is.na(logliks)
  errIdx = which(errs)
  errPeds = pedlist[errs]
  if(any(errs)) {
    pedlist = pedlist[!errs]
    logliks = logliks[!errs]
  }

  if(progbar) {
    close(pb) # Close progress bar
  }

  if(sortResults) {
    if(verbose)
      cat("Sorting by descending likelihood.\n")
    ord = order(logliks, decreasing = TRUE)
    pedlist = pedlist[ord]
    logliks = logliks[ord]
  }

  time = Sys.time() - st
  if(verbose)
    cat("Total time used: ", format(time, digits = 3), "\n")

  structure(list(pedlist = pedlist,
                 logliks = logliks,
                 kappa = kappa,
                 alleleMatrix = alleleMatrix,
                 loci = loci,
                 errPeds = errPeds,
                 errIdx = errIdx),
            class = "pedrec")
}


#' @export
`[.pedrec` = function(x, i) {
  structure(list(pedlist = x$pedlist[i],
                 logliks = x$logliks[i],
                 kappa = x$kappa,
                 alleleMatrix = x$alleleMatrix,
                 loci = x$loci),
            class = "pedrec")
}

#' @export
`[[.pedrec` = function(x, i) {
  y = x$pedlist[[i]]
  setMarkers(y, alleleMatrix = x$alleleMatrix, locusAttributes = x$loci)
}

#' @importFrom glue glue
#' @export
print.pedrec = function(x, ...) {
  am = x$alleleMatrix
  print(glue::glue("
    Pedigree reconstruction result.
    Input: {nrow(am)} individuals typed with {ncol(am)/2} markers.
    Ouput: {length(x$pedlist)} pedigrees sorted by likelihood ({length(x$errPeds)} failed)."
  ))
}


#' @importFrom graphics par plot text title
#' @export
plot.pedrec = function(x, top = NULL, nrow = NA, titles = "LR",
                       labs = rownames(x$alleleMatrix),
                       hatched = labs, col = list(red = labs), ...) {

  if(!is.null(top))
    x = x[seq_len(top)]

  L = length(x$pedlist)
  if(is.na(nrow))
    nrow = if(L<6) 1 else floor(sqrt(L))

  op = par(mfrow = c(nrow, ceiling(L/nrow)))
  on.exit(par(op))

  # Titles
  if(identical(titles, "LR")) {
    ll = x$loglik
    titles = bquote("Loglik" == .(round(ll[1], 2))) # bquote to avoid bold
    if(L > 1) {
      lr = sprintf("%.4g", exp(ll[1] - ll))
      titles2 = sapply(2:L, function(i) bquote(LR[1:.(i)] == .(lr[i])))
      titles = c(titles, titles2)
    }
  }

  for(i in seq_len(L)) {
    ped = x$pedlist[[i]]

    # Title
    tit = if(!is.null(titles)) titles[i] else NULL

    # Margin
    mar = if(is.null(tit)) c(1.5,1.5,1.5,1.5) else c(1.5,1.5,3,1.5)

    if(is.pedList(ped)) {
      plot(-1:1,-1:1, type="n", axes = F, xlab="", ylab="")
      mess = sprintf("Nr %d: disconnected.\n\nPlot separately with\n`plotPedList()`.", i)
      text(0, 0, mess, cex = 1.3)
      title(tit, line= -2)
      next
    }

    plot(ped, labs = labs, hatched = hatched, col = col,
         margin = mar, title = tit, ...)
  }
}
