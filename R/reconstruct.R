#' Pedigree reconstruction
#'
#' Reconstruct the most likely pedigree from genotype data
#'
#' @param x A `ped` object or a list of such.
#' @param ids A vector of ID labels from `x`. By default, the genotyped members of `x` are used.
#' @param alleleMatrix A matrix with two columns for each marker.
#' @param loci A list of marker attributes.
#' @param pedlist A list of pedigrees. If NULL, [buildPeds()] is used to
#'   generate a list.
#' @param pairwise A logical. If TRUE, an initial stage of pairwise IBD
#'   estimation is done, in order to infer certain parent-child pairs, as well
#'   as certain *non*-parent-child pairs. When this option is used, arguments to
#'   `knownPO` and `notPO` are ignored.
#' @param founderInb A number in the interval `[0,1]`, used as background
#'   inbreeding level in all founders
#' @param sortResults A logical. If TRUE, the output is sorted so that the most
#'   likely pedigree comes first.
#' @param verbose A logical; verbose output or not
#' @param ... Additional parameters passed on to [buildPeds()], e.g.,
#'   `sex`, `knownPO`, `notPO`, `connected`, `maxLinearInbreeding`, `genderSym`.
#'
#' @return An object of class `reconResult`, which is essentially list with the following entries:
#'
#'   * `pedlist` : A list of pedigrees, equal to the input argument `pedlist` if this was given. If `sortResults =
#'   TRUE` then the list is sorted so that the most likely pedigrees come first
#'
#'   * `logliks` : A numerical vector containing the pedigree log-likelihoods
#'
#'   * `alleleMatrix` : The input allele matrix
#'
#'   * `labels` : A character vector of ID labels
#'
#' @examples
#'
#' # Simulate genotype data for a trio family (increase N!)
#' x = forrel::markerSim(nuclearPed(1), N = 50, alleles = 1:3,
#'                       seed = 123, verbose = FALSE)
#'
#' # Extract marker data and reconstruct
#' res = reconstruct(x, pairwise = TRUE)
#'
#' # Plot most likely pedigrees
#' plot(res[1:6])
#'
#' # Alternative workflow: Extract data manually...
#' m = getAlleles(x)
#' loci = getLocusAttributes(x)
#' sex = getSex(x)
#'
#' # ...and then reconstruct
#' res2 = reconstruct(alleleMatrix = m, loci = loci, sex = sex,
#'                   pairwise = TRUE)
#'
#' #
#' y = forrel::markerSim(nuclearPed(children = c("s1", "s2")), ids = c("s1", "s2"),
#'                       N = 50, alleles = 1:2, seed = 123, verbose = FALSE)
#'
#' # Reconstruct
#' res3 = reconstruct(y, connected = FALSE)
#' plot(res3)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom forrel IBDestimate showInTriangle
#' @export
reconstruct = function(x, ids, alleleMatrix = NULL, loci = NULL,
                       pedlist = NULL, pairwise = FALSE, sortResults = TRUE,
                       founderInb = 0, verbose = TRUE, ...) {
  if(!missing(x)) {
    if(!is.ped(x) || is.pedList(x))
      stop2("Argument `x` must be a `ped` object, or a list of such")

    if(missing(ids))
      ids = typedMembers(x)

    if(is.null(alleleMatrix))
      alleleMatrix = getAlleles(x, ids)
    else
      stop2("When a ped/pedlist `x` is given as input, `alleleMatrix` must be NULL. (Data is extracted from `x`)")

    loci = getLocusAttributes(x, markers = loci)
  }
  else {
    ids = rownames(alleleMatrix)
    if(is.null(ids))
      ids = as.character(1:nrow(alleleMatrix))
  }

  # Change to numeric ids
  idsnum = seq_along(ids)
  rownames(alleleMatrix) = idsnum

  ### Build pedigree list
  if(is.null(pedlist)) {
    args = list(...)

    # Sex: If `x` is given, extract from this
    if(!missing(x)) {
      if(!is.null(args$sex))
        stop2("When a ped/pedlist `x` is given as input, `sex` must be NULL (it is extracted from `x`)")
      args$sex = getSex(x, ids)
    }

    # Do pairwise analysis first?
    if(pairwise) {
      if(length(args$knownPO) > 0 | length(args$knownPO) > 0)
        stop2("`knownPO` and `notPO` must be NULL when `pairwise = TRUE`")

      if(verbose) cat("Performing pairwise estimation in order to establish (unordered) parent-child pairs\n")
      POresult = inferPO(alleleMatrix, loci, list = TRUE)

      if(verbose) {
        forrel::showInTriangle(POresult$kappa, labels = TRUE, new = TRUE)
        cat("Certain parent-child:", toString(sapply(POresult$PO, paste, collapse = "-")), "\n")
        cat("Certain NON-parent-child:", toString(sapply(POresult$notPO, paste, collapse = "-")), "\n")
      }

      args$knownPO = POresult$PO
      args$notPO = POresult$notPO
    }

    args$ids = idsnum
    args$verbose = verbose
    pedlist = do.call(buildPeds, args)
  }

  npeds = length(pedlist)
  if(npeds == 0)
    stop2("Empty pedigree list")

  if(verbose) {
    cat("\nComputing likelihoods of", npeds, "pedigrees\n")
    pb = txtProgressBar(min = 0, max = npeds, style = 3) # Progress bar

  }

  # Compute likelihoods
  logliks = vapply(seq_len(npeds), function(i) {

    # Update progressbar
    if(verbose) setTxtProgressBar(pb, i)

    ped = pedlist[[i]]

    # Attach marker data
    if(is.ped(ped)) {
      x = setMarkers(ped, alleleMatrix = alleleMatrix, locusAttributes = loci)
      if(founderInb > 0) founderInbreeding(x, founders(x)) = founderInb
    }
    else {
      x = lapply(ped, function(comp) {
        y = setMarkers(comp, alleleMatrix = alleleMatrix, locusAttributes = loci)
        if(founderInb > 0) founderInbreeding(y, founders(y)) = founderInb
        y
      })
    }

    # Compute loglikelihood
    #tryCatch(loglikTotal(x), error = function(e) {plot(x); NA})
    loglikTotal(x)
  },
  FUN.VALUE = 0)

  if(verbose) {
    close(pb) # Close progress bar
  }

  if(sortResults) {
    ord = order(logliks, decreasing = TRUE)
    pedlist = pedlist[ord]
    logliks = logliks[ord]
  }

  structure(list(pedlist = pedlist,
                 logliks = logliks,
                 alleleMatrix = alleleMatrix,
                 labels = ids),
            class = "reconResult")
}

#' @export
`[.reconResult` = function(x, i) {
  structure(list(pedlist = x$pedlist[i],
                 logliks = x$logliks[i],
                 alleleMatrix = x$alleleMatrix,
                 labels = x$labels),
            class = "reconResult")
}

#' @importFrom graphics par plot text title
#' @export
plot.reconResult = function(x, titles = "logliks", nrow = NA, id.labels = x$labels,
                            shaded = x$labels, col = list(red = x$labels), ...) {

  L = length(x$pedlist)
  if(is.na(nrow))
    nrow = if(L<6) 1 else floor(sqrt(L))

  op = par(mfrow = c(nrow, ceiling(L/nrow)))
  on.exit(par(op))

  if(identical(titles, "logliks"))
    titles = paste("Loglik =", round(x$logliks, 2))

  labs = x$labels

  for(i in seq_along(x$pedlist)) {
    ped = x$pedlist[[i]]

    # Title
    tit = if(!is.null(titles)) titles[i] else NULL

    # Margin
    mar = if(is.null(tit)) c(1.5,1.5,1.5,1.5) else c(1.5,1.5,3,1.5)

    if(is.pedList(ped)) {
      plot(-1:1,-1:1, type="n", axes = F, xlab="", ylab="")
      mess ="Disconnected pedigree.\n\nPlot separately with\n`pedtools::plotPedList()`."
      text(0, 0, mess, cex = 1.3)
      title(tit, line= -2)
      next
    }

    ped = relabel(ped, old = seq_along(labs), labs)

    plot(ped, id.label = id.labels, shaded = shaded, col = col,
         margin = mar, title = tit, ...)
  }
}
