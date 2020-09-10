#' Pedigree reconstruction
#'
#' Reconstruct the most likely pedigree from genotype data
#'
#' @param x A `ped` object or a list of such.
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
#'   inbreeding level in all founders
#' @param sortResults A logical. If TRUE, the output is sorted so that the most
#'   likely pedigree comes first.
#' @param pairwise,genderSym Deprecated.
#' @param verbose A logical; verbose output or not
#' @inheritParams buildPeds
#'
#' @return An object of class `pedrec`, which is essentially list with the
#'   following entries:
#'
#'   * `pedlist`: A list of pedigrees, either built by [buildPeds()] or as
#'   supplied in the input argument `pedlist`. If `sortResults = TRUE`, the
#'   list is sorted so that the most likely pedigrees come first
#'
#'   * `logliks`: A numerical vector containing the pedigree log-likelihoods
#'
#'   * `alleleMatrix`: A matrix of marker alleles
#'
#'   * `loci`: A list of marker locus attributes
#'
#'   * `labels`: A character vector of ID labels
#'
#' @examples
#'
#' # Load `forrel` (for simulation)
#' library(forrel)
#'
#' ### Example 1: Family trio
#' x = nuclearPed()
#' x = markerSim(x, N = 50, alleles = 1:3, seed = 123, verbose = FALSE)
#'
#' # Reconstruct
#' res = reconstruct(x, inferPO = TRUE)
#'
#' # Plot most likely pedigrees
#' plot(res, top = 6)
#'
#' # Alternative workflow: Extract data manually...
#' m    = getAlleles(x)
#' loci = getLocusAttributes(x)
#' sex  = getSex(x)
#'
#' # ...and then reconstruct
#' res2 = reconstruct(alleleMatrix = m, loci = loci, sex = sex,
#'                   inferPO = TRUE)
#'
#' stopifnot(identical(res, res2))
#'
#' ### Example 2: Full siblings
#' ids = c("s1", "s2")
#' y = nuclearPed(children = ids)
#' y = markerSim(y, N = 50, ids = ids, seed = 123, verbose = FALSE)
#'
#' # Reconstruct and plot
#' res3 = reconstruct(y, connected = FALSE)
#' plot(res3)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom forrel showInTriangle
#' @export
reconstruct = function(x, ids, alleleMatrix = NULL, loci = NULL,
                       pedlist = NULL, inferPO = FALSE, sex = NULL,
                       age = NULL, knownPO = NULL, allKnown = FALSE,
                       notPO = NULL, noChildren = NULL, connected = TRUE,
                       maxLinearInb = Inf, sexSymmetry = TRUE,
                       sortResults = TRUE,
                       founderInb = 0, pairwise = NULL, genderSym = NULL,
                       verbose = TRUE) {
  if(!is.null(pairwise)) {
    message("Argument `pairwise` has been renamed to `inferPO` and will be removed in a future version")
    inferPO = pairwise
  }
  if(!is.null(genderSym)) {
    message("Argument `genderSym` has been renamed to `sexSymmetry` and will be removed in a future version")
    sexSymmetry = genderSym
  }
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
      ids = as.character(1:nrow(alleleMatrix))
  }

  # Change to numeric ids
  if(!is.null(knownPO))
    knownPO = lapply(knownPO, function(p) sort(match(p, ids)))
  if(!is.null(notPO))
    notPO = lapply(notPO, function(p) sort(match(p, ids)))

  # If age is named, convert to full vector
  if(is.numeric(age) && !is.null(names(age))) {
    if(!all(names(age) %in% ids))
      stop2("Unknown name in `age`: ", setdiff(names(age), ids))
    age.full = rep(NA_real_, length(ids))
    age.full[match(names(age), ids)] = age
    age = age.full
  }
  rownames(alleleMatrix) = seq_along(ids)

  kappa = NULL

  ### Build pedigree list
  if(is.null(pedlist)) {

    # (Optional) pairwise analysis to establish parent-child relationships
    if(inferPO) {

      if(length(knownPO) > 0 | length(knownPO) > 0)
        stop2("`knownPO` and `notPO` must be NULL when `inferPO = TRUE`")

      POresult = inferPO(alleleMatrix, loci, list = TRUE)
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

    pedlist = buildPeds(ids = ids, sex = sex, age = age,
                        knownPO = knownPO, allKnown = allKnown,
                        notPO = notPO, noChildren = noChildren,
                        connected = connected, maxLinearInb = maxLinearInb,
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
                 labels = ids,
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
                 loci = x$loci,
                 labels = x$labels),
            class = "pedrec")
}

#' @export
`[[.pedrec` = function(x, i) {
  y = x$pedlist[[i]]
  y = setMarkers(y, alleleMatrix = x$alleleMatrix, locusAttributes = x$loci)
  relabel(y, new = x$labels, old = seq_along(x$labels))
}

#' @importFrom glue glue
#' @export
print.pedrec = function(x, ...) {
  print(glue::glue("
    Pedigree reconstruction result.
    Input: {length(x$labels)} individuals typed with {ncol(x$alleleMatrix)/2} markers.
    Ouput: {length(x$pedlist)} pedigrees sorted by likelihood ({length(x$errPeds)} failed)."
  ))
}


#' @importFrom graphics par plot text title
#' @export
plot.pedrec = function(x, top = NULL, nrow = NA, titles = "LR", labs = x$labels,
                       hatched = x$labels, col = list(red = x$labels), ...) {

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

    plot(ped, labs = labs, hatched = hatched, col = col,
         margin = mar, title = tit, ...)
  }
}
