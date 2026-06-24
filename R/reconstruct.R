#' Pedigree reconstruction
#'
#' Reconstructs the most likely pedigree from genotype data.
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
#' ## Parallel likelihood computation
#' Likelihood calculations can be distributed across several R processes using
#' the **mirai** package. To use this, start mirai workers before calling
#' `reconstruct()`:
#'
#' ```r
#' mirai::daemons(8)  # Start 8 workers
#' res = reconstruct(x)
#' mirai::daemons(0)  # Stop workers
#' ```
#'
#' The call to `mirai::daemons(8)` starts 8 worker processes. These workers
#' remain available for later calls until stopped with `mirai::daemons(0)`.
#'
#' If no mirai workers are running, `reconstruct()` still works and runs the
#' likelihood calculations synchronously.
#'
#' The best number of workers depends on the number of pedigrees, marker data
#' size, and available memory. For large reconstruction runs, values such as 6,
#' 8, 12 or more may give substantial gains.
#'
#' @param x A `pedtools::ped` object or a list of such.
#' @param ids A vector of ID labels from `x`. By default, the genotyped members
#'   of `x` are used.
#' @param alleleMatrix A matrix with two columns for each marker. By default
#'   extracted from `x`
#' @param loci A list of marker attributes. By default extracted from `x`.
#' @param pedlist A list of pedigrees to optimise over. If NULL, [buildPeds()]
#'   is used to generate a list.
#' @param inferPO A logical. If TRUE, an initial stage of pairwise IBD
#'   estimation is done to infer high-confidence parent/child pairs, and also
#'   *non*-parent/child pairs. When this option is used, `knownPO` and `notPO`
#'   must be NULL.
#' @param founderInb A number in the interval `[0,1]`, used as background
#'   inbreeding level in all founders. Default: 0.
#' @param sortResults A logical. If TRUE (default), the output is sorted so that
#'   the most likely pedigree comes first.
#' @param numCores Deprecated. Use `mirai::daemons()` to set the number of
#'   workers.
#' @param verbose A logical.
#' @inheritParams buildPeds
#'
#' @return An object of class `pedrec`, which is essentially list with the
#'   following entries:
#'
#'   * `labs`: The individual labels as given in `ids`.
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
#' # Built-in dataset `trioData`
#' x = singletons(1:3) |>
#'   setMarkers(alleleMatrix = trioData, locusAttributes ="snp12")
#'
#' res = reconstruct(x, inferPO = TRUE, age = "1 > 2")
#'
#' # Plot most likely pedigrees
#' plot(res, top = 6)
#'
#'
#' #--------------------
#' # Example 2: Siblings
#' #--------------------
#' library(forrel)
#'
#' ids = c("s1", "s2")
#'
#' # Create pedigree and simulate profiles with 20 STR markers
#' y = nuclearPed(children = ids) |>
#'   profileSim(markers = NorwegianFrequencies[1:20], ids = ids, seed = 123)
#'
#' # Reconstruct allowing 2 extra individuals and any inbreeding
#' res2 = reconstruct(y, extra = 2, maxInbreeding = 1)
#' plot(res2, top = 6)
#'
#' # With mutation modelling
#' y = setMutmod(y, model = "equal", rate = 0.01)
#' res3 = reconstruct(y, extra = 2, maxInbreeding = 1)
#' plot(res3, top = 6)
#'
#' @importFrom mirai daemons daemons_set info mirai_map
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
reconstruct = function(x, ids, extra = "parents", alleleMatrix = NULL, loci = NULL,
                       pedlist = NULL, inferPO = FALSE, sex = NULL,
                       age = NULL, knownPO = NULL, knownSub = NULL, allKnown = FALSE,
                       notPO = NULL, noChildren = NULL, connected = TRUE,
                       maxInbreeding = 1/16, linearInb = FALSE, sexSymmetry = TRUE,
                       sortResults = TRUE, founderInb = 0, numCores = 1,
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
    loci = getLocusAttributes(x, markers = loci, checkComps = TRUE)
  }
  else {
    ids = rownames(alleleMatrix)
    if(is.null(ids))
      ids = rownames(alleleMatrix) = as.character(1:nrow(alleleMatrix))

    # Clean up unknown alleles (NA in `getAlleles()`)
    NAstrings = c("", "0", "-", "NA")
    alleleMatrix[alleleMatrix %in% NAstrings] = NA_character_

    loci = pedtools:::checkLocusAttribs(loci)
  }

  # Remove untyped rows of allele matrix
  typed = rowSums(is.na(alleleMatrix)) < ncol(alleleMatrix)
  if(sum(typed) < length(ids))
    alleleMatrix = alleleMatrix[typed, , drop = FALSE]

  ### Prep: allele lumping + prepare fast marker creation
  data = prepareData(alleleMatrix, loci)
  loci = data$loci
  amatList = data$amatList

  kappa = NULL

  ### Build pedigree list
  if(is.null(pedlist)) {

    # (Optional) pairwise analysis to establish parent-child relationships
    if(inferPO) {

      if(length(knownPO) > 0 | length(notPO) > 0)
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
                        knownPO = knownPO, knownSub = knownSub, allKnown = allKnown,
                        notPO = notPO, noChildren = noChildren,
                        connected = connected, maxInbreeding = maxInbreeding,
                        linearInb = linearInb, sexSymmetry = sexSymmetry, verbose = verbose)
  }

  npeds = length(pedlist)
  if(npeds == 0)
    stop2("Empty pedigree list")

  if(verbose)
    cat("\nComputing the likelihood of", npeds, "pedigrees.\n")

  if(!missing(numCores)) {
    s = "`numCores` has been deprecated; use `mirai::daemons()` for parallelisation. See `?reconstruct`."
    warning(s, call. = FALSE)
  }

  progbar = verbose && interactive()

  # Use mirai workers if available, otherwise run synchronously
  hasDaemons = mirai::daemons_set()
  nworkers = if(hasDaemons) mirai::info()[["connections"]] else 1L


  # Multiple workers
  if(nworkers > 1) {

    # Split pedigrees into batches of 8
    nchunks = min(npeds, 8L * nworkers)
    chunkSize = ceiling(npeds / nchunks)
    pedChunks = split(pedlist, ceiling(seq_len(npeds) / chunkSize))

    if(verbose)
      message(sprintf("Parallelisation: %d mirai workers; %d batches of ~%d pedigrees",
              nworkers, length(pedChunks), chunkSize))

    # Compute likelihoods in parallel
    mm = mirai::mirai_map(pedChunks, lapply,
      .args = list(FUN = loglikWrapper, amatList = amatList, loci = loci, founderInb = founderInb))

    loglikMirai = if(progbar) mm[mirai::.progress] else mm[]
    logliks = unlist(loglikMirai, recursive = TRUE, use.names = FALSE)

  }
  else {

    if(verbose)
      message("No parallelisation; use `mirai::daemons(n)` to set up multiple workers.")

    if(progbar)
      pb = txtProgressBar(min = 0, max = npeds, style = 3)

    logliks = vapply(seq_len(npeds), FUN.VALUE = numeric(1), function(i) {
      if(progbar) setTxtProgressBar(pb, i)
      loglikWrapper(pedlist[[i]], amatList, loci, founderInb)
    })

    if(progbar) close(pb)
  }

  # Deal with failed likelihoods
  errs = is.na(logliks)
  errIdx = which(errs)
  errPeds = pedlist[errs]
  if(any(errs)) {
    pedlist = pedlist[!errs]
    logliks = logliks[!errs]
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
    cat("Total time used:", format(time, digits = 3), "\n")

  structure(list(labs = ids,
                 pedlist = pedlist,
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
  structure(list(labs = x$labs,
                 pedlist = x$pedlist[i],
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
  typed = if(length(x$labs) == nrow(am)) "all" else nrow(am)

  print(glue::glue("
    Pedigree reconstruction result.
    Input: {length(x$labs)} individuals; {typed} typed with {ncol(am)/2} markers.
    Ouput: {length(x$pedlist)} pedigrees sorted by likelihood ({length(x$errPeds)} failed)."
  ))
}


#' @importFrom graphics par plot text title
#' @export
plot.pedrec = function(x, top = NULL, nrow = NA, titles = "LR",
                       labs = x$labs, highlight = x$labs,
                       hatched = rownames(x$alleleMatrix), ...) {

  if(!is.null(top))
    x = x[seq_len(top)]

  L = length(x$pedlist)
  if(is.na(nrow))
    nrow = if(L<6) 1 else floor(sqrt(L))

  op = par(mfrow = c(nrow, ceiling(L/nrow)))
  on.exit(par(op))

  # Titles
  if(identical(titles, "LR")) {
    ll = x$logliks
    titles = bquote("Loglik" == .(round(ll[1], 2))) # bquote to avoid bold
    if(L > 1) {
      lr = sprintf("%.4g", exp(ll[1] - ll))
      titles2 = sapply(2:L, function(i) bquote(LR[1:.(i)] == .(lr[i])))
      titles = c(titles, titles2)
    }
  }

  # Highlight typed members
  if(length(highlight) > 0)
    col = list(red = highlight)
  else
    col = NULL

  for(i in seq_len(L)) {
    ped = x$pedlist[[i]]
    if(is.pedList(ped))
      ped = unclass(ped)

    # Title
    tit = if(!is.null(titles)) titles[i] else NULL

    # Margin
    mar = if(is.null(tit)) c(1,2,1.5,2) else c(1,2,3,2)

    plot(ped, labs = labs, hatched = hatched, col = col,
         margin = mar, title = tit, ...)
  }
}
