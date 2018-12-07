#' Pedigree reconstruction
#'
#' Reconstruct the most likely pedigree from genotype data
#'
#' @param alleleMatrix A matrix with two columns for each marker
#' @param loci A list of marker attributes
#' @param pedlist A list of pedigrees. If NULL, [buildPeds()] is used to
#'   generate a list
#' @param pairwise A logical. If TRUE, an initial stage of pairwise IBD
#'   estimation is done, in order to infer certain parent-child pairs, as well
#'   as certain *non*-parent-child pairs. When this option is used, arguments to
#'   `knownPO` and `notPO` are ignored
#' @param founderInb A number in the interval `[0,1]`, used as background
#'   inbreeding level in all founders
#' @param sortResults A logical. If TRUE, the output is sorted so that the most
#'   likely pedigree comes first
#' @param verbose A logical; verbose output or not
#' @param ... Additional parameters passed on to [buildPeds()]
#'
#' @return A list with two entries:
#'
#'   * `pedlist` : A list of pedigrees. This equals (except possibly for
#'   sorting) the input argument `pedlist` if this was given. If `sortResults =
#'   TRUE` then the list is sorted so that the most likely pedigrees come first
#'
#'   * `logliks` : A numerical vector containing the pedigree log-likelihoods
#'
#'   * `alleleMatrix` The input allele matrix
#'
#' @examples
#'
#' # Simulate genotype data for a trio family (increase N!)
#' x = forrel::markerSim(nuclearPed(1), N = 50, alleles = 1:3, seed = 123)
#'
#' # Extract allele matrix and locus attributes (frequencies a.s.o.)
#' m = getAlleles(x)
#' loci = lapply(x$markerdata, attributes)
#'
#' # Reconstruct the most likely pedigree from the data
#' res = reconstruct(m, loci, sex = c(1, 2, 1), connected = TRUE,
#'                   genderSym = TRUE, pair=TRUE)
#'
#' # Plot the best pedigrees
#' plotBestPeds(res)
#'
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom forrel IBDestimate showInTriangle
#' @export
reconstruct = function(alleleMatrix, loci, pedlist = NULL, pairwise = F, founderInb = 0, sortResults = T, verbose = T, ...) {
  ids = rownames(alleleMatrix)
  if(is.null(ids))
    ids = idsnum = 1:nrow(alleleMatrix)
  else {
    idsnum = suppressWarnings(as.integer(ids))
    if(anyNA(idsnum))
      stop2("Non-integer ID labels detected in `alleleMatrix`: ",
            ids[is.na(idsnum)])
  }
  if(verbose)
    cat("ID labels inferred from allele matrix:", toString(ids[order(idsnum)]), "\n")

  if(is.null(pedlist)) {
    if(!identical(idsnum, 1:max(idsnum)))
      stop2("\nFor now, only basic ID labels are supported. When `pedlist` is not supplied,\n",
            "the rownames of `alleleMatrix` must be either NULL or 1,2,...,N")

    args = list(...)

    if(pairwise) {
      if(length(args$knownPO) > 0 | length(args$knownPO) > 0)
        stop2("`knownPO` and `notPO` must be NULL when `pairwise = TRUE`")

      if(verbose) cat("Performing pairwise estimation in order to establish (unordered) parent-child pairs\n")
      POresult = inferPO(alleleMatrix, loci, list = TRUE)

      if(verbose) {
        forrel::showInTriangle(POresult$kappa, labels = T, new=T)
        cat("Certain parent-child:\n")
        print(POresult$PO, row.names = F)
        cat("\nCertain NON-parent-child:\n")
        print(POresult$UN, row.names = F)
        cat("\n")
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
    return(NULL)

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
      x = setMarkers(ped, allele_matrix = alleleMatrix, locus_annotations = loci)
      if(founderInb > 0) founder_inbreeding(x, founders(x)) = founderInb
    }
    else {
      x = lapply(ped, function(comp) {
        y = setMarkers(comp, allele_matrix = alleleMatrix, locus_annotations = loci)
        if(founderInb > 0) founder_inbreeding(y, founders(y)) = founderInb
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
    ord = order(logliks, decreasing = T)
    pedlist = pedlist[ord]
    logliks = logliks[ord]
  }

  list(pedlist = pedlist, logliks = logliks, alleleMatrix = alleleMatrix)
}


