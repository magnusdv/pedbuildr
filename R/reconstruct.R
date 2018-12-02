#' Pedigree reconstruction
#'
#' Reconstruct the most likely pedigree from genotype data
#'
#' @param alleleMatrix A matrix with two columns for each marker
#' @param loci A list of marker attributes
#' @param pedlist A list of pedigrees. If NULL, [buildPeds()] is used to
#'   generate a list
#' @param founderInb A number in the interval `[0,1]`, used as background
#'   inbreeding level in all founders
#' @param sortResults A logical. If TRUE, the output is sorted so that the most
#'   likely pedigree comes first
#' @param verbose A logical; verbose output or not
#' @param ... Additional parameters passed on to [buildPeds()]
#'
#' @return A list with two entries:
#'
#'  * `pedlist` : A list of pedigrees. This is the same as the input `pedlist`, but possibly sorted (if `sortResults = TRUE`)
#'
#'  * `logliks` : A numerical vector containing the pedigree log-likelihoods
#'
#' @examples
#' \dontrun{ # Requires the `forrel` package
#'
#' # Simulate genotype data for a trio family (increase N!)
#' x = forrel::markerSim(nuclearPed(1), N = 5, alleles = 1:3, seed = 123)
#'
#' # Extract allele matrix and locus attributes (frequencies a.s.o.)
#' m = getAlleles(x)
#' loci = lapply(x$markerdata, attributes)
#'
#' # Reconstruct the most likely pedigree from the data
#' res = reconstruct(m, loci, sex = c(1, 2, 1), connected = TRUE, genderSym = T)
#'
#' # Plot the best pedigrees
#' plotBestPeds(res)
#'
#' }
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
reconstruct = function(alleleMatrix, loci, pedlist = NULL, founderInb = 0, sortResults = T, verbose = T, ...) {
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
    pedlist = buildPeds(idsnum, verbose = verbose, ...)
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
    }
    else {
      x = lapply(ped, function(comp)
        setMarkers(comp, allele_matrix = alleleMatrix, locus_annotations = loci))
    }

    # Founder inbreeding
    if(founderInb > 0)
      founder_inbreeding(x, founders(x)) = founderInb

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

  list(pedlist = pedlist, logliks = logliks)
}


