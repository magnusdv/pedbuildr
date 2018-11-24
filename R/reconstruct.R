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
#'
#' @examples
#' x = forrel::markerSim(nuclearPed(1), N = 10, alleles = 1:3)
#' m = getAlleles(x)
#' loci = lapply(x$markerdata, attributes)
#'
#' res = reconstruct(m, loci, sex = c(1, 2, 1), connected = TRUE)
#'
#'
#' @export
reconstruct = function(alleleMatrix, loci, pedlist = NULL, founderInb = 0, sortResults = T, verbose = T, ...) {
  rownms = rownames(alleleMatrix)
  if(!is.null(rownms) && !all.equal(rownms, as.character(seq_along(rownms))))
    stop2("Pedigree labels other than 1,2, ... are not supported yet.\n",
          "(Hence rownames of `alleleMatrix` must be 1,2, ... or NULL.)")

  if(is.null(pedlist)) {
    if(verbose) cat("Building pedigree list\n")
    pedlist = buildPeds(1:nrow(alleleMatrix), verbose = verbose, ...)
  }
  else {
    if(verbose) cat("Pedigree list:", length(pedlist), "pedigrees\n")
  }

  if(verbose) cat("\nComputing likelihoods of", length(pedlist), "pedigrees...")
  logliks = vapply(pedlist, function(ped) {

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

  if(verbose) cat("done!\n")

  if(sortResults) {
    ord = order(logliks, decreasing = T)
    pedlist = pedlist[ord]
    logliks = logliks[ord]
  }

  list(pedlist = pedlist, logliks = logliks)
}

# Plot top list
#plotPeds(tutlist[top], nrow=3, cex = 1.2, margins = c(3,3,3,3),
#         titles = sprintf("lnlik = %.2f; Fmax = %.2f", liks[top], inb[top]))


