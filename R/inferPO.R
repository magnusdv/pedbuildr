# Estimate pairwise IBD coefficients and return
# list of certain parent-child and *not* parent-child.

#' @importFrom utils packageVersion
inferPO = function(alleleMatrix, loci, list = FALSE, verbose = TRUE) {
  ids = rownames(alleleMatrix)
  if(is.null(ids))
    ids = 1:nrow(alleleMatrix)

  slist = lapply(ids, function(i)
    setMarkers(singleton(i, sex = 0), alleleMatrix = alleleMatrix, locusAttributes = loci))

  pairs = .comb2(ids)

  if(packageVersion("forrel") <= "1.2") {
    kappa = forrel::IBDestimate(slist, pairs)
    PO = kappa[kappa$k0 < 0.01 & kappa$k2 < 0.5, , drop = FALSE]
    notPO = kappa[kappa$N != 0 & kappa$k0 > 0.5, , drop = FALSE]
  }
  else {
    kappa = forrel::ibdEstimate(slist, pairs, verbose = verbose)
    PO = kappa[kappa$k0 < 0.01 & kappa$k2 < 0.5, ]
    notPO = kappa[kappa$N != 0 & kappa$k0 > 0.5, ]
  }

  if(list) {
    PO = lapply(seq_len(nrow(PO)), function(r) PO[r, 1:2])
    notPO = lapply(seq_len(nrow(notPO)), function(r) notPO[r, 1:2])
  }

  list(PO = PO, notPO = notPO, kappa = kappa)
}
