# Estimate pairwise IBD coefficients and return
# list of certain parent-child and *not* parent-child.
inferPO = function(alleleMatrix, loci, list = FALSE) {
  ids = rownames(alleleMatrix)
  if(is.null(ids))
    ids = idsnum = 1:nrow(alleleMatrix)

  slist = lapply(ids, function(i)
    setMarkers(singleton(i, sex = 0), alleleMatrix = alleleMatrix, locusAttributes = loci))

  pairs = .comb2(ids)
  kappa = forrel::IBDestimate(slist, pairs)

  PO = kappa[kappa$k0 < 0.01 & kappa$k2 < 0.5, , drop = F]
  notPO = kappa[kappa$N != 0 & kappa$k0 > 0.5, , drop = F]

  if(list) {
    PO = lapply(seq_len(nrow(PO)), function(r) as.numeric(PO[r, 1:2]))
    notPO = lapply(seq_len(nrow(notPO)), function(r) as.numeric(notPO[r, 1:2]))
  }

  list(PO = PO, notPO = notPO, kappa = kappa)
}
