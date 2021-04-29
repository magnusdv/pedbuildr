# Lump alleles and prepare alleleMatrix for fast marker creation

#' @importFrom pedmut isLumpable lumpedMatrix mutationModel
prepareData = function(alleleMatrix, loci, verbose = FALSE) {

  # Number of markers
  nMark = ncol(alleleMatrix)/2

  # Marker names in matrix (or NULL)
  if(!is.null(matNames <- colnames(alleleMatrix)))
    matNames = sub("\\.[^.]*$", "", matNames[2 * seq_len(nMark)])
  if(any(matNames == "" | matNames %in% seq_len(nMark)))
    matNames = NULL

  hasNames = !is.null(matNames)

  # Reduce and sort `loci`
  if(hasNames) {

    # Loci names: Either list names or 'name' attrib
    locNames = names(loci) %||% vapply(loci, function(a) a$name %||% "", "")
    idx = match(matNames, locNames, 0L)

    if(any(idx == 0L))
      stop2("Unknown marker name: ", matNames[idx == 0])

    loci = loci[idx]
  }
  else if(length(loci) != nMark)
    stop2("Allele matrix incompatible with list of loci")

  # Empty list of internal 2-column allele matrices
  amatList = vector("list", nMark)

  # Lumping starts here!
  for(k in seq_along(loci)) {

    attrs = loci[[k]]
    origAlleles = attrs$alleles

    # Observed alleles
    obsAll = alleleMatrix[, c(2*k-1, 2*k)]
    obs = unique.default(obsAll[!is.na(obsAll)])

    doLump = TRUE

    # No lumping if all, or all but one, are observed
    if (length(obs) >= length(origAlleles) - 1) {
      mess = sprintf("Lumping not needed: %d of %d alleles observed", length(obs), length(origAlleles))
      doLump = FALSE
    }

    if(doLump) {
      obsIdx = match(obs, origAlleles, 0L)
      lump = if(!length(obs)) origAlleles else origAlleles[-obsIdx]

      mut = attrs$mutmod

      # No lumping if mutation model is present and not lumpable for this lump
      if(!is.null(mut) && !isLumpable(mut, lump)) {
        mess = "Mutation model is not lumpable"
        doLump = FALSE
      }
    }

    if(doLump) {
      # Update alleles and freqs
      obsIdxSorted = sort.int(obsIdx, method = "shell")
      keepAls = origAlleles[obsIdxSorted]
      loci[[k]]$alleles = c(keepAls, "lump")

      obsFreq = attrs$afreq[obsIdxSorted]
      loci[[k]]$afreq = c(obsFreq, 1 - sum(obsFreq))

      # Update mutation model
      if (!is.null(mut)) {
        lumpedMale = pedmut::lumpedMatrix(mut$male, lump)
        lumpedFemale = pedmut::lumpedMatrix(mut$female, lump)
        loci[[k]]$mutmod = pedmut::mutationModel(list(female = lumpedFemale, male = lumpedMale))
      }

      mess = sprintf("Lumping: %d -> %d alleles", length(origAlleles), length(keepAls) + 1)
    }
    else {
      keepAls = origAlleles
    }

    if(verbose)
      message("Marker ", if(hasNames) matNames[k] else k, ": ", mess)

    # Update internal marker matrix
    amat = match(obsAll, keepAls, nomatch = 0L, incomparables = NA_character_)
    dim(amat) = dim(obsAll)
    rownames(amat) = rownames(obsAll)
    amatList[[k]] = amat
  }

  list(loci = loci, amatList = amatList)
}


setMarkersFAST = function(x, amatList, loci) {
  # If pedlist, recurse over components
  if(is.pedList(x)) {
    y = lapply(x, function(comp) setMarkersFAST(comp, amatList, loci))
    return(y)
  }

  ids = x$ID
  sex = x$SEX

  # genotyped indivs (extract from first amat)
  idsG = .myintersect(ids, rownames(amatList[[1]]))

  # template allele matrix with `ids` as row names
  mTmp = integer(2 * length(ids))
  dim(mTmp) = c(length(ids), 2)
  rownames(mTmp) = ids

  # Create marker objects
  mlist = lapply(seq_along(loci), function(k) {
    m = mTmp
    m[idsG, ] = amatList[[k]][idsG, ]

    loc = loci[[k]]

    newMarker(m, alleles = loc$alleles, afreq = loc$afreq, name = loc$name, chrom = loc$chrom,
              posMb = loc$posMb, mutmod = loc$mutmod, pedmembers = ids, sex = sex)
  })

  class(mlist) = "markerList"
  x$MARKERS = mlist
  x
}
