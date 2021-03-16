
# Build pedigrees between a set of individuals and a specified number of "extras".
# Missing parents are added in a minimal way.
# Note: knownPO, notPO, noChildren must be internal numerics
buildPedsExtra = function(labs, sex, extra = 0, ageMat = NULL, knownPO = NULL, allKnown = FALSE,
                      notPO = NULL, noChildren = NULL, connected = TRUE,
                      maxLinearInb = Inf, sexSymmetry = TRUE, verbose = FALSE) {

  st = Sys.time()
  N = length(labs)
  Ntot = N + extra
  sex = c(sex, rep(0L, extra))

  # Convert to list with known/imp relations *with lesser ID* for each
  knownPO_long = lapply(1:Ntot, function(i)
    unlist(lapply(knownPO, function(p) if(max(p) == i) min(p))))
  notPO_long = lapply(1:Ntot, function(i)
    unlist(lapply(notPO, function(p) if(max(p) == i) min(p))))

  # Convert age matrix to list format (works with NULL)
  y = unname(ageMat[, "younger"])
  o = unname(ageMat[, "older"])
  ageList = lapply(1:N, function(i) list(younger = y[o == i], older = o[y == i]))

  if(verbose) cat("\nBuilding pedigree list:\n")

  # Initial step: Create adjacency matrix for the first indiv.
  a1 = adjMatrix(sex = sex[1], connected = TRUE)
  attr(a1, "anc") = list(1L)

  Mlist = list(a1)

  # Wrapper for adding an individual
  addFUN = function(a, k) {
    addIndividual(a, addedSex = sex[k], origN = N, final = k == Ntot,
                  connected = connected, knownPO = knownPO_long[[k]], notPO = notPO_long[[k]],
                  ageList = ageList[1:k], noChildren = noChildren, maxLinearInb = maxLinearInb)
  }

  # Add main individuals one by one
  for(k in 1 + seq_len(N-1)) { # safer than 2:N
    Mlist = unlist(lapply(Mlist, function(a) addFUN(a, k)), recursive = FALSE)
    if(verbose && k < N)
      cat(sprintf("  First %d: %d candidates (%s)\n", k, length(Mlist), ftime(st)))
  }

  ALLSOLS = validSolutions(Mlist, connected = connected) # not checkLonely or sexSymmetry here

  if(verbose)
    cat(sprintf("  All %d + 0 extra: %d solutions | %d candidates (%s)\n",
                N, length(ALLSOLS), length(Mlist), ftime(st)))


  if(extra) {
    intermediate = vector("list", length = extra)

    for(ex in seq_len(extra)) {
      k = N + ex
      Mlist = unlist(lapply(Mlist, function(a) addFUN(a, k)), recursive = FALSE)
      nCand0 = length(Mlist)
      Mlist = removeDuplicates(Mlist, sexSymmetry)
      intermediate[[ex]] = validSolutions(Mlist, N, connected)

      if(verbose) {
        nCand = length(Mlist)
        cat(sprintf("  All %d + %d extra: %d solutions | %d candidates | %d duplicates removed (%s)\n",
                    N, ex, length(intermediate[[ex]]), nCand, nCand0 - nCand, ftime(st)))
      }
    }

    ALLSOLS = c(ALLSOLS, unlist(intermediate, recursive = FALSE))
  }

  if(verbose)
    cat("  Total solutions:", length(ALLSOLS), "\n  Converting to ped\n")

  peds = lapply(ALLSOLS, function(a) adj2ped(addMissingParents1(a), labs))
  if(verbose)
    cat("  Time used:", ftime(st), "\n")

  invisible(peds)
}


# Main function for adding 1 individual
addIndividual = function(a, addedSex, origN, final, connected = FALSE,
                         knownPO = NULL, notPO = NULL, ageList = NULL,
                         noChildren = NULL, maxLinearInb = Inf) {

  # If unknown sex, add male and female in two steps
  if(addedSex == 0L) {
    # Male
    addedMale = addIndividual(a, 1L, origN = origN, final = final, connected = connected,
                              knownPO = knownPO, notPO = notPO, ageList = ageList,
                              noChildren = noChildren, maxLinearInb = maxLinearInb)
    # Female
    addedFemale = addIndividual(a, 2L, origN = origN, final = final, connected = connected,
                                knownPO = knownPO, notPO = notPO, ageList = ageList,
                                noChildren = noChildren, maxLinearInb = maxLinearInb)

    return(c(addedMale, addedFemale))
  }

  n = dim(a)[1]
  k = n + 1
  extra = n >= origN

  checkAge = length(ageList) > 0
  thisAge = ageList[[k]]
  checkInb = maxLinearInb < Inf

  prev = seq_len(n)
  sex = attr(a, "sex")

  # Potential parents
  mal = which(sex == 1)
  fem = which(sex == 2)

  excludeAsParent = c(notPO, noChildren, thisAge$younger)
  potentialFa = c(0, .mysetdiff(mal, excludeAsParent, makeUnique = FALSE))
  potentialMo = c(0, .mysetdiff(fem, excludeAsParent, makeUnique = FALSE))

  # Potential children
  if(k %in% noChildren)
    potentialCh = integer(0)
  else {
    excludeAsChild = c(notPO, thisAge$older)
    missingParent = switch(addedSex,
                           which(.colSums(a[mal, , drop = F], length(mal), n) == 0),
                           which(.colSums(a[fem, , drop = F], length(fem), n) == 0))
    potentialCh = .mysetdiff(missingParent, excludeAsChild, makeUnique = FALSE)
  }

  # Prepare main loop
  len = length(potentialFa) * length(potentialMo) * 2^length(potentialCh)
  RES = vector("list", length = len)

  aExt = c(rbind(a, rep(FALSE, n)), rep(FALSE, n + 1))
  dim(aExt) = c(k, k)
  sexExt = as.integer(c(sex, addedSex))
  tmpA = newAdjMatrix(aExt, sex = sexExt)

  tmpAnc = attr(a, "anc")
  length(tmpAnc) = k # extend by 1

  # Main loop!
  i = 0
  for(fa in potentialFa) for(mo in potentialMo) {

    # Vector of actual parents
    par = c(fa, mo)
    par = par[par > 0]
    npar = length(par)

    # Quick check of parent-child inbreeding
    if(maxLinearInb == 0 && npar == 2 && any(a[par, par]))
      next

    # Ancestors of parents (including parents themselves)
    if(npar == 0)
      parAnc = integer(0)
    else if(npar == 1)
      parAnc = tmpAnc[[par]]
    else
      parAnc = unique.default(c(tmpAnc[[par[1]]], tmpAnc[[par[2]]]))

    # Potential children: Exclude ancestors
    potCh = .mysetdiff(potentialCh, parAnc, makeUnique = FALSE)

    # Force known PO
    known = .mysetdiff(knownPO, par, makeUnique = FALSE)
    if(!all(known %in% potCh))
      next

    # Loop over subsets of children
    for(chs in powerset(potCh, force = known)) {
      # message("Fa = ", fa, "; mo = ", mo, "; ch = ", toString(chs))

      # If final extra, skip if uninformative: a) leaf  b) founder with 1 child
      if(final && extra && (length(chs) == 0 || (length(chs) == 1 && npar == 0)))
        next

      # Another quick check of parent-child inbreeding
      if(maxLinearInb == 0 && any(a[c(par, chs), c(par, chs)]))
        next

      # Add relations
      newA = tmpA
      newA[par, k] = TRUE
      newA[k, chs] = TRUE

      # More general check for linear inbreeding
      if(checkInb && linearInb2(newA, minDist = maxLinearInb + 1))
        next

      # Check for superfluous extras
      if(extra) {
        hasLonely = hasLonelyExtras(newA, origN)
        if(final && hasLonely)
          next
        attr(newA, "hasLonely") = hasLonely
      }

      # Check connectedness (used by adj2ped even if `connected = F`)
      if(final || k >= origN) {
        isConn = isConnected(newA)
        if(final && connected && !isConn)
          next
        attr(newA, "connected") = isConn
      }

      # Update ancestors
      if(!final || checkAge) {
        newAnc = tmpAnc
        newAnc[[k]] = c(k, parAnc)

        # Identify descendants
        isDesc = vapply(prev, function(d)
          !(d %in% par) && (d %in% chs || any(chs %in% newAnc[[d]])), logical(1))
        desc = which(isDesc)

        # Add new ancestors to all descendants
        if(length(parAnc))
          newAnc[desc] = lapply(newAnc[desc], function(v) unique.default(c(v, k, parAnc)))
        else
          newAnc[desc] = lapply(newAnc[desc], function(v) c(v, k))

        attr(newA, "anc") = newAnc
      }

      # Check age data
      if(checkAge && ageViol(newA, newAnc, ageList))
        next

      # Canonical ordering of extras
      if(extra && k > origN + 1)
        newA = canonicalPerm(newA, origN)

      # Add to output list
      RES[[i+1]] = newA
      i = i+1
    }
  }

  # Remove NULLs
  length(RES) = i

  # Return
  RES
}


# Detect uninformative extras: i) leaves and ii) founders with 1 child
hasLonelyExtras = function(a, N) {
  Ntot = dim(a)[1]
  Nex = Ntot - N
  if(Nex == 0)
    return(FALSE)

  ### Leaves

  # Number of children of each extra
  extraCh = .rowSums(a[N + seq_len(Nex), , drop = FALSE], Nex, Ntot)
  if(any(extraCh == 0))
    return(TRUE)

  ### Uninformative founders

  # Number of parents of each extra
  extraPar = .colSums(a[, N + seq_len(Nex), drop = FALSE], Ntot, Nex)
  if(any(extraPar == 0 & extraCh == 1))
    return(TRUE)

  FALSE
}

# Ad hoc standard permutation of the extras
# NB: Not perfect
canonicalPerm = function(a, N) {
  n = dim(a)[1]
  # If at most 1 extra, return
  if(n <= N + 1)
    return(a)

  Nseq = seq_len(N)
  exSeq = N + seq_len(n - N)
  sex = attr(a, "sex")

  # Order according to children, parents and sex
  block1 = lapply(Nseq, function(i) a[exSeq, i])
  block2 = lapply(Nseq, function(i) a[i, exSeq])
  args = c(block1, block2, list(sex[exSeq], decreasing = TRUE, method = "shell"))
  ord = do.call(order, args)
  p = c(Nseq, N + ord)

  # Permute
  a[] = a[p, p]
  attr(a, "sex") = sex[p]
  if(!is.null(anc <- attr(a, "anc"))) {
    pAnc = lapply(anc[p], function(v) match(v, p))
    attr(a, "anc") = pAnc
  }

  a
}

removeDuplicates = function(DA, sexSymmetry) {
  if(sexSymmetry)
    dups = duplicated.default(lapply(DA, as.integer))
  else
    dups = duplicated.default(lapply(DA, function(da) c(as.integer(da), attr(da, "sex"))))

  DA[!dups]
}

validSolutions = function(DA, checkLonely = FALSE, connected = FALSE) {

  if(length(DA) == 0)
    return(DA)

  if(checkLonely) {
    DA = DA[!sapply(DA, attr, "hasLonely")]
    # print(length(DA))

    if(length(DA) == 0)
      return(DA)
  }

  if(connected) {
    DA = DA[sapply(DA, attr, "connected")]
    # print(length(DA))

    if(length(DA) == 0)
      return(DA)
  }

  DA
}

ageViol = function(a, anc, ageList) {
  for(i in 1:dim(a)[1]) {
    y = ageList[[i]]$younger
    if(length(y) && any(y %in% anc[[i]]))
      return(TRUE)
  }
  return(FALSE)
}
