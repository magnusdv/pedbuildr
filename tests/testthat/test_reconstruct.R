.data = function(nMark = 6L) {
  am = trioData[, seq_len(nMark), drop = FALSE]
  singletons(1:3, sex = c(1, 2, 1)) |>
    setMarkers(alleleMatrix = am, locusAttributes = "snp12")
}

.hyps = function() {
  list(nuclearPed(), singletons(1:3, sex = c(1, 2, 1)))
}

.trueLL = function(ped, am, loci) {
  ped |>
    setMarkers(alleleMatrix = am, locusAttributes = loci) |>
    pedprobr::likelihood(logbase = exp(1)) |>
    sum()
}

.recon = function(...) reconstruct(..., verbose = FALSE)

test_that("reconstruct gives correct likelihoods", {
  x = .data()
  peds = .hyps()
  am = getAlleles(x)
  loci = getLocusAttributes(x, checkComps = TRUE)

  res = .recon(x, pedlist = peds, sortResults = FALSE)

  expected = vapply(peds, \(p) .trueLL(p, am, loci), numeric(1))

  expect_equal(res$logliks, expected, tolerance = 1e-12)
})

test_that("ped input and explicit marker input are equivalent", {
  x = .data()
  peds = .hyps()
  am = getAlleles(x)
  loci = getLocusAttributes(x, checkComps = TRUE)

  from_x = .recon(x, pedlist = peds, sortResults = FALSE)

  explicit = .recon(alleleMatrix = am, loci = loci, sex = c(1, 2, 1),
                         pedlist = peds, sortResults = FALSE)

  expect_identical(explicit$labs, from_x$labs)
  expect_equal(explicit$logliks, from_x$logliks, tolerance = 1e-12)
})

test_that("sorting changes pedigree and likelihood order together", {
  x = .data()
  peds = .hyps()

  raw   = .recon(x, pedlist = peds, sortResults = FALSE)
  sorted = .recon(x, pedlist = peds, sortResults = TRUE)

  ord = order(raw$logliks, decreasing = TRUE)
  expect_identical(sorted$pedlist, raw$pedlist[ord])
  expect_identical(sorted$logliks, raw$logliks[ord])
})

test_that("regular mutation-model lumping preserves likelihood", {
  ped = nuclearPed()

  x = ped |>
    addMarker(geno = c("1/2", "2/3", "2/2"), alleles = 1:8,
              afreq = rep(1/8, 8), mutmod = "equal", rate = 0.01)

  res = .recon(x, pedlist = list(ped), sortResults = FALSE)

  expected = sum(pedprobr::likelihood(x, logbase = exp(1)))
  expect_equal(res$logliks, expected, tolerance = 1e-12)
})

test_that("inferPO rejects explicit notPO input", {
  x = .data(2L)

  expect_error(.recon(x, inferPO = TRUE, notPO = list(c("1", "2"))),
               "must be NULL")
})

