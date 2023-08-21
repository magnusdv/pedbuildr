## code to prepare `trioData` dataset goes here

library(forrel)

# Simulate data for three individuals, whose true pedigree is the following:
x = nuclearPed(fa = "X1", mother = "X2", children = 1:2)
x = addChildren(x, father = 2, mother = "X3", nch = 1, id = 3)

# Plot
plot(x, hatched = 1:3, col = list(red = 1:3))

# Simulate 100 SNPs
x = markerSim(x, N = 100, ids = 1:3, alleles = 1:2, seed = 1729)

# Genotype matrix
trioData = getGenotypes(x, ids = 1:3)

# Save
usethis::use_data(trioData, overwrite = TRUE)

#########################

# Test reconstruction (see also ?trioData)
y = singletons(1:3) |> setMarkers(alleleMatrix = trioData, locusAttributes = "snp12")
summary(y)

res = reconstruct(y)
plot(res, top = 6)
