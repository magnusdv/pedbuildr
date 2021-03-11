## code to prepare `trioData` dataset goes here

library(forrel)

# Simulate data for three individuals, whose true pedigree is the following:
x = nuclearPed(fa = "X1", mother = "X2", children = 1:2)
x = addChildren(x, father = 2, mother = "X3", nch = 1, id = 3)

# Plot
plot(x, hatched = 1:3, col = list(red = 1:3))

# Simulate 100 SNPs
x = markerSim(x, N = 100, ids = 1:3, alleles = 1:2, seed = 1729)

# Convert to data frame and remove ped info
trioData = as.data.frame(x)[internalID(x, 1:3), ]
trioData$fid = trioData$mid = 0
rownames(trioData) = NULL

usethis::use_data(trioData, overwrite = TRUE)

#########################

# Test reconstruction (see also ?trioData)
y = as.ped(trioData, locusAttributes = "snp-12")
summary(y)

res = reconstruct(y)
plot(res, top = 6)

plot(reconstruct(y, extra = 3, age = "1 > 2,3", linearInb = FALSE), 6)
