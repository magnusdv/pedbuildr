# This code creates an object called `partitions` and stores it in "R/sysdata.rda".
# This is a lookup table of set partitions used in the `addMissingParents()` function.

library(partitions)

# List of all set partitions with n elements
setPartitions = function(n) {
  m = setparts(n)
  lapply(1:ncol(m), function(i) m[, i])
}


# Create list of lists
partitions = lapply(1:9, setPartitions)

# Store in data/ folder
usethis::use_data(partitions, internal = TRUE)
