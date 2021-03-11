#' Reconstruction example with three individuals
#'
#' This dataset contains simulated genotypes for 3 males at 100 SNP markers.
#'
#' @format A data frame with 3 rows and 104 columns. The first 4 columns contain
#'   pedigree info:
#'
#'   * `id`: Individual ID
#'
#'   * `fid`: Paternal ID, where 0 means missing father
#'
#'   * `mid`: Maternal ID, where 0 means missing mother
#'
#'   * `sex`: 1 = male; 2 = female
#'
#'   The remaining 100 columns contain the genotypes, in the form a/b, where a
#'   and b are the observed alleles (labelled 1 and 2).
#'
"trioData"
