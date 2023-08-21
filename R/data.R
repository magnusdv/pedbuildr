#' Reconstruction example with three individuals
#'
#' This dataset contains simulated genotypes for 3 individuals at 100 SNP
#' markers.
#'
#' @format A matrix with 3 rows and 100 columns. Each entry contains a genotype
#'   in the form a/b.
#'
#' @examples
#' trioData[, 1:10]
#'
#' x = list(singleton(1), singleton(2), singleton(3)) |>
#'   setMarkers(alleleMatrix = trioData, locusAttributes ="snp12")
#'
#' x
#'
"trioData"


#' Pedigree of Tutankhamun
#'
#' A reconstructed pedigree of the Egyptian Pharaoh Tutankhamun, with genotypes
#' for 8 STR markers, as published by Hawass et al.
#'
#' @format A data frame with 7 rows and 12 columns:
#'
#'   * `id`,`fid`,`mid`,`sex`: Pedigree columns in standard format
#'   * `D13S317`, ...: Genotype columns for 8 markers
#'
#' @source Hawass et al. *Ancestry and pathology in King Tutankhamun's family*.
#'   Jama (2010).
#'
#' @examples
#' # Pedigree as published
#' plot(Tutankhamun)
#'
#' # Simple reconstruction, assuming all directly related
#' res = reconstruct(Tutankhamun, extra = 0, inferPO = TRUE, maxInbreeding = 1)
#' plot(res, top = 4)
#'
#' # Published ped is most likely (with these assumptions)
#' identical(res[[1]], Tutankhamun)
#'
"Tutankhamun"
