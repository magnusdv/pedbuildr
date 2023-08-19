#' Reconstruction example with three individuals
#'
#' This dataset contains simulated genotypes for 3 males at 100 SNP markers.
#'
#' @format A data frame with 3 rows and 104 columns. The first 4 columns contain
#'   pedigree info in standard format:
#'
#'   * `id`: Individual ID
#'   * `fid`: Paternal ID, where 0 means missing father
#'   * `mid`: Maternal ID, where 0 means missing mother
#'   * `sex`: 1 = male; 2 = female
#'
#'   The remaining 100 columns contain the genotypes, in the form a/b, where a
#'   and b are the observed alleles (labelled 1 and 2).
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
#'   See help page for [trioData] for further details on format.
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
