# pedbuildr 0.3.0

## New features

* New S3 class `pedCollection` with for output of `buildPeds()`, with separate `print()`, `plot()` and subsetting methods 

* `buildPeds()` has a new argument `maxInbreeding` indicating the highest tolerated inbreeding coefficient in the constructed pedigrees. The default value is set to 1/16 (as in children of first-cousin parents), which is appropriate for typical applications in human pedigrees. Set to 0 to disallow inbreeding completely, or to 1 to allow any inbreeding.

* New dataset `Tutankhamun` based on the pedigree and genotype data published by Hawass et al. *Ancestry and pathology in King Tutankhamun's family*. Jama (2010).

* General overhaul of code and organisation; the main functions run significantly faster now.


# pedbuildr 0.2.0

* Initial CRAN version
