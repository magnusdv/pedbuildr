# pedbuildr 0.4.0

* **pedbuildr** now uses [**mirai**](https://mirai.r-lib.org/) for parallel likelihood computations in `reconstruct()`, giving substantial speedups in many cases. Users can enable parallelisation by starting workers with `mirai::daemons()` before calling `reconstruct()`. See `?reconstruct` for details.

* The `numCores` argument has been deprecated, as parallelisation is now controlled through the mirai framework.

* More robust likelihood calculations, due to more flexible loop breaking in **pedtools** and **pedprobr**. As a result, the likelihood function should now in principle "never" fail, removing the need for separate handling of such failures. For back compatibility, the `errPed` and `errIdx` elements of the output of `reconstruct()` are still present.

* `buildPeds()` is significantly faster than before in the default case with `linearInb = 0`, i.e., when mating between lineally related individuals is disallowed.


# pedbuildr 0.3.0

* New S3 class `pedCollection` with for output of `buildPeds()`, with separate `print()`, `plot()` and subsetting methods 

* `buildPeds()` has a new argument `maxInbreeding` indicating the highest tolerated inbreeding coefficient in the constructed pedigrees. The default value is set to 1/16 (as in children of first-cousin parents), which is appropriate for typical applications in human pedigrees. Set to 0 to disallow inbreeding completely, or to 1 to allow any inbreeding.

* New dataset `Tutankhamun` based on the pedigree and genotype data published by Hawass et al. *Ancestry and pathology in King Tutankhamun's family*. Jama (2010).

* General overhaul of code and organisation; the main functions run significantly faster now.


# pedbuildr 0.2.0

* Initial CRAN version
