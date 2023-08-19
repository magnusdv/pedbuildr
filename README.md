
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pedbuildr

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/pedbuildr)](https://CRAN.R-project.org/package=pedbuildr)
[![](https://cranlogs.r-pkg.org/badges/grand-total/pedbuildr?color=yellow)](https://cran.r-project.org/package=pedbuildr)
[![](https://cranlogs.r-pkg.org/badges/last-month/pedbuildr?color=yellow)](https://cran.r-project.org/package=pedbuildr)
<!-- badges: end -->

The goal of **pedbuildr** is to reconstruct pedigrees from genotype
data. This is done by optimising the likelihood over all possible
pedigrees subject to given restrictions. As part of the **ped suite**
ecosystem of R packages for pedigree analysis, it uses on
[pedtools](https://CRAN.R-project.org/package=pedtools) for handling
pedigrees, and imports
[pedprobr](https://CRAN.R-project.org/package=pedprobr) for calculating
pedigree likelihoods.

See also the `ibdEstimate()` function of
[forrel](https://CRAN.R-project.org/package=forrel), which does pairwise
relatedness estimation.

## Installation

The **pedbuildr** package can be installed from CRAN as follows:

``` r
install.packages("pedbuildr")
```

Alternatively, the latest development version is available from GitHub.

``` r
# install.packages("devtools")
devtools::install_github("magnusdv/pedbuildr")
```

## A reconstruction example

The built-in dataset `trioData` contains simulated genotypes for three
individuals at 100 SNP markers. As a simple demonstration we will try to
reconstruct the pedigree connecting these individuals.

To get started, load **pedbuildr**.

``` r
library(pedbuildr)
#> Loading required package: pedtools
```

The `trioData` is a data frame in so-called *ped format*. Here are the
first few columns:

``` r
trioData[, 1:10]
#>   id fid mid sex <1> <2> <3> <4> <5> <6>
#> 1  1   0   0   1 1/2 1/2 1/1 1/2 1/2 1/2
#> 2  2   0   0   1 1/1 1/2 1/2 2/2 1/2 2/2
#> 3  3   0   0   1 1/1 1/2 1/2 2/2 1/1 1/2
```

We convert the data into a pedigree object, using the `as.ped()`
function from **pedtools**.

``` r
x = as.ped(trioData, locusAttributes = "snp12")
summary(x)
#> List of 3 singletons.
#> Labels: 1 (male), 2 (male), 3 (male).
#> 100 attached markers.
```

The `locusAttributes` argument tells R that all the markers are
diallelic with alleles 1 and 2. By default, the alleles have equal
frequencies.

To reconstruct the pedigree, simply run `reconstruct()`:

``` r
res = reconstruct(x)
#> Pedigree parameters:
#>   ID labels: 1, 2, 3
#>   Sex: 1, 1, 1
#>   Extra: parents
#>   Age info: -
#>   Known PO: -
#>   Known non-PO: -
#>   No children: -
#>   Connected only: TRUE
#>   Symmetry filter: TRUE
#>   Max inbreeding: 0.0625
#>   Linear inbreeding: FALSE
#> 
#> Building pedigree list:
#>   Undirected adjacency matrices: 8 
#>   Directed adjacency matrices: 16 
#>   After adding parents: 57 
#>   Connected solutions: 44 
#>   Excessive inbreeding: 0 
#> 
#> Computing the likelihood of 44 pedigrees.
#> Sorting by descending likelihood.
#> Total time used:  2.87 secs
```

A tailor-made `plot` function makes it easy to visualise the most likely
pedigrees:

``` r
plot(res, top = 6)
```

<img src="man/figures/README-reconstruct-trio-1.png" style="display: block; margin: auto;" />

## Further options

*A priori* there are infinitely many possible pedigrees connecting a set
of individuals. (For example, two individuals may be *k*’th cousins for
any *k* = 1,2,… .) In order to obtain a manageable search space,
`reconstruct()` offers a range of restriction parameters:

- `extra`: The number of extra individuals allowed to connect the
  original individuals. (See further explanations below.)
- `age`: A numerical age vector, or a character vector describing age
  inequalities. For example, `age = c("A > B,C", "B > D")` excludes B
  and C as ancestors of A, and D as an ancestor of B (but no assumption
  is made for C vs. D).
- `inferPO`: If TRUE, an initial stage of pairwise IBD estimation is
  done, in order to infer certain parent-child pairs, as well as certain
  non-parent-child pairs.
- `knownPO`: Known parent–offspring pairs.
- `notPO`: Pairs known not to be parent–offspring.
- `allKnown`: If TRUE, then `knownPO` is the complete list of
  parent–offspring pairs.
- `noChildren`: Individuals known to have no children.
- `maxInbreeding`: The highest permitted inbreeding coefficient in the
  pedigree. By default this is set to 1/16, which is appropriate for
  most projects involving human pedigrees.
- `linearInb`: The allowed level of inbreeding between linear
  descendants. For example, `linearInb = 1` allows mating between
  parent–child, but not grandparent–grandchild. The default value
  `FALSE` disallows all inbreeding of this type.
- `connected`: If TRUE (default), only connected pedigrees are
  considered.
- `sexSymmetry`: If TRUE (default) pedigrees are considered equal if
  they differ only in the sexes of the added parents, e.g., paternal
  versus maternal half-siblings.

Let us re-run the reconstruction of `trioData` adding a few of these
restrictions. We allow 3 extra individuals and indicate that individual
1 is older than the others. Furthermore, we ask the program to infer
parent-child relationships automatically, and disallow linear
inbreeding.

``` r
res2 = reconstruct(x, extra = 3, age = "1 > 2,3", inferPO = TRUE)
#> Pairwise estimation:
#>   PO: 2-3 
#>   non-PO: 1-3 
#> 
#> Pedigree parameters:
#>   ID labels: 1, 2, 3
#>   Sex: 1, 1, 1
#>   Extra: 3
#>   Age info: 1>2, 1>3
#>   Known PO: 2-3
#>   Known non-PO: 1-3
#>   No children: -
#>   Connected only: TRUE
#>   Symmetry filter: TRUE
#>   Max inbreeding: 0.0625
#>   Linear inbreeding: FALSE
#> 
#> Building pedigree list:
#>   First 2: 2 candidates (0.0038 secs)
#>   All 3 + 0 extra: 1 solutions | 3 candidates (0.00572 secs)
#>   All 3 + 1 extra: 9 solutions | 30 candidates | 21 duplicates removed (0.0163 secs)
#>   All 3 + 2 extra: 35 solutions | 266 candidates | 501 duplicates removed (0.349 secs)
#>   All 3 + 3 extra: 183 solutions | 183 candidates | 459 duplicates removed (1.7 secs)
#>   Total solutions: 228 
#>   Converting to ped
#>   Excessive inbreeding: 114 
#>   Time used: 2.33 secs 
#> 
#> Computing the likelihood of 114 pedigrees.
#> Sorting by descending likelihood.
#> Total time used:  22.6 secs
```

The most likely results this time are shown below:

``` r
plot(res2, top = 6)
```

<img src="man/figures/README-reconstruct-trio-2-1.png" style="display: block; margin: auto;" />

We see that the same pedigree “wins”, but some esoteric alternatives
have appeared among the runners-up.

## More about `extra`

Arguably the most important parameter to `reconstruct()` is `extra`,
which controls the size of the pedigrees to consider. It can be either a
nonnegative integer, or the word “parents”. If an integer, it sets the
maximum number of extra members used to connect the original
individuals. (The final pedigrees may contain further extras still,
since missing parents are added at the end.)

If `extra = "parents"`, a special algorithm is invoked. First all
directed acyclic graphs between the original individuals are generated,
and then parents are added in all possible ways. This is (currently) the
default behaviour, since it avoids setting an *ad hoc* number of
“extras”. However, it only works well in relatively small cases.
