% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reconstruct.R
\name{reconstruct}
\alias{reconstruct}
\title{Pedigree reconstruction}
\usage{
reconstruct(
  x,
  ids,
  extra = "parents",
  alleleMatrix = NULL,
  loci = NULL,
  pedlist = NULL,
  inferPO = FALSE,
  sex = NULL,
  age = NULL,
  knownPO = NULL,
  knownSub = NULL,
  allKnown = FALSE,
  notPO = NULL,
  noChildren = NULL,
  connected = TRUE,
  maxInbreeding = 1/16,
  linearInb = FALSE,
  sexSymmetry = TRUE,
  sortResults = TRUE,
  founderInb = 0,
  numCores = 1,
  verbose = TRUE
)
}
\arguments{
\item{x}{A \code{pedtools::ped} object or a list of such.}

\item{ids}{A vector of ID labels from \code{x}. By default, the genotyped members
of \code{x} are used.}

\item{extra}{Either the word "parents" (default), or a non-negative integer.
See Details.}

\item{alleleMatrix}{A matrix with two columns for each marker. By default
extracted from \code{x}}

\item{loci}{A list of marker attributes. By default extracted from \code{x}.}

\item{pedlist}{A list of pedigrees to optimise over. If NULL, \code{\link[=buildPeds]{buildPeds()}}
is used to generate a list.}

\item{inferPO}{A logical. If TRUE, an initial stage of pairwise IBD
estimation is done to infer high-confidence parent/child pairs, and also
\emph{non}-parent/child pairs. When this option is used, arguments to \code{knownPO}
and \code{notPO} are ignored.}

\item{sex}{A vector of the same length as \code{labs}, with entries 1 (male) or 2
(female).}

\item{age}{A numeric or character vector. If numeric, and \code{age[i] < age[j]},
then individual \code{i} will not be an ancestor of individual \code{j}. The numbers
themselves are irrelevant, only the partial ordering. (No inference is made
about individuals of equal age.) Alternatively, for finer control, \code{age}
may be a character vector of inequalities, e.g., \code{age = c("1>2", "1>3")}.}

\item{knownPO}{A list of vectors of length 2, containing the ID labels of
pairs known to be parent-offspring. By default, both directions are
considered; use \code{age} to force a specific direction.}

\item{knownSub}{A \code{ped} object involving a subset of the \code{labs} individuals.}

\item{allKnown}{A logical. If TRUE, no other pairs than \code{knownPO} will be
assigned as parent-offspring. If FALSE (default), all pairs except those in
\code{notPO} are treated as potential parent-offspring.}

\item{notPO}{A list of vectors of length 2, containing the ID labels of pairs
known \emph{not} to be parent-offspring.}

\item{noChildren}{A vector of ID labels, indicating individuals without
children of their own.}

\item{connected}{A logical. If TRUE (default), only connected pedigrees are
returned.}

\item{maxInbreeding}{A single numeric indicating the highest permitted
inbreeding coefficient. Default: 1/16 (as with first-cousin parents.)}

\item{linearInb}{A parameter controlling the maximum separation of linearly
related spouses. Either TRUE (allow all linear inbreeding), FALSE (disallow
all) or a non-negative integer. For example, \code{linearInb = 1} allows
parent/child mating, but not grandparent/grandchild or more distant linear
relatives. Default: FALSE.}

\item{sexSymmetry}{A logical. If TRUE (default), pedigrees which are equal
except for the gender distribution of the \emph{added} parents, are regarded as
equivalent, and only one of each equivalence class is returned. Example:
paternal vs. maternal half sibs.}

\item{sortResults}{A logical. If TRUE (default), the output is sorted so that
the most likely pedigree comes first.}

\item{founderInb}{A number in the interval \verb{[0,1]}, used as background
inbreeding level in all founders. Default: 0.}

\item{numCores}{A positive integer. The number of cores used in
parallelisation. Default: 1.}

\item{verbose}{A logical.}
}
\value{
An object of class \code{pedrec}, which is essentially list with the
following entries:
\itemize{
\item \code{labs}: The individual labels as given in \code{ids}.
\item \code{pedlist}: A list of pedigrees, either built by \code{\link[=buildPeds]{buildPeds()}} or as
supplied in the input argument \code{pedlist}. If \code{sortResults = TRUE}, the list
is sorted so that the most likely pedigrees come first
\item \code{logliks}: A numerical vector of pedigree log-likelihoods
\item \code{kappa}: A data frame with pairwise estimates (if \code{inferPO = TRUE})
\item \code{alleleMatrix}: A matrix of marker alleles
\item \code{loci}: A list of marker locus attributes
\item \code{errPeds}: A list of pedigrees for which the likelihood calculation
failed
\item \code{errIdx}: The indices of pedigrees in \code{errPeds} as elements of \code{pedlist}
}
}
\description{
Reconstructs the most likely pedigree from genotype data.
}
\details{
The parameter \code{extra} controls which of two algorithms are used to create the
pedigree list.

If \code{extra} is a nonnegative integer, it determines the number of extra
individuals allowed in the iterative pedigree construction. These extras
start off with undetermined sex, meaning that both males and females are
used. It should be noted that the final pedigrees may contain additional
extras, since missing parents are added at the end.

If \code{extra} is the word "parents", the algorithm is not iterative. It first
generates all directed acyclic graphs between the original individuals. Then
their parents are added and merged in all possible ways. This option has the
advantage of not requiring an explicit/ad hoc number of "extras", but works
best in smaller cases.
}
\examples{
#-----------------
# Example 1: Trio
#-----------------

# Built-in dataset `trioData`
x = singletons(1:3) |>
  setMarkers(alleleMatrix = trioData, locusAttributes ="snp12")

res = reconstruct(x, inferPO = TRUE, age = "1 > 2")

# Plot most likely pedigrees
plot(res, top = 6)


#--------------------
# Example 2: Siblings
#--------------------
library(forrel)

ids = c("s1", "s2")

# Create pedigree and simulate profiles with 20 STR markers
y = nuclearPed(children = ids) |>
  profileSim(markers = NorwegianFrequencies[1:20], ids = ids, seed = 123)

# Reconstruct allowing 2 extra individuals and any inbreeding
res2 = reconstruct(y, extra = 2, maxInb = 1)
plot(res2)

# With mutation modelling
y = setMutmod(y, model = "equal", rate = 0.01)
res3 = reconstruct(y, extra = 2, maxInb = 1)
plot(res3)

}
