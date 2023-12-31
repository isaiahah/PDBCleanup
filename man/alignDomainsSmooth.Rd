% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alignDomains.R
\name{alignDomainsSmooth}
\alias{alignDomainsSmooth}
\title{Independently align multiple structure domains smoothly.}
\usage{
alignDomainsSmooth(fixed, mobile, domains, smoothing = 2)
}
\arguments{
\item{fixed}{A protein structure of class "pdb" (from bio3d) to
act as the fixed template.}

\item{mobile}{A protein structure of class "pdb" (from bio3d) which is
smoothly aligned to the fixed structure.}

\item{domains}{The locations of domains to align. Provided as a list of
lists, where each outer item describes one domain and each inner item is
a two-item list containing an integer vector of the domain's indices in
the fixed protein and an integer vector of the domain's indices in the
mobile protein. Each pair must contain the same number of residues across
structures.}

\item{smoothing}{An optional numeric value used in smoothing. A lower value
causes less smoothing, a larger value causes greater smoothing.
The default value is 2.}
}
\value{
A protein structure of class "pdb" with the residues from the mobile
structure smoothly aligned to the fixed structure on each domain.
}
\description{
A function to align a mobile protein structure against a fixed template
protein structure by independently aligning each pair of provided domains.
The alignment is smoothened by performing a rigid alignment then calculating
the final displacement of each atom as a weighted average of the displacement
of nearby atoms. For atom y, the weight assigned to atom x with displacement
d_x is w = exp(- ||y - x||^2 / u^2) where u is the smoothing factor which
the user can optionally set (default 2 A).
The locations of atoms in the mobile protein structure which are not part of
a domain do not move.
}
\details{
WARNING: this function is slow (~10 seconds on 7000 atoms).
}
\examples{
\dontrun{
# library(bio3d)
# Open the predicted structure 6ofs_predicted.pdb provided by the package.
# This is the predicted structure of an E coli zinc protease.
predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                 package = "PDBCleanup")
predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
# Open the experimental structure 6ofs_experimental.pdb provided by the
# package. This is the experimentally determined structure of the same E coli
# zinc protease.
experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
                                 package = "PDBCleanup")
experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)

# Define the domains. Here, they are the same because we use alternate
# structures from one species rather than structures from different species
domains6ofs <- list(list(36:486, 36:486), list(511:930, 511:930))

# The domains are placed differently between structures. This alignment
# smoothly moves the predicted domains into the places predicted by the
# experimental domains, which is beneficial if the predicted structure
# has unlikely atom placements (eg bad bond distances)
alignedPredicted6ofs <- alignDomainsSmooth(experimental6ofs, predicted6ofs,
                                           domains6ofs)
}

}
\references{
Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.

Terwilliger, T.C., Liebschner, D., Croll, T.I. et al.
AlphaFold predictions are valuable hypotheses and accelerate but do not
replace experimental structure determination. Nat Methods (2023).
https://doi.org./10.1038/s41592-023-02087-4
}
