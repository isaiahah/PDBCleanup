#' Independently align multiple structure domains rigidly.
#'
#' A function to align a mobile protein structure against a fixed template
#' protein structure by independently aligning each pair of provided domains.
#' The alignment is rigid, eg the domain is moved as one unit. Further,
#' the locations of atoms in the mobile protein structure which are not part of
#' a domain do not move.
#'
#' @param fixed A protein structure of class "pdb" (from bio3d) to
#'    act as the fixed template.
#' @param mobile A protein structure of class "pdb" (from bio3d) which is
#'    aligned to the fixed structure.
#' @param domains The locations of domains to align. Provided as a list of
#'    lists, where each outer item describes one domain and each inner item is
#'    a two-item list containing an integer vector of the domain's indices in
#'    the fixed protein and an integer vector of the domain's indices in the
#'    mobile protein. Each pair must contain the same number of residues across
#'    structures.
#'
#' @returns A protein structure of class "pdb" with the residues from the mobile
#'    structure rigidly aligned to the fixed structure on each domain.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' @examples
#' # library(bio3d)
#' # Open the predicted structure 6ofs_predicted.pdb provided by the package.
#' # This is the predicted structure of an E coli zinc protease.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' # Open the experimental structure 6ofs_experimental.pdb provided by the
#' # package. This is the experimentally determined structure of the same E coli
#' # zinc protease.
#' experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                  package = "PDBCleanup")
#' experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
#'
#' # Define the domains. Here, they are the same because we use alternate
#' # structures from one species rather than structures from different species
#' domains6ofs <- list(list(36:486, 36:486), list(511:930, 511:930))
#'
#' # The domains are placed differently between structures. This alignment moves
#' # the predicted domains into the places predicted by the experimental
#' # domains, which is beneficial if the experimental domain locations is better.
#' alignedPredicted6ofs <- alignDomainsRigid(experimental6ofs, predicted6ofs,
#'                                           domains6ofs)
#'
#' @export
#' @import bio3d
alignDomainsRigid <- function(fixed, mobile, domains) {
  # Check input types
  if (!(inherits(fixed, "pdb"))) {
    stop("fixed structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }
  if (!(inherits(mobile, "pdb"))) {
    stop("mobile structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }
  if (!("list" %in% class(domains))) {
    stop("Domains should be a list: See alignDomains documentation")
  } else {
    ; # Valid type
  }
  if (length(domains) == 0) {
    stop("Must provide at least 1 domain to align")
  } else {
    ; # Valid length
  }
  for (i in seq_along(domains)) {
    if (!("list" %in% class(domains[[i]])) || length(domains[[i]]) != 2) {
      stop("All domains must be length 2 lists: See alignDomains documentation")
    } else {
      ; # Valid domains
    }
  }
  for (i in seq_along(domains)) {
    if (length(domains[[i]][1]) != length(domains[[i]][2])) {
      stop("Each domain pair should have an equal length in both proteins")
    } else {
      ; # Valid domain lengths
    }
  }

  # Copy the mobile protein structure by trimming to all atoms
  allMobileSelection <- bio3d::atom.select(pdb = mobile,
                                           resno = unique(mobile$atoms$resno))
  aligned <- bio3d::trim(pdb = mobile, allMobileSelection)

  # Loop over domains, performing alignment and storing aligned coordinates
  for (i in seq_along(domains)) {
    # Store the boolean vector of which atoms are in the domain
    fixedResidues <- domains[[i]][[1]]
    fixedInDomain <- (fixed$atom$resno %in% fixedResidues)
    mobileResidues <- domains[[i]][[2]]
    mobileInDomain <- (mobile$atom$resno %in% mobileResidues)
    mobileDomainRows <- which(mobileInDomain)

    # Extract the xyz coordinates of atoms in the domain
    fixedDomainAtom <- fixed$atom[fixedInDomain, ]
    fixedDomainxyz <- as.matrix(fixed$atom[fixedInDomain, c("x", "y", "z")])
    fixedDomainxyz <- matrix(t(fixedDomainxyz), nrow = 1) # Reshape
    mobileDomainxyz <- as.matrix(mobile$atom[mobileInDomain, c("x", "y", "z")])
    mobileDomainxyz <- matrix(t(mobileDomainxyz), nrow = 1) # Reshape

    # Fit along alpha carbons by finding CA atom indices then converting to
    # xyz indices (atom i corresponds to xyz indices 3i-2, 3i-1, and 3i)
    fixedDomainCA <- which(fixed$atom[fixedInDomain, "elety"] == "CA")
    fixedDomainFit <- c(3 * fixedDomainCA - 2, 3 * fixedDomainCA - 1,
                        3 * fixedDomainCA)
    mobileDomainCA <- which(mobile$atom[mobileInDomain, "elety"] == "CA")
    mobileDomainFit <- c(3 * mobileDomainCA - 2, 3 * mobileDomainCA - 1,
                        3 * mobileDomainCA)

    # Align domains
    alignedDomainxyz <- bio3d::fit.xyz(fixedDomainxyz, mobileDomainxyz,
                                       fixedDomainFit, mobileDomainFit)

    # Save the aligned domain in atoms and xyz attributes
    for (i in seq_along(mobileDomainRows)) {
      row <- mobileDomainRows[i]
      # Set atom coordinates
      aligned$atom[row, "x"] <- alignedDomainxyz[3 * (i - 1) + 1]
      aligned$atom[row, "y"] <- alignedDomainxyz[3 * (i - 1) + 2]
      aligned$atom[row, "z"] <- alignedDomainxyz[3 * (i - 1) + 3]
    }
  }
  # Set final xyz
  alignedxyz <- as.matrix(aligned$atom[c("x", "y", "z")])
  alignedxyz <- matrix(t(alignedxyz), nrow = 1) # Reshape
  class(alignedxyz) <- class(aligned$xyz)
  aligned$xyz <- alignedxyz

  return(aligned)
}


#' Independently align multiple structure domains smoothly.
#'
#' A function to align a mobile protein structure against a fixed template
#' protein structure by independently aligning each pair of provided domains.
#' The alignment is smoothened by performing a rigid alignment then calculating
#' the final displacement of each atom as a weighted average of the displacement
#' of nearby atoms. For atom y, the weight assigned to atom x with displacement
#' d_x is w = exp(- ||y - x||^2 / u^2) where u is the smoothing factor which
#' the user can optionally set (default 2 A).
#' The locations of atoms in the mobile protein structure which are not part of
#' a domain do not move.
#'
#' WARNING: this function is slow (~10 seconds on 7000 atoms).
#'
#' @param fixed A protein structure of class "pdb" (from bio3d) to
#'    act as the fixed template.
#' @param mobile A protein structure of class "pdb" (from bio3d) which is
#'    smoothly aligned to the fixed structure.
#' @param domains The locations of domains to align. Provided as a list of
#'    lists, where each outer item describes one domain and each inner item is
#'    a two-item list containing an integer vector of the domain's indices in
#'    the fixed protein and an integer vector of the domain's indices in the
#'    mobile protein. Each pair must contain the same number of residues across
#'    structures.
#' @param smoothing An optional numeric value used in smoothing. A lower value
#'    causes less smoothing, a larger value causes greater smoothing.
#'    The default value is 2.
#'
#' @returns A protein structure of class "pdb" with the residues from the mobile
#'    structure smoothly aligned to the fixed structure on each domain.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' Terwilliger, T.C., Liebschner, D., Croll, T.I. et al.
#' AlphaFold predictions are valuable hypotheses and accelerate but do not
#' replace experimental structure determination. Nat Methods (2023).
#' https://doi.org./10.1038/s41592-023-02087-4
#'
#' @examples
#' \dontrun{
#' # library(bio3d)
#' # Open the predicted structure 6ofs_predicted.pdb provided by the package.
#' # This is the predicted structure of an E coli zinc protease.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' # Open the experimental structure 6ofs_experimental.pdb provided by the
#' # package. This is the experimentally determined structure of the same E coli
#' # zinc protease.
#' experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                  package = "PDBCleanup")
#' experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
#'
#' # Define the domains. Here, they are the same because we use alternate
#' # structures from one species rather than structures from different species
#' domains6ofs <- list(list(36:486, 36:486), list(511:930, 511:930))
#'
#' # The domains are placed differently between structures. This alignment
#' # smoothly moves the predicted domains into the places predicted by the
#' # experimental domains, which is beneficial if the predicted structure
#' # has unlikely atom placements (eg bad bond distances)
#' alignedPredicted6ofs <- alignDomainsSmooth(experimental6ofs, predicted6ofs,
#'                                            domains6ofs)
#' }
#'
#' @export
#' @import bio3d
alignDomainsSmooth <- function(fixed, mobile, domains, smoothing = 2.0) {
  # Check fixed type
  if (!(inherits(fixed, "pdb"))) {
    stop("fixed structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }
  # Check mobile type
  if (!(inherits(mobile, "pdb"))) {
    stop("mobile structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }
  # Check domain type and values
  if (!("list" %in% class(domains))) {
    stop("Domains should be a list: See alignDomains documentation")
  } else {
    ; # Valid type
  }
  if (length(domains) == 0) {
    stop("Must provide at least 1 domain to align")
  } else {
    ; # Valid length
  }
  for (i in seq_along(domains)) {
    if (!("list" %in% class(domains[[i]])) || length(domains[[i]]) != 2) {
      stop("All domains must be length 2 lists: See alignDomains documentation")
    } else {
      ; # Valid domains
    }
  }
  for (i in seq_along(domains)) {
    if (length(domains[[i]][1]) != length(domains[[i]][2])) {
      stop("Each domain pair should have an equal length in both proteins")
    } else {
      ; # Valid domain lengths
    }
  }
  # Check type of smoothing argument
  if (!(is.numeric(smoothing)) || (length(smoothing) != 1)) {
    stop("smoothing must be a single numeric value")
  } else {
    ; # Valid type and length
  }

  # Copy the mobile protein structure by trimming to all atoms
  allMobileSelection <- bio3d::atom.select(pdb = mobile,
                                           resno = unique(mobile$atoms$resno))
  aligned <- bio3d::trim(pdb = mobile, allMobileSelection)

  # Loop over domains, performing rigid, alignment, smoothing, and
  # storing aligned coordinates
  for (i in seq_along(domains)) {
    # Store the boolean vector of which atoms are in the domain
    fixedResidues <- domains[[i]][[1]]
    fixedInDomain <- (fixed$atom$resno %in% fixedResidues)
    mobileResidues <- domains[[i]][[2]]
    mobileInDomain <- (mobile$atom$resno %in% mobileResidues)
    mobileDomainRows <- which(mobileInDomain)

    # Extract the xyz coordinates of atoms in the domain
    fixedDomainAtom <- fixed$atom[fixedInDomain, ]
    fixedDomainxyz <- as.matrix(fixed$atom[fixedInDomain, c("x", "y", "z")])
    fixedDomainxyz <- matrix(t(fixedDomainxyz), nrow = 1) # Reshape
    mobileDomainxyz <- as.matrix(mobile$atom[mobileInDomain, c("x", "y", "z")])
    mobileDomainxyz <- matrix(t(mobileDomainxyz), nrow = 1) # Reshape

    # Fit along alpha carbons by finding CA atom indices then converting to
    # xyz indices (atom i corresponds to xyz indices 3i-2, 3i-1, and 3i)
    fixedDomainCA <- which(fixed$atom[fixedInDomain, "elety"] == "CA")
    fixedDomainFit <- c(3 * fixedDomainCA - 2, 3 * fixedDomainCA - 1,
                        3 * fixedDomainCA)
    mobileDomainCA <- which(mobile$atom[mobileInDomain, "elety"] == "CA")
    mobileDomainFit <- c(3 * mobileDomainCA - 2, 3 * mobileDomainCA - 1,
                         3 * mobileDomainCA)

    # Rigidly align domains
    alignedDomainxyz <- bio3d::fit.xyz(fixedDomainxyz, mobileDomainxyz,
                                       fixedDomainFit, mobileDomainFit)

    # Calculate displacements
    displacements <- alignedDomainxyz - mobileDomainxyz
    # Compute square distances
    distances <- matrix(NA, nrow = length(mobileDomainRows),
                        ncol = length(mobileDomainRows))
    for (i in seq_along(mobileDomainRows)) {
      for (j in seq_along(mobileDomainRows)) {
        if (i > j) { # Avoid duplicate calls of (i, j) and (j, i)
          xDist <- (mobileDomainxyz[3 * (i - 1) + 1] -
                      mobileDomainxyz[3 * (j - 1) + 1]) ** 2
          yDist <- (mobileDomainxyz[3 * (i - 1) + 2] -
                      mobileDomainxyz[3 * (j - 1) + 2]) ** 2
          zDist <- (mobileDomainxyz[3 * (i - 1) + 3] -
                      mobileDomainxyz[3 * (j - 1) + 3]) ** 2
          dist <- (xDist + yDist + zDist)
          distances[i, j] <- dist
          distances[j, i] <- dist
        } else if (i == j) {
          distances[i, j] <- 0
        } else {
          ; # Handled in i > j
        }
      }
    }
    # Perform smoothing and save to output structure
    for (i in seq_along(mobileDomainRows)) {
      # Calculate weights: w = e^(-d^2/u^2) where d^2 is square distance and
      # u is the smoothing parameter
      weights <- 2.7182 ** (-distances[i,] / smoothing ** 2)
      weights <- weights / sum(weights)
      weights <- rep(weights, each = 3)
      weightedDisplacements <- matrix(weights * displacements, nrow = 3)
      row <- mobileDomainRows[i]
      # Update atom coordinates
      aligned$atom[row, "x"] <- (aligned$atom[row, "x"] +
                                   sum(weightedDisplacements[1,], na.rm = TRUE))
      aligned$atom[row, "y"] <- (aligned$atom[row, "y"] +
                                   sum(weightedDisplacements[2,], na.rm = TRUE))
      aligned$atom[row, "z"] <- (aligned$atom[row, "z"] +
                                   sum(weightedDisplacements[3,], na.rm = TRUE))
    }
  }
  # Set final xyz
  alignedxyz <- as.matrix(aligned$atom[c("x", "y", "z")])
  alignedxyz <- matrix(t(alignedxyz), nrow = 1) # Reshape
  class(alignedxyz) <- class(aligned$xyz)
  aligned$xyz <- alignedxyz

  return(aligned)
}

# [END]
