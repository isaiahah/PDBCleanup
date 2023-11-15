#' Independently align multiple domains in a protein structure.
#'
#' A function to align a mobile protein structure against a fixed template
#' protein structure by independently aligning each pair of provided domains.
#' The locations of atoms in the mobile protein structure which are not part of
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
#'    mobile protein.
#'
#' @returns A protein structure of class "pdb" with the residues from the mobile
#'    structure aligned to the fixed structure on each domain.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' @examples
#' # library(bio3d)
#' # Load an AlphaFold prediction.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' # Load the experimental prediction.
#' experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                  package = "PDBCleanup")
#' experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
#'
#' # Define the domains. Here, they are the same because we use alternate
#' # structures from one species rather than structures from different species
#' domains6ofs <- list(list(36:486, 36:486), list(511:930, 511:930))
#'
#' alignedPredicted6ofs <- alignDomains(experimental6ofs, predicted6ofs,
#'                                      domains6ofs)
#'
#' @export
#' @import bio3d
alignDomains <- function(fixed, mobile, domains) {
  # Check input types
  if (!("pdb" %in% class(fixed))) {
    stop("Fixed structure must be a pdb object from bio3d")
  }
  if (!("pdb" %in% class(mobile))) {
    stop("Mobile structure must be a pdb object from bio3d")
  }
  if (!("list" %in% class(domains))) {
    stop("Domains should be a list: See alignDomains documentation")
  }
  if (length(domains) == 0) {
    stop("Must provide at least 1 domain to align")
  }
  for (i in seq_along(domains)) {
    if (!("list" %in% class(domains[[i]])) || length(domains[[i]]) != 2) {
      stop("All domains must be length 2 lists: See alignDomains documentation")
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
    fixedDomainxyz <- as.matrix(fixed$atom[fixedInDomain, c("x", "y", "z")])
    fixedDomainxyz <- matrix(t(fixedDomainxyz), nrow = 1) # Reshape
    mobileDomainxyz <- as.matrix(mobile$atom[mobileInDomain, c("x", "y", "z")])
    mobileDomainxyz <- matrix(t(mobileDomainxyz), nrow = 1) # Reshape

    # Align domains
    alignedDomainxyz <- bio3d::fit.xyz(fixedDomainxyz, mobileDomainxyz)

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

# [END]