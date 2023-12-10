#' Select the residues of a protein structure at the provided indices.
#'
#' A function to select the specified residues of a given protein structure and
#' return a modified structure with all other residues removed.
#' This is primarily a helper function for the functions selectHighResPredicted
#' and selectHighResExperimental. However, it is exported if users have more
#' general applications.
#'
#' @param structure A bio3d pdb protein structure which the function will
#'     select the specified residues from.
#' @param resno A vector of integers which are the residues to select.
#'
#' @return A bio3d pdb protein structure containing only the specified residues.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' @examples
#' # library(bio3d)
#' # Open a PDB structure. The package provides one at 6ofs_predicted.pdb, which
#' # is the experimentally determined structure of an E coli zinc protease.
#' experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                     package = "PDBCleanup")
#' experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
#' # Define the residues of interest. These are the two structured domains.
#' domainsResno6ofs <- c(36:486, 511:930)
#' # Select the domains from the full 6ofs structure
#' domains6ofs <- selectResno(structure = experimental6ofs,
#'                            resno = domainsResno6ofs)
#'
#' @export
#' @importFrom bio3d atom.select trim
selectResno <- function(structure, resno) {
  # Check type of structure argument
  if (!(inherits(structure, "pdb"))) {
    stop("Provided structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }

  # Check type of provided residues
  if (!(inherits(resno, "integer"))) {
    stop("Provided residue positions must be a vector of integers")
  } else {
    ; # Valid type
  }

  # Trim the structure to the high resolution residues
  resnoSelection <- bio3d::atom.select(pdb = structure,
                                       resno = resno, type="ATOM")
  selectedStructure <- bio3d::trim(pdb = structure, resnoSelection)
  # Trim removes sequence regions, so copy full sequence
  selectedStructure$seqres <- structure$seqres

  return(selectedStructure)
}


#' Select high resolution regions of a protein structure.
#'
#' A function to select high resolution regions in a protein structure.
#'
#' In an experimental structure, these regions are determined using the
#' B-factor, which measures the average uncertainty in atom location using the
#' formula B = 4*pi^2 u^2 / 3 where u is the uncertainty in angstroms. A low
#' B-factor indicates a higher quality region, and the function selects all
#' residues containing atoms with a B-factor below a threshold.
#'
#' In a predicted structure, these regions are determined using the pLDDT
#' stored in the B-factor column, which is AlphaFold's prediction for that
#' residue's score on the local distance difference test. A high pLDDT
#' indicates the predicted residue is likely close to the true position, while
#' a low pLDDT incidates less confidence in the predicted residue location.
#' The function selects all residues with pLDDT above a threshold.
#'
#' @param structure A bio3d pdb protein structure which the function will
#'     select the high quality residues from.
#' @param predicted A logical value indicating if the structure is predicted or
#'     experimental. Predicted structures are interpreted using B-factor (lower
#'     B-factor is high resolution) and experimental structures are interpreted
#'     using pLDDT (higher pLDDT is high resolution).
#' @param threshold An optional numeric value specifying the threshold for high
#'     resolution regions. In predicted structures, the default is 70 and all
#'     residues with higher pLDDT are high resolution. In experimental
#'     structures, the default is 100 (~2.75 A uncertainty) and all residues
#'     with lower B-factor are marked high resolution.
#'
#' @return A bio3d pdb protein structure containing only high-quality residues
#'     satisfying the specified threshold.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' Zhoutong Sun, Qian Liu, Ge Qu, Yan Feng, and Manfred T. Reetz.
#'     Utility of B-Factors in Protein Science: Interpreting Rigidity,
#'     Flexibility, and Internal Motion and Engineering Thermostability.
#'     Chemical Reviews 2019 119 (3), 1626-1665.
#'     https://doi.org/10.1021/acs.chemrev.8b00290
#'
#' @examples
#' # library(bio3d)
#' # Open the predicted structure 6ofs_predicted.pdb provided by the package.
#' # This is the predicted structure of an E coli zinc protease.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                     package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' # Use default threshold of pLDDT >= 70
#' predicted6ofsHighRes <- selectHighResolution(structure = predicted6ofs,
#'                                              predicted = TRUE)
#'
#' # Open the experimental structure 6ofs_experimental.pdb provided by the
#' # package. This is the experimentally determined structure of the same E coli
#' # zinc protease.
#' experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                     package = "PDBCleanup")
#' experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
#' # Use custom threshold of bfactor <= 80, or under ~2.5 A mean deviation
#' experimental6ofsHighRes <- selectHighResolution(structure = experimental6ofs,
#'                                                 predicted = FALSE,
#'                                                 threshold = 80)
#'
#' @export
selectHighResolution <- function(structure, predicted, threshold = NA) {
  # Check type of structure argument
  if (!(inherits(structure, "pdb"))) {
    stop("Provided structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }

  # Check type of predicted argument
  if (!(is.logical(predicted)) || (length(predicted) != 1)) {
    stop("Predicted must be a single logical value")
  } else {
    ; # Valid type and length
  }

  # Check type and value of threshold argument
  if (!(is.na(threshold))) {
    if (!(is.numeric(threshold)) || (length(threshold) != 1)) {
      stop("Threshold must be a single numeric value")
    } else {
      ; # Valid type and length
    }
    if (threshold <= 0) {
      warning("Provided threshold is negative")
    } else {
      ; # Valid threshold value
    }
  }

  atoms <- structure$atom
  # Determine the high resolution residues from input specification
  if (predicted == TRUE) { # Predicted structure, select high B-factor
    if (is.na(threshold)) { # Default threshold is 70
      threshold <- 70
    } else {
      ; # Threshold already set
    }

    # Select residues above B-factor threshold
    highResolutionResno <- unique(atoms[atoms$"b" >= threshold, "resno"])
  } else { # Experimental structure, select low B-factor
    if (is.na(threshold)) { # Default threshold is 100
      threshold <- 100
    } else {
      ; # Threshold already set
    }

    # Select regions below B-factor threshold
    highResolutionResno <- unique(atoms[atoms$"b" <= threshold, "resno"])
  }

  # Select high resolution structure from selected residue positions
  highResStructure <- selectResno(structure, highResolutionResno)
  return(highResStructure)
}


#' Select ordered regions of a protein structure.
#'
#' A function to select all intrinsically ordered residues of a protein
#' structure and remove all intrinsically disordered residues. Residues are
#' classified using the FoldIndex score (implemented in the idpr package rather
#' than re-implementing), which predicts disorder at each residue with the
#' formula D = 2.785 H - | R | - 1.151 where H is the scaled mean hydropathy
#' and R is the scaled mean charge within a sequence window around the residue.
#' It ranges from -1 to 1 where positive scores indicate ordered residues while
#' negative scores indicate disordered residues.
#' The user optionally provides a threshold (default 0) to denote order and
#' optionally provides a windowSize (default 21) for this calculation.
#'
#' The first and last windowSize // 2 residues are removed, as they lack a
#' sufficient window size to calculate FoldIndex score.
#'
#' @param structure A bio3d pdb protein structure which the function will
#'     select the ordered residues from.
#' @param threshold An optional numerical value specifying the threshold for
#'     ordered regions. Defaults to scores greater than 0 denoting order.
#' @param windowSize An optional odd integer value specifying the window size
#'     used when calculating disorder score. Defaults to 21.
#'
#' @return A bio3d pdb protein structure containing only ordered residues with
#'     disorder scores satisfying the specified threshold.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' McFadden, W. M., and Yanowitz, J. L. (2022). idpr: A package for
#' profiling and analyzing Intrinsically Disordered Proteins in R.
#' PLOS ONE, 17(4), e0266929. doi:10.1371/journal.pone.0266929.
#'
#' Prilusky, J., Felder, C. E., et al. (2005). FoldIndex: a simple tool to
#' predict whether a given protein sequence is intrinsically unfolded.
#' Bioinformatics, 21(16), 3435-3438.
#'
#' @examples
#' # Open the experimental structure 6ofs_experimental.pdb provided by the
#' # package. This is the experimentally determined structure of an E coli
#' # zinc protease.
#' experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                     package = "PDBCleanup")
#' experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
#' # Use custom threshold 0.1 and window size 51
#' experimental6ofsOrdered <- selectOrdered(experimental6ofs,
#'                                          threshold = 0.1,
#'                                          windowSize = 51)
#'
#' @importFrom idpr foldIndexR
#' @export
selectOrdered <- function(structure, threshold = 0, windowSize = 21) {
  # Check type of structure argument
  if (!(inherits(structure, "pdb"))) {
    stop("Provided structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }

  # Check type of threshold argument
  if (!(is.numeric(threshold)) || (length(threshold) != 1)) {
    stop("threshold must be a single numeric value")
  } else {
    ; # Valid type and length
  }

  # Check type of windowSize argument
  if (!(is.numeric(windowSize)) || (length(windowSize) != 1)) {
    stop("windowSize must be a single numeric value")
  } else {
    ; # Valid type and length
  }
  # Check windowSize is odd
  if (windowSize %% 2 != 1) {
    stop("windowSize must be odd")
  } else {
    ; # windowSize is odd
  }

  # Convert sequence stored in structure from 3 letter codes to 1 letter codes
  AAcode3to1 <- c("ALA" = "A", "ARG" = "R", "ASN" = "N", "ASP" = "D",
                 "CYS" = "C", "GLU" = "E", "GLN" = "Q", "GLY" = "G",
                 "HIS" = "H", "ILE" = "I", "LEU" = "L", "LYS" = "K",
                 "MET" = "M", "PHE" = "F", "PRO" = "P", "SER" = "S",
                 "THR" = "T", "TRP" = "W", "TYR" = "Y", "VAL" = "V")
  sequence <- rep("", length(structure$seqres))
  for (i in seq_along(structure$seqres)) {
    sequence[i] <- AAcode3to1[structure$seqres[i]]
  }

  # Compute disorder scores and select ordered regions
  disorderScores <- idpr::foldIndexR(sequence, window = windowSize,
                                     plotResults = FALSE)
  orderedResno <- disorderScores[disorderScores$foldIndex > threshold,
                                 "Position"]

  # Select the ordered structure from the ordered residue positions
  orderedStructure <- selectResno(structure, orderedResno)
  return(orderedStructure)
}

# [END]
