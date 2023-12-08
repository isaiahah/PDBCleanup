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
#'     select the specified residues from.
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

# [END]
