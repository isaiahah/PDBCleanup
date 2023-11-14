#' Select High Resolution Regions of a Protein Structure
#'
#' A function to select high resolution regions in a protein structure.
#' High resolution regions are selected using B-factor.
#' The user provides whether the structure is an AlphaFold predicted structure,
#' meaning high B-factor indicates high resolution, or an experimental
#' structure, meaning low B-factor indicates high resolution. They can
#' optionally provide a threshold value, otherwise the function uses a default
#' threshold.
#'
#' @param structure A protein structure of class "pdb" (Bio3d format) to extract
#'     high resolution regions from.
#' @param predicted A boolean indicating whether the provided structure is a
#'     predicted structure. If TRUE, select regions with high B-factor as
#'     high resolution. If FALSE, select regions with low B-factor as high
#'     resolution.
#'     Defaults to TRUE.
#' @param threshold An optional parameter specifying the threshold
#'    value to identify high resolution regions. If predicted is TRUE, defaults
#'    to 70 per AlphaFold guidelines. If predicted is FALSE, defaults to
#'    60, which corresponds to 1.5 A mean deviation.
#'
#' @return A protein structure of class "pdb" containing only high resolution
#'    regions which satisfy the specified B-factor requirement.
#'
#' @references
#' Grant, B.J. et al. (2006) Bioinformatics 22, 2695--2696.
#'
#' Jumper, J., Evans, R., Pritzel, A. et al.
#'     Highly accurate protein structure prediction with AlphaFold.
#'     Nature 596, 583–589 (2021).
#'     https://doi.org/10.1038/s41586-021-03819-2
#'
#' Zhoutong Sun, Qian Liu, Ge Qu, Yan Feng, and Manfred T. Reetz.
#'     Utility of B-Factors in Protein Science: Interpreting Rigidity,
#'     Flexibility, and Internal Motion and Engineering Thermostability.
#'     Chemical Reviews 2019 119 (3), 1626-1665.
#'     https://doi.org/10.1021/acs.chemrev.8b00290
#'
#' @example
#' Example
#'
#' @export
#' @importFrom bio3d atom.select trim
selectHighResolution <- function(structure,
                                 predicted = TRUE,
                                 threshold = NA) {
  # Create a copy of the structure to select high resolution regions
  selectionAll <- bio3d::atom.select(pdb = structure,
                                     resno = unique(structure$atom$resno))
  highResStructure <- bio3d::trim(pdb = structure, selectionAll)
  atoms <- highResStructure$atom

  if (predicted) { # Predicted structure, select high B-factor
    if (is.na(threshold)) { # Default threshold is 70
      threshold = 70
    } else {
      ;
    }

    # Select regions above B-factor threshold
    highResAtoms <- atoms[atoms$"b" >= threshold, ]
    highResStructure$atom <- highResAtoms
  } else { # Not predicted structure, select low B-factor
    if (is.na(threshold)) { # Default threshold is 60
      threshold = 60
    } else {
      ;
    }

    # Select regions below B-factor threshold
    highResAtoms <- atoms[atoms$"b" <= threshold, ]
  }

  # Update the structure's atom attribute
  highResStructure$atom <- highResAtoms
  # Update the structure's xyz attribute
  highResxyz <- as.matrix(highResAtoms[c("x", "y", "z")])
  highResxyz <- matrix(highResxyz, nrow = 1, byrow = TRUE) # Reshape
  class(highResxyz) <- class(structure$xyz)
  highResStructure$xyz <- highResxyz
  # Update the structure's Calpha attribute
  highResCalpha <- (highResAtoms$elety == "CA")

  return(highResStructure)
}

# [END]
