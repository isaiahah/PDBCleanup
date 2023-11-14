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
#' @examples
#' # library(bio3d)
#' # Select high quality regions in an AlphaFold prediction.
#' 6ofsPredictedFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' 6ofsPredicted <- bio3d::read.pdb(6ofsPredictedFile)
#' # Use default threshold of >= 70 for predicted structures
#' 6ofsPredictedHighRes <- selectHighResolution(structure = 6ofs_predicted,
#'                                              predicted = TRUE)
#'
#' # Select high quality regions in an experimental prediction.
#' 6ofsExperimentalFile <- system.file("extdata", "6ofs_experimental.pdb",
#'                                     package = "PDBCleanup")
#' 6ofsExperimental <- bio3d::read.pdb(6ofsExperimentalFile)
#' # Use custom threshold of <= 80, or < 1.75 A mean deviation
#' 6ofsExperimentalHighRes <- selectHighResolution(structure = 6ofsExperimental,
#'                                                 predicted = FALSE,
#'                                                 threshold = 80)
#'
#' @export
#' @importFrom bio3d atom.select trim
selectHighResolution <- function(structure,
                                 predicted,
                                 threshold = NA) {
  atoms <- structure$atom
  if (predicted) { # Predicted structure, select high B-factor
    if (is.na(threshold)) { # Default threshold is 70
      threshold = 70
    } else {
      ;
    }

    # Select residues above B-factor threshold
    highResolutionIndices <- unique(atoms[atoms$"b" >= threshold, "resno"])
  } else { # Not predicted structure, select low B-factor
    if (is.na(threshold)) { # Default threshold is 100
      threshold = 100
    } else {
      ;
    }

    # Select regions below B-factor threshold
    highResolutionIndices <- unique(atoms[atoms$"b" <= threshold, "resno"])
  }

  # Trim the structure to the high resolution residues
  highResolutionSelection <- bio3d::atom.select(pdb = structure,
                                                resno = highResolutionIndices)
  highResStructure <- bio3d::trim(pdb = structure, highResolutionSelection)
  # Trim removes sequence regions, so copy full sequence
  highResStructure$seqres <- structure$seqres

  return(highResStructure)
}

# [END]
