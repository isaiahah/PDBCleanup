#' Visualize the 3D structure of a protein.
#'
#' A function to visualize the provided structure with a 3D cartoon model.
#' This function opens a browser window with the interactive 3D model.
#'
#' @param structure A bio3d pdb protein structure to visualize as a 3D model.
#' @param color An optional character string specifying the color for the 3D
#'    model. Default value is "#00cccc".
#'
#' @returns The interactive viewer showing the 3D structure model. Shown in the
#'    viewer in RStudio, otherwise shown in a browser window.
#'
#' @references
#' Su W, Johnston B (2021). r3dmol: Create Interactive 3D Visualizations of
#'     Molecular Data. R package version 0.1.2,
#'     <https://CRAN.R-project.org/package=r3dmol>.
#'
#' @examples
#' # library(bio3d)
#' # Load and visualize an AlphaFold prediction.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' viewStructure(predicted6ofs)
#'
#' @export
#' @import r3dmol
viewStructure <- function(structure, color = "#00cccc") {
  # Check type of structure argument
  if (!(inherits(structure, "pdb"))) {
    stop("Provided structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }

  if (!(is.character(color)) || (length(color) != 1)) {
    stop("Provided color must be a character string with one item")
  } else {
    ; # Valid type and length
  }

  # Create the 3D model from the structure
  model3d <- r3dmol::m_bio3d(structure)
  # Visualize the 3D model
  structureVisual <- r3dmol::r3dmol() %>%
    r3dmol::m_add_model(data = model3d) %>%
    r3dmol::m_zoom_to() %>%
    r3dmol::m_set_style(style = r3dmol::m_style_cartoon(color = color))

  return(structureVisual)
}

# [END]
