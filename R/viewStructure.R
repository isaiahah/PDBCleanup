#' Visualize the 3D structure of a protein.
#'
#' A function to visualize the provided structure with a 3D cartoon model.
#' This function opens a browser window with the interactive 3D model.
#'
#' @param structure A protein structure of class "pdb" (from Bio3d) to
#'    visualize as a 3D model.
#'
#' @returns The interactive viewer showing the 3D structure model. Shown in the
#'    viewer in RStudio, otherwise shown in a browser window.
#'
#' @references
#' Su W, Johnston B (2021). _r3dmol: Create Interactive 3D Visualizations of
#'     Molecular Data_. R package version 0.1.2,
#'     <https://CRAN.R-project.org/package=r3dmol>.
#'
#' @example
#' Example
#'
#' @export
#' @import r3dmol
viewStructure <- function(structure) {
  # Create the 3D model from the structure
  model3d <- m_bio3d(structure)
  # Visualize the 3D model
  structureVisual <- r3dmol() %>%
    m_add_model(data = model3d) %>%
    m_zoom_to() %>%
    m_set_style(style = m_style_cartoon(color = "#00cccc"))

  return(structureVisual)
}

# [END]
