#' Plot the quality of a protein structure along its sequence.
#'
#' A function to visualize the B-factor of a protein structure as a bar plot
#' along the sequence. The user specifies whether the structure is predicted
#' to specify the y-axis label and gives a custom plot title.
#'
#' @param structure A protein structure of class "pdb" (from bio3d) to
#'    visualize B-values from.
#' @param title The title for the generated plot.
#' @param predicted Whether the protein structure is predicted, so the y axis
#'    label is "Predicted Structure Quality," or experimental, so the y axis
#'    label is "Experimental B-factor."
#'
#' @returns The histogram plot of the B-factor along the protein sequence.
#'
#' @references
#' H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag
#'     New York, 2016.
#'
#' @examples
#' # library(bio3d)
#' # library(ggplot2)
#' # Load and visualize the quality along an AlphaFold prediction.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' plotProteinQuality(predicted6ofs,
#'                    title = "Quality of Predicted 6ofs structure",
#'                    predicted = TRUE)
#' @export
#' @import ggplot2
plotProteinQuality <- function(structure, title, predicted) {
  # Check structure argument
  if (!("pdb" %in% class(structure))) {
    stop("Provided structure must be a pdb object from bio3d")
  }

  if (predicted) {
    ytitle <- "Predicted Structure Quality"
  } else {
    ytitle <- "Experimental B-factor"
  }

  proteinQuality <- unique(structure$atom[ , c("resno", "b")])
  proteinQualityPlot <- ggplot2::ggplot(data = proteinQuality) +
    ggplot2::aes(x = resno, y = b) +
    ggplot2::geom_col() +
    ggplot2::labs(title = title, x = "Sequence Position", y = ytitle) +
    ggplot2::theme_minimal()
  return(proteinQualityPlot)
}

# [END]
