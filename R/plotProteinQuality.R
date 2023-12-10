#' Plot the quality of a protein structure along its sequence.
#'
#' A function to visualize the data in the B-factor column of a protein
#' structure as a bar plot along the sequence. The user specifies the plot title
#' and optionally the x-axsis and y-axis title.
#'
#' @param structure A bio3d pdb protein structure to plot quality values from.
#' @param title A character string containing the title for the generated plot.
#' @param xtitle An optional character string containing the x-axis title.
#'    Default value is "Sequence Position"
#' @param ytitle An optional character string containing the y-axis title.
#'    Default value is "Predicted Structure Quality"
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
#' # Open the predicted structure 6ofs_predicted.pdb provided by the package.
#' # This is the predicted structure of an E coli zinc protease.
#' predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
#'                                  package = "PDBCleanup")
#' predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
#' plotProteinQuality(predicted6ofs,
#'                    title = "Quality of Predicted 6ofs structure")
#' @export
#' @import ggplot2
plotProteinQuality <- function(structure, title,
                               xtitle = "Sequence Position",
                               ytitle = "Predicted Sequence Quality") {
  # Check type of structure argument
  if (!(inherits(structure, "pdb"))) {
    stop("Provided structure must be a pdb object from bio3d")
  } else {
    ; # Valid type
  }

  # Check type of title
  if (!(is.character(title)) || (length(title) != 1)) {
    stop("Provided title must be a character string with one item")
  } else {
    ; # Valid type and length
  }

  # Check type of xtitle
  if (!(is.character(xtitle)) || (length(xtitle) != 1)) {
    stop("Provided xtitle must be a character string with one item")
  } else {
    ; # Valid type and length
  }

  # Check type of ytitle
  if (!(is.character(ytitle)) || (length(ytitle) != 1)) {
    stop("Provided ytitle must be a character string with one item")
  } else {
    ; # Valid type and length
  }

  proteinQuality <- unique(structure$atom[ , c("resno", "b")])
  proteinQualityPlot <- ggplot2::ggplot(data = proteinQuality) +
    ggplot2::aes(x = resno, y = b) +
    ggplot2::geom_col() +
    ggplot2::labs(title = title, x = xtitle, y = ytitle) +
    ggplot2::theme_minimal()
  return(proteinQualityPlot)
}

# [END]
