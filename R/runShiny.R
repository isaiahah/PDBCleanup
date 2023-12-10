#' Launch Shiny App for PDBCleanup
#'
#' A function that launches the Shiny app for PDBCleanup.
#' This app provides a web interface for the package functions
#' selectHighResolution, selectOrdered, and alignDomains from a user-uploaded
#' file. The code has been placed in \code{./inst/shiny-scripts}.
#'
#' @return No return value. Instead, open a Shiny page.
#'
#' @examples
#' \dontrun{
#' PDBCleanup:runShiny()
#' }
#'
#' @references
#' Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, Xie Y, Allen J,
#' McPherson J, Dipert A, Borges B (2023).
#' shiny: Web Application Framework for R.
#' R package version 1.8.0, https://CRAN.R-project.org/package=shiny>
#'
#' Based on:
#' Silva, A. (2022) TestingPackage: An Example R Package For BCB410H.
#' Unpublished. URL https://github.com/anjalisilva/TestingPackage.
#'
#' @export
#' @importFrom shiny runApp

runShiny <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "PDBCleanup")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}

# [END]
