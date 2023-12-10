# This example is adapted from:
# Silva, A. (2022) TestingPackage: An Example R Package For BCB410H.
# Unpublished. URL https://github.com/anjalisilva/TestingPackage.

library(shiny)
library(bio3d)

# Define UI
ui <- fluidPage(
  # Set title
  titlePanel("Preprocess and Cleanup PDB Protein Structure Files"),

  sidebarLayout(
    # Set sidebar, with input fields
    sidebarPanel(
      # Description and instructions
      tags$p("PDBCleanup is a package offering utilities to improve the
             quality of PDB protein structure files. Its features include
             selecting the high quality and ordered sections of structures, or
             aligning individual domains of a structure against a different
             reference structure. For more information on these features, their
             inputs, and how to use them in R, see the PDBCleanup R man page"),
      br(),
      tags$p("Instructions: Enter the PDB protein structure file, choose the
             desired operation tab, and enter the remaining parameters. Default
             values are shown below. Press 'Run' to perform the operation.
             The output structure will be shown on the main panel and available
             for download as a PDB file."),
      br(),
      fileInput(inputId = "PDB1",
                label = "Provide a protein structure for analysis. File should be a PDB.",
                accept = ".pdb"),

      # Input panels
      tabsetPanel(type="tabs",
                  tabPanel(
                    "Select High Resolution Regions",
                    radioButtons(inputId = "predicted",
                                 label = "Is the structure predicted? High resolution means high pLDDT in predicted structures and low B-factor in experimental structures.",
                                 choices = c("Predicted",
                                             "Experimental")),
                    textInput(inputId = "threshold",
                              label = "Provide a threshold to determine high resolution region. Leave NA to use default of 70 pLDDT or 80 B-factor.",
                              "NA"),
                    br(),
                    actionButton(inputId = "selectHighRes",
                                 label = "Run (select high resolution)")
                    ),
                  tabPanel(
                    "Select Ordered Regions",
                    textInput(inputId = "threshold",
                                 label = "Provide a FoldIndex threshold to determine ordered regions. FoldIndex ranges from -1 to 1, with values above 0 indicating ordered.",
                                 "0"),
                    textInput(inputId = "windowSize",
                              label = "Provide the window size used to calculate FoldIndex. Must be odd.",
                              "21"),
                    br(),
                    actionButton(inputId = "selectOrdered",
                                 label = "Run (select ordered)")
                    ),
                  tabPanel(
                    "Align Domains",
                    fileInput(inputId = "PDB2",
                              label = "Provide a reference protein structure to align against. File should be a PDB.",
                              accept = ".pdb"),
                    textInput(inputId = "domains1",
                              label = "Provide the domains in the first structure as a semicolon-separated list of ranges without spaces, eg 1-5;10-15",
                              "36-486;511-930"),
                    textInput(inputId = "domains2",
                              label = "Provide the domains in the second structure in the same format.",
                              "36-486;511-930"),
                    br(),
                    actionButton(inputId = "alignDomains",
                                 label = "Run (align domains)")
                  )
                )
    ),

    # Set main panel, with output
    mainPanel(
      tabsetPanel(
        tabPanel("View and Download Processed Protein",
                 h3("Choose the desired function, input data, and click 'Run'"),
                 h3("3D model of processed protein structure:"),
                 br(),
                 r3dmolOutput(outputId = "model"),
                 br(),
                 downloadButton(outputId = "downloadPDB", label = "Download .pdb")
        ),
        tabPanel("View Quality Along Protein",
                 h3("Choose the desired function, input data, and click 'Run'"),
                 h3("Quality Score along the processed protein structure:"),
                 br(),
                 plotOutput(outputId = "qualityPlot")
        )
      )
    ),

  )
)


# Define server
server <- function(input, output) {
  functionOutputs <- reactiveValues(structure = NULL)

  # Reactive for selectHighRes
  observeEvent(eventExpr = input$selectHighRes, {
    if (!is.null(input$PDB1$datapath)) {
      # Load structure
      structure <- bio3d::read.pdb(input$PDB1$datapath)
      # Interpret predicted variable
      if (input$predicted == "Predicted") {
        predicted <- TRUE
      } else {
        predicted <- FALSE
      }
      # Save output values
      functionOutputs$structure <- PDBCleanup::selectHighResolution(structure,
                                                                    predicted = predicted,
                                                                    threshold = as.numeric(input$threshold))
    }
  })

  # Reactive for selectOrdered
  observeEvent(eventExpr = input$selectOrdered, {
    if (!is.null(input$PDB1$datapath)) {
      # Load structure
      structure <- bio3d::read.pdb(input$PDB1$datapath)
      # Save output values
     functionOutputs$structure <- PDBCleanup::selectHighResolution(structure,
                                                                   threshold = as.numeric(input$threshold),
                                                                   windowSize = as.numeric(input$windowSize))
    }
  })

  # Reactive for alignDomains
  observeEvent(eventExpr = input$alignDomains, {
    if (!is.null(input$PDB1$datapath) && !is.null(input$PDB2$datapath)) {
      # Load structures
      structure1 <- bio3d::read.pdb(input$PDB1$datapath)
      structure2 <- bio3d::read.pdb(input$PDB2$datapath)
      # Load domains
      ranges1 <- strsplit(strsplit(input$domains1, ";")[[1]], "-")
      domains1 <- vector(mode = "list", length = length(ranges1)) # Preallocate
      for (i in seq_along(ranges1)) {
        domains1[[i]] = ranges1[[i]][1]:ranges1[[i]][2]
      }
      ranges2 <- strsplit(strsplit(input$domains2, ";")[[1]], "-")
      domains2 <- vector(mode = "list", length = length(ranges2)) # Preallocate
      for (i in seq_along(ranges1)) {
        domains2[[i]] = ranges2[[i]][1]:ranges1[[i]][2]
      }
      allDomains <- list(domains1, domains2)
      # Save output values
      functionOutputs$structure <- PDBCleanup::alignDomains(fixed = structure2,
                                                            mobile = structure1,
                                                            domains = allDomains)
    }
  })

  # Define outputs
  output$model <- renderR3dmol({
    if (!is.null(functionOutputs$structure)) {
      PDBCleanup::viewStructure(functionOutputs$structure)
    }
  })
  output$downloadPDB <- downloadHandler(
    filename = function() {
      return("PDBCleanup_output.pdb")
    },
    content = function(file) {
      bio3d::write.pdb(functionOutputs$structure, file)
    }
  )
  output$qualityPlot <- renderPlot(
    if (!is.null(functionOutputs$structure)) {
      PDBCleanup::plotProteinQuality(functionOutputs$structure,
                                     "Protein Quality Plot")
    }
  )
}


# Create the Shiny app
shiny::shinyApp(ui, server)

# [END]
