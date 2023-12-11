# This example is adapted from:
# Silva, A. (2022) TestingPackage: An Example R Package For BCB410H.
# Unpublished. URL https://github.com/anjalisilva/TestingPackage.

library(shiny)
library(shinyalert)
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
      tags$p("Note: Smooth domain alignment is not supported in Shiny due to
             the function's long runtime."),
      br(),
      shinyalert::useShinyalert(force = TRUE),
      uiOutput("PDBdata1"),
      actionButton(inputId = "PDBDataDetails1", label = "Demo Dataset Details"),
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
                    uiOutput("PDBdata2"),
                    actionButton(inputId = "PDBDataDetails2", label = "Demo Dataset Details"),
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

  # Define outputs: structure model, downloadable PDB, and quality plot
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

  # Define PDB data download elements
  PDBurl1 <- a("Predicted 6ofs PDB",
               href="https://alphafold.ebi.ac.uk/files/AF-P31828-F1-model_v4.pdb")
  output$PDBdata1 <- renderUI({
    tagList("Download Demo Protein Structure:", PDBurl1)
  })
  PDBurl2 <- a("Experimental 6ofs PDB",
               href="https://files.rcsb.org/download/6OFS.pdb")
  output$PDBdata2 <- renderUI({
    tagList("Download Demo Protein Structure:", PDBurl2)
  })

  # Define popup with PDB data information
  observeEvent(input$PDBDataDetails1, {
    shinyalert(title = "Predicted Protein Structure",
               text = "A computationally predicted structure of E coli protein Pqql, a probable zinc protease. It was predicted by AlphaFold2 based on the protein sequence for UniProt entry P31828.
               Citation: Varadi M., Anyango S., Deshpande M., et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research, Volume 50, Issue D1, 7 January 2022, Pages D439–D444. https://doi.org/10.1093/nar/gkab1061",
               type = "info")
  })
  observeEvent(input$PDBDataDetails2, {
    shinyalert(title = "Experimental Protein Structure",
               text = "An experimental structure of E coli protein Pqql, a probable zinc protease. The protein has Uniprot ID P31828 and was uploaded to the PDB under ID 60fs.
               Citation: Grinter R, Leung PM, Wijeyewickrema LC, Littler D, Beckham S, Pike RN, et al. (2019) Protease-associated import systems are widespread in Gram-negative bacteria. PLoS Genet 15(10): e1008435. https://doi.org/10.1371/journal.pgen.1008435",
               type = "info")
  })
}


# Create the Shiny app
shiny::shinyApp(ui, server)

# [END]
