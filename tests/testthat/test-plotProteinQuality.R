library(PDBCleanup)

# Minimal tests, as hard to automatically check graph.
# Recommended to check visually instead.

test_that("Fails on invalid input", {
  expect_error(PDBCleanup::plotProteinQuality(NA, NA, NA, NA))
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  expect_error(PDBCleanup::plotProteinQuality(NA, "Title", "xtitle",
                                              "ytitle"))
  expect_error(PDBCleanup::plotProteinQuality(predicted6ofs, NA, "xtitle",
                                              "ytitle"))
  expect_error(PDBCleanup::plotProteinQuality(predicted6ofs, "Title", NA,
                                              "ytitle"))
  expect_error(PDBCleanup::plotProteinQuality(predicted6ofs, "Title", "xtitle",
                                              NA))
})

test_that("Passes on valid input", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  qualityPlot <- PDBCleanup::plotProteinQuality(predicted6ofs, "6ofs Quality",
                                                "Residue Number", "pLLDT")
  expect_true(TRUE) # Assert code reached here without error
})
