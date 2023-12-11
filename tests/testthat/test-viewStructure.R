library(PDBCleanup)

# Minimal tests, as hard to automatically check structure representation.
# Recommended to check visually instead.

test_that("Fails on invalid structure", {
  expect_error(PDBCleanup::viewStructure(NA))
})

test_that("Passes on valid structure", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  structureModel <- PDBCleanup::viewStructure(predicted6ofs)
  expect_true(TRUE) # Assert code reached here without error
})
