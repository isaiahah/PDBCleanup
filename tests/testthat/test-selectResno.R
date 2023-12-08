library(PDBCleanup)

test_that("Error on invalid structure object", {
  badStructure <- NA
  expect_error(PDBCleanup::selectResno(badStructure, c()))
})

test_that("Test on valid set of residue indices", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  domainsResno <- c(36:486, 511:930)

  highResolution <- PDBCleanup::selectResno(structure = predicted6ofs,
                                            resno = domainsResno)
  # Predicted structures have no HETATOM, so do not need to filter to ATOM
  expect_equal(nrow(highResolution$atom),
               sum(predicted6ofs$atom$"resno" %in% domainsResno))
  expect_equal(all(highResolution$atom$"resno" %in% domainsResno), TRUE)
})

test_that("Test on fully invalid set of residue indices", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  domainsResno <- c(10000:15000)

  highResolution <- PDBCleanup::selectResno(structure = predicted6ofs,
                                            resno = domainsResno)
  # Predicted structures have no HETATOM, so do not need to filter to ATOM
  expect_equal(nrow(highResolution$atom), 0)
})

test_that("Test on partially invalid set of residue indices", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  domainsResno <- c(511:930, 10000:15000)

  highResolution <- PDBCleanup::selectResno(structure = predicted6ofs,
                                            resno = domainsResno)
  # Predicted structures have no HETATOM, so do not need to filter to ATOM
  # Predicted structures have no HETATOM, so do not need to filter to ATOM
  expect_equal(nrow(highResolution$atom),
               sum(predicted6ofs$atom$"resno" %in% c(511:930)))
  expect_equal(all(highResolution$atom$"resno" %in% c(511:930)), TRUE)
})
