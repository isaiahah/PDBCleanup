library(PDBCleanup)
library(idpr)

test_that("Error on invalid types of inputs", {
  badStructure <- NA
  expect_error(PDBCleanup::selectOrdered(badStructure))

  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)

  expect_error(PDBCleanup::selectOrdered(predicted6ofs, threshold = "30"))
  expect_error(PDBCleanup::selectOrdered(predicted6ofs, windowSize = "51"))
  expect_error(PDBCleanup::selectOrdered(predicted6ofs, windowSize = 30))
})

test_that("Error on even window size", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)

  expect_error(PDBCleanup::selectOrdered(predicted6ofs, windowSize = 30))
})

# Defined to avoid repeated code in following tests
sequence3to1 <- function(sequence3) {
  AAcode3to1 <- c("ALA" = "A", "ARG" = "R", "ASN" = "N", "ASP" = "D",
                  "CYS" = "C", "GLU" = "E", "GLN" = "Q", "GLY" = "G",
                  "HIS" = "H", "ILE" = "I", "LEU" = "L", "LYS" = "K",
                  "MET" = "M", "PHE" = "F", "PRO" = "P", "SER" = "S",
                  "THR" = "T", "TRP" = "W", "TYR" = "Y", "VAL" = "V")
  sequence1 <- rep("", length(sequence3))
  for (i in seq_along(sequence3)) {
    sequence1[i] <- AAcode3to1[sequence3[i]]
  }
  return(sequence1)
}

test_that("Test with default values", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)

  orderedStructure <- PDBCleanup::selectOrdered(predicted6ofs)
  disorderScores <- idpr::foldIndexR(sequence3to1(predicted6ofs$seqres),
                                     window = 21, plotResults = FALSE)
  expect_equal(length(unique(orderedStructure$atom$resno)),
               sum(disorderScores$foldIndex > 0))
  expect_equal(all(disorderScores[disorderScores$Position %in% orderedStructure$atom$resno, "foldIndex"] > 0),
               TRUE)
})

test_that("Test with custom threshold", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)

  threshold <- 0.1
  orderedStructure <- PDBCleanup::selectOrdered(predicted6ofs,
                                                threshold = threshold)
  disorderScores <- idpr::foldIndexR(sequence3to1(predicted6ofs$seqres),
                                     window = 21, plotResults = FALSE)
  expect_equal(length(unique(orderedStructure$atom$resno)),
               sum(disorderScores$foldIndex > threshold))
  expect_equal(all(disorderScores[disorderScores$Position %in% orderedStructure$atom$resno, "foldIndex"] > threshold),
               TRUE)
})

test_that("Test with custom window size", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)

  windowSize <- 51
  orderedStructure <- PDBCleanup::selectOrdered(predicted6ofs,
                                                windowSize = windowSize)
  disorderScores <- idpr::foldIndexR(sequence3to1(predicted6ofs$seqres),
                                     window = windowSize, plotResults = FALSE)
  expect_equal(length(unique(orderedStructure$atom$resno)),
               sum(disorderScores$foldIndex > 0))
  expect_equal(all(disorderScores[disorderScores$Position %in% orderedStructure$atom$resno, "foldIndex"] > 0),
               TRUE)
})
