library(PDBCleanup)

test_that("Error on invalid structure object", {
  badStructure <- NA
  expect_error(PDBCleanup::selectHighResolution(badStructure,
                                                predicted = TRUE))
})

test_that("Test on predicted structure with default threshold", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  highResolution <- PDBCleanup::selectHighResolution(structure = predicted6ofs,
                                                     predicted = TRUE)

  expect_s3_class(highResolution, "pdb")
  expect_equal(highResolution$seqres, predicted6ofs$seqres)
  expect_equal(3 * nrow(highResolution$atom), length(highResolution$xyz))
  expect_equal(nrow(highResolution$atom),
               sum(predicted6ofs$atom$"b" >= 70))
  expect_equal(all(highResolution$atom$"b" >= 70), TRUE)
})

test_that("Test on predicted structure with custom threshold", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  threshold <- 40
  highResolution <- PDBCleanup::selectHighResolution(structure = predicted6ofs,
                                                     predicted = TRUE,
                                                     threshold = threshold)

  expect_s3_class(highResolution, "pdb")
  expect_equal(highResolution$seqres, predicted6ofs$seqres)
  expect_equal(3 * nrow(highResolution$atom), length(highResolution$xyz))
  expect_equal(nrow(highResolution$atom),
               sum(predicted6ofs$atom$"b" >= threshold))
  expect_equal(all(highResolution$atom$"b" >= threshold), TRUE)
})

test_that("Test on experimental structure with default threshold", {
  experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
                                   package = "PDBCleanup")
  experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
  highResolution <- PDBCleanup::selectHighResolution(
    structure = experimental6ofs,
    predicted = FALSE
  )

  expect_s3_class(highResolution, "pdb")
  expect_equal(highResolution$seqres, experimental6ofs$seqres)
  expect_equal(3 * nrow(highResolution$atom), length(highResolution$xyz))
  belowThreshold <- (experimental6ofs$atom$"b" <= 100)
  expect_equal(unique(highResolution$atom$resno),
               unique(experimental6ofs$atom[belowThreshold, "resno"]))
})

test_that("Test on experimental structure with custom threshold", {
  experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
                                      package = "PDBCleanup")
  experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
  threshold <- 80
  highResolution <- PDBCleanup::selectHighResolution(
    structure = experimental6ofs,
    predicted = FALSE,
    threshold = threshold
  )

  expect_s3_class(highResolution, "pdb")
  expect_equal(highResolution$seqres, experimental6ofs$seqres)
  expect_equal(3 * nrow(highResolution$atom), length(highResolution$xyz))
  belowThreshold <- (experimental6ofs$atom$"b" <= threshold)
  expect_equal(unique(highResolution$atom$resno),
               unique(experimental6ofs$atom[belowThreshold, "resno"]))
})
