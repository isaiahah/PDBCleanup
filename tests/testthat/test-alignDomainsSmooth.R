library(PDBCleanup)

# Tests adapted from alignDomainsRigid tests, as similar goals
test_that("Test on invalid inputs", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
                                      package = "PDBCleanup")
  experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
  domains6ofs <- list(list(36:486, 36:486), list(511:930, 511:930))

  expect_error(PDBCleanup::alignDomainsSmooth(NA, NA, NA))
  expect_error(PDBCleanup::alignDomainsSmooth(predicted6ofs, NA, domains6ofs))
  expect_error(PDBCleanup::alignDomainsSmooth(NA, experimental6ofs, domains6ofs))
  expect_error(PDBCleanup::alignDomainsSmooth(predicted6ofs, experimental6ofs, NA))
})

test_that("Reduced RMSD after alignment", {
  predicted6ofsFile <- system.file("extdata", "6ofs_predicted.pdb",
                                   package = "PDBCleanup")
  predicted6ofs <- bio3d::read.pdb(predicted6ofsFile)
  experimental6ofsFile <- system.file("extdata", "6ofs_experimental.pdb",
                                      package = "PDBCleanup")
  experimental6ofs <- bio3d::read.pdb(experimental6ofsFile)
  domains6ofs <- list(list(36:486, 36:486), list(511:930, 511:930))

  predictedAligned <- PDBCleanup::alignDomainsSmooth(experimental6ofs, predicted6ofs,
                                                    domains6ofs)
  # Aligned should diverge from unaligned
  expect_gt(rmsd(predictedAligned, predicted6ofs,
                 a.inds = 1:nrow(predictedAligned$atom),
                 b.inds = 1:nrow(predicted6ofs$atom)), 0)
  # Aligned should be better fit than unaligned
  expect_gt(rmsd(predicted6ofs, experimental6ofs,
                 a.inds = 1:nrow(predicted6ofs$atom),
                 b.inds = 1:nrow(predicted6ofs$atom)),
            rmsd(predictedAligned, experimental6ofs,
                 a.inds=1:nrow(predictedAligned$atom),
                 b.inds=1:nrow(predictedAligned$atom)))
  # No other guarantees, due to unmoved regions outside domains
})

