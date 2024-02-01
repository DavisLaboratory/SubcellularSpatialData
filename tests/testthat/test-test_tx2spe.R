test_that("Transcript aggregation works", {
  library(SpatialExperiment)
  library(dplyr)

  data(tx_small)
  required_cols = c("sample_id", "cell", "gene", "genetype", "x", "y", "counts", "technology")

  # check required columns throw errors
  for (i in required_cols) {
    expect_error(tx2spe(tx_small[, setdiff(colnames(tx_small), i)]), "missing")
  }
  expect_error(tx2spe(tx_small[, setdiff(colnames(tx_small), "region")], "region"), "missing")

  # bins > 0
  expect_error(tx2spe(tx_small, "hex", 0), "greater")
  expect_error(tx2spe(tx_small, "square", 0), "greater")

  # test types
  expect_error(tx2spe(mutate(tx_small, x = as.character(x))), "numeric")
  expect_error(tx2spe(mutate(tx_small, y = as.character(y))), "numeric")
  expect_error(tx2spe(mutate(tx_small, counts = as.character(counts))), "numeric")
})
