# Regenerate tests/testthat/fixtures/v1_reference.rds.
#
# Run from the package root, against the version whose numbers should be
# pinned:
#   Rscript --vanilla tests/testthat/fixtures/make_v1_reference.R
#
# Originally generated from tseLCA 1.1.1 (commit b65d12d). Regenerate only for
# deliberate numerical changes, and document each one in NEWS.md.

pkgload::load_all(".", quiet = TRUE)
source("tests/testthat/helper-v1-reference.R")

data_list <- v1_data()
ref <- lapply(names(v1_configs), function(nm) {
  message("Fitting ", nm)
  v1_extract(v1_fit(v1_configs[[nm]], data_list))
})
names(ref) <- names(v1_configs)
attr(ref, "tseLCA_version") <- as.character(utils::packageVersion("tseLCA"))
attr(ref, "created") <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

saveRDS(ref, "tests/testthat/fixtures/v1_reference.rds", version = 2)
message("Wrote ", length(ref), " reference fits")
