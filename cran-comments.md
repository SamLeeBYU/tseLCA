## Test environments
* local machine, R 4.6.1
* win-builder (release and devel)
* Ubuntu Linux (R-devel)
* macOS arm64 (R-devel)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Downstream dependencies
There are no downstream dependencies for this package.

## Resubmission

This is a resubmission following reviewer feedback from Konstanze Lauseker.
Changes made:

- DESCRIPTION: Removed single quotes from acronyms (BCH, ML, LCA, etc.) and 
  added explanations of all acronyms in the description text.
- R/* : Replaced T/F with TRUE/FALSE throughout, including in 
  man/parse_rebase.Rd example.
- Added \value tags to all exported functions missing them, including 
  bk2018_params.Rd.
- inst/examples: Removed all default file writes to the home filespace; 
  examples now write only to tempdir().
- inst/examples: Removed rm(list = ls()) calls.
- inst/examples: Removed any install.packages() calls.
- Title: corrected hyphenated compound capitalization to "Three-Step" per
  reviewer (Uwe Ligges) request.