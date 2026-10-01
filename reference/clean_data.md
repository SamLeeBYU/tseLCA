# Prepare and validate data for tseLCA estimation

Prepare and validate data for tseLCA estimation

## Usage

``` r
clean_data(
  data,
  Y.names,
  Zp.names = NULL,
  Zo.name = NULL,
  incomplete = FALSE,
  include.intercept = TRUE,
  verbose = FALSE,
  Zp.formula = NULL,
  Y.levels = NULL
)
```

## Arguments

- data:

  A data.frame.

- Y.names:

  Character vector of item column names.

- Zp.names:

  Character vector of covariate column names, or `NULL`. Ignored when
  `Zp.formula` is given.

- Zo.name:

  Single distal outcome column name, or `NULL`.

- incomplete:

  Logical. If `TRUE`, use FIML for partially-observed Y.

- include.intercept:

  Logical. Include an intercept in the covariate design built from
  `Zp.names`.

- verbose:

  Logical. Print row-drop messages.

- Zp.formula:

  One-sided formula for the covariate design (e.g.
  `~ age + factor(region)`), or `NULL` to build one from `Zp.names`.

- Y.levels:

  Named list of indicator categories, as returned by
  `.recode_indicators()`, when `data` already holds 0-based codes;
  `NULL` recodes the indicators here.

## Value

A named list with:

- Y.obs:

  N_Y x K expanded one-hot indicator matrix for Steps 1 & 2.

- mDesign:

  N_Y x K design/mask matrix (NULL when incomplete = FALSE).

- ivItemcat:

  Integer vector of category counts per item.

- Y.levels:

  Named list of the categories of each item.

- keep_Y:

  Integer indices of rows kept for Steps 1 & 2 (into original n).

- Z_mat:

  n_Z x (Q+1) covariate design matrix, or NULL.

- Zp.formula, Z_terms, Z_xlevels:

  The covariate formula, its terms, and the factor levels used, or NULL.

- keep_step3_Z_in_Y:

  Positions of Z-complete rows within keep_Y.

- Zo_mat:

  N_Zo x 1 distal outcome matrix, or NULL.

- keep_step3_Zo_in_Y:

  Positions of Zo-complete rows within keep_Y.

- keep_step3_Zo:

  Indices of Zo-complete rows (into original n).

- keep_step3_Zo_in_Z:

  With covariates, positions of the distal rows within the covariate
  rows (distal rows then also need complete covariates); otherwise NULL.
