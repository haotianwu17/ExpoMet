#' Example exposome dataset for ExpoNet
#'
#' A de-identified example dataset used to demonstrate
#' the main ExpoNet workflow.
#'
#' @format A data frame with 395 rows and 5,254 columns. The first column,
#' `exp_id`, contains participant identifiers. The remaining 5,253 columns
#' (`C0001` through `C5252`, plus `C9999`) contain simulated exposomic
#' feature measurements.
"exposome_example"

#' Example phenotype dataset for ExpoNet
#'
#' A simulated phenotype dataset used to demonstrate
#' the main ExpoNet workflow.
#'
#' @format A data frame with rows as participants and columns including:
#' \describe{
#'   \item{exp_id}{The ID for the participants.}
#'   \item{sex}{Participant sex, coded as 1 or 2.}
#'   \item{age}{Participant age in years, ranging from 18 to 80.}
#'   \item{smoking}{Smoking status, coded as 1, 2, or 3.}
#'   \item{Disease}{Disease status, coded as "Yes" or "No".}
#' }
"pheno_example"

#' Simulated feature annotation key
#'
#' A simulated annotation table corresponding to features in the
#' `exposome_example` dataset.
#'
#' @format A data frame containing simulated exposomic feature identifiers
#' and associated feature annotation information.
"simulated_feature_key_final"
