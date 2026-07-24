#' Example exposome dataset for ExpoMet
#'
#' A de-identified example dataset used to demonstrate
#' the main ExpoMet workflow.
#'
#' @format A data frame with rows as participants and columns including:
#' \describe{
#'   \item{exp_id}{The ID for the participants.}
#'   \item{C0001}{The values for exposomic feature 'C0001'.}
#' }
"exposome_example"

#' Example phenotype dataset for ExpoMet
#'
#' A simulated phenotype dataset used to demonstrate
#' the main ExpoMet workflow.
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

#' Example dataset of the feature info for ExpoMet, associated with the "exposome_example" file
#'
#'
#' @format A data frame with rows as exposomic features, and columns are the associated info:
#' \describe{
#'   \item{feature_name}{Correspond to the column names in "exposome_example"}
#' }
"simulated_feature_key_final"
