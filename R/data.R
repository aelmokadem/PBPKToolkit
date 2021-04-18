#' NHANES dataset (nhanesData)
#'
#' Dataset containing anthropometric data for individuals 2 years of age and above extracted from NHANES database.
#'
#' @format A data frame with 25769 rows and 8 variables:
#' \describe{
#'   \item{SEX}{Sex; 1=male and 2=female}
#'   \item{AGE_YR}{Age in years}
#'   \item{AGE_MN}{Age in months}
#'   \item{RACE}{Race}
#'   \item{BW}{Body weight, in kg}
#'   \item{HT}{Height, in cm}
#'   \item{BMI}{Body mass index, in kg/m^2}
#'   \item{BSA}{Body surface area, in m^2}
#' }
#' @source \url{https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm}
NULL

#' Percentages of blood content (BC)
#'
#' A dataset containing the percentages of blood content in 20 human organs.
#'
#' @format A data frame with 20 rows and 2 variables:
#' \describe{
#'   \item{organ}{Human body organs}
#'   \item{bloodPerc}{Percentage of blood content}
#' }
#' @source \url{https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4}
NULL

#' Coefficient of variation (normSD)
#'
#' A dataset containing the coefficient of variation for normally distributed organ volumes.
#'
#' @format A data frame with 17 rows and 2 variables:
#' \describe{
#'   \item{organ}{Human body organs}
#'   \item{cv}{Coefficient of variation}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/17431751/}
NULL

#' Geometric standard deviation (lnormSD)
#'
#' A dataset containing the geometric standard deviation for log normally distributed organs.
#'
#' @format A data frame with 4 rows and 2 variables:
#' \describe{
#'   \item{organ}{Human body organs}
#'   \item{geomSD}{Geometric standard deviation}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/17431751/}
NULL

#' Scaling factors (SF)
#'
#' A dataset containing the scaling factors for allometric scaling of organs with height.
#'
#' @format A data frame with 4 rows and 2 variables:
#' \describe{
#'   \item{organ}{Human body organs}
#'   \item{SF}{Scaling factors}
#' }
#' @source \url{https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm}
NULL


