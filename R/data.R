#' Unified tissue composition data (TCData)
#'
#' Unified tissue composition data used for the partition coefficient calculation.
#'
#' @format A data frame with 13 rows and 13 variables:
#' \describe{
#'   \item{tissue}{Tissue}
#'   \item{f_water}{Fraction of water}
#'   \item{f_lipids}{Fraction of lipids}
#'   \item{f_proteins}{Fraction of proteins}
#'   \item{f_pl}{Fraction of phospholipids}
#'   \item{f_n_l}{Fraction of neutral lipids}
#'   \item{f_n_pl}{Fraction of neutral phospholipids}
#'   \item{f_a_pl}{Fraction of acidic phospholipids}
#'   \item{pH}{pH}
#'   \item{f_ew}{Fraction of extracellular water}
#'   \item{f_iw}{Fraction of intracellular water}
#'   \item{AR}{Albumin ratio}
#'   \item{LR}{Lipoprotein ratio}
#' }
#' @source \url{https://dmd.aspetjournals.org/content/48/10/903}
NULL

#' NHANES dataset (nhanesData)
#'
#' Dataset containing anthropometric data for individuals 2 years of age and above extracted from NHANES database.
#'
#' @format A data frame with 25769 rows and 8 variables:
#' \describe{
#'   \item{SEX}{Sex; 1=male and 2=female}
#'   \item{AGE_YR}{Age in years}
#'   \item{AGE_MN}{Age in months}
#'   \item{RACE}{Race; 1=non-Hispanic white, 2=non-Hispanic black, 3=Mexican American, 4=other race, 5=other Hispanic}
#'   \item{BW}{Body weight, in kg}
#'   \item{HT}{Height, in cm}
#'   \item{BMI}{Body mass index, in kg/m^2}
#'   \item{BSA}{Body surface area, in m^2}
#' }
#' @source \url{https://github.com/Open-Systems-Pharmacology/OSPSuite.Documentation/wiki/Create-Individual-Algorithm}
NULL

#' ICRP dataset (icrpData)
#'
#' Dataset containing physiologic parameters for typical individuals.
#'
#' @format A data frame with 6 rows and 47 variables:
#' \describe{
#'   \item{age}{Age, in years}
#'   \item{bw_m}{Male body weight, in kg}
#'   \item{bw_f}{Female body weight, in kg}
#'   \item{ht_m}{Male height, in cm}
#'   \item{ht_f}{Female height, in cm}
#'   \item{bsa_m}{Male body surface area, in m^2}
#'   \item{bsa_f}{Female body surface area, in m^2}
#'   \item{ad_m}{Male adipose volume, in L}
#'   \item{ad_f}{Female adipose volume, in L}
#'   \item{bo_m}{Male bone volume, in L}
#'   \item{bo_f}{Female bone volume, in L}
#'   \item{br_m}{Male brain volume, in L}
#'   \item{br_f}{Female brain volume, in L}
#'   \item{go_m}{Male gonads volume, in L}
#'   \item{go_f}{Female gonads volume, in L}
#'   \item{he_m}{Male heart volume, in L}
#'   \item{he_f}{Female heart volume, in L}
#'   \item{ki_m}{Male kidneys volume, in L}
#'   \item{ki_f}{Female kidneys volume, in L}
#'   \item{la_int_m}{Male large intestines volume, in L}
#'   \item{la_int_f}{Female large intestines volume, in L}
#'   \item{li_m}{Male liver volume, in L}
#'   \item{li_f}{Female liver volume, in L}
#'   \item{lu_m}{Male lungs volume, in L}
#'   \item{lu_f}{Female lungs volume, in L}
#'   \item{mu_m}{Male muscles volume, in L}
#'   \item{mu_f}{Female muscles volume, in L}
#'   \item{pa_m}{Male pancreas volume, in L}
#'   \item{pa_f}{Female pancreas volume, in L}
#'   \item{sk_m}{Male skin volume, in L}
#'   \item{sk_f}{Female skin volume, in L}
#'   \item{sm_int_m}{Male small intestines volume, in L}
#'   \item{sm_int_f}{Female small intestines volume, in L}
#'   \item{sp_m}{Male spleen volume, in L}
#'   \item{sp_f}{Female spleen volume, in L}
#'   \item{st_m}{Male stomach volume, in L}
#'   \item{st_f}{Female stomach volume, in L}
#'   \item{bl_m}{Male blood volume, in L}
#'   \item{bl_f}{Female blood volume, in L}
#'   \item{co_m}{Male cardiac output flow, in L/min}
#'   \item{co_f}{Female cardiac output flow, in L/min}
#'   \item{th_m}{Male thymus volume, in L}
#'   \item{th_f}{Female thymus volume, in L}
#'   \item{ln_m}{Male lymph nodes volume, in L}
#'   \item{ln_f}{Female lymph nodes volume, in L}
#'   \item{ot_m}{Male other organs volume, in L}
#'   \item{ot_f}{Female other organs volume, in L}
#' }
#' @source \url{https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4}
NULL

#' Percentages of blood flow (flow)
#'
#' A dataset containing the percentages of blood flows from cardiac output for human organs.
#'
#' @format A data frame with 21 rows and 2 variables:
#' \describe{
#'   \item{organ}{Human body organs}
#'   \item{flowPerc}{Percentage of blood flow}
#' }
#' @source \url{https://journals.sagepub.com/doi/pdf/10.1177/ANIB_32_3-4}
NULL

#' Percentages of blood content (BC)
#'
#' A dataset containing the percentages of blood content in human organs.
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

#' Physiological data for mAb PBPK model (mabShahData)
#'
#' Physiological data to be used for the monoclonal antibody PBPK model.
#'
#' @format A data frame with 18 rows and 9 variables:
#' \describe{
#'   \item{Organ}{Organ}
#'   \item{Vtot}{Total volume}
#'   \item{Vplas}{Plasma volume}
#'   \item{Vbc}{Blood cell volume}
#'   \item{Vi}{Interstitial volume}
#'   \item{Ve}{Endosomal volume}
#'   \item{Vc}{Cellular volume}
#'   \item{Qplas}{Plasma flow}
#'   \item{Qbc}{Blood cell flow}
#' }
#' @source \url{https://pubmed.ncbi.nlm.nih.gov/22143261/}
NULL

