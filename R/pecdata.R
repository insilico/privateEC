# pecdata.R - Bill White - 4/27/17

#' An MRI data set used used in the paper referenced below.
#'
#' The fMRI data includes 80 unmedicated MDD (52 females, age±sd. 33±11) and 80
#' healthy controls (HCs) (41 females, age±sd. 31±10). We used AFNI(Cox, 1996)
#' to process the rs-fMRI data and extract 3003 z-transformed correlation
#' coefficients between 78 brain regions identified by a functional region of
#' interest atlas (Shirer, et al., 2012). 3 bad subject removed for final
#' dimension: 157 x 3004 (3003 + pheno)
#'
#' @docType data
#' @keywords datasets
#' @name rsfMRIcorrMDD
#' @usage data(rsfMRIcorrMDD)
#' @references
#' Trang Le, W. K. Simmons, M. Misaki, B.C. White, J. Savitz, J. Bodurka,
#' and B. A. McKinney. “Privacy preserving evaporative cooling feature
#' selection and classification with Relief-F and Random Forests,”
#' Bioinformatics. Accepted. https://doi.org/10.1093/bioinformatics/btx298. 2017
#' @format
#' Data frame with 157 rows, 3003 all-pairs correlation variables and a pheno variable +/- 1.
NULL
