#' Data frame with amino acid names
#'
#' Data frame with full amino acid names, shorthand and one letter
#'
#' @docType data
#'
#' @usage data(AminoAcids)
#'
#' @format An object of class \code{\link[base:data.frame]{"data.frame"}}.
#'
#' @keywords datasets
"AminoAcids"

#' A mock SBA list
#'
#' A list with a mock autoimmunity SBA data set.
#'
#' @docType data
#'
#' @usage data(MockSBA)
#'
#' @format An object of class \code{\link[base:list]{"list"}} with 6 elements:
#' \describe{
#'   \item{MFI}{A data.frame with raw MFI values, samples as rows and beads as columns.}
#'   \item{COUNT}{A data.frame with bead count, samples as rows and beads as columns.}
#'   \item{SAMPLES}{A data.frame with sample information, samples as rows.}
#'   \item{BEADS}{A data.frame with bead information, beads as rows.}
#'   \item{CT}{A data.frame with raw MFI values from the coupling efficiency test, replicates as rows and beads as columns.}
#'   \item{FILTERINFO}{A string vector, starting will NULL but is expanded depending on performed filtering steps.}
#' }
#' @keywords datasets
"MockSBA"
