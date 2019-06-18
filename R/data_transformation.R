#' MAD normalization
#'
#' Sample based normalization to number of Median Absolute Deviations from the median (MADs)
#'
#' @details The input values will be normalized per sample to the number of MADs from the median
#' using the algorithm MADs = (MFI - median )/MAD, where MAD is calculated using mad(constant=1)
#'
#' The input values should be MFI values, and structured as a list, even if only one data set is used, see examples.
#'
#' @param x List of MFI values with two levels per element: level one = assay data sets ; level two =  bead subsets (e.g. wih and w/o controls)
#' @param constant Constant for mad() function, default is 1 (compared to 1.4826 in base function).
#' @return List of MADs, with same structure as input list.
#' @examples Input structure examples:
#' # One assay data set with one subset
#' list(Assay=list(All=SBA@X))
#' # One assay data set but with two different subsets
#' list(Assay=list(All=SBA@X,
#' WithoutControls=SBA@X[,SBA@Beads$Type != "Control"]))
#' # Two assay data sets with two different subsets
#' list(CSF=list(All=SBA_csf@X,
#' WithoutControls=SBA_csf@X[,SBA_csf@Beads$Type != "Control"]),
#' Plasma=list(All=SBA_plasma@X,
#' WithoutControls=SBA_plasma@X[,SBA_plasma@Beads$Type != "Control"]))
#' @export

mads <- function(x, constant=1, ...) {
  lapply(x,
         function(assay) lapply(assay,
                                function(selection) {
                                  (selection-apply(selection, 1, function(x) median(x, ...)))/
                                    apply(selection, 1, function(x) mad(x, constant=constant, ...)) } ))
}
