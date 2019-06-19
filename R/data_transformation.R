#' MAD normalization
#'
#' Sample based normalization to number of Median Absolute Deviations (MADs)
#' from the median for Autoimmunity profiling data.
#'
#' @details The input values will be normalized per sample to the number of MADs from the median
#' using the algorithm MADs = (MFI - median )/MAD, where MAD is calculated using mad(constant=1)
#'
#' The input values should be MFI values, and structured as a list,
#' even if only one data set is used (see examples).
#'
#' @param x List of MFI values with two levels per element: level one = assay data sets ;
#' level two =  bead subsets (e.g. wih and w/o controls)
#' @param constant Constant for mad() function, default is 1 (compared to 1.4826 in base function).
#' @param ... Further arguments passed do \link[stats]{median} and \link[stats]{mad}
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

ap_mads <- function(x, constant=1, ...) {
  lapply(x, function(assay)
    lapply(assay,function(selection) {
      (selection-apply(selection, 1, function(x) median(x, ...)))/
        apply(selection, 1, function(x) mad(x, constant=constant, ...)) } ))
}


#' Cutoff key
#'
#' Create a cutoff key for scoring of Autoimmunity Profiling data.
#'
#' @details The input values will be binned into discrete bins (scores).
#'
#' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs).
#' @return data.frame with three columns:
#'
#'    [,1] MADs cutoff value
#'
#'    [,2] Corresponding score value
#'
#'    [,3] Corresponding color using the Zissou1 palette in \link[wesanderson]{wes_palette}
#' @export

ap_cutoffs <- function(MADlimits=seq(0,70,5)){

  xmad_score <- data.frame(xmad=c(NA, MADlimits),
                           score=c(0, 1:length(MADlimits)/10),
                           color=as.character(wes_palette(name = "Zissou1",
                                                          n = length(MADlimits)+1, type = "continuous")))
  rownames(xmad_score) <- c("Below0xMAD",paste0(MADlimits,rep("xMAD",length(MADlimits))))
  return(xmad_score)
}


#' Scoring
#'
#' Binning of MADs values in Autoimmunity Profiling.
#'
#' @details The input values will be binned into discrete bins (scores).
#'
#' The input values should be MADs values, and structured as a list
#' (preferably the output from function ap_mads()), even if only one data set is used
#' (see examples in \link[rappp]{ap_mads}(.
#'
#' @param x List of MADs values with two levels per element: level one = assay data sets ;
#' level two =  bead subsets (e.g. wih and w/o controls).
#' It is recommended is to use the the output from \link[rappp]{ap_mads}).
#' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs).
#' @param rightmost.closed,left.open logical, see \link[base]{findInterval} for details.
#' @param check.names logical, see \link[base]{data.frame} for details
#' @param ... Further arguments passed do \link[base]{findInterval}
#' @return List with two main elements
#'     [[1]] Cutoff key as data.frame with cutoff values, scores and colors
#'     [[2]] scored data, with same structure as input list.
#' @export

ap_scoring <- function(x, MADlimits=seq(0,70,5),
                       rightmost.closed=FALSE, left.open=FALSE,
                       check.names=FALSE, ...) {

  xmad_score <- ap_cutoffs(MADlimits)

  scores <- lapply(x, function(assay)
    lapply(assay, function(selection)
      data.frame(matrix(t(apply(selection, 1, function(row)
        findInterval(row, xmad_score$xmad[-1],
                     rightmost.closed = rightmost.closed, left.open = left.open, ...))),
        ncol=dim(selection)[2],
        dimnames=list(rownames(selection),
                      colnames(selection)) )/10, check.names = check.names)))

  output <- list(Cutoff_key=xmad_score,
                 Scoring=scores)
}

#' Binary
#'
#' Create binary matrices based on scored Autoimmunity profiling data.
#'
#' @details The input values will be binned into binary data, consisting of 0 and 1.
#'
#' The input values should be scoring values, and structured as a list
#' (preferably the output from function ap_scoring()), even if only one data set is used
#' (see examples in \link[rappp]{ap_mads}).
#'
#' @param x List of scoring values with two levels per element: level one = assay data sets ;
#' level two =  bead subsets (e.g. wih and w/o controls).
#' It is recommended to use to element Scoring in the output from \link[rappp]{ap_scoring}).
#' @param cutoffs data.frame with at least one column named score with the desired cutoffs to use,
#' and rownames you want to have as identifier for each cutoff.
#' It is recommended is to use the Cutoff_key element in the output from \link[rappp]{ap_scoring}).
#'
#' @return List with binary data.frames
#' @export

ap_binary <- function(x, cutoffs) {

  binary_list <- lapply(x, function(assay)
    lapply(assay, function(selection)
      lapply(cutoffs$score, function(cutoff)
        data.frame(ifelse(selection >= cutoff, 1, 0), check.names = FALSE))))

  binary_list <- lapply(binary_list, function(assay)
    lapply(assay, function(selection) { names(selection) <- rownames(cutoffs) ; selection }))

  return(binary_list)
}
