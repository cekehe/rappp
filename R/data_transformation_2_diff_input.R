#' MAD normalization
#'
#' Sample based normalization to number of Median Absolute Deviations (MADs)
#' from the median for Autoimmunity profiling data.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @param constant Constant for \code{\link[stats:mad]{mad()}} function, default is 1 (compared to 1.4826 in base function).
#' @param ... Further arguments passed to \code{\link[stats:median]{median()}} and \code{\link[stats:mad]{mad()}}
#' @details The input values will be normalized per sample to the number of MADs from the median
#' using the algorithm MADs = (MFI - median )/MAD, where MAD is calculated using \code{mad(constant=1)}.
#'
#' The x list needs to include at least the elements:
#'
#'     MFI = assay mfi,
#'
#'     BEADS = Beads info (Filtered column with information about filtering),
#'
#' @return Updated input x with the new list element MADS.
#' @export

ap_mads2 <- function(x, constant=1, ...) {

  tmp_data <- x$MFI[, which(x$BEADS$Filtered == "")]

  mads <- (tmp_data - apply(tmp_data, 1, function(i) median(i, ...)))/
        apply(tmp_data, 1, function(i) mad(i, constant=constant, ...))

  mads <- data.frame(mads, NA, check.names=F)[, match(colnames(x$MFI), colnames(mads), nomatch=dim(mads)[2]+1)]
  colnames(mads) <- colnames(x$MFI)

  x <- append(x, list(MADS=mads))

  return(x)
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
#'    [,3] Corresponding color using the Zissou1 palette in \code{\link[wesanderson]{wes_palette}}
#' @export

ap_cutoffs2 <- function(MADlimits=seq(0,70,5)){

  xmad_score <- data.frame(xmad=c(NA, MADlimits),
                           score=c(0, 1:length(MADlimits)/10),
                           color=as.character(wes_palette(name = "Zissou1",
                                                          n = length(MADlimits)+1, type = "continuous")))
  rownames(xmad_score) <- c("Below0xMAD",paste0(MADlimits,rep("xMAD",length(MADlimits))))
  return(xmad_score)
}


#' #' Scoring
#' #'
#' #' Binning of MADs values in Autoimmunity Profiling.
#' #'
#' #' @details The input values will be binned into discrete bins (scores).
#' #'
#' #' The input values should be MADs values, and structured as a list
#' #' (preferably the output from function ap_mads()), even if only one data set is used
#' #' (see examples in \link[rappp]{ap_mads}(.
#' #'
#' #' @param x List of MADs values with two levels per element: level one = assay data sets ;
#' #' level two =  bead subsets (e.g. wih and w/o controls).
#' #' It is recommended is to use the the output from \link[rappp]{ap_mads}).
#' #' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs).
#' #' @param rightmost.closed,left.open logical, see \link[base]{findInterval} for details.
#' #' @param check.names logical, see \link[base]{data.frame} for details
#' #' @param ... Further arguments passed do \link[base]{findInterval}
#' #' @return List with two main elements
#' #'     [[1]] Cutoff key as data.frame with cutoff values, scores and colors
#' #'     [[2]] scored data, with same structure as input list.
#' #' @export
#'
#' ap_scoring2 <- function(x, MADlimits=seq(0,70,5),
#'                        rightmost.closed=FALSE, left.open=FALSE,
#'                        check.names=FALSE, ...) {
#'
#'   xmad_score <- ap_cutoffs(MADlimits)
#'
#'   scores <- lapply(x, function(assay)
#'     lapply(assay, function(selection)
#'       data.frame(matrix(t(apply(selection, 1, function(row)
#'         findInterval(row, xmad_score$xmad[-1],
#'                      rightmost.closed = rightmost.closed, left.open = left.open, ...))),
#'         ncol=dim(selection)[2],
#'         dimnames=list(rownames(selection),
#'                       colnames(selection)) )/10, check.names = check.names)))
#'
#'   output <- list(Cutoff_key=xmad_score,
#'                  Scoring=scores)
#' }
#'
#' #' Binary
#' #'
#' #' Create binary matrices based on scored Autoimmunity profiling data.
#' #'
#' #' @details The input values will be binned into binary data, consisting of 0 and 1.
#' #'
#' #' The input values should be scoring values, and structured as a list
#' #' (preferably the output from function ap_scoring()), even if only one data set is used
#' #' (see examples in \link[rappp]{ap_mads}).
#' #'
#' #' @param x List of scoring values with two levels per element: level one = assay data sets ;
#' #' level two =  bead subsets (e.g. wih and w/o controls).
#' #' It is recommended to use to element Scoring in the output from \link[rappp]{ap_scoring}).
#' #' @param cutoffs data.frame with at least one column named score with the desired cutoffs to use,
#' #' and rownames you want to have as identifier for each cutoff.
#' #' It is recommended is to use the Cutoff_key element in the output from \link[rappp]{ap_scoring}).
#' #'
#' #' @return List with binary data.frames
#' #' @export
#'
#' ap_binary2 <- function(x, cutoffs) {
#'
#'   binary_list <- lapply(x, function(assay)
#'     lapply(assay, function(selection)
#'       lapply(cutoffs$score, function(cutoff)
#'         data.frame(ifelse(selection >= cutoff, 1, 0), check.names = FALSE))))
#'
#'   binary_list <- lapply(binary_list, function(assay)
#'     lapply(assay, function(selection) { names(selection) <- rownames(cutoffs) ; selection }))
#'
#'   return(binary_list)
#' }
#'
#' #' Cutoff selection
#' #'
#' #' Select cutoff based on the slope of the density of scores per antigen.
#' #'
#' #' @details A cutoff will be selected for each antigen based on the
#' #' distribution of the scores for the antigen.
#' #' The algorithm will search for a local min nearest the highest
#' #' peak in a density plot using bandwidth=0.1.
#' #'
#' #' The input values should be scoring values, and structured as a list
#' #' (preferably the output from function ap_scoring()), even if only one data set is used
#' #' (see examples in \link[rappp]{ap_mads}).
#' #'
#' #' @param x List of scoring values with two levels per element: level one = assay data sets ;
#' #' level two =  bead subsets (e.g. wih and w/o controls).
#' #' It is recommended to use to element Scoring in the output from \link[rappp]{ap_scoring}).
#' #' @param cutoffs data.frame with at least two columns:
#' #'
#' #'     - One column named score with the desired cutoffs to use
#' #'
#' #'     - One column named xmad with the corresponding MAD cutoff values
#' #'
#' #' It is recommended is to use the Cutoff_key element in the output from \link[rappp]{ap_scoring}).
#' #' @param slope_cutoff Arbitrary slope cutoff value. Can be chosen freely.
#' #' @param offset Offset used to prevent script from finding the peak (as slope = 0 there).
#' #' @param bw Bandwidth for density funciton, default set to 0.1.
#' #' @return List with two main elements
#' #'
#' #'     [[1]] Density output used for cutoff selection, with same structure as input list.
#' #'
#' #'     [[2]] Calculated antigen specific cutoffs, with same structure as input list.
#' #' @export
#'
#' ap_cutoff_selection2 <- function(x,
#'                                 cutoffs,
#'                                 slope_cutoff=-0.5,
#'                                 offset=0.1,
#'                                 bw=0.1) {
#'
#'   dens <- lapply(x, function(y) rep(list(NULL), length(y)))
#'
#'   slope_cutoff_scores <- rep(list(NULL), length(x))
#'   for(assay in seq_along(x)){
#'     slope_cutoff_scores[[assay]] <- rep(list(NULL), length(x[[assay]]))
#'
#'     for(selection in seq_along(x[[assay]])){
#'       inputdata <- x[[assay]][[selection]]
#'
#'       ## Apply density function on the scores to get x and y values for density-plots
#'       dens[[assay]][[selection]] <- apply(inputdata, 2, function(y) density(y, bw=bw))
#'
#'       # Calculate the slope in all points (except the last) by looking ahead one point. Do this for all beads.
#'       slope <- lapply(dens[[assay]][[selection]], function(y) diff(y$y)/diff(y$x))
#'
#'       slope_cutoff_indices <- rep(NA, length.out = length(slope))
#'       names(slope_cutoff_indices) <- names(slope)
#'       for (i in 1:length(slope)) { #For each bead
#'         # AJF code is based on starting from the median, but I changed to highest peak on left or right side of plot.
#'         if (dens[[assay]][[selection]][[i]]$x[which.max(dens[[assay]][[selection]][[i]]$y)] <= max(cutoffs$score)/2) { #If the highest peak is on the left hand side of the plot (most cases).
#'           lookupStartIdx <- which.min(abs(dens[[assay]][[selection]][[i]]$x-(dens[[assay]][[selection]][[i]]$x[which.max(dens[[assay]][[selection]][[i]]$y)]+offset))) #Find start index just to the right of highest peak.
#'           lookupVector <- (lookupStartIdx):length(slope[[i]]) #Start looking from the start index to the end of the plot.
#'           CO.fun <- function(i,j,slope_cutoff) { #And set the cutoff function to look for the first point where the slope is above the chosen CO value
#'             slope[[i]][j] > slope_cutoff #Removed from AJF code (incorporated in start index instead): & dens[[assay]][[selection]][[i]]$x[j] >= score.median + minCOdistanceFromMedian
#'           }
#'         } else if (dens[[assay]][[selection]][[i]]$x[which.max(dens[[assay]][[selection]][[i]]$y)] > max(cutoffs$score)/2){ #Else, if the highest peak is on the right hand side of the plot
#'           lookupStartIdx <- which.min(abs(dens[[assay]][[selection]][[i]]$x-(dens[[assay]][[selection]][[i]]$x[which.max(dens[[assay]][[selection]][[i]]$y)]-offset)))#Find start index just to the left of highest peak.
#'           lookupVector <- (lookupStartIdx):1 #Start looking from the start index to the beginning of the plot.
#'           CO.fun <- function(i,j,slope_cutoff) { #And set the cutoff function to look for the first point where the slope is below the negative chosen CO value
#'             slope[[i]][j] < -slope_cutoff #Removed from AJF code (incorporated in start index instead): & dens[[assay]][[selection]][[i]]$x[j] <= score.median - minCOdistanceFromMedian
#'           }
#'         } else {
#'           stop("Error in slope-based cutoff assignment")
#'         }
#'         for (j in lookupVector) { #Now, for each slope step
#'           if (CO.fun(i,j,slope_cutoff)) { #If the slope at this step fulfils the cutoff function
#'             slope_cutoff_indices[i] <- j #Then store that slope step for that bead
#'             break() #And exit the slope step loop to continue with next bead
#'           }
#'         }
#'       } # End loop i for each bead (slope)
#'
#'       for (i in seq_along(dens[[assay]][[selection]])) {
#'         if (!is.na(slope_cutoff_indices[i])) {
#'           slope_cutoff_scores[[assay]][[selection]][i] <- dens[[assay]][[selection]][[i]]$x[slope_cutoff_indices[i]] #Find the score at the index where the slope is above cutoff. These are the cutoffs that are red lines in the density plot!
#'         } else {
#'           slope_cutoff_scores[[assay]][[selection]][i] <- max(dens[[assay]][[selection]][[i]]$x)+0.1 #If no cutoff was found (i.e. no samples were reactive), set the cutoff to 0.1 above the maximum score for which there is density.
#'         }
#'       } # End loop i for each bead (dens[[assay]][[selection]])
#'       names(slope_cutoff_scores[[assay]][[selection]]) <- colnames(x[[assay]][[selection]])
#'     } #End loop for each selection
#'     names(slope_cutoff_scores[[assay]]) <- names(x[[assay]])
#'     names(dens[[assay]]) <- names(x[[assay]])
#'   } #End loop for each assay
#'   names(slope_cutoff_scores) <- names(x)
#'
#'   # Translate the continuous slope cutoff values to the above discrete score value
#'   ag_score_cutoffs <- lapply(slope_cutoff_scores, function(assay)
#'     lapply(assay, function(selection) tibble(bead=names(selection),
#'                                              score=ceiling(selection*10)/10,
#'                                              xmad=cutoffs$xmad[match(score, cutoffs$score)])))
#'
#'   output <- list(dens=dens,
#'                  Slope_cutoff_discrete=ag_score_cutoffs,
#'                  Slope_cutoff=slope_cutoff_scores)
#'   return(output)
#' }
#'
#'
#' #' Full AP data transformation
#' #'
#' #' Wrapper function for full Autoimmunity Profiling data transformations.
#' #'
#' #' @details The input values should be MFI values, and structured as a list,
#' #' even if only one data set is used (see examples).
#' #'
#' #' @param x List of MFI values with two levels per element: level one = assay data sets ;
#' #' level two =  bead subsets (e.g. wih and w/o controls)
#' #' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs).
#' #' @param ... See respective functions for details: \link[rappp]{ap_mads}, \link[rappp]{ap_scoring},
#' #' \link[rappp]{ap_binary}, \link[rappp]{ap_cutoff_selection}
#' #' @return List with 7 main elements
#' #'
#' #'     [[1]] Input MFI
#' #'
#' #'     [[2]] MADs
#' #'
#' #'     [[3]] Scores
#' #'
#' #'     [[4]] Binary
#' #'
#' #'     [[5]] Antigen specific cutoffs selected based on density slope.
#' #'
#' #'     [[6]] Density information
#' #'
#' #'     [[7]] Cutoff translation key
#' #' @export
#'
#' ap_norm2 <- function(x, MADlimits=seq(0,70,5), ...){
#'
#'   tmp_mads <- ap_mads(x, ...)
#'
#'   tmp_score <- ap_scoring(tmp_mads, MADlimits=MADlimits, ...)
#'
#'   tmp_binary <- ap_binary(tmp_score$Scoring, cutoffs=tmp_score$Cutoff_key, ...)
#'
#'   tmp_slope <- ap_cutoff_selection(tmp_score$Scoring, cutoffs=tmp_score$Cutoff_key, ...)
#'
#'   output <- list(MFI=x,
#'                  MADs=tmp_mads,
#'                  Scoring=tmp_score$Scoring,
#'                  Binary=tmp_binary,
#'                  Slope_cutoff=tmp_slope$Slope_cutoff_discrete,
#'                  Score_density=tmp_slope$dens,
#'                  Cutoff_key=tmp_score$Cutoff_key)
#'
#'   return(output)
#' }
