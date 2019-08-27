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
#'     BEADS = Beads info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any beads with no text (ie. "") or "NegControl" in such column will be included in the transformation.
#'
#' @return Updated input x with the new list element
#'
#'     MADs = assay MADs.
#' @export

ap_mads2 <- function(x, constant=1, ...) {

  if("Filtered" %in% colnames(x$BEADS)){
  tmp_data <- x$MFI[, which(x$BEADS$Filtered == "" | grepl("NegControl", x$BEADS$Filtered))]
  } else {
    tmp_data <- x$MFI
  }

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


#' Scoring
#'
#' Binning of MADs values in Autoimmunity Profiling.
#'
#' @param x List with at least one element, see Deatils for naming and content.
#' It is recommended to use the the output from \code{\link[rappp:ap_mads2]{ap_mads2()}}.
#' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs).
#' @param rightmost.closed,left.open logical, see \code{\link[base:findInterval]{findInterval()}} for details.
#' @param check.names logical, see \code{\link[base:data.frame]{data.frame()}} for details
#' @param ... Further arguments passed do \code{\link[base:findInterval]{findInterval()}}
#' @details The input values will be binned into discrete bins (scores).
#'
#' The x list needs to include at least the element:
#'
#'     MADs = assay MADs,
#'
#' @return Updated input x with the new list elements
#'
#'     COKEY = Cutoff key as data.frame with cutoff values, scores and colors
#'
#'     SCORE = scored data
#' @export

ap_scoring2 <- function(x, MADlimits=seq(0,70,5),
                       rightmost.closed=FALSE, left.open=FALSE,
                       check.names=FALSE, ...) {

  xmad_score <- ap_cutoffs2(MADlimits)

  tmp_data <- x$MADS

  scores <- data.frame(matrix(t(apply(tmp_data, 1,
                                      function(i) findInterval(i, xmad_score$xmad[-1],
                                                               rightmost.closed = rightmost.closed,
                                                               left.open = left.open, ...))),
                              ncol=dim(tmp_data)[2],
                              dimnames=list(rownames(tmp_data),
                                            colnames(tmp_data)) )/10, check.names = check.names)

  x <- append(x, list(COKEY=xmad_score,
                           SCORE=scores))

  return(x)
}

#' Binary
#'
#' Create binary matrices based on scored Autoimmunity profiling data.
#'
#' @param x List with at least one element, see Deatils for naming and content.
#' It is recommended to use to element SCORE in the output from \code{\link[rappp:ap_scoring2]{ap_scoring2()}}.
#' @param cutoffs data.frame with at least one column named score with the desired cutoffs to use,
#' and rownames you want to have as identifier for each cutoff.
#' It is recommended to use the COKEY element in the output from \code{\link[rappp:ap_scoring2]{ap_scoring2()}}.
#' @details The input values will be binned into binary data, consisting of 0 and 1.
#'
#' The x list needs to include at least the element:
#'
#'     SCORE = scored data,
#'
#' @return Updated input x with the new list element
#'
#'     BINARY = list with one data.frame per cutoff
#'
#' @export

ap_binary2 <- function(x, cutoffs) {

  tmp_data <- x$SCORE

  binary_list <- lapply(cutoffs$score, function(cutoff)
        data.frame(ifelse(tmp_data >= cutoff, 1, 0), check.names = FALSE))

  names(binary_list) <- rownames(cutoffs)

  x <- append(x, list(BINARY=binary_list))

  return(x)
}

#' Cutoff selection
#'
#' Select cutoff based on the slope of the density of scores per antigen.
#'
#' @param x List with at least one element, see Deatils for naming and content.
#' It is recommended to use to element Scoring in the output from \code{\link[rappp:ap_scoring2]{ap_scoring2()}}.
#' @param cutoffs data.frame with at least one column named score with the desired cutoffs to use,
#' and rownames you want to have as identifier for each cutoff.
#' It is recommended to use the COKEY element in the output from \code{\link[rappp:ap_scoring2]{ap_scoring2()}}.
#' @param slope_cutoff Arbitrary slope cutoff value. Can be chosen freely.
#' @param offset Offset used to prevent script from finding the peak (as slope = 0 there).
#' @param bw Bandwidth for density funciton, default set to 0.1.
#' @details A cutoff will be selected for each antigen based on the
#' distribution of the scores for the antigen. The algorithm will search for a
#' local min nearest the highest peak in a density plot using bandwidth=0.1.
#'
#' The x list needs to include at least the element:
#'
#'     SCORE = scored data,
#'
#' @return Updated input x with the new list elements
#'
#'     DENS = Density output used for cutoff selection,
#'
#'     AGCO_CONT = Calculated antigen specific cutoffs, continues values,
#'
#'     AGCO = Calculated antigen specific cutoffs, translated into the descrete cutoff steps,
#'
#'     BINARY_CO = Binary table based on the antigen specific cutoffs.
#'
#' @export

ap_cutoff_selection2 <- function(x,
                                cutoffs,
                                slope_cutoff=-0.5,
                                offset=0.1,
                                bw=0.1) {

      inputdata <- x$SCORE
      inputdata <- inputdata[,-which(apply(inputdata, 2, function(i) sum(is.na(i))) == dim(inputdata)[1])]

      ## Apply density function on the scores to get x and y values for density-plots
      dens <- apply(inputdata, 2, function(y) density(y, bw=bw))

      # Calculate the slope in all points (except the last) by looking ahead one point. Do this for all beads.
      slope <- lapply(dens, function(y) diff(y$y)/diff(y$x))

      slope_cutoff_indices <- rep(NA, length.out = length(slope))
      names(slope_cutoff_indices) <- names(slope)
      for (i in 1:length(slope)) { #For each bead
        # AJF code is based on starting from the median, but I changed to highest peak on left or right side of plot.
        if (dens[[i]]$x[which.max(dens[[i]]$y)] <= max(cutoffs$score)/2) { #If the highest peak is on the left hand side of the plot (most cases).
          lookupStartIdx <- which.min(abs(dens[[i]]$x-(dens[[i]]$x[which.max(dens[[i]]$y)]+offset))) #Find start index just to the right of highest peak.
          lookupVector <- (lookupStartIdx):length(slope[[i]]) #Start looking from the start index to the end of the plot.
          CO.fun <- function(i,j,slope_cutoff) { #And set the cutoff function to look for the first point where the slope is above the chosen CO value
            slope[[i]][j] > slope_cutoff #Removed from AJF code (incorporated in start index instead): & dens[[i]]$x[j] >= score.median + minCOdistanceFromMedian
          }
        } else if (dens[[i]]$x[which.max(dens[[i]]$y)] > max(cutoffs$score)/2){ #Else, if the highest peak is on the right hand side of the plot
          lookupStartIdx <- which.min(abs(dens[[i]]$x-(dens[[i]]$x[which.max(dens[[i]]$y)]-offset)))#Find start index just to the left of highest peak.
          lookupVector <- (lookupStartIdx):1 #Start looking from the start index to the beginning of the plot.
          CO.fun <- function(i,j,slope_cutoff) { #And set the cutoff function to look for the first point where the slope is below the negative chosen CO value
            slope[[i]][j] < -slope_cutoff #Removed from AJF code (incorporated in start index instead): & dens[[i]]$x[j] <= score.median - minCOdistanceFromMedian
          }
        } else {
          stop("Error in slope-based cutoff assignment")
        }
        for (j in lookupVector) { #Now, for each slope step
          if (CO.fun(i,j,slope_cutoff)) { #If the slope at this step fulfils the cutoff function
            slope_cutoff_indices[i] <- j #Then store that slope step for that bead
            break() #And exit the slope step loop to continue with next bead
          }
        }
      } # End loop i for each bead (slope)

      slope_cutoff_scores <- rep(NA, length.out = length(dens))
      for (i in seq_along(dens)) {
        if (!is.na(slope_cutoff_indices[i])) {
          slope_cutoff_scores[i] <- dens[[i]]$x[slope_cutoff_indices[i]] #Find the score at the index where the slope is above cutoff. These are the cutoffs that are red lines in the density plot!
        } else {
          slope_cutoff_scores[i] <- max(dens[[i]]$x)+0.1 #If no cutoff was found (i.e. no samples were reactive), set the cutoff to 0.1 above the maximum score for which there is density.
        }
      } # End loop i for each bead (dens)
      names(slope_cutoff_scores) <- colnames(inputdata)

  # Translate the continuous slope cutoff values to the above discrete score value
      ag_score_cutoffs <- tibble(bead=names(slope_cutoff_scores),
                                 score=ceiling(slope_cutoff_scores*10)/10,
                                 xmad=cutoffs$xmad[match(score, cutoffs$score)])

      # Add NA elements for non-included beads to match original data
      dens <- dens[match(colnames(x$SCORE), names(dens))]
      names(dens) <- colnames(x$SCORE)

      binary_cutoff <- data.frame(do.call(cbind, lapply(1:dim(inputdata)[2], function(i)
        ifelse(inputdata[,i] >= ag_score_cutoffs$score[i], 1, 0))))
      colnames(binary_cutoff) <- colnames(inputdata)

      ag_score_cutoffs <- ag_score_cutoffs[match(colnames(x$SCORE), ag_score_cutoffs$bead),]
      ag_score_cutoffs$bead <- colnames(x$SCORE)

      slope_cutoff_scores <- slope_cutoff_scores[match(colnames(x$SCORE), names(slope_cutoff_scores))]
      names(slope_cutoff_scores) <- colnames(x$SCORE)

      binary_cutoff <- data.frame(binary_cutoff, NA, check.names=F)[, match(colnames(x$SCORE), colnames(binary_cutoff),
                                                                            nomatch=dim(binary_cutoff)[2]+1)]
      rownames(binary_cutoff) <- rownames(x$SCORE)
      colnames(binary_cutoff) <- paste0(ag_score_cutoffs$bead, "_co", ag_score_cutoffs$xmad, "xMAD")

      x <- append(x, list(DENS=dens,
                          AGCO_CONT=slope_cutoff_scores,
                          AGCO=ag_score_cutoffs,
                          BINARY_CO=binary_cutoff))
  return(x)
}

#' Calculate reactivity frequencies
#'
#' Create binary matrices based on scored Autoimmunity profiling data.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @param samplegroups factor vector of groupings. Only samples with an assigned level are included in plots.
#'     If left as \code{NULL} (default), the all non-filtered, if filetring done otherwise all, will be assigned "Sample".
#' @details
#'
#' The x list needs to include at least the elements:
#'
#'     BINARY = list with one data.frame per cutoff
#'
#'     BINARY_CO = Binary table based on the antigen specific cutoffs.
#'
#' @return Updated input x with the new list element
#'
#'     REACTSUM_AG = number of reactive samples per antigen and sample group,
#'
#'     REACTFREQ_AG = reactivity frequency per antigen and sample group,
#'
#'     REACTSUM_SAMP = number of reactive antigens per sample,
#'
#'     REACTFREQ_SAMP = reactivity frequency per sample,
#'
#' @export

ap_reactsummary2 <- function(x, samplegroups=NULL) {

  data_bin <- append(x$BINARY, list(Ag_selected=x$BINARY_CO))

  print("set samplegroups")
  if(is.null(samplegroups)){
    if("Filtered" %in% colnames(x$SAMPLES)){
      samplegroups <- factor(ifelse(x$SAMPLES$'Filtered' == "", "Sample", NA))
    } else {
      samplegroups <- factor(rep("Sample", dim(data_cont)[1]))
    }
  }
  data_size <- table(samplegroups)
  print(data_size)
  n_ag <- lapply(data_bin, function(i) apply(i, 1, function(l) sum(!is.na(l))))

  print("Per antigen")
  # Calculate per antigen
  data_sum_ag <- lapply(data_bin, function(i) apply(i, 2, function(l) aggregate(l, by=list(samplegroups), FUN=sum)))
  names(data_sum_ag) <- names(data_bin)

  data_freq_ag <- lapply(1:length(data_sum_ag),
                         function(cutoff) lapply(data_sum_ag[[cutoff]],
                                                 function(antigen) round(antigen$x/data_size*100,1)))
  names(data_freq_ag) <- names(data_sum_ag)


  data_sum_ag <- lapply(data_sum_ag, function(cutoff) do.call(cbind, lapply(cutoff, function(antigen) antigen$x)))
  rownames(data_sum_ag) <- levels(samplegroups)
  data_freq_ag <- lapply(data_freq_ag, function(cutoff) do.call(cbind, cutoff))

  print("Per sample")
  # Calculate per sample
  data_sum_samp <- lapply(data_bin,
                          function(i) apply(i, 1,
                                            function(l) sum(l, na.rm=T)))
  names(data_sum_samp) <- names(data_bin)

  data_freq_samp <- lapply(1:length(data_sum_samp), function(cutoff) round(data_sum_samp[[cutoff]]/n_ag[[cutoff]]*100,1))
  names(data_freq_samp) <- names(data_sum_samp)

  # Add to input
  x <- append(x,
              list(REACTSUM_AG=data_sum_ag,
                   REACTFREQ_AG=data_freq_ag,
                   REACTSUM_SAMP=data_sum_samp,
                   REACTFREQ_SAMP=data_freq_samp))

  return(x)
}

#' Full AP data transformation
#'
#' Wrapper function for full Autoimmunity Profiling data transformations.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs).
#' @param ... See respective functions for details:
#'     \code{\link[rappp:ap_mads2]{ap_mads2()}}, \code{\link[rappp:ap_scoring2]{ap_scoring2()}},
#'     \code{\link[rappp:ap_binary2]{ap_binary2()}}, \code{\link[rappp:ap_cutoff_selection2]{ap_cutoff_selection2()}}.
#' @details The x list needs to include at least the elements:
#'
#'     MFI = assay mfi,
#'
#'     BEADS = Beads info (Filtered column with information about filtering),
#'
#' @return Updated input x with the new list elements
#'
#'     MADs = assay MADs,
#'
#'     COKEY = Cutoff key as data.frame with cutoff values, scores and colors,
#'
#'     SCORE = scored data,
#'
#'     BINARY = list with one data.frame per cutoff,
#'
#'     DENS = Density output used for cutoff selection,
#'
#'     AGCO = Calculated antigen specific cutoffs, translated into the descrete cutoff steps,
#'
#'     AGCO_CONT = Calculated antigen specific cutoffs, continues values.
#'
#' @export

ap_norm2 <- function(x, MADlimits=seq(0,70,5), ...){

  print("Doing MADs transformation")
  x <- ap_mads2(x, ...)

  print("Doing Scoring")
  x <- ap_scoring2(x, MADlimits=MADlimits, ...)

  print("Doing Binary transformation")
  x <- ap_binary2(x, cutoffs=x$COKEY, ...)

  print("Finding cutoffs")
  x <- ap_cutoff_selection2(x, cutoffs=x$COKEY, ...)

  print("Summarize reactivities")
  x <- ap_reactsummary2(x, ...)

  return(x)
}
