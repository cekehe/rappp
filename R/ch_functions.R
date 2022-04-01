#' Coupling efficency test
#'
#' Flags beads with signal similar to empty bead in coupling test, produces plot if wanted.
#'
#'
#' @param x List with at least three elements, see Deatils for naming and content.
#' @param empty_bead Column index for empty bead.
#' @param empty_co_multiple Number of sd above empty for cutoff.
#' @param types Which types of beads should be included in flagging? See details.
#' @param shouldplot Logical, should a plot be made?
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @param ... Further arguments passed to \code{\link[base:mean]{mean()}} and \code{\link[stats:sd]{sd()}}.
#' @details The x list needs to include at least the elements
#'     CT = coupling test mfi,
#'
#'     BEADS = Beads info. See below for required columns.
#'
#'     FILTERINFO = Vector with info on which filter steps has been done.
#'
#' The BEADS element needs at least the columns:
#'
#'     "Type" with info about type of content on bead,
#'     should include at least what is set in argument types (exact match).
#'     Eg. "PrEST" for PrESTs, or Full_length for full length representations,
#'
#'     "Plate" with numerical coupling plate number(s).
#'
#' Note: The function plots to a layout containing two areas.
#'
#' @return Updated input x with relevant filtering info and a pdf
#'    with plot (if \code{shouldplot=TRUE} and \code{shouldpdf=TRUE}).
#' @export

ap_ct <- function(x, empty_bead, empty_co_multiple=3, types="PrEST",
                  shouldplot=TRUE, shouldpdf=TRUE, filename="coupling_efficiency.pdf",
                  width=25, height=6, useDingbats=FALSE, ...) {

    empty_co <- mean(x$CT[,empty_bead], ...) + empty_co_multiple*sd(x$CT[,empty_bead], ...)

    types <- paste0("^", paste0(types, collapse="$|^"), "$")

    if(shouldplot){
      if(shouldpdf){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)
      }
      par(mar=c(12,4,2,8), cex.axis=0.8)
      layout(matrix(c(1,1,1,1,2), nrow=1))

      bs=beeswarm(x$CT, pch=16, las=2, corral="gutter", xaxt="n",
                  main="Copupling efficiency test", ylab="Signal intensity [MFI]",
                  pwcol=rep(ifelse(grepl("empty|bare|blank|neutravidin", colnames(x$CT),ignore.case=T), "orange",
                                   ifelse(grepl("his6abp|hisabp", colnames(x$CT),ignore.case=T), "darkgreen",
                                          ifelse(grepl("hig|anti-human", colnames(x$CT),ignore.case=T), "blue",
                                                 ifelse(grepl("ebna", colnames(x$CT),ignore.case=T), "purple",
                                                        ifelse(!grepl(types, x$BEADS$Type), "lightgrey",
                                                               ifelse(apply(x$CT, 2, median, na.rm=T) < empty_co, "red",
                                                                      "darkgrey")))))),
                            each=dim(x$CT)[1]))

          vert_lines <- seq(min(bs$x)-0.5, max(bs$x)+0.5, 1)
      abline(v=vert_lines, lty=2,
             col=ifelse(c(duplicated(x$BEADS$Plate),F), "lightgrey", "black"))

      axis(1, at=rollmean(vert_lines, 2), labels=unique(bs$x.orig), cex.axis=0.5, las=2, tick=F)

      tmp <- vert_lines[which(!c(duplicated(x$BEADS$Plate),F))]
      mtext(text=paste("Plate ",unique(x$BEADS$Plate)), at=diff(tmp)/2+tmp[-length(tmp)],
            side=1, line=0, cex=0.7, font=2)

      legend(par("usr")[2], par("usr")[4],
             legend=c("Empty", "His6ABP", "ahIgX", "EBNA1", "Passed", "Flagged", "Not relevant"),
             col=c("orange", "darkgreen", "blue", "purple", "darkgrey", "red", "lightgrey"),
             pch=16, xpd=NA)
      abline(h=empty_co, lty=2)
      textxy(X=par("usr")[2], Y=empty_co, offset=0.55, cex=1,
             labs=paste0("mean(empty)+", empty_co_multiple,"*sd(empty)=", round(empty_co, 0)), xpd=NA)

      tmp_text <- data.frame(Name=colnames(x$CT), Type=x$BEADS$Type)
      if(length(which(apply(x$CT, 2, median, na.rm=T) < empty_co &
                      grepl(types, tmp_text$Type))) > 0){
        tmp_text <- matrix(tmp_text$Name[which(apply(x$CT, 2, min) < empty_co &
                                          grepl(types, tmp_text$Type))], ncol=1)
        ap_textplot(tmp_text, mar=c(2,2,1,2),
                 show.rownames=F, show.colnames=F, hadj=0, valign="top", cex=0.8)
        mtext("Beads with low coupling efficiency signal", font=2, cex=0.9, xpd=NA)
      } else {
        frame()
        mtext("No beads displayed \n low coupling efficiency signal.", font=2, cex=0.7, line=-3)
      }
      if(shouldpdf){
      dev.off()
      }
    }

    # Annotate filtering in BEADS
    if(length(which(colnames(x$BEADS) == "Flagged")) == 0){
      x$BEADS <- data.frame(Flagged="", x$BEADS, stringsAsFactors=F)
    }

    x$BEADS$Flagged <- ifelse(apply(x$CT, 2, mean) < empty_co &
                                grepl(types, x$BEADS$Type, ignore.case=T),
                              paste0(x$BEADS$Flagged,", Coupling"),
                              paste(x$BEADS$Flagged))
    x$BEADS$Flagged <- gsub("^, ", "", x$BEADS$Flagged)

    x$FILTERINFO <- c(x$FILTERINFO, "CouplingEfficiency")

    return(x)
}


#' Loading control
#'
#' Filter samples with low MFI for the anti-human IgX bead.
#'
#' @param x List with at least three elements, see Deatils for naming and content.
#' @param IgX_bead Column index for anti-IgX bead.
#' @param IgType Which Imunoglobulin is measured, default is G.
#' @param IgX_cutoff MFI cutoff value for filtering.
#' @param cosfac Median absolute deviation multipliers in vector c(upper, lower),
#'     for drawing lines and detecting potential outliers.
#' @param internal_sampID Column name in SAMPLES with internal sample IDs, such as from LIMS.
#'     Replicates (named with one of pool|rep|mix|commercial, not case sensitive)
#'     and blanks (named with one of empty|blank|buffer, not case sensitive) must be stated in this column.
#' @param external_sampID Column name in SAMPLES with externaö sample IDs, such as given from the collaborator or user.
#' @param shouldplot Logical, should a plot be made?
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @param ... Arguments are passed to base plot function.
#' @details The x list needs to include at least the elements:
#'
#'     MFI = assay mfi,
#'
#'     SAMPLES = Sample info. See below for required columns.
#'
#'     FILTERINFO = Vector with info on which filter steps has been done.
#'
#' The SAMPLES element needs at least the columns:
#'
#'     "AssayWell" with Well IDs in the assay plate, e.g A01, B01 etc.,
#'
#' Note: The function plots to a layout containing seven areas.
#'
#' @return Updated input x with relevant filtering info and a pdf
#'     with plot (if \code{shouldplot=TRUE} and \code{shouldpdf=TRUE}).
#' @export

ap_igx <- function(x, IgX_bead, IgType="G", IgX_cutoff=5000, cosfac=c(3, -3),
                   internal_sampID="sample_name", external_sampID="tube_label",
                   shouldplot=TRUE, shouldpdf=TRUE, filename="anti-humanIgX.pdf",
                   width=12, height=6, useDingbats=FALSE, ...) {

  plotdata <- unlist(x$MFI[,IgX_bead])
  sampledata <- x$SAMPLES
  SamplesNames <- sampledata[,internal_sampID]

  tmp <- subset(plotdata, !grepl("empty|blank|buffer", SamplesNames, ignore.case=T))
  cosIgG <- median(tmp, na.rm=T)+cosfac*mad(tmp, constant = 1, na.rm=T)
  cosIgG <- c(cosIgG, IgX_cutoff)

  which_lowIgG <- which(plotdata<cosIgG[3])

  if(shouldplot){
    if(shouldpdf){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)
    }
    layout(matrix(c(1,1,1,2,3,
                    1,1,1,4,5,
                    1,1,1,6,7), nrow=3, byrow=T))
    par(mar=c(5,5,4,4))

    plot(1:length(plotdata), plotdata, cex=0.6, ...,
         xlab="Samples in analysis order",ylab="Signal intensity (MFI)",main=paste0("Anti-hIg", IgType, ""),
         col=ifelse(grepl("empty|blank|buffer", SamplesNames, ignore.case=T),2,
                    ifelse(grepl("pool|rep|mix|commercial", SamplesNames, ignore.case=T),5, 4)))

    legend(par("usr")[1], par("usr")[4], horiz=T, yjust=0.1, bty="n",
           legend=c("Sample","Replicate","Buffer"),
           cex=0.8,
           fill=c(4, 5, 2), border=NA,
           xpd=NA)

    abline(h=cosIgG, lty=2)
    calibrate::textxy(X=rep(par("usr")[2], 3),Y=cosIgG, labs=c(paste0(cosfac,"xMAD+median (samples)"), paste0("Filter cutoff (",IgX_cutoff,")")),
                      offset=0.6, xpd=NA)

    # Display samples outside boundaries
    par(mar=c(5,1,4,1))
    plottext_all <- data.frame(AssayWell=sampledata$AssayWell,
                               InternalID=sampledata[,internal_sampID],
                               Subject=sampledata[,external_sampID],
                               MFI=plotdata)
    plottext_all$MFIgroup <- ifelse(plottext_all$MFI > cosIgG[1], paste0("above ", cosIgG[1]),
                                    ifelse(plottext_all$MFI < cosIgG[3], paste0("below ", cosIgG[3]),
                                           ifelse(plottext_all$MFI < cosIgG[2], paste0("between ", cosIgG[3]," & ", cosIgG[2]), NA)))
    mfi_groups <- c("above", "between", "below")
    for(i in mfi_groups){
      if(sum(grepl(i, plottext_all$MFIgroup)) > 0){
        plottext <- plottext_all[grep(i, plottext_all$MFIgroup),]
        tmp_name <- unique(plottext$MFIgroup)
        plottext <- plottext[order(plottext$MFI, decreasing=T), -which(colnames(plottext) == "MFIgroup")]

        if(dim(plottext)[1] > 20){
          ap_textplot(plottext[1:20,],
                      halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top", xpd=NA)
          ap_textplot(plottext[21:dim(plottext)[1],],
                      halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top", xpd=NA)
          mtext(paste0("anti-hIg", IgType, " MFI ", tmp_name), font=2, cex=0.6, xpd=NA, at=-0.5)
        } else {
          ap_textplot(plottext,
                      halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top", xpd=NA)
          mtext(paste0("anti-hIg", IgType, " MFI ", tmp_name), font=2, cex=0.6)
          frame()
        }

      } else {
        frame()
        frame()
      }

    }
    if(shouldpdf){
      dev.off()
    }
  }

  # Annotate filtering in SAMPLES
  if(!("Filtered" %in% colnames(x$SAMPLES))){
    x$SAMPLES <- data.frame(Filtered="", x$SAMPLES, stringsAsFactors=F, check.names=F)
  }
  if(length(which_lowIgG) > 0) {
    tmp_remove <- rownames(sampledata)[which_lowIgG]
    if(length(grep("empty|blank|buffer", tmp_remove, ignore.case=T)) > 0){
      tmp_remove <- tmp_remove[-grep("empty|blank|buffer", tmp_remove, ignore.case=T)]
    }
    if(length(tmp_remove) > 0){
      x$SAMPLES$Filtered <- ifelse(rownames(x$SAMPLES) %in% tmp_remove,
                                   paste0(x$SAMPLES$Filtered, ", hIg", IgType),
                                   paste(x$SAMPLES$Filtered))
      x$SAMPLES$Filtered <- gsub("^, ", "", x$SAMPLES$Filtered)
    }
  }
  x$FILTERINFO <- c(x$FILTERINFO, paste0("anti-hIg", IgType))

  return(x)
}


#' Bead count
#'
#' Filter and/or flag samples and beads with low bead count.
#'
#' @param x List with at least four elements, see Deatils for naming and content.
#' @param internal_sampID Column name in SAMPLES with internal sample IDs, such as from LIMS.
#'     Replicates (named with one of pool|rep|mix|commercial, not case sensitive)
#'     and blanks (named with one of empty|blank|buffer, not case sensitive) must be stated in this column.
#' @param external_sampID Column name in SAMPLES with externaö sample IDs, such as given from the collaborator or user.
#' @param Aglabels Column name in BEADS with antigen names to be used in pdf.
#' @param protein Column name in BEADS with short protein or gene name.
#' @param agID Column name in BEADS with antigen identifier, eg. PrEST ID or product number.
#' @param samp_co Cutoff for filtering samples with low median count.
#' @param bead_flag Cutoff for flagging beads with low counts.
#' @param bead_filter Cutoff for filtering beads with low counts.
#' @param N_filter Accepted number of samples with low count per bead ID.
#' @param bead_dispense How many wells are beads dispensed in per aspiration?
#'     If all in one dispense, then input is NA, 0 , or 1.
#'     If evenly divided, then input is one value with number of wells per dispense, e.g. 96.
#'     If unevenly distributed, then input is a vector of how many wells per dispense,
#'     e.g. \code{c(rep(96, 3), rep(48, 2))} for first three dispenses being 96 at a time,
#'     and the last 96 being divided in two dispenses.
#' @param luminex_wash After how many wells are there washes in the Luminex?
#' @param presampfilter Logical, should samples with annotation under Filtered be removed prior to count evaluation?
#' @param shouldplot Logical, should a plot be made?
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @param cex_axis_samp Text size on x-axis in sample-based plots.
#' @param cex_axis_ag Text size on x-axis in bead-based plots.
#' @param ... Other arguments passed to ap_textplot.
#' @details The x list needs to include at least the elements:
#'
#'     COUNT = bead count,
#'
#'     BEADS = Beads info. See below for required columns.
#'
#'     SAMPLES = Sample info. See below for required columns.
#'
#'     FILTERINFO = Vector with info on which filter steps has been done.
#'
#' The BEADS element needs at least the columns:
#'
#'     "BeadID" with bead ID number,
#'
#'     separate columns with full desired labels, protein/gene names and product ID numbers.
#'
#' The SAMPLES element needs at least the columns:
#'
#'     "AssayWell" with Well IDs in the assay plate, e.g A01, B01 etc.,
#'
#'     "Filtered" if \code{presampfilter=TRUE}
#'
#' Note: The function plots to a layout containing up to four areas.
#'
#' @return Updated input x with relevant filtering and/or flagging info and a pdf
#'     with plots (if \code{shouldplot=TRUE} and \code{shouldpdf=TRUE}).
#' @export

ap_count <- function(x, internal_sampID="sample_name", external_sampID="tube_label",
                     Aglabels="Gene_agID", protein="GeneShort", agID="PrEST",
                     samp_co=32, bead_flag=32, bead_filter=16, N_filter=0,
                     bead_dispense=32, luminex_wash=96, presampfilter=FALSE,
                     shouldplot=TRUE, shouldpdf=TRUE, filename="bead_count.pdf",
                     width=12, height=10, useDingbats=FALSE,
                     cex_axis_samp = 0.1, cex_axis_ag = 0.1, ...) {

  plotdata <- t(x$COUNT)
  sampledata <- x$SAMPLES
  beaddata <- x$BEADS

  if(presampfilter){
    plotdata <- plotdata[, which(sampledata$Filtered == "")]
    sampledata <- sampledata[which(sampledata$Filtered == ""), ]
  }

  which_lowSB <- which(apply(plotdata, 2, function(x) median(x, na.rm=T)) < samp_co) # Checks which samples have low count in general, e.g. due to faulty bead dispens in well

  if(shouldplot & shouldpdf){
    pdf(filename, width=width, height=height, useDingbats=useDingbats)
  }

  if(length(bead_dispense) == 1){
    if(is.na(bead_dispense) | bead_dispense == 0){
      bead_dispense <- "burlywood1"
    } else {
      bead_dispense <- c("burlywood1", "burlywood4")[c(rep(1, bead_dispense),
                                                       rep(2, bead_dispense))]
    }
  } else if (sum(bead_dispense) == ncol(plotdata)){
    tmp <- NULL
    for(i in seq_along(bead_dispense)){
      if(i %% 2 != 0){
        tmp <- c(tmp, rep(1, bead_dispense[i]))
      } else {
        tmp <- c(tmp, rep(2, bead_dispense[i]))
      }
    }
    bead_dispense <- c("burlywood1", "burlywood4")[tmp]
  } else {
    stop("bead_dispense not correctly entered")
  }

  if(is.na(luminex_wash) | luminex_wash == 0){
    luminex_wash <- -100
  } else if(length(luminex_wash) == 1){
    luminex_wash <- c(seq(luminex_wash, dim(plotdata)[2], luminex_wash), dim(plotdata)[2])+0.5
  } else {
    luminex_wash <- luminex_wash+0.5
  }

  for(state in c("before", "after")){

    if(shouldplot){

      if(state == "before"){
        layout(matrix(c(1,1,1,2,
                        3,3,3,4), nrow=2, byrow=T))
        par(mar=c(5, 4, 3, 10))
      } else {
        layout(matrix(c(1,
                        2), nrow=2, byrow=T))
        par(mar=c(5, 4, 3, 10))
      }

      # Per sample ALL DATA
      boxplot(plotdata, pch=16, cex=0.6, ylim=c(0, max(plotdata, na.rm=T)), las=1, names=F, xaxt="n",
              col=bead_dispense,
              border=bead_dispense,
              ylab="Bead count per sample")
      text(1:dim(plotdata)[2],par("usr")[3]-1, labels = rownames(sampledata),
           srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=cex_axis_samp)
      abline(h=c(16,32, median(plotdata, na.rm=T)), lty=2, col=c("grey","red", "cornflowerblue"))
      abline(v=luminex_wash)
      legend(par("usr")[2], par("usr")[4],
             legend=c(bead_filter,
                      paste0("Failed (", samp_co, ")"),
                      paste0("Median (", median(plotdata, na.rm=T), ")"),
                      rep("Bead dispensing batch", length(unique(bead_dispense))),
                      ifelse(luminex_wash[1] != -100, "Luminex wash", "")),
             lty=c(rep(2,3), rep(0, length(unique(bead_dispense))+1)),
             pch=c(rep(NA, 3+length(unique(bead_dispense))), ifelse(luminex_wash[1] != -100, 73, 0)),
             col=c("grey","red", "cornflowerblue", rep(0, length(unique(bead_dispense))), ifelse(luminex_wash[1] != -100, "black", 0)),
             fill=c(rep(0, 3), unique(bead_dispense), 0),
             border=c(rep(0,3), unique(bead_dispense), 0),
             xpd=NA, cex=0.7, bty="n")
      mtext("Sample wells, in order of analysis", side=1, cex=0.7, line=1)
      mtext("Sample bead count", side=3, cex=1, font=2, line=1)
      if(state == "after"){
        mtext(paste0("Filtered: Samples with median bead count < ", samp_co,
                     " and analytes with >", N_filter, " samples with bead count < ", bead_filter, " removed"),
              side=3, cex=0.6, line=0)
      }
    }

    if(state == "before"){
      if(length(which_lowSB) > 0 & shouldplot){
        lowSB <- sampledata[which_lowSB, ]
        ap_textplot(data.frame(
          AssayWell=lowSB$AssayWell,
          InternalID=lowSB[,internal_sampID],
          ExternalID=lowSB[,external_sampID],
          MedianCount=apply(plotdata, 2, function(x) median(x, na.rm=T))[which_lowSB],
          LowestCount=apply(plotdata, 2, function(x) min(x, na.rm=T))[which_lowSB],
          HighestCount=apply(plotdata, 2, function(x) max(x, na.rm=T))[which_lowSB]),
          show.rownames=F, valign="top", halign="left",
          hadj=0, vadj=0, mar=c(1, 2, 3, 8), xpd=NA, ...)#, cex=0.4)
        title(paste0("Samples with median bead count < ", samp_co, ", (N=", dim(lowSB)[1], ")"), xpd=NA)
      } else if(shouldplot){
        ap_textplot(matrix("No samples filtered based on bead count."), show.rownames=F, show.colnames=F)
      }

      # REMOVE LOW SAMPLES
      if(length(which_lowSB) > 0) { # Removes samples from plotdata with low median bead count
        plotdata <- t(plotdata[, -which_lowSB])
        sampledata <- sampledata[-which_lowSB, ]
      } else {
        plotdata <- t(plotdata)
      }

      # Per antigen
      which_failAB <- which(apply(plotdata, 2, function(x) length(which(x < bead_filter))) > N_filter)
      which_lowAB <- which(apply(plotdata, 2, function(x) length(which(x < bead_flag))) > 0)
    } else {
      plotdata <- t(plotdata)
    }

    if(shouldplot){
      boxplot(plotdata, pch=16, cex=0.6, ylim=c(0, max(plotdata, na.rm=T)), las=1, names=F, xaxt="n",
              ylab="Bead count per analyte")
      text(1:dim(plotdata)[2],par("usr")[3]-1, labels = beaddata[,Aglabels],
           srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=cex_axis_ag)
      abline(h=c(bead_filter,
                 bead_flag,
                 median(unlist(plotdata), na.rm=T)), lty=2, col=c("red", "orange", "cornflowerblue"))
      legend(par("usr")[2], par("usr")[4],
             legend=c(paste0("Failed if N>", N_filter, " (", bead_filter,")"),
                      paste0("Flagged (", bead_flag,")"),
                      paste0("Median (", median(unlist(plotdata), na.rm=T), ")")),
             lty=c(rep(2,3)),
             col=c("red", "orange", "cornflowerblue"),
             fill=c(rep(0,3)),
             border=c(rep(0,3)),
             xpd=NA, cex=0.7, bty="n")
      mtext("Analyte bead count", side=3, cex=1, font=2, line=1)
      if(state == "before"){
        mtext(paste0("Samples with median bead count < ", samp_co, " removed (see above plot)"),
              side=3, cex=0.6, line=0)
      } else {
        mtext(paste0("Filtered: Samples with median bead count < ", samp_co,
                     " and analytes with >", N_filter, " samples with bead count < ", bead_filter, " removed"),
              side=3, cex=0.6, line=0)
      }
    }

    if(state == "before"){
      if(length(which_lowAB) > 0){
        lowAB <- beaddata[which_lowAB, ]
        lowAB <- data.frame(
          BeadID=lowAB$BeadID,
          Analyte=lowAB[,protein],
          ID=lowAB[,agID],
          MedianCount=apply(plotdata, 2, function(x) median(x, na.rm=T))[which_lowAB],
          LowestCount=apply(plotdata, 2, function(x) min(x, na.rm=T))[which_lowAB],
          Nbelow16=apply(plotdata, 2, function(x) length(which(x < 16)))[which_lowAB],
          HighestCount=apply(plotdata, 2, function(x) max(x, na.rm=T))[which_lowAB])
        lowAB <- data.frame(lowAB, Action=ifelse(lowAB$LowestCount > bead_filter |
                                                   lowAB$Nbelow16 <= N_filter, "Flagged", "Filtered"))
        lowAB <- lowAB[order(lowAB$Action, lowAB$LowestCount),]

        if(shouldplot){
          ap_textplot(lowAB, show.rownames=F, valign="top", halign="left",
                      hadj=0, vadj=0, mar=c(1, 1, 3, 8), xpd=NA, ...)#, cex=0.3)
          title(paste0("Analytes with any bead count < ", bead_flag, " (N=",dim(lowAB)[1],")"), xpd=NA)
        }
      } else {
        if(shouldplot){
          ap_textplot(matrix("No analytes filtered or flagged based on bead count."), show.rownames=F, show.colnames=F)
        }
      }

      # REMOVE LOW ANTIGENS
      if(length(which_failAB) > 0) {
        plotdata <- t(plotdata[, -which_failAB])
        beaddata <- beaddata[-which_failAB, ]
      } else {
        plotdata <- t(plotdata)
      }
    }

  }

  if(shouldplot & shouldpdf){
    dev.off()
  }

  # Annotate filtering in SAMPLES and BEADS
  # SAMPLES
  if(!("Filtered" %in% colnames(x$SAMPLES))){
    x$SAMPLES <- data.frame(Filtered="", x$SAMPLES, stringsAsFactors=F, check.names=F)
  }
  if(length(which_lowSB) > 0){
    x$SAMPLES$Filtered <- ifelse(rownames(x$SAMPLES) %in% names(which_lowSB),
                                 paste0(x$SAMPLES$Filtered,", Count"),
                                 paste(x$SAMPLES$Filtered))
    x$SAMPLES$Filtered <- gsub("^, ", "", x$SAMPLES$Filtered)
  }

  # BEADS filtered
  if(!("Filtered" %in% colnames(x$BEADS))){
    x$BEADS <- data.frame(Filtered="", x$BEADS, stringsAsFactors=F, check.names=F)
  }

  if(length(which_lowAB) > 0){
    if(length(which(lowAB$Action == "Filtered")) > 0){
      x$BEADS$Filtered <- ifelse(rownames(x$BEADS) %in% rownames(lowAB)[which(lowAB$Action == "Filtered")],
                                 paste0(x$BEADS$Filtered,", Count"),
                                 paste(x$BEADS$Filtered))
      x$BEADS$Filtered <- gsub("^, ", "", x$BEADS$Filtered)
    }
  }

  # BEADS flagged
  if(!("Flagged" %in% colnames(x$BEADS))){
    x$BEADS <- data.frame(Flagged="", x$BEADS, stringsAsFactors=F, check.names=F)
  }

  if(length(which_lowAB) > 0){
    if(length(which(lowAB$Action == "Flagged")) > 0){
      x$BEADS$Flagged <- ifelse(rownames(x$BEADS) %in% rownames(lowAB)[which(lowAB$Action == "Flagged")],
                                paste0(x$BEADS$Flagged,", Count"),
                                paste(x$BEADS$Flagged))
      x$BEADS$Flagged <- gsub("^, ", "", x$BEADS$Flagged)
    }
  }

  x$FILTERINFO <- c(x$FILTERINFO, "beadcount")

  return(x)
}

#' Signal overview
#'
#' Boxplots of signals per antigen and sample in different orders.
#'
#' @param x List with at least one elements, see Deatils for naming and content.
#' @param grepneg Regular expression for color of negative control beads.
#' @param greppos Regular expression for color of positive control beads.
#' @param includeFilter Logical, should antigen boxes be colored based on filtering?
#' @param Aglabels Column in x$BEADS with matching names as columns in x$MFI. Only used if includeFilter is TRUE.
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @param ... Further arguments passed to \code{\link[graphics:boxplot]{boxplot()}}.
#' @details The x list needs to include at least the element
#'     MFI = assay mfi.
#'
#'     BEADS = Beads info. Only needed if includeFilter is TRUE.
#'
#'     Data points with the value NA or 0 will be set to 1 for the plotting to allow for
#'     logarithmic scale without filtering any beads or samples.
#'
#' Note: The function plots to a layout containing three areas.
#'
#' @export

ap_overview <- function(x,
                        grepneg = "his6abp|hisabp|empty|bare|biotin|neutravidin|neg_",
                        greppos = "anti-h|hIg|ebna",
                        includeFilter = TRUE,
                        Aglabels = "Gene_agID",
                        shouldpdf=TRUE,
                        filename="Signal_overview.pdf",
                        width=25, height=15, useDingbats=FALSE, ...){

  if(shouldpdf){
    pdf(filename,
        width=width, height=height, useDingbats=useDingbats)
  }
  par(mfcol=c(3, 1), mar=c(15,5,4,2))

  tmp_data <- x$MFI
  tmp_data[which(is.na(tmp_data) | tmp_data == 0, arr.ind=T)] <- 1

  ## Antigens
  plotdata <- list('Bead ID'=tmp_data,
                   median=tmp_data[,order(apply(tmp_data, 2, median, na.rm=T))],
                   max=tmp_data[,order(apply(tmp_data, 2, max, na.rm=T))])

  for(i in 1:length(plotdata)){

    colors <- ifelse(grepl(grepneg, colnames(plotdata[[i]]), ignore.case=T),
                     as.color("brown", 0.7),
                     ifelse(grepl(greppos, colnames(plotdata[[i]]), ignore.case=T),
                            as.color("darkolivegreen", 0.7),
                            0))

    if(includeFilter & "Filtered" %in% colnames(x$BEADS)){
      tmp_beads <- x$BEADS[match(colnames(plotdata[[i]]), x$BEADS[, Aglabels]),]

      colors <- ifelse(tmp_beads$Filtered != "", "grey", colors)
    }

    boxplot(plotdata[[i]], pch=16, cex=0.5, log="y", las=2, xaxt="n", ...,
            main=paste0("Antigens, sorted by ", names(plotdata)[i]), ylab="log(MFI) [AU]",
            outcol=as.color(ifelse(colors == 0, "black", colors), 0.6),
            col=colors)

    cex_xaxis <- c(1,1,0.75, 0.35, 0.25, 0.1)[findInterval(dim(plotdata[[i]])[2], c(1, seq(96, 96*5, 96)))]
    axis(1, at=1:dim(plotdata[[i]])[2], labels=colnames(plotdata[[i]]), cex.axis=cex_xaxis, las=2)
  }

  ## Samples
  tmp <- data.frame(t(tmp_data), check.names=F)
  plotdata <- list('analysis order'=tmp,
                   median=tmp[,order(apply(tmp, 2, median, na.rm=T))],
                   max=tmp[,order(apply(tmp, 2, max, na.rm=T))])

  for(i in 1:length(plotdata)){
    boxplot(plotdata[[i]], pch=16, cex=0.5, log="y", las=2,  xaxt="n", ...,
            main=paste0("Samples, sorted by ", names(plotdata)[i]), ylab="log(MFI) [AU]",
            outcol=ifelse(grepl("empty|buffer|blank", colnames(plotdata[[i]]), ignore.case=T), as.color("brown", 0.7),
                          ifelse(grepl("rep|pool|mix|commercial", colnames(plotdata[[i]]), ignore.case=T), as.color("cornflowerblue", 0.7), as.color("black", 0.5))),
            col=ifelse(grepl("empty|buffer|blank", colnames(plotdata[[i]]), ignore.case=T), as.color("brown", 0.7),
                       ifelse(grepl("rep|pool|mix|commercial", colnames(plotdata[[i]]), ignore.case=T), as.color("cornflowerblue", 0.7), as.color("black", 0.5))))

    cex_xaxis <- c(1,1,0.75, 0.35, 0.25, 0.1)[findInterval(dim(plotdata[[i]])[2], c(1, seq(96, 96*5, 96)))]
    axis(1, at=1:dim(plotdata[[i]])[2], labels=colnames(plotdata[[i]]), cex.axis=cex_xaxis, las=2)
  }
  if(shouldpdf){
    dev.off()
  }
}


#' Check replicate samples
#'
#' various plots in one PDF to assess the replicates and overall reproducibility.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @param iter How many times  random samples should be iterated.
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @details The x list needs to include at least the element
#'     MFI = assay mfi,
#'
#'     SAMPLES = Sample info. See below for required columns.
#'
#' The SAMPLES element needs at least the columns:
#'
#'     "sample_name" with sample names, preferably LIMS-IDs, where
#'     replicates (named with one of pool|rep|mix|commercial)
#'
#'     "AssayNum" with assay number (vector with 1s if only one assay),
#'
#' Note: The function plots to a layout, number of areas depending on number
#'     of different replicates and assays.
#'
#' @export

ap_rep <- function(x, iter=500, shouldpdf=TRUE, filename="replicates.pdf",
                   width=12, height=12, useDingbats=FALSE){

  ## INPUT
  if(sum(grepl("pool|rep|mix|commercial", x$SAMPLES$sample_name, ignore.case=T)) == 0){
    warning("No replicates found (substrings pool, rep, mix or commercial missing in sample_name)")
  } else {
      if("AssayNum" %in% colnames(x$SAMPLES)){
      data <- append(split(x$MFI, x$SAMPLES$AssayNum), list(x$MFI))
      names(data) <- c(paste0("Assay_", names(data)[-length(data)]), "All_combined")

      samples <- append(split(x$SAMPLES, x$SAMPLES$AssayNum), list(x$SAMPLES))
      names(samples) <- c(paste0("Assay_", names(samples)[-length(samples)]), "All_combined")
    } else {
      data <- list(Assay=x$MFI)
      samples <-  list(Assay=x$SAMPLES)
    }

    ## CALCULATIONS
    # True replicates
    {
      replicates <- lapply(1:length(data), function(i) split(data[[i]][grep("pool|rep|mix|commercial", samples[[i]]$sample_name, ignore.case=T),],
                                                             samples[[i]]$sample_name[grep("pool|rep|mix|commercial", samples[[i]]$sample_name, ignore.case=T)]))
      replicates <- lapply(replicates, function(i) i[which(lapply(i, function(l) dim(l)[1]) > 1)])
      nrreplicates <- lapply(replicates, function(i) lapply(i, function(l) dim(l)[1]))

      CVs_rep <- lapply(replicates, function(i) lapply(i, function(l) apply(l, 2, cv, digits=5, na.rm=T)))

      cor_rep_s <- lapply(replicates, function(i) lapply(i, function(l) cor(t(l), method="spearman", use="pairwise.complete.obs")))
      cor_rep_s <- lapply(cor_rep_s, function(i) lapply(i, function(l) l[upper.tri(l)]))

      cor_rep_p <- lapply(replicates, function(i) lapply(i, function(l) cor(t(l), method="pearson", use="pairwise.complete.obs")^2))
      cor_rep_p <- lapply(cor_rep_p, function(i) lapply(i, function(l) l[upper.tri(l)]))
    }

    # Iterate over random sets of samples
    {
      tmp_data <- lapply(1:length(data), function(i) data[[i]][-grep("pool|rep|mix|commercial",
                                                                     samples[[i]]$sample_name, ignore.case=T),])

      rand_samp <- lapply(1:length(replicates), function(i)
        lapply(1:length(replicates[[i]]), function(l)
          lapply(1:iter, function(m) tmp_data[[i]][sample(1:dim(tmp_data[[i]])[1], nrreplicates[[i]][[l]], replace=F),])))

      CVs_rand <- lapply(rand_samp, function(i) lapply(i, function(l) lapply(l, function(m)
        apply(m, 2, cv, digits=5, na.rm=T))))

      cor_rand_s_tmp <- lapply(rand_samp, function(i) lapply(i, function(l) lapply(l, function(m)
        cor(t(m), method="spearman", use="pairwise.complete.obs"))))
      cor_rand_s <- lapply(cor_rand_s_tmp, function(i) lapply(i, function(l) lapply(l, function(m)
        m[upper.tri(m)])))

      cor_rand_p_tmp <- lapply(rand_samp, function(i) lapply(i, function(l) lapply(l, function(m)
        cor(t(m), method="pearson", use="pairwise.complete.obs")^2)))
      cor_rand_p <- lapply(cor_rand_p_tmp, function(i) lapply(i, function(l) lapply(l, function(m)
        m[upper.tri(m)])))

      CVs_rand_m <- lapply(CVs_rand, function(i) lapply(i, function(l)
        apply(do.call(cbind, l), 1, median, na.rm=T)))
      cor_rand_s_m <- lapply(cor_rand_s, function(i) lapply(i, function(l)
        apply(do.call(cbind, l), 1, median, na.rm=T)))
      cor_rand_p_m <- lapply(cor_rand_p, function(i) lapply(i, function(l)
        apply(do.call(cbind, l), 1, median, na.rm=T)))

      CVs_rand_m <- lapply(1:length(replicates), function(i) {
        names(CVs_rand_m[[i]]) <- names(replicates[[i]]) ; CVs_rand_m[[i]]})
      cor_rand_s_m <- lapply(1:length(replicates), function(i) {
        names(cor_rand_s_m[[i]]) <- names(replicates[[i]]) ; cor_rand_s_m[[i]]})
      cor_rand_p_m <- lapply(1:length(replicates), function(i) {
        names(cor_rand_p_m[[i]]) <- names(replicates[[i]]) ; cor_rand_p_m[[i]]})
    }

    names(CVs_rep) <- names(data)
    names(cor_rep_s) <- names(data)
    names(cor_rep_p) <- names(data)

    names(CVs_rand_m) <- paste0(names(data), "_random")
    names(cor_rand_s_m) <- paste0(names(data), "_random")
    names(cor_rand_p_m) <- paste0(names(data), "_random")

    dens_rep_s <- lapply(cor_rep_s, function(i) lapply(i, function(l) density(l)))
    dens_rep_p <- lapply(cor_rep_p, function(i) lapply(i, function(l) density(l)))
    dens_rand_s <- lapply(cor_rand_s_m, function(i) lapply(i, function(l) density(l)))
    dens_rand_p <- lapply(cor_rand_p_m, function(i) lapply(i, function(l) density(l)))

    assay_cv <- lapply(data, function(y) apply(y, 2, function(x) cv(x, na.rm=T, digits=5)))
    assay_max <- lapply(data, function(y) apply(y, 2, function(x) max(x, na.rm=T)))
    assay_mean <- lapply(data, function(y) apply(y, 2, function(x) mean(x, na.rm=T)))
    assay_median <- lapply(data, function(y) apply(y, 2, function(x) median(x, na.rm=T)))

    ## PLOTS
    if(shouldpdf){
    pdf(filename,
        width=12, height=12, useDingbats=F)
    }
    par(mar=c(10,4,3,1))

    # CV boxplots
    data_melt <- rbind(melt(CVs_rep), melt(CVs_rand_m))
    boxplot(value~L1+L2, data=data_melt, outcol=0, las=1, col=grey.colors(2), names=NA, ylim=c(0, max(data_melt$value, na.rm=T)),
            ylab="CVs [%]", xlab=NA, main="CVs between replicates \n one point = one antigen")
    beeswarm(value~L1+L2, data=data_melt, pch=16, cex=0.5, corral="gutter", add=T)
    legend("topleft",
           legend=c("True replicates",
                    paste0("False replicates (", iter," iterations)"),
                    "CV=10%"),
           fill=c(grey.colors(2), NA),
           border=c(rep("black", 2), 0),
           lty=c(rep(NA, 2), 2),
           col=c(rep(NA, 2), "red"))
    abline(h=10, col="red", lty=2)
    tmp_text <- paste0(rep(unlist(nrreplicates), each=2), " samples")
    mtext(tmp_text, at=seq_along(tmp_text), side=1, line=0.5, cex=0.8)
    text(1:length(tmp_text),par("usr")[3]-3, labels = levels(interaction(data_melt$L1, data_melt$L2)),
         srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=0.9)

  # Correlations
    par(mfrow=c(sum(unlist(lapply(cor_rep_s, length))),2),
        mar=c(4,4,3,1))
    for(l in 1:length(dens_rep_s)){
      for(m in 1:length(dens_rep_s[[l]])){

        plot_rep <- dens_rep_s[[l]][[m]]
        plot_rand <- dens_rand_s[[l]][[m]]

        plot(range(0, 1), range(plot_rep$y, plot_rand$y), type = "n", xlab = "Spearman's rho",
             ylab = "Density", main=paste0(names(dens_rep_s)[l],", ", names(dens_rep_s[[l]])[m], ": ",
                                           nrreplicates[[l]][[m]], " samples"))
        lines(plot_rep, col = "black")
        lines(plot_rand, col = "cornflowerblue")
        legend("topleft", legend=c("True replicates", paste0("False replicates (", iter," iterations)")),
               lty=1, col=c("black","cornflowerblue"), cex=0.7)

        plot_rep <- dens_rep_p[[l]][[m]]
        plot_rand <- dens_rand_p[[l]][[m]]

        plot(range(0, 1), range(plot_rep$y, plot_rand$y), type = "n", xlab = bquote("Pearson's R"^"2"),
             ylab = "Density", main=paste0(names(dens_rep_s)[l],", ", names(dens_rep_s[[l]])[m], ": ",
                                           nrreplicates[[l]][[m]], " samples"))
        lines(plot_rep, col = "black")
        lines(plot_rand, col = "cornflowerblue")
      }
    }

  # Assay CV vs replicate CV
  par(mfrow=c(2,3))
  for(l in 1:length(CVs_rep)){
    for(i in 1:length(CVs_rep[[l]])){
    plot(CVs_rep[[l]][[i]], assay_cv[[l]], pch=16, cex=0.7, xlim=c(0,50),
         xlab="CVs of replicates, per antigen [%]", ylab="CVs of assay, per antigen [%]",
         main=paste0(names(CVs_rep)[l], ": ", names(CVs_rep[[l]])[i]))
    textxy(CVs_rep[[l]][[i]][which(CVs_rep[[l]][[i]] > 10)],
         assay_cv[[l]][which(CVs_rep[[l]][[i]] > 10)], cex=0.4, offset=0.55,
         labs=names(assay_cv[[l]])[which(CVs_rep[[l]][[i]] > 10)])
    abline(v=10, lty=2, col="grey")
    abline(h=10, lty=2, col="grey")
    }
  }

  # Assay CV vs assay max
  par(mfrow=c(2,3))
  for(l in 1:length(assay_max)){
    plot(assay_max[[l]], assay_cv[[l]], pch=16, cex=0.7,
         xlim=c(0, max(assay_max[[l]], na.rm=T)), ylim=c(0, max(assay_cv[[l]], na.rm=T)),
         xlab="Max signal intensity per antigen [MFI]", ylab="CVs of assay, per antigen [%]", main=names(assay_max)[l])
    abline(h=10, lty=2, col="grey")
  }

  # Assay CV vs assay median
  par(mfrow=c(2,3))
  for(l in 1:length(assay_median)){
    plot(assay_median[[l]], assay_cv[[l]], pch=16, cex=0.7,
         xlim=c(0, max(assay_median[[l]], na.rm=T)), ylim=c(0, max(assay_cv[[l]], na.rm=T)),
         xlab="Median signal intensity per antigen [MFI]", ylab="CVs of assay, per antigen [%]", main=names(assay_median)[l])
    abline(h=10, lty=2, col="grey")
  }

  # Assay CV vs assay mean
  par(mfrow=c(2,3))
  for(l in 1:length(assay_mean)){
    plot(assay_mean[[l]], assay_cv[[l]], pch=16, cex=0.7,
         xlim=c(0, max(assay_mean[[l]], na.rm=T)), ylim=c(0, max(assay_cv[[l]], na.rm=T)),
         xlab="Mean signal intensity per antigen [MFI]", ylab="CVs of assay, per antigen [%]", main=names(assay_mean)[l])
    abline(h=10, lty=2, col="grey")
  }

  if(shouldpdf){
  dev.off()
  }
  }
}


#' tSNE plots with different perplexities
#'
#' tSNE plots with different perplexities and designs.
#'
#' @param z the data as a data.frame
#' @param perp a vector of values to test as perplexities, eg. c(2,5,10,50).
#' @param sqrt logical, should sqrt(Nsamples) be added to perplexities?
#'     Potentially the lowest needed perplexity.
#' @param iterations number of iterations per tSNE
#' @param groups grouping for colors and lines
#' @param names text to display as point labels
#' @param legend TRUE or FALSE
#' @param legendname Title of legend
#' @param main Title of page
#' @param filename String with filename and desired path, end with .pdf
#' @param height height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#'
#' @export

tsne_perp <- function(z, perp=c(2,5,10,50), sqrt=TRUE, iterations=1000, groups, names,
                      legend=T, legendname=NULL, main=NULL,
                      filename="t-SNE_perplexities.pdf", height=16, useDingbats=F) {

  n=1
    if(sqrt){
    perp <- sort(c(perp, round(sqrt(dim(z)[1]))))
    }

    g <- rep(list(NULL), (length(perp)))

    for(p in 1:length(perp)){
      tsne_tmp = Rtsne(z, check_duplicates=F, perplezity=perp[p], maz_iter=iterations)

      x=tsne_tmp$Y[,1]
      y=tsne_tmp$Y[,2]

      data.df<-data.frame(x=x, y=y, groups=groups, names=names)

      # With lines & points
      g[[n]]<-ggplot(data.df, aes(x=x, y=y, group=groups, color=groups))+
        geom_point(size=1, alpha=0.7, show.legend = legend)+ geom_line(alpha=0.5, show.legend = legend)+
        theme_light()+ scale_colour_manual(values=hue_pal()(length(levels(groups))),
                                           name=legendname)+
        labs(title=main, subtitle=paste0("Perplexity ", perp[p], ", ", iterations, " iterations"), x="t-SNE1", y="t-SNE2")
      n=n+1

      # With lines & labels
      g[[n]]<-ggplot(data.df, aes(x=x, y=y, group=groups, color=groups, label=names))+
        geom_text(size=2, alpha=0.7, show.legend = F)+ geom_line(alpha=0.5, show.legend = legend)+
        geom_point(alpha=0, show.legend = legend)+
        theme_light()+ scale_colour_manual(values=hue_pal()(length(levels(groups))),
                                           name=legendname)+
        guides(colour = guide_legend(override.aes = list(shape=19, size = 1, alpha=0.7)))+
        labs(title=main, subtitle=paste0("Perplexity ", perp[p], ", ", iterations, " iterations"), x="t-SNE1", y="t-SNE2")
      g[[n]] <- ggplot_gtable(ggplot_build(g[[n]]))
      g[[n]]$layout$clip[g[[n]]$layout$name == "panel"] <- "off"
      n=n+1

      # Only points
      g[[n]]<-ggplot(data.df, aes(x=x, y=y, group=groups, color=groups))+
        geom_point(size=1, alpha=0.7, show.legend = legend)+
        theme_light()+ scale_colour_manual(values=hue_pal()(length(levels(groups))),
                                           name=legendname)+
        labs(title=main, subtitle=paste0("Perplexity ", perp[p], ", ", iterations, " iterations"), x="t-SNE1", y="t-SNE2")
      n=n+1

      # Only labels
      g[[n]]<-ggplot(data.df, aes(x=x, y=y, group=groups, color=groups, label=names))+
        geom_text(size=2, alpha=0.7, show.legend = F)+
        geom_point(alpha=0, show.legend = legend)+
        theme_light()+ scale_colour_manual(values=hue_pal()(length(levels(groups))),
                                           name=legendname)+
        guides(colour = guide_legend(override.aes = list(shape=19, size = 1, alpha=0.7)))+
        labs(title=main, subtitle=paste0("Perplexity ", perp[p], ", ", iterations, " iterations"), x="t-SNE1", y="t-SNE2")
      g[[n]] <- ggplot_gtable(ggplot_build(g[[n]]))
      g[[n]]$layout$clip[g[[n]]$layout$name == "panel"] <- "off"
      n=n+1
    }

  lay <- matrix(1:(4*length(perp)), nrow=4)

  ggsave(file=filename,
         grid.arrange(grobs=g, layout_matrix = lay),
         width=4.5*length(perp)+ifelse(legend, length(perp)*2*length(levels(groups))/10, 0),
         height=height, useDingbats=useDingbats)
}

#' Overview of signals in relation to neg control beads
#'
#' Plot overview of singals in relation to the neg control beads (Empty, His6ABP, Neutravidin or Biotin).
#' Based on output from Autoimmunity Profiling scoring function \code{\link[rappp:ap_scoring2]{ap_scoring2()}}.
#'
#' @param x List with at least four elements, see Deatils for naming and content.
#' @param center_prob value in [0,1] passed to prob in \code{\link[stats:quantile]{quantile()}}, defaults to the median.
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @details The x list needs to include at least the element
#'
#'     MFI = assay mfi, column names of negative control columns should include empty|bare|blank|his6abp|hisabp|neutravidin,
#'
#'     SCORE = scored data, column names of negative control columns should include empty|bare|blank|his6abp|hisabp|neutravidin,
#'
#'     CUTOFF_KEY = Cutoff key as data.frame with cutoff values, scores and colors.
#'
#'     BEADS = Beads info, including a column called "Filtered" with the value "NegControl" for the negative control beads.
#'     Any beads with no text (ie. "" or NA) or "NegControl" in such column will be included.
#'
#'     SAMPLES = Sample info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any samples with no text (ie. "" or NA) in such column will be included.
#'
#' Note: The function plots to a layout containing three areas.
#'
#' @export

ap_negbeads <- function(x,
                        center_prob = 0.5,
                        shouldpdf = TRUE,
                        filename = "neg-control-beads.pdf",
                        width = 15,
                        height = 10,
                        useDingbats = FALSE){

  negctrl <- list(Empty="empty|bare|blank|buffer",
                  His6ABP="his6abp|hisabp",
                  Neutravidin="neutravidin",
                  Biotin="Biotin",
                  "Other neg ag"="^neg_M")

  if(shouldpdf){
  pdf(filename, width=width, height=height, useDingbats=useDingbats)
  }

  layout(matrix(c(1,1,1,2,3,4), ncol=3, byrow=T))
  par(mar=c(4,4,4,5))

    plotdata <- x$MFI
    plotdata_score <- x$SCORE

  if("Filtered" %in% colnames(x$BEADS)){
    plotdata <- plotdata[, which(is.na(x$BEADS$Filtered) |
                                x$BEADS$Filtered == "" |
                                grepl("NegControl", x$BEADS$Filtered))]
    plotdata_score <- plotdata_score[, which(is.na(x$BEADS$Filtered) |
                                        x$BEADS$Filtered == "" |
                                        grepl("NegControl", x$BEADS$Filtered))]
  }

  if("Filtered" %in% colnames(x$SAMPLES)){
    plotdata <- plotdata[which(is.na(x$SAMPLES$Filtered) |
                                x$SAMPLES$Filtered == "" |
                                grepl("NegControl", x$SAMPLES$Filtered)),]
    plotdata_score <- plotdata_score[which(is.na(x$SAMPLES$Filtered) |
                                        x$SAMPLES$Filtered == "" |
                                        grepl("NegControl", x$SAMPLES$Filtered)),]
  }

  plotcolor <- x$CUTOFF_KEY
  legend_col <- c(Empty="magenta", His6ABP="chartreuse1",
                  Neutravidin="cyan", Biotin="aquamarine1",
                  "Other neg ag"="black")

  beeswarm(data.frame(t(plotdata)), log=T, corral="gutter", cex=0.5, las=2,
           pwcol=ifelse(grepl(negctrl$Empty, rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T), legend_col["Empty"],
                        ifelse(grepl(negctrl$His6ABP,rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T), legend_col["His6ABP"],
                               ifelse(grepl(negctrl$Neutravidin,rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T), legend_col["Neutravidin"],
                                      ifelse(grepl(negctrl$Biotin,rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T), legend_col["Biotin"],
                                             ifelse(grepl(negctrl$"Other neg ag",rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T), legend_col["Other neg ag"],
                                             as.color(paste(plotcolor$color[t(plotdata_score)*10+1]), 0.4)))))),
           pwpch=rep(ifelse(grepl(paste0(unlist(negctrl), collapse="|"), colnames(plotdata), ignore.case=T), 16, 1), dim(plotdata)[1]),
           ylab="Signal intensity [AU]", cex.axis=0.5)

  legend_text <- names(negctrl)[which(unlist(lapply(negctrl, function(i) sum(grepl(i, colnames(plotdata), ignore.case=T)))) > 0)]

  legend(par("usr")[2], 10^par("usr")[4], xpd=T, cex=0.7,
         legend=c(legend_text, rev(rownames(plotcolor))),
         col=c(legend_col[legend_text], paste(rev(plotcolor$color))),
         pch=c(rep(16, length(legend_text)), rep(1, length(plotcolor$color))))

  par(pty="s")
  for(i in seq_along(negctrl)){
    if(sum(grepl(negctrl[[i]], colnames(plotdata), ignore.case=T)) > 0){
      plot(apply(plotdata, 1, quantile, prob = center_prob, na.rm = T),
           plotdata[,grep(negctrl[[i]], colnames(plotdata), ignore.case=T)],
           las=1,
           xlab=paste0(ifelse(center_prob == 0.5, "Median", paste0(round(center_prob*100, 0), " percentile")), " signal per sample"),
           ylab=paste0(names(negctrl)[[i]]," bead signal"),
           main=paste0(names(negctrl)[i], " bead"),
           col=paste(plotcolor$color[plotdata_score[,grep(negctrl[[i]], colnames(plotdata), ignore.case=T)]*10+1]))
    }
  }

  if(shouldpdf){
  dev.off()
  }
}

#' Calculate reactivity frequencies
#'
#' Calculate number of reactivities and the corresponding frequencies in Autoimmunity Profiling data.
#'
#' @param x List with at least two elements, see Details for naming and content.
#' @param samplegroups factor vector of groupings. Only samples with an assigned level are included in plots.
#'     If left as \code{NULL} (default), the all non-filtered, if filetring done otherwise all, will be assigned "Sample".
#' @param percdec integer indicating the number of decimal places in percentage value.
#' @param check.names logical, altered default from \code{\link[base:data.frame]{data.frame()}}.
#' @details
#'
#' The x list needs to include at least the elements:
#'
#'     SAMPLES = Sample info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any samples with no text (ie. "") in such column will be included.
#'     Column "sample_name" with sample names needed.
#'
#'     BINARY = list with one data.frame per cutoff
#'
#'     Optional element:
#'
#'     BINARY_CO = Binary table based on antigen specific cutoffs
#'     from \code{\link[rappp:ap_cutoff_selection2]{ap_cutoff_selection2()}}.
#'
#' @return A list with the elements
#'
#'     SAMPLEGROUPS = annnotation of which group each sample has been assigned,
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

ap_reactsummary2 <- function(x,
                             samplegroups = NULL,
                             percdec = 1,
                             check.names = FALSE) {

  if("BINARY_CO" %in%  names(x)){
    data_bin <- append(x$BINARY, list(Selected_co=x$BINARY_CO))
  } else {
    data_bin <- x$BINARY
  }

  if(is.null(samplegroups)){

    if("Filtered" %in% colnames(x$SAMPLES)){
      samplegroups <- factor(ifelse(is.na(x$SAMPLES$Filtered) | x$SAMPLES$Filtered == "", "Sample", NA), exclude=c(paste(NA), NA))
    } else {
      samplegroups <- factor(rep("Sample", dim(x$SAMPLES)[1]))
    }

  } else {

    if("Filtered" %in% colnames(x$SAMPLES)){
      samplegroups[which(!(is.na(x$SAMPLES$Filtered) | x$SAMPLES$Filtered == ""))] <- NA
      samplegroups <- factor(samplegroups)
    } # else: samplegroups is already defined, so no need for an else action.

  }

  data_size <- table(samplegroups)
  n_ag <- lapply(data_bin, function(i) apply(i, 1, function(l) sum(!is.na(l))))

  # Calculate per antigen
  data_sum_ag <- lapply(data_bin, function(i) apply(i, 2, function(l) aggregate(l, by=list(samplegroups), FUN=sum)))

  data_freq_ag <- lapply(data_sum_ag,
                         function(cutoff) lapply(cutoff,
                                                 function(antigen) round(antigen$x/data_size*100, digits=percdec)))

  data_sum_ag <- lapply(data_sum_ag,
                        function(cutoff) data.frame(do.call(cbind,
                                                            lapply(cutoff,
                                                                   function(antigen) antigen$x)),
                                                    check.names = check.names))
  data_sum_ag <- lapply(data_sum_ag, function(i) {
    rownames(i) <- levels(samplegroups) ;
    colnames(i) <- colnames(data_sum_ag[[length(data_sum_ag)]]) ; i } )

  data_freq_ag <- lapply(data_freq_ag, function(cutoff) data.frame(do.call(cbind, cutoff),
                                                                   check.names = check.names))
  data_freq_ag <- lapply(data_freq_ag, function(i) {
    colnames(i) <- colnames(data_freq_ag[[length(data_freq_ag)]]) ; i } )

  if(length(levels(samplegroups)) > 1){
    comparisons <-  combn(levels(samplegroups), 2)
    fisher_p <- rep(list(NULL), dim(comparisons)[2])
    freq_diff <- rep(list(NULL), dim(comparisons)[2])

    for(t in 1:dim(comparisons)[2]){
      test_groups <- factor(ifelse(paste(samplegroups) %in%
                                     comparisons[,t], paste(samplegroups), NA))

      tmp_fisher <- matrix(NA, nrow=length(data_bin), ncol=dim(data_bin[[1]])[2])
      tmp_diff <- matrix(NA, nrow=length(data_bin), ncol=dim(data_bin[[1]])[2])

      # Binary based tests
      for(b in 1:length(data_bin)){
        testdata <- data_bin[[b]]

        tmp_fisher[b,] <- apply(testdata,2,
                                function(x) fisher.test(test_groups,
                                                        factor(x,levels=0:1))$p.value)

        tmp_diff[b,] <- unlist(data_freq_ag[[b]][comparisons[1,t],] -
                                 data_freq_ag[[b]][comparisons[2,t],])

      }
      colnames(tmp_fisher) <- colnames(data_freq_ag[[length(data_freq_ag)]])
      rownames(tmp_fisher) <- names(data_bin)
      fisher_p[[t]] <- data.frame(tmp_fisher, check.names=F)


      colnames(tmp_diff) <- colnames(data_freq_ag[[length(data_freq_ag)]])
      rownames(tmp_diff) <- names(data_bin)
      freq_diff[[t]] <- data.frame(tmp_diff, check.names=F)

    }
    names(fisher_p) <- paste0(casefold(comparisons[1,], upper=T), "vs",
                              casefold(comparisons[2,], upper=T))

    names(freq_diff) <- names(fisher_p)
  }

  data_sum_ag <- data.frame(do.call(rbind, data_sum_ag), check.names=F)
  data_freq_ag <- data.frame(do.call(rbind, data_freq_ag), check.names=F)

  # Calculate per sample
  data_sum_samp <- lapply(data_bin,
                          function(i) apply(i, 1,
                                            function(l) sum(l, na.rm=T)))
  names(data_sum_samp) <- names(data_bin)

  data_freq_samp <- lapply(1:length(data_sum_samp),
                           function(cutoff) round(data_sum_samp[[cutoff]]/n_ag[[cutoff]]*100,1))
  names(data_freq_samp) <- names(data_sum_samp)

  data_sum_samp <- data.frame(do.call(cbind, data_sum_samp), check.names=F)
  data_freq_samp <- data.frame(do.call(cbind, data_freq_samp), check.names=F)

  # Output
  if(length(levels(samplegroups)) > 1){
  output <- list(SAMPLEGROUPS=data.frame(Sample=x$SAMPLES$sample_name,
                                         Grouping=paste(samplegroups)),
                 REACTSUM_AG=data_sum_ag,
                 REACTFREQ_AG=data_freq_ag,
                 REACTSUM_SAMP=data_sum_samp,
                 REACTFREQ_SAMP=data_freq_samp,
                 COMPARISONS=comparisons,
                 FISHER_P=fisher_p,
                 FREQ_DIFF=freq_diff)
  } else {
    output <- list(SAMPLEGROUPS=data.frame(Sample=x$SAMPLES$sample_name,
                                           Grouping=paste(samplegroups)),
                   REACTSUM_AG=data_sum_ag,
                   REACTFREQ_AG=data_freq_ag,
                   REACTSUM_SAMP=data_sum_samp,
                   REACTFREQ_SAMP=data_freq_samp)
  }

  return(output)
}

#' Reactivity result plots
#'
#' Plot beeswarm, density and frequency plots for each antigen.
#' Based on output from Autoimmunity Profiling wrapper function \code{\link[rappp:ap_norm2]{ap_norm2()}}.
#'
#' @param x list with at least nine elements, see Deatils for naming and content.
#' @param samplegroups factor vector of groupings. NB! Do not include . (period) in names.
#'     Only samples with an assigned level are included in plots.
#'     If left as \code{NULL} (default), then all non-filtered if filtering has been done,
#'     otherwise all, will be assigned "Sample".
#'     Passed to \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} to calculate frequencies.
#' @param groupcolors A matrix or data.frame with a column named "group" with group names and
#'     a column named "color" with a color for each group (one row per group, ie. factor level).
#'     Alternatively, a vector with colors (will be assigned to the factor levels in order).
#' @param agtoplot indices for which antigens to plot, default is all.
#'     Character vector with column names of what to plot also ok.
#' @param percdec integer indicating the number of decimal places in percentage value.
#' @param cofisher Cutoff in fisher plot.
#' @param shouldpdf Logical, should it plot to pdf?
#' @param filename string with filename and desired path, end with .pdf
#' @param height height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats logical, altered default from \code{\link[grDevices:pdf]{pdf()}}.
#' @param check.names logical, altered default from \code{\link[base:data.frame]{data.frame()}}.
#' @details the x list needs to include at least the element
#'
#'     MADs = assay MADs,
#'
#'     SCORE = scored data
#'
#'     BINARY = list with one data.frame per cutoff,
#'
#'     BINARY_CO = binary table based on antigen specific cutoffs.
#'
#'     ANTIGEN_CUTOFFS = calculated antigen specific cutoffs, translated into the descrete cutoff steps,
#'
#'     CUTOFF_KEY = cutoff key as data.frame with cutoff values, scores and colors.
#'
#'     SAMPLES = sample info. Including column "sample_name" with sample names, preferably LIMS-IDs, where
#'     replicates (named with one of pool|rep|mix|commercial)
#'     and blanks (named with one of empty|blank|buffer) are also stated,
#'     If any wells should be excluded then these should be annotated in a column called "Filtered".
#'     Any samples with no text (ie. "" or NA) in such column will be included.
#'
#'     BEADS = beads info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any beads with no text (ie. "" or NA) will be included in the transformation.
#'
#'     DENS = Density output used for cutoff selection,
#'
#' Note: The function plots to a layout containing up to 16 areas.
#'
#' @return A list with the elements (output from \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}})
#'
#'     SAMPLEGROUPS = annnotation of which group each sample has been assigned,
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

ap_agresults <- function(x,
                         samplegroups=NULL,
                         groupcolors=2:6,
                         agtoplot=NULL,
                         percdec=1,
                         cofisher=0.05,
                         shouldpdf=TRUE,
                         filename="AntigenResults.pdf",
                         height=18,
                         useDingbats=FALSE,
                         check.names=FALSE) {

  print("Calculating frequencies")
  react_summary <- ap_reactsummary2(x,
                                    samplegroups = samplegroups,
                                    percdec = percdec,
                                    check.names = check.names)
  samplegroups <- factor(react_summary$SAMPLEGROUPS$Grouping, exclude=c(paste(NA), NA))
  n_groups <- length(levels(samplegroups))
  data_size <- table(samplegroups)
  if(n_groups > 1){
    n_comparisons <- dim(react_summary$COMPARISONS)[2]
  }

  if(!(class(groupcolors) %in% c("matrix","data.frame"))){
    groupcolors <- data.frame(group=levels(samplegroups),
                              color=groupcolors[seq_along(levels(samplegroups))])
  } else {
    groupcolors <- data.frame(groupcolors, check.names=F)
    groupcolors <- groupcolors[match(paste(levels(samplegroups)), paste(groupcolors$group)), ]
  }

    print("Extract data")
    if("Filtered" %in% colnames(x$BEADS)){
      data_cont <- x$MADS[, which(is.na(x$BEADS$Filtered) | x$BEADS$Filtered == "")]
      data_score <- x$SCORE[, which(is.na(x$BEADS$Filtered) | x$BEADS$Filtered == "")]
      cutoffs <- x$ANTIGEN_CUTOFFS[which(is.na(x$BEADS$Filtered) | x$BEADS$Filtered == ""), ]
    } else {
      data_cont <- x$MADS
      data_score <- x$SCORE
      cutoffs <- x$ANTIGEN_CUTOFFS
    }

    cokey <- x$CUTOFF_KEY
    if_selected_co <- sum(grepl("Selected_co", rownames(react_summary$REACTSUM_AG))) > 0

    if(if_selected_co){
      data_sum <- react_summary$REACTSUM_AG[grep("Selected_co", rownames(react_summary$REACTSUM_AG)), ]
      data_freq <- react_summary$REACTFREQ_AG[grep("Selected_co", rownames(react_summary$REACTFREQ_AG)), ]
      data_freq_all <- react_summary$REACTFREQ_AG[-grep("Selected_co", rownames(react_summary$REACTFREQ_AG)), ]
    } else {
      data_freq_all <- react_summary$REACTFREQ_AG
    }

    print("set agtoplot")
    if(is.null(agtoplot)){
      agtoplot <- 1:dim(data_cont)[2]
    } else if(is.character(agtoplot)){
      agtoplot <- match(agtoplot, colnames(data_cont))
    }

    print("initiate pdf")
      # Create PDF
    if(shouldpdf){
     pdf(filename,
          width=ifelse(n_groups > 1, 20+n_groups*1, 15), height=height, useDingbats=useDingbats)
    }
      mar_top <- ifelse(n_groups > 1, ceiling((n_groups+1)/3)+3, 5)
      mtext_sub_line <- ifelse(n_groups > 1, ceiling((n_groups+1)/3), 2)
      par(mgp=c(3,1,0))

      if(n_groups > 1){
        layout(matrix(1:ifelse(length(agtoplot) > 4, 16, length(agtoplot)*4),
                      nrow=ifelse(length(agtoplot) > 4, 4, length(agtoplot)),
                      byrow=T))
      } else {
        layout(rbind(c(1,2,2,3,3),
                     t(sapply(seq(3,9,3), function(x)
                       c(1,2,2,3,3)+x)))[1:ifelse(length(agtoplot) > 4, 4, length(agtoplot)),])

      }

      n=1
    for(a in agtoplot){
      tmp_ag <- colnames(data_cont)[a]
      print(paste("Plotting ag", n, "of", length(agtoplot),"(",tmp_ag,")"))
      n=n+1

      if(if_selected_co){
      dens <- x$DENS[[tmp_ag]]
      tmp_which_co <- cutoffs$score[which(cutoffs$bead == tmp_ag)]*10+1
      tmp_cutoff <- cokey$xmad[tmp_which_co]
      }

      # MADs Beeswarm, antigen score coloring
      par(mar=c(6,5,mar_top,3))
      plotdata <- data_cont[,tmp_ag]
      boxplot(plotdata~samplegroups, col="lightgrey", outcol=0, las=2,
              ylab="MADs [AU]", xaxt="n", xlab=NA,
              ylim=c(min(plotdata, na.rm=T), ifelse(max(plotdata, na.rm=T) > 50, max(plotdata, na.rm=T), 50)))
      if(if_selected_co){ abline(h=tmp_cutoff, lty=2) }

      pwcols <- paste(cokey$color[data_score[,tmp_ag]*10+1])
      pwcols <- gsub("NA","#ffffff",pwcols)
      beeswarm(plotdata~samplegroups, pch=16, corral="gutter", corralWidth=0.5, cex=0.8, add=T,
               pwcol=as.color(pwcols, 0.8))

      legend(par("usr")[2], par("usr")[4], legend=rev(c("<0",cokey$xmad[-1])),
             title=expression(bold("MADs cutoff")),
             pch=16, cex=0.6, bty="n", xjust=0.2, title.adj=4,
             col=rev(paste(cokey$color)), xpd=NA)
      mtext("Visualization of signals.",
            line=mtext_sub_line,
            cex=0.65)

      if(if_selected_co){
      mtext(paste0("Above dashed line:\n(",tmp_cutoff," MADs)"), adj=0.7,
            side=1, at=par("usr")[1], line=0.6, cex=0.5)

        axis_text <- paste0(data_sum[, grep(paste0("^",tmp_ag,"(?=_co)"), colnames(data_sum), perl = T)], "/", data_size,
                            " (", round(data_freq[, grep(paste0("^",tmp_ag,"(?=_co)"), colnames(data_freq), perl = T)],
                                        digits=percdec), "%)\n",
                            levels(samplegroups))
        mtext(axis_text,
              side=1, at=1:n_groups, line=0.7+0.7*max(str_count(axis_text, "\\\n")), cex=0.7)

      } else {
        mtext(levels(samplegroups),
              side=1, at=1:n_groups, line=1.1, cex=0.7)
      }
      mtext(tmp_ag,
            line=mtext_sub_line+1,
            font=2)

      # Histrogram & Density
      par(mar=c(6,5,mar_top,1))
      h <- hist(data_score[,tmp_ag], breaks=seq(min(cokey$score)-0.1,max(cokey$score)+0.1, 0.1), prob=T, right=F,
                main=NA, xlim=c(-0.1, max(cokey$score)+0.1), xlab="MADs cutoff\nDensity bandwidth = 0.1", xaxt="n")
      axis(1, labels=c("<0",cokey$xmad[-1]), at=h$breaks[-c(1, length(h$breaks))], cex.axis=0.8)

      if(if_selected_co){
        mtext("Distribution of binned values, algorithm assigned cutoff at dashed line.",
              line=mtext_sub_line,
              cex=0.65)
        abline(v=(tmp_which_co-1)/10, lty=2)
        lines(dens,
              col="maroon")
      } else {
        mtext("Distribution of binned values", line=0.1, cex=0.65)
      }
        mtext(tmp_ag,
              line=mtext_sub_line+1,
              font=2)

      # Frequency
        par(mar=c(6,5,mar_top,1))
        plotdata <- data_freq_all[,grep(gsub("\\*", "\\\\*", paste0("^",tmp_ag,"(?=_co)")),
                                        colnames(data_freq_all), perl = T), drop=F]
      if(n_groups > 1){
         plotdata <- split(plotdata, do.call(rbind, strsplit(rownames(plotdata), "\\."))[,2])
         plotdata <- do.call(cbind, plotdata)
         colnames(plotdata) <- names(data_size)
      }
      plot(NULL, xlim=c(0,dim(cokey)[1]),
           ylim=c(0,100), xaxt="n", yaxt="n",
           ylab="Reactivity frequency [%]", xlab="MADs cutoff")
      axis(2, at=seq(0,100,10), labels=seq(0,100,10), las=1)
      axis(1, at=1:dim(cokey)[1], labels=c("<0",cokey$xmad[-1]), cex.axis=0.8)
      abline(h=seq(0,100,10), col="lightgrey", lty=2)
      abline(v=1:dim(cokey)[1], col="lightgrey", lty=2)
      if(if_selected_co){ abline(v=tmp_which_co, lty=2) }

      matplot(plotdata, type="l", lty=1:5, lwd=2, add=T,
              col=paste(groupcolors$color))

      if(n_groups > 1){
        mtext("Percentage of reactive samples per group at each exemplified cutoff.",
              line=mtext_sub_line, cex=0.65)
        legend(par("usr")[1], par("usr")[4], yjust=0, xpd=NA, bty="n", cex=0.8, ncol=3,
               lty=1:5, lwd=2, col=paste(groupcolors$color),
               legend=gsub("\\\n", "_", groupcolors$group), seg.len=5)
        mtext(tmp_ag, line=mtext_sub_line+1, font=2)
      } else {
        mtext("Percentage of reactive samples at each exemplified cutoff.", line=0.1, cex=0.65)
        mtext(tmp_ag, line=2, font=mtext_sub_line)
      }

      # Fisher's exact test
      if(n_groups > 1){
        par(mar=c(6,2,
                  mar_top,
                  ceiling(0.03572519*max(unlist(lapply(names(react_summary$FISHER_P), nchar)), na.rm=T)*19.3) )) # strwidth("M") = 0.03572519

      plotdata <- melt(map(react_summary$FISHER_P, as.matrix)) #melt(react_summary$FISHER_P)
      plotdata <- plotdata[which(unlist(lapply(strsplit(paste(plotdata$Var2), "_co"), function(i) i[[1]])) == tmp_ag),]
      if(sum(grepl("Selected", plotdata$Var1)) > 0){
        plotdata <- plotdata[-grep("Selected", plotdata$Var1), ]
      }

      plotdata <- matrix(plotdata$value, nrow=dim(x$CUTOFF_KEY)[1], dimnames=list(unique(plotdata$Var1), unique(plotdata$L1)))
      if(sum(plotdata > 1) > 0){
      plotdata[which(plotdata > 1, arr.ind=T)] <- 1 # Fix bug with floating aritmetics, some 1s will not be recognized as 1s in image (white field).
      }
      plotdata <- plotdata[,dim(plotdata)[2]:1, drop=F]

      breaks <- sort(c(1, 0.1, 0.05, 10^-seq(2,4,1), 0))
      break_col <- hcl.colors(length(breaks)-1, "BrwnYl", rev = F)
      image(x=1:dim(plotdata)[1], y=1:dim(plotdata)[2], z=plotdata, xaxt="n", yaxt="n", bty="n",
            col=break_col, ylab=NA, xlab="MADs cutoff", breaks=breaks)
      axis(1, at=1:dim(cokey)[1], labels=c("<0",cokey$xmad[-1]), cex.axis=0.8, tick=F)
      axis(4, at=1:dim(plotdata)[2], gsub("\\\n", "_", colnames(plotdata)), las=2, tick=F)
      legend(par("usr")[1], par("usr")[4], pt.bg=rev(break_col), pch=22, horiz=T, pt.cex=2.5,
             legend=paste0("<=",#expression("\u2264"),
                           format(rev(breaks[-1]), scientific=F, drop0trailing=T)),
             bty="n", yjust=0.2, xpd=NA, cex=1)
      abline(h=c(0.5, (1:dim(plotdata)[2])+0.5), lty=2, col="lightgrey")
      abline(v=c(0.5, (1:dim(plotdata)[1])+0.5), lty=2, col="lightgrey")
      if(if_selected_co){ abline(v=tmp_which_co+c(0.5,-0.5), lty=2) }

        mtext("Fisher's exact test p-values per pariwise group comparison at each exemplified cutoff.",
              line=mtext_sub_line, cex=0.65, adj=0)
          mtext(tmp_ag, line=mtext_sub_line+1, font=2)
      }

    }
      if(shouldpdf){
    dev.off()
      }
    return(react_summary)
  }

#' Analysis summary
#'
#' Summarizes number of filtered/flagged beads/samples, amino acid lenghts and protein representation.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @param filter Logical, should samples and beads be filtered away before summary?
#' @details The x list needs to include at least the elements
#'
#'     SAMPLES = Sample info. Including column "sample_name" with LIMS-IDs, and "Filtered" if argument filter is TRUE.
#'
#'     BEADS = Beads info. Including columns:
#'
#'     "Type" with info about type of content on bead, at least including "PrEST" for PrESTs,
#'
#'     "PrEST.seq..aa." with amino acid sequences,
#'
#'     "Filtered" with filtering annotation, e.g. from other ap_-functions,
#'
#'     "Flagged" with filtering annotation, e.g. from other ap_-functions.
#'
#' @export

ap_summary <- function(x, filter=TRUE) {

  if(filter & "Filtered" %in% colnames(x$SAMPLES) & "Filtered" %in% colnames(x$BEADS)){
    samples <- x$SAMPLES[which(is.na(x$SAMPLES$Filtered) | x$SAMPLES$Filtered == ""),]
    beads <- x$BEADS[which(is.na(x$BEADS$Filtered) | x$BEADS$Filtered == ""),]
  } else {
    beads <- x$BEADS
    samples <- x$SAMPLES
  }

  summarylist <- list(
    # Cohorts
    "Samples per cohort" = table(matrix(unlist(strsplit(as.character(samples$sample_name),"-")),
                                        ncol=2, byrow=T)[,1]),
    # Uniprot IDs
    "Table of Uniprot IDs" = sort(table(unlist(strsplit(as.character(beads$Uniprot[which(beads$Type == "PrEST")]),
                                                        ";|,")))),
    "Unique Uniprot IDs" = length(table(unlist(strsplit(as.character(beads$Uniprot[which(beads$Type == "PrEST")]),
                                                        ";|,")))),
    # Ensembl IDs
    "Table of Ensembl IDs" = sort(table(unlist(strsplit(as.character(beads$ENSG.ID[which(beads$Type == "PrEST")]),
                                                        ";|,")))),
    "Unique Ensembl IDs" = length(table(unlist(strsplit(as.character(beads$ENSG.ID[which(beads$Type == "PrEST")]),
                                                        ";|,")))),
    # Gene names
    "Table of Genes" = sort(table(unlist(strsplit(as.character(beads$Gene.name[which(beads$Type == "PrEST")]),
                                                  ";|,")))),
    "Unique Genes" = length(table(unlist(strsplit(as.character(beads$Gene.name[which(beads$Type == "PrEST")]),
                                                  ";|,")))),
    # Aminoacids
    "Table of sequence lenghts" = sort(apply(as.matrix(beads$PrEST.seq..aa.[
      grep("PrEST",beads$Type, ignore.case=T)],ncol=1), 1, nchar)), # Check aa-sequence length range
    "min sequence lenghts" = min(apply(as.matrix(beads$PrEST.seq..aa.[
      grep("PrEST",beads$Type, ignore.case=T)],ncol=1), 1, nchar)), # Check aa-sequence length range
    "max sequence lenghts" = max(apply(as.matrix(beads$PrEST.seq..aa.[
      grep("PrEST",beads$Type, ignore.case=T)],ncol=1), 1, nchar)), # Check aa-sequence length range
    "median sequence lenghts" = median(apply(as.matrix(beads$PrEST.seq..aa.[
      grep("PrEST",beads$Type, ignore.case=T)],ncol=1), 1, nchar)) # Check aa-sequence median length
  )

    # Filtered & Flagged summary
  beads <- x$BEADS
  samples <- x$SAMPLES
  if("Filtered" %in% colnames(beads) & "Flagged" %in% colnames(beads) & "Filtered" %in% colnames(samples)){
    summarylist <- append(summarylist,
                          list("Filtered antigens" = table(beads$Filtered),
                               "Flagged antigens" = table(beads$Flagged),
                               "Filtered samples" = table(samples$Filtered)))
  }

  return(summarylist)
}

#' Export data from Autoimmunity Profiling analysis to excel
#'
#' Exports chosen list elemnts to sheets in an Excel-file.
#' Default elements are based on output from \code{\link[rappp:ap_norm2]{ap_norm2()}}, and
#'  \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} or \code{\link[rappp:ap_agresults]{ap_agresults()}}
#'
#' @param x list with at least the elements to export, see Deatils for more information and exceptions to the rule.
#' @param elements character vector of names of the list elements to export.
#'     The sheets will be ordered in the same order as the vector.
#'     Specific MAD-cutoffs are not included in the default output but can be added, eg. "25xMAD".
#' @param filename string with filename and desired path, end with .xlsx.
#' @param shouldround logical, if TRUE MFI values are rounded to integers and MADs to two decimals.
#' @param row.names logical. If TRUE, the row names of the data frames are included in the Excel file worksheets.
#'     Deafult altered from \code{\link[WriteXLS:WriteXLS]{WriteXLS()}}.
#' @param save_rdata logical. If TRUE, an RData-file is saved, containg a list matching the Excel-file content.
#' @param ... arguments passed to \code{\link[WriteXLS:WriteXLS]{WriteXLS()}}.
#' @details The x list needs to include at least the elements specified under \code{elements}.
#'   It is recommended to append the output from \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} or
#'   \code{\link[rappp:ap_agresults]{ap_agresults()}} to the output from \code{\link[rappp:ap_norm2]{ap_norm2()}}
#'   and use the combined list as function input.
#'
#'   Exceptions for input element names:\cr
#'   - If an element is named BEADS in the input data, its name will be changed to ANTIGENS.
#'   Therefore, ANTIGENS is a default element to export while BEADS is not.\cr
#'   - Sums and frequencies at the cutoffs assigned by \code{\link[rappp:ap_cutoff_selection2]{ap_cutoff_selection2()}}
#'   will be combined into a new element called REACTIVITIES
#'  if the output from \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} or
#'  \code{\link[rappp:ap_agresults]{ap_agresults()}} is included in the output
#'  Therefore, REACTIVITIES is a default element to export although it is not present in the input data.\cr
#'  - If samplegorups were provided in \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} the
#'  resulting Fisher's exact test p-values and frequency differences will be extracted for the
#'  the cutoffs assigned by \code{\link[rappp:ap_cutoff_selection2]{ap_cutoff_selection2()}} and formatted
#'  into the new elements FISHER and DIFFERENCES.
#'  Therefore, FISHER and DIFFERENCES are default elements to export although they are not present in the input data.\cr
#'
#'  Line breaks (backslash n) in row and column names will be changed to underscore ("_").
#'
#'   If other list structures are used, it is most likely more convenient to just use
#'   \code{\link[WriteXLS:WriteXLS]{WriteXLS()}}, which this function is built on.
#'
#' @export

ap_excel <- function(x,
                     elements = c("MFI", "MADS", "SCORE", "BINARY_CO",
                                "REACTIVITIES", "FISHER", "DIFFERENCES",
                                "ANTIGEN_CUTOFFS", "CUTOFF_KEY",
                                "ANTIGENS", "SAMPLES", "COUNT"),
                     filename = "DataOutput.xlsx",
                     shouldround = TRUE,
                     row.names = TRUE,
                     save_rdata = TRUE,
                     ...) {
  excel <- x

  ## Handling basic elements in ap SBA structure
  if("BEADS" %in% names(excel)){
  names(excel)[which(names(excel) == "BEADS")] <- "ANTIGENS"
  }

  if(shouldround){
    if("MFI" %in% names(excel)){
      excel$MFI <- data.frame(apply(excel$MFI, 2, function(x) {class(x) <- "numeric" ; x } ), check.names=F)
      excel$MFI <- round(excel$MFI, 0)
    }
    if("MADS" %in% names(excel)){
      excel$MADS <- data.frame(apply(excel$MADS, 2, function(x) {class(x) <- "numeric" ; x } ), check.names=F)
      excel$MADS <- round(excel$MADS, 2)
    }
  }

  ## Extract specific binary matrices
  if(sum(grepl("xMAD", elements)) > 0 & "BINARY" %in% names(excel)){
    for(i in grep("xMAD", elements)){
      excel <- append(excel, excel$BINARY[elements[i]])
    }
  }

 ## Handling elements from ap_reactsummary2() if inluded
  if("REACTSUM_AG" %in% names(excel)){
  tmp_sum <- excel$REACTSUM_AG[grep("Selected", rownames(excel$REACTSUM_AG)),]
  rownames(tmp_sum) <- gsub("Selected_co","Sum", rownames(tmp_sum))
  }

  if("REACTSUM_AG" %in% names(excel)){
  tmp_freq <- excel$REACTFREQ_AG[grep("Selected", rownames(excel$REACTFREQ_AG)),]
  rownames(tmp_freq) <- gsub("Selected_co","Frequency", rownames(tmp_freq))
  }

  if("FISHER_P" %in% names(excel)){
    tmp_fisher <- melt(map(excel$FISHER_P, as.matrix))
    tmp_fisher <- tmp_fisher[grep("Selected", tmp_fisher$Var1),]
    tmp_fisher <- matrix(tmp_fisher$value, ncol=length(levels(tmp_fisher$Var2)), byrow=T,
                         dimnames=list(paste(unique(tmp_fisher$L1)), paste(unique(tmp_fisher$Var2))))
    rownames(tmp_fisher) <- gsub("Selected_co","", rownames(tmp_fisher))
  }

  if("FREQ_DIFF" %in% names(excel)){
    tmp_diff <- melt(map(excel$FREQ_DIFF, as.matrix))
    tmp_diff <- tmp_diff[grep("Selected", tmp_diff$Var1),]
    tmp_diff <- matrix(tmp_diff$value, ncol=length(levels(tmp_diff$Var2)), byrow=T,
                         dimnames=list(paste(unique(tmp_diff$L1)), paste(unique(tmp_diff$Var2))))
    rownames(tmp_diff) <- gsub("Selected_co","", rownames(tmp_diff))
  }

  if(exists("tmp_sum") & exists("tmp_freq")){
  excel <- append(excel,
                  list(REACTIVITIES=data.frame(t(tmp_sum),t(tmp_freq))))
  }

  if(exists("tmp_fisher") & exists("tmp_diff")){
    excel <- append(excel,
                    list(FISHER=data.frame(tmp_fisher, check.names=F),
                         DIFFERENCES=data.frame(tmp_diff, check.names=F)))
  }

  ## Final formatting
  excel <- excel[match(elements,
                       names(excel))]

  for(i in seq_along(excel)){
    if(class(excel[[i]])[1] %in% c("data.frame", "matrix")){
      rownames(excel[[i]]) <- gsub("\\\n", "_", rownames(excel[[i]]))
      colnames(excel[[i]]) <- gsub("\\\n", "_", colnames(excel[[i]]))
    }
  }

  ## Write to file
  if(save_rdata){
    save(excel, file=gsub("\\.xls|\\.xlsx", "\\.RData", filename))
  }

  openxlsx::write.xlsx(x=excel,
                       file = filename,
                       rowNames = row.names, ...)
}

#' Cutoff key image
#'
#' Create a cutoff key for scoring of Autoimmunity Profiling data and
#' also produce an image.
#' Uses \code{\link[rappp:ap_cutoffs2]{ap_cutoffs2()}}
#'
#' @param cutoffkey table matching output from \code{\link[rappp:ap_cutoffs2]{ap_cutoffs2()}},
#'     recommended input if \code{\link[rappp:ap_norm2]{ap_norm2()}} has been used.
#' @param MADlimits vector of MADs values used as boundaries for binning (≥MADs), eg. seq(0,70,5).
#'     Not used if cutoffkey is provided.
#' @param shouldpdf Logical, should it plot to png?
#' @param filename string with filename and desired path, end with .png
#' @return If MADlimits is provided a data.frame with three columns will be returned:
#'
#'    [,1] MADs cutoff value
#'
#'    [,2] Corresponding score value
#'
#'    [,3] Corresponding color using the Zissou1 palette in \code{\link[wesanderson]{wes_palette}}
#' @export

ap_cutoffs2image <- function(cutoffkey = NULL,
                             MADlimits = NULL,
                             shouldpng = TRUE,
                             filename = "CutoffColorKey.png") {

  if(!is.null(cutoffkey)){
    xmad_score <- cutoffkey
  } else if(!is.null(MADlimits)){
    xmad_score <- ap_cutoffs2(MADlimits = MADlimits)
  } else {
    warning("Either cutoffkey or MADlimits has to be provided.")
  }

  if(shouldpng){
    png(filename=filename)
  }

  par(mar=c(4,6,4,1))
  plot(x=xmad_score$score, y=rep(1,dim(xmad_score)[1]),
       ylim=c(0.9,1.1), xlim=c(min(xmad_score$score)-0.2, max(xmad_score$score)),
       xlab=NA, ylab=NA, xaxt="n", yaxt="n", frame.plot=F,
       pch=22, bg=paste(xmad_score$color), col="black", cex=4)
  text(x=min(xmad_score$score)-0.2, y=1, xpd=NA,
       labels="Color", font=2, offset=0, cex=1, adj=1)
  # Add score boxes
  points(x=xmad_score$score, y=rep(1.015,dim(xmad_score)[1]),
         pch=22, col="black", cex=4)
  textxy(X=xmad_score$score,
         Y=rep(1.015,dim(xmad_score)[1]),
         labs=sprintf("%.1f",xmad_score$score), offset=0, cex=0.8)
  text(x=min(xmad_score$score)-0.2, y=1.015, xpd=NA,
       labels="Score", font=2, offset=0, cex=1, adj=1)
  # Add xMAD boxes
  points(x=xmad_score$score, y=rep(1.03,dim(xmad_score)[1]),
         pch=22, col="black", cex=4)
  textxy(X=xmad_score$score,
         Y=rep(1.03,dim(xmad_score)[1]),
         labs=c("<0", xmad_score$xmad[-1]), offset=0, cex=0.8)
  text(x=min(xmad_score$score)-0.2, y=1.03, xpd=NA,
       labels="MADs cutoff (\u2265)", font=2, offset=0, cex=1, adj=1)

  if(shouldpng){
    dev.off()
  }

  if(is.null(cutoffkey)){
    return(xmad_score)
  }
}

#' Make peptides
#'
#' Get possible peptide designs across different lengths and overlaps.
#'
#' @param sequence string with amino acid sequence to make peptides from,
#'     see Examples for how to use with the Clipboard.
#' @param lengths number of amino acids per peptide, can be a vector of lenghts.
#' @param min_overlap the shortest allowed overlap between peptides
#' @param min_shift the shortest allowed shift between peptides
#' @return A list with one element per design
#' @examples
#' ## The input sequence can easily be copied from wherever,
#'     like Excel, LIMS etc and used directly in the function:
#' ## on PC
#' make_peptides(readClipboard())
#'
#' ## on MAC
#' make_peptides(t(read.table(pipe("pbpaste"), stringsAsFactors = FALSE)))
#'
#' @export

make_peptides <- function(sequence,
                          lengths = 15:20,
                          min_overlap = 2,
                          min_shift = 1) {

  if(length(sequence) > 1){
  sequence <- paste0(sequence,collapse="")
  }
  sequence <- substring(sequence, seq(1,nchar(sequence),1), seq(1,nchar(sequence),1))

  lengths <- lengths
  min_overlap <- min_overlap
  min_shift <- min_shift

  ## Search for possible length & overlap with no NAs
  {
    overlaps <- lapply(lengths, function(i) min_overlap:(i-min_shift))
    peptides_list <- rep(list(NA), length(unlist(overlaps)))
    m=1
    for(l in seq_along(lengths)){
      aa_length <- lengths[l]

      for(o in overlaps[[l]]){
        aa_overlap <- o

        start_aas <- seq(1, length(sequence), aa_length-aa_overlap)
        if(start_aas[which.max(start_aas)-1] >= (length(sequence)-aa_length+1)){
          start_aas <- start_aas[-which(start_aas %in% (length(sequence)-aa_length+1):length(sequence))[-1]]
        }

        peptides <- matrix(NA, nrow=length(start_aas), ncol=aa_length)
        n=1
        for(i in start_aas){
          peptides[n,] <- sequence[i:(i+aa_length-1)]
        n=n+1
          }

          peptides_list[[m]] <- peptides
          names(peptides_list)[m] <- paste0("Length",aa_length,"_Overlap",aa_overlap,"_N", dim(peptides)[1])
          m=m+1

      }
    }

    peptides_list <- peptides_list[-which(unlist(lapply(peptides_list, function(i) sum(is.na(i)))) > 0)]
    peptides_collapsed_noNA <- lapply(peptides_list, function(x) apply(x, 1, function(y) paste0(y, collapse="")))
  }
  return(peptides_collapsed_noNA)
}

#' Check peptides
#'
#' Get a summary of risk factors for reduced peptide synthesis yield.
#'
#' @param to_check a vector of sequences to check.
#' @param map_to an optional sequence used to extract amino acid indices of peptides.
#' @return A data.frame with various summaries.
#' @details Information from Sigma Aldrich (Merck) about peptide design:\cr
#' https://www.sigmaaldrich.com/technical-documents/articles/biology/designing-peptides.html
#'
#' Amino Acid Classifications:\cr
#' Hydrophobic (non-polar): Ala, Ile, Leu, Met, Phe, Trp, Val\cr
#' Uncharged (polar): Asn, Cys, Gly, Gln, Pro, Ser, Thr, Tyr\cr
#' Acidic (polar): Asp, Glu\cr
#' Basic (polar): His, Lys, Arg\cr
#' TIP – Keep hydrophobic amino acid content below 50% of the total sequence length and
#' to include at least one charged amino acid for every five amino acids.
#' At a physiological pH, Asp, Glu, Lys and Arg will contain charged side chains.
#' A single conservative replacement, such as replacing Ala with Gly or adding polar amino acids
#' to the N- or C-terminus may improve solubility.
#'
#' There are several strategies for improving peptide stability, which will lead to higher
#' purity and optimal solubility. Amino acid composition of the peptide sequence impacts the
#' overall stability and considerations should be made for the following scenarios:
#'
#' 1. Multiple Cys, Met or Trp amino acids may be difficult to obtain in high purity partly
#' due to the susceptibility of oxidation and/or side reactions.\cr
#' TIP – Choose sequences which minimize these residues or choose conservative replacements
#' for these amino acids.  Norleucine can substitute for Met and Ser can  be  a less reactive
#' replacement for Cys.  If overlapping peptides from a protein sequence are being designed,
#' shifting the starting point of each peptide may also create a better balance between
#' hydrophobic and hydrophilic amino acid residues.
#'
#' 2. N-terminal Gln (Q) is unstable and may cyclize to pyroglutamate when exposed to the
#' acidic conditions of cleavage.\cr
#' TIP – Amidate the N-terminus of the sequence or substitute this amino acid.
#'
#' 3. Asparagine (N) has a protecting group that is difficult to remove when placed at the
#' N-terminus of a peptide sequence.\cr
#' TIP – Remove the Asn at this location, substitute with another amino acid or lengthen the
#' peptide by one amino acid residue.
#'
#' 4. Multiple prolines (P) or adjacent serines (S)  in a sequence can result in a product that
#' is lower in purity or contains many
#' deletion products. Multiple prolines can also undergo a cis/trans isomerization, resulting in
#' an apparent lower purity product.
#'
#' 5. Beta sheet formation is a concern as it causes incomplete solvation of the growing peptide
#' chain and will result in a higher incidence of deletion sequences in the final product. \cr
#' TIP – Avoid sequences that contain multiple or adjacent Val, Ile, Tyr, Phe, Trp, Leu, Gln and Thr.
#' Break the pattern by making conservative replacements, for example, inserting a Gly or Pro at every
#' third residue, replacing Gln with Asn, or replacing Thr with Ser.
#'
#' @export

check_peptides <- function(to_check, map_to=NULL){

  to_check_list <- lapply(as.character(to_check), function(x)
    substring(x, seq(1,nchar(x),1), seq(1,nchar(x),1)))

  length_aa <- as.vector(sapply(to_check, nchar))

  N_terminal <- unlist(lapply(to_check_list, function(x) ifelse(x[1] %in% c("Q","N"), "Yes", "No")))

  beta_sheet <- unlist(lapply(to_check_list, function(x) sum(x %in% c("V","I","Y","F","W","L","Q","T"))))

  P_and_S <- unlist(lapply(to_check_list, function(x) sum(x %in% c("P","S"))))

  oxidation_risk <- unlist(lapply(to_check_list, function(x) sum(x %in% c("M","C","W"))))

  hydrophobic <- unlist(lapply(to_check_list, function(x) sum(x %in% c("A", "I", "L", "M", "F", "W", "V"))))

  acidic <- unlist(lapply(to_check_list, function(x) sum(x %in% c("D", "E"))))

  charged <- unlist(lapply(to_check_list, function(x) sum(x %in% c("D", "E","L","R"))))

  check_summary <- data.frame(Sequences=to_check,
                              Length=length_aa,
                              Risky_N_terminal=N_terminal,
                              Hydrophobic_sum=hydrophobic,
                              Hydrophobic_perc=round(hydrophobic/length_aa*100),
                              # Hydrophobic_risky=ifelse(hydrophobic/length_aa > 0.45, "Yes", "No"),
                              Acidic_sum=acidic,
                              Acidic_perc=round(acidic/length_aa*100),
                              # Acidic_risky=ifelse(acidic/length_aa > 0.5, "Yes", "No"),
                              Charged_sum=charged,
                              Charged_perc=round(charged/length_aa*100),
                              Beta_sheet_aa_sum=beta_sheet,
                              Beta_sheet_aa_perc=round(beta_sheet/length_aa*100),
                              P_and_S_sum=P_and_S,
                              P_and_S_perc=round(P_and_S/length_aa*100),
                              Ox_risk_sum=oxidation_risk,
                              Ox_risk_perc=round(oxidation_risk/length_aa*100))

  if(!is.null(map_to)){
    aa_index <- unlist(lapply(to_check_list, function(x)
      paste0("aa",which(rollapply(map_to, length(x), identical, x)), "-",
             (which(rollapply(map_to, length(x), identical, x))+length(x)-1))))

    check_summary <- add_column(check_summary, Mapping = aa_index, .after = 2)
  }

  return(check_summary)

}


#' Align sequences
#'
#' Visualize how sequences align to a full length protein (exact matches).
#'
#' @param inputfile file with sequence information, see Details for needed formatting and columns.
#' @param shouldpdf Logical, should it plot to pdf?
#' @param shouldpdf Logical, should it plot to pdf?
#' @param outputfile filename for the output file.
#' @param gene string with genename (or other name-identifier) for the alignment.
#' @param uniprot string with Uniprot ID for the alignment.
#' @param writesequence logical, should aa sequences for peptides be stated on the right axis?
#' @param shouldtextplot logical, should a table with sequence info be plotted below the alignment?
#' @param textplotcolumns columns to include if printing a table below the alignment.
#' @details NB! The sequences have to be identical matches to consequtive stretches
#'     in the full length sequence you are aligning to.
#'
#'     The input file needs to be tab-delimited, have one row per peptide/protein and the columns\cr
#'     - Name (identifier to be written above the sequence)\cr
#'     - Type ("Full_length" for sequence to align against,
#'             "Peptide" for peptides (sequence will be printed on the side),
#'             any other label can be used but will not be considered in the function)
#'     - Sequence1 (1 in column name also if only one sequence per row)\cr
#'     - Sequence2 (only if mosaic, not necessary if only one sequence per row)\cr
#'     - Color (color for the sequence bar)\cr
#'     - Include (values 1 or 0, which sequences in input file should be included in plot?)
#'
#' If printing table underneath extra columns may be needed in the input file
#'     depending on what information you want to print.
#'
#' Note: The function plots to a layout containing up to two areas.
#'
#' @export

align_sequences <- function(inputfile,
                            shouldpdf = TRUE,
                            outputfile = "SequenceAlignment.pdf",
                            gene = NULL,
                            uniprot = NULL,
                            writesequence = TRUE,
                            shouldtextplot = FALSE,
                            textplotcolumns = NULL) {

  sequences <- read.delim(inputfile, na.strings="")
  all <- lapply(as.character(sequences$Sequence1), function(x) substring(x, seq(1,nchar(x),1), seq(1,nchar(x),1)))
  names(all) <- sequences$Name
  FLlength <- length(all[[which(sequences$Type == "Full_length")]])

  if(sum(grepl("Sequence2", colnames(sequences))) > 0){
    for(i in 1:dim(sequences)[1]){
      if(!is.na(sequences$Sequence2[i])){
        tmp <- as.character(sequences$Sequence2[i])
        all[[i]] <- list(all[[i]], substring(tmp, seq(1,nchar(tmp),1), seq(1,nchar(tmp),1)))
      }
    }
    rm(tmp)
  }
  all <- all[which(sequences$Include == 1)]
  sequences <- sequences[which(sequences$Include == 1),]

  if(shouldpdf){
  pdf(file=outputfile,
      height=ifelse(length(all)*0.5 < 3, 5, length(all)*0.5),
      width=25, useDingbats=F)
  }

  if(shouldtextplot){
    layout(matrix(c(1,1,2), ncol=1))
  } else {
    layout(matrix(1, ncol=1))
  }

  if(writesequence){
  par(lend="butt", mar=c(2,2,2,25))
  } else {
    par(lend="butt", mar=c(2,2,2,2))
  }
  y_at <- seq(0, length(all)*0.5-0.5, 0.5)
  aasteps <- ifelse(FLlength < 200, 10, ifelse(FLlength < 500, 50, 100))
  plot(0,0, col=0, xlim=c(1,FLlength), ylim=c(0,sum(sequences$Include)*0.5+0.5), xaxt="n",xlab=NA, las=1, yaxt="n", ylab=NA)
  axis(1, at=seq(1, FLlength, aasteps),
       labels=as.character(seq(1, FLlength, aasteps)),
       las=1)
  abline(v=seq(1, FLlength, aasteps), lty=2, col="lightgrey")

  if(writesequence){
    ylabs <- ifelse(sequences$Type == "Peptide", paste0(sequences$Sequence1), "")
    axis(4, at=y_at, labels=ylabs, las=2, cex=0.8)
    abline(h=y_at, lty=2, col="lightgrey")
  }

  n=0
  for(i in 1:length(all)){
    if(length(all[[i]]) != 2){
      x0=which(rollapply(all[[1]], length(all[[i]]), identical, all[[i]]))
      x1=(which(rollapply(all[[1]], length(all[[i]]), identical, all[[i]]))+length(all[[i]])-1)
      segments(x0=x0, y0=n,
               x1=x1, y1=n,
               lwd=5, col=as.character(sequences$Color[i]))

      label=paste0(names(all)[i]," (",paste0(paste0(x0,"-",x1), collapse=";"),", ",length(all[[i]])," aa)")
      text(x=ifelse(x0[1] < FLlength*0.75, x0[1], x1[1]), y=(n+0.1), adj=ifelse(x0[1] < FLlength*0.75, 0, 1),
           labels=label,
           cex=0.9, offset=0, xpd=NA)
    } else {
      x0_1=which(rollapply(all[[1]], length(all[[i]][[1]]), identical, all[[i]][[1]]))
      x1_1=(which(rollapply(all[[1]], length(all[[i]][[1]]), identical, all[[i]][[1]]))+length(all[[i]][[1]])-1)
      segments(x0=x0_1, y0=n,
               x1=x1_1, y1=n,
               lwd=5, col=as.character(sequences$Color[i]))
      x0_2=which(rollapply(all[[1]], length(all[[i]][[2]]), identical, all[[i]][[2]]))
      x1_2=(which(rollapply(all[[1]], length(all[[i]][[2]]), identical, all[[i]][[2]]))+length(all[[i]][[2]])-1)
      segments(x0=x0_2, y0=n,
               x1=x1_2, y1=n,
               lwd=5, col=as.character(sequences$Color[i]))

      label=paste0(names(all)[i]," (",
                   paste0(paste0(x0_1,"-",x1_1), collapse=";"),", ",
                   paste0(paste0(x0_2,"-",x1_2), collapse=";"),
                   ", ",length(all[[i]][[1]]),"+",length(all[[i]][[2]])," aa)")
      text(x=ifelse(x0_1[1] < FLlength*0.75, x0_1[1], x1_1[1]), y=(n+0.1),
           adj=ifelse(x0_1[1] < FLlength*0.75, 0, 1),
           labels=label,
           cex=0.9, offset=0, xpd=NA)
    }
    n=n+0.5
  }
  mtext(paste0(gene," (Uniprot ",uniprot,")"), side=3, font=2)

  if(shouldtextplot){
    ap_textplot(sequences[,textplotcolumns][-1,], show.rownames=F)
  }

  if(shouldpdf){
    dev.off()
  }
}

#' Find amino acid indices
#'
#' Align query sequnces to full length sequences to find start and stopp indices of sequences.
#'
#' @param x List with at least one element, see Deatils for naming and content.
#' @param fasta_list A list of sequences in fasta format. Either this or sequence_list should be provided.
#' @param sequence_list A list of only sequences. Either this or fasta_list should be provided.
#' @details NB! The sequences have to be identical matches to consequtive stretches
#'     in the full length sequence you are aligning to. A future update will include partial matches.
#'
#'     The x list needs to include at least the element
#'
#'     BEADS = Beads info with identifiable rownames and with at least the columns:
#'
#'     "Sequence" with one letter abriviation for the aminoa acids,
#'
#'     "Uniprot" withthe uniprot ID the sequence originates from
#'
#'     If a list of fasta sequences is provided (\code{fasta_list}), it should be the output of
#'     \code{\link[seqinr:read.fasta]{read.fasta()}} in the package \code{seqinr}.
#'
#'     If a list of only sequenes is provided (\code{sequence_list}),
#'     each element should correspond to a Uniprot ID, isoform specific if applicable,
#'      and contain a consequitive string of one letter amino acid abbreviations.
#'
#' @return Updated input x with new columns in the BEADS element for the start and stopp indices of each tested isoform.
#' @export

ap_aaindex <- function(x, fasta_list=NULL, sequence_list=NULL){

  if(!is.null(fasta_list)){
    align_to <- fasta_list
  } else if (!is.null(sequence_list)){
    align_to <- sequence_list
  } else {
    stop("No list of sequence(s) to align to provided.")
  }

  for(p in seq_along(align_to)){
    tmp_isoform <- unlist(map(strsplit(names(align_to)[p], "\\|"), 2))
    tmp_uniprot <- unlist(strsplit(tmp_isoform, "-"))[1]

    if(!is.null(fasta_list)){
      FLseq <- getSequence(align_to[[p]])
    } else {
      FLseq <- substring(align_to[[p]], seq(1,nchar(align_to[[p]]),1), seq(1,nchar(align_to[[p]]),1))
    }

    tmp_seq <- x$BEADS[which(x$BEADS$Sequence != "" & grepl(tmp_uniprot, x$BEADS$Uniprot)),
                       "Sequence",drop=F]
    tmp_seq$Sequence <- as.character(tmp_seq$Sequence)
    tmp_seq <- apply(tmp_seq, 1, function(x) substring(x, seq(1,nchar(x),1), seq(1,nchar(x),1)))

    if(length(tmp_seq) > 1){
      ## Find alignments
      all <- append(list(FullLength=FLseq), tmp_seq)
      tmp_start <- rep(NA, length(all))
      tmp_stopp <- rep(NA, length(all))
      for(i in 1:length(all)){ # Find alignment
        x0=which(rollapply(all[[1]], length(all[[i]]), identical, all[[i]]))
        x1=(which(rollapply(all[[1]], length(all[[i]]), identical, all[[i]]))+length(all[[i]])-1)
        if(length(x0) == 0){ # Test if sequence is in alignment
          all[[i]] <- NA ; print(tmp_isoform) ; print(names(all)[i]) ; names(all)[i] <- "NoAlign" # Removes the sequence not in alignment
        } else {
          tmp_start[i] <- x0
          tmp_stopp[i] <- x1
        }
      }
      names(tmp_start) <- names(all)
      names(tmp_stopp) <- names(all)
      if("NoAlign" %in% names(all)){
        tmp_start <- tmp_start[-grep("NoAlign", names(all))]
        tmp_stopp <- tmp_stopp[-grep("NoAlign", names(all))]
        all <- all[-grep("NoAlign", names(all))]
      }
    }

    x$BEADS[, paste0("Start_",tmp_isoform)] <- tmp_start[match(rownames(x$BEADS), names(tmp_start))]
    x$BEADS[, paste0("Stopp_",tmp_isoform)] <- tmp_stopp[match(rownames(x$BEADS), names(tmp_stopp))]
  }
  return(x)
}
