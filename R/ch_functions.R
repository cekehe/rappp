#' Coupling efficency test
#'
#' Flags beads with signal similar to empty bead in coupling test, produces plot if wanted.
#'
#'
#' @param x List with at least three elements, see Deatils for naming and content.
#' @param empty_bead Column index for empty bead.
#' @param empty_co_multiple Number of sd above empty for cutoff.
#' @param shouldplot Logical, should a plot be made?
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
#'     "Type" with info about type of content on bead, at least including "PrEST" for PrESTs,
#'
#'     "Plate" with numerical coupling plate number(s).
#'
#' @return Updated input x with relevant filtering info and a pdf with plot (if \code{shouldplot=T}).
#' @export

ap_ct <- function(x, empty_bead, empty_co_multiple=3,
                  shouldplot=T, filename="coupling_efficiency.pdf", width=25, height=6, useDingbats=F, ...) {

    empty_co <- mean(x$CT[,empty_bead], ...) + empty_co_multiple*sd(x$CT[,empty_bead], ...)

    if(shouldplot){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)
      par(mar=c(12,4,2,8), cex.axis=0.8)
      layout(matrix(c(1,1,1,1,2), nrow=1))
      bs=beeswarm(x$CT, pch=16, las=2, corral="gutter", xaxt="n",
                  main="Copupling efficiency test", ylab="Signal intensity [MFI]",
                  pwcol=rep(ifelse(grepl("empty|bare|blank", colnames(x$CT),ignore.case=T), "orange",
                                   ifelse(grepl("his6abp|hisabp", colnames(x$CT),ignore.case=T), "darkgreen",
                                          ifelse(grepl("hig|anti-human", colnames(x$CT),ignore.case=T), "blue",
                                                 ifelse(grepl("ebna", colnames(x$CT),ignore.case=T), "purple",
                                                        ifelse(apply(x$CT, 2, mean) < empty_co, "red", "darkgrey"))))),
                            each=dim(x$CT)[1]))

      axis(1, at=bs$x[seq(1, dim(bs)[1], 3)], labels=unique(bs$x.orig), cex.axis=0.5, las=2, tick=F)

      vert_lines <- seq(min(bs$x)-0.5, max(bs$x)+0.5, 1)
      abline(v=vert_lines, lty=2,
             col=ifelse(c(duplicated(x$BEADS$Plate),F), "lightgrey", "black"))
      tmp <- vert_lines[which(!c(duplicated(x$BEADS$Plate),F))]
      mtext(text=paste("Plate ",unique(x$BEADS$Plate)), at=diff(tmp)/2+tmp[-length(tmp)],
            side=1, line=0, cex=0.7, font=2)

      legend(par("usr")[2], par("usr")[4],
             legend=c("Empty", "His6ABP", "ahIgX", "EBNA1", "Flagged"),
             col=c("orange", "darkgreen", "blue", "purple", "red"),
             pch=16, xpd=NA)
      abline(h=empty_co, lty=2)
      textxy(X=par("usr")[2], Y=empty_co, offset=0.55, cex=1,
             labs=paste0("mean(empty)+", empty_co_multiple,"*sd(empty)=", round(empty_co, 0)), xpd=NA)

      tmp_text <- matrix(colnames(x$CT)[which(apply(x$CT, 2, mean) < empty_co)], ncol=1)
      if(length(grep("HPRR", tmp_text, ignore.case=T)) > 0){
        tmp_text <- matrix(tmp_text[grep("HPRR", tmp_text, ignore.case=T)], ncol=1)
        ap_textplot(tmp_text, mar=c(2,2,1,2),
                 show.rownames=F, show.colnames=F, hadj=0, valign="top", cex=0.8)
        mtext("Protein fragments with low coupling efficiency signal", font=2, cex=0.9, xpd=NA)
      } else {
        frame()
        mtext("No protein fragments displayed \n low coupling efficiency signal.", font=2, cex=0.7, line=-3)
      }

      dev.off()
    }

    # Annotate filtering in BEADS
    if(length(which(colnames(x$BEADS) == "Flagged")) == 0){
      x$BEADS <- data.frame(Flagged="", x$BEADS, stringsAsFactors=F)
    }

    x$BEADS$Flagged <- ifelse(apply(x$CT, 2, mean) < empty_co &
                                grepl("PrEST", x$BEADS$Type, ignore.case=T),
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
#' @param IgX_bead Column index for empty bead.
#' @param IgType Which Imunoglobulin is measured, default is G.
#' @param IgX_cutoff MFI cutoff value for filtering.
#' @param cosfac Median absolute deviation multipliers in vector c(upper, lower),
#'     for drawing lines and detecting potential outliers.
#' @param shouldplot Logical, should a plot be made?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
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
#'     "sample_name" with sample names, preferably LIMS-IDs, where
#'     replicates (named with one of pool|rep|mix|commercial)
#'     and blanks (named with one of empty|blank|buffer) are also stated,
#'
#'     "AssayNum" with assay number (vector with 1s if only one assay, support for up to 5 assys in one plot),
#'
#'     "AssayWell" with Well IDs in the assay plate, e.g A01, B01 etc.,
#'
#'     "tube_label" with alternative sample names, eg. from collaborator,
#'
#' @return Updated input x with relevant filtering info and a pdf with plot (if \code{shouldplot=T}).
#' @export

ap_igx <- function(x, IgX_bead, IgType="G", IgX_cutoff=5000, cosfac=c(3, -3),
                   shouldplot=T, filename="anti-humanIgX.pdf", width=12, height=6, useDingbats=F) {

    plotdata <- unlist(x$MFI[,IgX_bead])
    sampledata <- x$SAMPLES
    SamplesNames <- sampledata$sample_name
    AssayNum <- sampledata$AssayNum

    cosIgG <- median(plotdata, na.rm=T)+cosfac*mad(plotdata, constant = 1, na.rm=T)
    tmp <- plotdata[grepl("empty|blank|buffer", SamplesNames, ignore.case=T)]
    cosIgG <- c(cosIgG, IgX_cutoff)

    which_lowIgG <- which(plotdata<cosIgG[3])

    if(shouldplot){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)
      layout(matrix(c(1,1,1,2,3,
                      1,1,1,4,5,
                      1,1,1,6,7), nrow=3, byrow=T))
      par(mar=c(5,5,4,4))

      plot(1:length(plotdata), plotdata, cex=0.6, pch=c(16:18,6,8)[AssayNum],
           xlab="Samples in analysis order",ylab="Signal intensity (MFI)",main=paste0("Total hIg", IgType, ""),
           col=ifelse(grepl("empty|blank|buffer", SamplesNames, ignore.case=T),2,
                      ifelse(grepl("pool|rep|mix|commercial", SamplesNames, ignore.case=T),5, 4)))
      if(length(unique(AssayNum)) > 1){
        legend(par("usr")[1], par("usr")[4], horiz=T, yjust=0.1, bty="n",
               legend=c("Sample","Replicate","Buffer", paste("Assay", unique(AssayNum))),
               cex=0.8,
               pch=c(rep(NA, 3), c(16:18,6,8)[unique(AssayNum)]),
               fill=c(4, 5, 2, rep(NA, length(unique(AssayNum)))), border=NA,
               col=c(rep(NA, 3), rep("grey", length(unique(AssayNum)))),
               xpd=NA)
      } else {
        legend(par("usr")[1], par("usr")[4], horiz=T, yjust=0.1, bty="n",
               legend=c("Sample","Replicate","Buffer"),
               cex=0.8,
               fill=c(4, 5, 2), border=NA,
               xpd=NA)
      }
      abline(h=cosIgG, lty=2)
      textxy(X=rep(par("usr")[2], 3),Y=cosIgG, labs=c(paste0(cosfac,"xMAD+median (all)"), "Filter cutoff"),
             offset=0.6, xpd=NA)

      par(mar=c(5,1,4,1))
      frame()
      frame()
      # Display samples with high but still outliers
      if(length(which(plotdata<cosIgG[2] & plotdata>cosIgG[3])) > 0) {
        plottext <- data.frame(AssayWell=sampledata$AssayWell,
                               InternalID=sampledata$sample_name,
                               Subject=sampledata$tube_label,
                               MFI=plotdata)[which(plotdata<cosIgG[2] & plotdata>cosIgG[3]),]
        plottext <- plottext[order(plottext$MFI, decreasing=T),]

        if(dim(plottext)[1] > 20){
          ap_textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top")
          ap_textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top")
          mtext(paste0("anti-hIg", IgType, " MFI between ",cosIgG[3]," & ", cosIgG[2]), font=2, cex=0.6, xpd=NA, at=-0.5)
        } else {
          ap_textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top")
          mtext(paste0("anti-hIg", IgType, " MFI between ",cosIgG[3]," & ", cosIgG[2]), font=2, cex=0.6)
          frame()
        }

      } else {
        frame()
        frame()
      }

      # Display and remove samples with low total IgG signal
      if(length(which_lowIgG) > 0) {
        plottext <- data.frame(AssayWell=sampledata$AssayWell,
                               InternalID=sampledata$sample_name,
                               Subject=sampledata$tube_label,
                               MFI=plotdata)[which_lowIgG,]
        plottext <- plottext[order(plottext$MFI, decreasing=T),]

        if(dim(plottext)[1] > 20){
          ap_textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top")
          ap_textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top")
          mtext(paste0("anti-hIg", IgType, " MFI below ",cosIgG[3]), font=2, cex=0.6, xpd=NA, at=-0.5)
        } else {
          ap_textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=0.7, valign="top")
          mtext(paste0("anti-hIg", IgType, " MFI below ",cosIgG[3]), font=2, cex=0.6)
          frame()
        }

      } else {
        frame()
        frame()
      }
      dev.off()
    }

    # Annotate filtering in SAMPLES
    if(length(which(colnames(x$SAMPLES) == "Filtered")) == 0){
      x$SAMPLES <- data.frame(Filtered="", x$SAMPLES, stringsAsFactors=F)
    }
    if(length(which_lowIgG) > 0) {
      tmp_remove <- rownames(sampledata)[which_lowIgG]
      tmp_remove <- tmp_remove[-grep("empty", tmp_remove, ignore.case=T)]
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
#' @param labels Column name in BEADS with antigen names to be used in pdf.
#' @param protein Column name in BEADS with short protein or gene name.
#' @param agID Column name in BEADS with antigen identifier, eg. PrEST ID or product number.
#' @param samp_co Cutoff for filtering samples with low median count.
#' @param bead_flag Cutoff for flagging beads with low counts.
#' @param bead_filter Cutoff for filtering beads with low counts.
#' @param N_filter Accepted number of samples with low count per bead ID.
#' @param shouldplot Logical, should a plot be made?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
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
#'     "Plate" with coupling plate number(s),
#'
#'     "BeadID" with bead ID number,
#'
#'     separate columns with full desired labels, protein/gene names and product ID numbers.
#'
#' The SAMPLES element needs at least the columns:
#'
#'     "sample_name" with sample names, preferably LIMS-IDs, where
#'     replicates (named with one of pool|rep|mix|commercial)
#'     and blanks (named with one of empty|blank|buffer) are also stated,
#'
#'     "AssayWell" with Well IDs in the assay plate, e.g A01, B01 etc.,
#'
#'     "tube_label" with alternative sample names, eg. from collaborator,
#'
#'     "WashCol", a numerical vector with groupings for during-run-washes in Luminex.
#'     Eg. c(rep(1, 96), rep(2, 96), rep(1, 96), rep(2, 96)) if washes every 96th well in a 384-well plate.
#'     Colors boxes accordingly. There are 6 colors available.
#'
#' @return Updated input x with relevant filtering and/or flagging info and a pdf with plots (if \code{shouldplot=T}).
#' @export

ap_count <- function(x, labels="Gene_HPRR", protein="GeneShort", agID="PrEST",
                     samp_co=32, bead_flag=32, bead_filter=16, N_filter=0,
                   shouldplot=T, filename="bead_count.pdf", width=12, height=10, useDingbats=F) {

    plotdata <- t(x$COUNT)
    sampledata <- x$SAMPLES
    beaddata <- x$BEADS

    which_lowSB <- which(apply(plotdata, 2, function(x) median(x, na.rm=T)) < samp_co) # Checks which samples have low count in general, e.g. due to faulty bead dispens in well

    if(shouldplot){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)

      layout(matrix(c(1,1,1,2,
                      3,3,3,4), nrow=2, byrow=T))
      par(mar=c(4, 4, 2, 6))

      # Per sample ALL DATA
      boxplot(plotdata, pch=16, cex=0.6, ylim=c(0, max(plotdata, na.rm=T)), las=1, names=F, xaxt="n",
              col=c("burlywood1", "burlywood4",
                    "darkorchid2","darkorchid4",
                    "darkolivegreen2","darkolivegreen4")[sampledata$WashCol],
              border=c("burlywood1", "burlywood4",
                       "darkorchid2","darkorchid4",
                       "darkolivegreen2","darkolivegreen4")[sampledata$WashCol],
              main="Sample bead count", ylab="Bead count per sample")
      text(1:dim(plotdata)[2],par("usr")[3]-1, labels = rownames(sampledata),
           srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=0.1)
      abline(h=c(16,32, median(plotdata, na.rm=T)), lty=2, col=c("grey","red", "cornflowerblue"))
      legend(par("usr")[2], par("usr")[4],
             legend=c(bead_filter,
                      paste0("Failed (", samp_co, ")"),
                      paste0("Median (", median(plotdata, na.rm=T), ")"),
                      rep("Wash cycle grouping", length(unique(sampledata$WashCol)))),
             lty=c(rep(2,3), rep(0,length(unique(sampledata$WashCol)))),
             col=c("grey","red", "cornflowerblue", rep(0, length(unique(sampledata$WashCol)))),
             fill=c(rep(0, 3), c("burlywood1", "burlywood4",
                                 "darkorchid2","darkorchid4",
                                 "darkolivegreen2","darkolivegreen4")[unique(sampledata$WashCol)]),
             border=c(rep(0,3), rep("black",length(unique(sampledata$WashCol)))),
             xpd=NA, cex=0.7)
      mtext("Sample wells, in order of analysis", side=1, cex=0.7, line=1)

      if(length(which_lowSB) > 0){
        lowSB <- sampledata[which_lowSB, ]
        tp <- ap_textplot(data.frame( # function in package: gplots
          AssayWell=lowSB$AssayWell,
          InternalID=lowSB$sample_name,
          Subject=lowSB$tube_label,
          MedianCount=apply(plotdata, 2, function(x) median(x, na.rm=T))[which_lowSB],
          LowestCount=apply(plotdata, 2, function(x) min(x, na.rm=T))[which_lowSB],
          HighestCount=apply(plotdata, 2, function(x) max(x, na.rm=T))[which_lowSB]),
          cex=0.4, cmar=1.5, show.rownames=F, valign="top")
        title(paste0("Samples with median bead count < ,", samp_co, ", (N=", dim(lowSB)[1], ")"), xpd=NA)

      } else {
        ap_textplot(matrix("No samples filtered based on bead count."), show.rownames=F, show.colnames=F)
      }
    }

    # REMOVE LOW SAMPLES
    if(length(which_lowSB) > 0) { # Removes samples from plotdata with low median bead count
      plotdata <- t(plotdata[, -which_lowSB])
      sampledata <- sampledata[-which_lowSB, ]
    } else {
      plotdata <- t(plotdata)
    }

    # Per antigen
    which_failAB <- which(apply(plotdata, 2, function(x) length(which(x < bead_filter))) > 0)
    which_lowAB <- which(apply(plotdata, 2, function(x) length(which(x < bead_flag))) > 0)

    if(shouldplot){
      bp=boxplot(plotdata, pch=16, cex=0.6, ylim=c(0, max(plotdata, na.rm=T)), las=1, names=F, xaxt="n",
                 ylab="Bead count per analyte",
                 col=c("blueviolet","burlywood4","chartreuse4","coral4")[beaddata$Plate], # color boxes based on which coupling plate they come from
                 border=c("blueviolet","burlywood4","chartreuse4","coral4")[beaddata$Plate])
      text(1:dim(plotdata)[2],par("usr")[3]-1, labels = beaddata[,labels],
           srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=0.1)
      abline(h=c(bead_filter,
                 bead_flag,
                 median(unlist(plotdata), na.rm=T)), lty=2, col=c("red", "orange", "cornflowerblue"))
      legend(par("usr")[2], par("usr")[4],
             legend=c(paste0("Failed if N>", N_filter, " (", bead_filter,")"),
                      paste0("Flagged (", bead_flag,")"),
                      paste0("Median (", median(unlist(plotdata), na.rm=T), ")"),
                      paste0("Coupling plate ", unique(beaddata$Plate))),
             lty=c(rep(2,3), rep(0, length(unique(beaddata$Plate)))),
             col=c("red", "orange", "cornflowerblue", rep(0, length(unique(beaddata$Plate)))),
             fill=c(rep(0,3), unique(c("blueviolet","burlywood4","chartreuse4","coral4")[beaddata$Plate])),
             border=c(rep(0,3), rep(1, length(unique(beaddata$Plate)))),
             xpd=NA, cex=0.7)
      mtext("Analyte bead count", side=3, cex=1, font=2, line=1)
      mtext(paste0("Samples with median bead count < ", samp_co, " removed (see above plot)"),
            side=3, cex=0.6, line=0)
    }

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
      lowAB <- lowAB[order(lowAB$LowestCount),]
      lowAB <- data.frame(lowAB, Action=ifelse(lowAB$LowestCount > bead_filter | lowAB$Nbelow16 <= N_filter, "Flagged", "Filtered"))

      if(shouldplot){
        tp <- ap_textplot(lowAB, cex=0.3, cmar=1.5, show.rownames=F, valign="top", xpd=NA) # function in package: gplots
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

    ## FILTERED PLOTS
    if(shouldplot){
      layout(matrix(c(1,
                      2), nrow=2, byrow=T))
      par(mar=c(4, 4, 2, 8))

      # Per sample LOW SAMPLES AND ANALYTES REMOVED
      boxplot(plotdata, pch=16, cex=0.6, ylim=c(0, max(plotdata, na.rm=T)), las=1, names=F, xaxt="n",
              col=c("burlywood1", "burlywood4","darkorchid2","darkolivegreen2","darkolivegreen4")[sampledata$WashCol],
              border=c("burlywood1", "burlywood4","darkorchid2","darkolivegreen2","darkolivegreen4")[sampledata$WashCol],
              ylab="Bead count per sample")
      text(1:dim(plotdata)[2], par("usr")[3]-1, labels = rownames(sampledata),
           srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=0.1)
      abline(h=c(16,32, median(plotdata, na.rm=T)), lty=2, col=c("grey","red", "cornflowerblue"))
      legend(par("usr")[2], par("usr")[4],
             legend=c(bead_filter,
                      paste0("Failed (", samp_co, ")"),
                      paste0("Median (", median(plotdata, na.rm=T), ")"),
                      rep("Wash cycle grouping", length(unique(sampledata$WashCol)))),
             lty=c(rep(2,3), rep(0,length(unique(sampledata$WashCol)))),
             col=c("grey","red", "cornflowerblue", rep(0, length(unique(sampledata$WashCol)))),
             fill=c(rep(0, 3), c("burlywood1", "burlywood4","darkorchid2","darkolivegreen2","darkolivegreen4")[unique(sampledata$WashCol)]),
             border=c(rep(0,3), rep("black",length(unique(sampledata$WashCol)))),
             xpd=NA, cex=0.7)
      mtext("Sample wells, in order of analysis", side=1, cex=0.7, line=1)
      mtext("Samples bead count", side=3, cex=1, font=2, line=1)
      mtext(paste0("Filtered: Samples with median bead count < ", samp_co,
            " and analytes with >", N_filter, " samples with bead count < ", bead_filter, " removed"),
            side=3, cex=0.6, line=0)

      # Per analyte LOW SAMPLES AND ANALYTES REMOVED
      plotdata <- t(plotdata) ; Info_plotdata <- "analyte_count"
      bp=boxplot(plotdata, pch=16, cex=0.6, ylim=c(0, max(plotdata, na.rm=T)), las=1, names=F, xaxt="n",
                 ylab="Bead count per analyte",
                 col=c("blueviolet","burlywood4","chartreuse4","coral4")[beaddata$Plate],  # color boxes based on which coupling plate they come from
                 border=c("blueviolet","burlywood4","chartreuse4","coral4")[beaddata$Plate])
      text(1:dim(plotdata)[2],par("usr")[3]-1, labels = beaddata$Gene_HPRR,
           srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=0.1)
      abline(h=c(16, 32, median(unlist(plotdata), na.rm=T)), lty=2, col=c("red", "orange", "cornflowerblue"))
      legend(par("usr")[2], par("usr")[4],
             legend=c(paste0("Failed if N>", N_filter, " (", bead_filter,")"),
                      paste0("Flagged (", bead_flag,")"),
                      paste0("Median (", median(unlist(plotdata), na.rm=T), ")"),
                      paste0("Coupling plate ", unique(beaddata$Plate))),
             lty=c(rep(2,3), rep(0, length(unique(beaddata$Plate)))),
             col=c("red", "orange", "cornflowerblue", rep(0, length(unique(beaddata$Plate)))),
             fill=c(rep(0,3), unique(c("blueviolet","burlywood4","chartreuse4","coral4")[beaddata$Plate])),
             border=c(rep(0,3), rep(1, length(unique(beaddata$Plate)))),
             xpd=NA, cex=0.7)
      mtext("Analyte bead count", side=3, cex=1, font=2, line=1)
      mtext(paste0("Filtered: Samples with median bead count < ", samp_co,
                   " and analytes with >", N_filter, " samples with bead count < ", bead_filter, " removed"),
            side=3, cex=0.6, line=0)
      dev.off()
    }

    # Annotate filtering in SAMPLES and BEADS
    # SAMPLES
    if(length(which(colnames(x$SAMPLES) == "Filtered")) == 0){
      x$SAMPLES <- data.frame(Filtered="", x$SAMPLES, stringsAsFactors=F)
    }
    if(length(which_lowSB) > 0){
     x$SAMPLES$Filtered <- ifelse(rownames(x$SAMPLES) %in% names(which_lowSB),
                                  paste0(x$SAMPLES$Filtered,", Count"),
                                  paste(x$SAMPLES$Filtered))
     x$SAMPLES$Filtered <- gsub("^, ", "", x$SAMPLES$Filtered)
    }

    # BEADS filtered
    if(length(which(colnames(x$BEADS) == "Filtered")) == 0){
      x$BEADS <- data.frame(Filtered="", x$BEADS, stringsAsFactors=F)
    }

    if(length(which(lowAB$Action == "Filtered")) > 0){
      x$BEADS$Filtered <- ifelse(x$BEADS$Gene_HPRR %in% rownames(lowAB)[which(lowAB$Action == "Filtered")],
                                 paste0(x$BEADS$Filtered,", Count"),
                                 paste(x$BEADS$Filtered))
      x$BEADS$Filtered <- gsub("^, ", "", x$BEADS$Filtered)
    }

    # BEADS flagged
    if(length(which(colnames(x$BEADS) == "Flagged")) == 0){
      x$BEADS <- data.frame(Flagged="", x$BEADS, stringsAsFactors=F)
    }

    if(length(which(lowAB$Action == "Flagged")) > 0){
    x$BEADS$Flagged <- ifelse(x$BEADS$Gene_HPRR %in% rownames(lowAB)[which(lowAB$Action == "Flagged")],
                              paste0(x$BEADS$Flagged,", Count"),
                              paste(x$BEADS$Flagged))
    x$BEADS$Flagged <- gsub("^, ", "", x$BEADS$Flagged)
}

  x$FILTERINFO <- c(x$FILTERINFO, "beadcount")

  return(x)
}

#' Signal overview
#'
#' Boxplots of signals per antigen and sample in different orders.
#'
#' @param x List with at least one elements, see Deatils for naming and content.
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @param ... Further arguments passed to \code{\link[graphics:boxplot]{boxplot()}}.
#' @details The x list needs to include at least the element
#'     MFI = assay mfi,
#' @export

ap_overview <- function(x,
                  filename="Signal_overview.pdf",
                  width=20, height=15, useDingbats=F, ...){

  pdf(filename,
      width=width, height=height, useDingbats=useDingbats)
  par(mfcol=c(3, 1), mar=c(10,5,4,2))

    ## Antigens
    plotdata <- list('Bead ID'=x$MFI,
                     median=x$MFI[,order(apply(x$MFI, 2, median, na.rm=T))],
                     max=x$MFI[,order(apply(x$MFI, 2, max, na.rm=T))])

    for(i in 1:length(plotdata)){
    boxplot(plotdata[[i]], pch=16, cex=0.5, log="y", las=2, ...,
            main=paste0("Antigens, sorted by ", names(plotdata)[i]), ylab="log(MFI) [AU]",
            outcol=ifelse(grepl("his6abp|hisabp|empty|bare", colnames(plotdata[[i]]), ignore.case=T), as.color("brown", 0.7),
                          ifelse(grepl("anti-h|hIg|ebna", colnames(plotdata[[i]]), ignore.case=T), as.color("darkolivegreen", 0.7), as.color("black", 0.5))),
            col=ifelse(grepl("his6abp|hisabp|empty|bare", colnames(plotdata[[i]]), ignore.case=T), as.color("brown", 0.7),
                       ifelse(grepl("anti-h|hIg|ebna", colnames(plotdata[[i]]), ignore.case=T), as.color("darkolivegreen", 0.7), 0)))
    }

    ## Samples
    tmp <- data.frame(t(x$MFI), check.names=F)
    plotdata <- list('analysis order'=tmp,
                     median=tmp[,order(apply(tmp, 2, median, na.rm=T))],
                     max=tmp[,order(apply(tmp, 2, max, na.rm=T))])

    for(i in 1:length(plotdata)){
    boxplot(plotdata[[i]], pch=16, cex=0.5, log="y", las=2,  ...,
            main=paste0("Samples, sorted by ", names(plotdata)[i]), ylab="log(MFI) [AU]",
            outcol=ifelse(grepl("empty|buffer|blank", colnames(plotdata[[i]]), ignore.case=T), as.color("brown", 0.7),
                          ifelse(grepl("rep|pool|mix|commercial", colnames(plotdata[[i]]), ignore.case=T), as.color("cornflowerblue", 0.7), as.color("black", 0.5))),
            col=ifelse(grepl("empty|buffer|blank", colnames(plotdata[[i]]), ignore.case=T), as.color("brown", 0.7),
                       ifelse(grepl("rep|pool|mix|commercial", colnames(plotdata[[i]]), ignore.case=T), as.color("cornflowerblue", 0.7), as.color("black", 0.5))))
    }
  dev.off()
}


#' Check replicate samples
#'
#' various plots in one PDF to assess the replicates and overall reproducibility.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @param iter How many times  random samples should be iterated.
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
#' @export

ap_rep <- function(x, iter=500, filename="replicates.pdf", width=12, height=12, useDingbats=F){
  ## INPUT
  data <- append(split(x$MFI, x$SAMPLES$AssayNum), list(x$MFI))
  names(data) <- c(paste0("Assay_", names(data)[-length(data)]), "All_combined")

  samples <- append(split(x$SAMPLES, x$SAMPLES$AssayNum), list(x$SAMPLES))
  names(samples) <- c(paste0("Assay_", names(samples)[-length(samples)]), "All_combined")

  ## CALCULATIONS
  CVs <- matrix(NA, ncol=length(data), nrow=dim(x$MFI)[2])
  cor_samp_s <- rep(list(NULL), length(data))
  cor_samp_p <- rep(list(NULL), length(data))

  CVs_r_m <- matrix(NA, ncol=length(data), nrow=dim(x$MFI)[2])
  cor_samp_s_r_m <- rep(list(NULL), length(data))
  cor_samp_p_r_m <- rep(list(NULL), length(data))

  nrreplicates <- rep(NA, length(data))

  for(l in 1:length(data)){
    # Calculate for replicates
    replicates <- data[[l]][grep("pool|rep|mix|commercial", samples[[l]]$sample_name, ignore.case=T),]
    nrreplicates[l] <- dim(replicates)[1]
    CVs[,l] <- apply(replicates, 2, cv, digits=5, na.rm=T)

    cor_samp_s[[l]] <- cor(t(replicates), method="spearman", use="pairwise.complete.obs")
    cor_samp_s[[l]] <- cor_samp_s[[l]][upper.tri(cor_samp_s[[l]])]

    cor_samp_p[[l]] <- cor(t(log(replicates)), method="pearson", use="pairwise.complete.obs")^2
    cor_samp_p[[l]] <- cor_samp_p[[l]][upper.tri(cor_samp_p[[l]])]

    # Iterate over random sets of samples
    CVs_r <- matrix(NA, ncol=iter, nrow=dim(x$MFI)[2])
    tmp_data <- data[[l]][-grep("pool|rep|mix|commercial", samples[[l]]$sample_name, ignore.case=T),]
    cor_samp_s_r <- matrix(NA, ncol=iter, nrow=(dim(tmp_data)[1]^2-dim(tmp_data)[1])/2)
    cor_samp_p_r <- matrix(NA, ncol=iter, nrow=(dim(tmp_data)[1]^2-dim(tmp_data)[1])/2)
    for(j in 1:iter){
      rand_samp <- tmp_data[sample(1:dim(tmp_data)[1], nrreplicates[l], replace=F),]
      CVs_r[,j] <- apply(rand_samp, 2, cv, digits=5, na.rm=T)

      cor_samp_s_tmp <- cor(t(rand_samp), method="spearman", use="pairwise.complete.obs")
      cor_samp_s_r[,j] <- cor_samp_s_tmp[upper.tri(cor_samp_s_tmp)]

      cor_samp_p_tmp <- cor(t(log(rand_samp)), method="pearson", use="pairwise.complete.obs")^2
      cor_samp_p_r[,j] <- cor_samp_p_tmp[upper.tri(cor_samp_p_tmp)]
    }
    CVs_r_m[,l] <- apply(CVs_r, 1, median, na.rm=T)
    cor_samp_s_r_m[[l]] <- apply(cor_samp_s_r, 1, median, na.rm=T)
    cor_samp_p_r_m[[l]] <- apply(cor_samp_p_r, 1, median, na.rm=T)
  }
  colnames(CVs) <- names(data)
  colnames(CVs_r_m) <- paste0(names(data), "_Random")
  names(cor_samp_s) <- names(data)
  names(cor_samp_p) <- names(data)
  names(cor_samp_s_r_m) <- paste0(names(data), "_Random")
  names(cor_samp_p_r_m) <- paste0(names(data), "_Random")

  assay_cv <- lapply(data, function(y) apply(y, 2, function(x) cv(x, na.rm=T, digits=5)))
  assay_max <- lapply(data, function(y) apply(y, 2, function(x) max(x, na.rm=T)))
    assay_mean <- lapply(data, function(y) apply(y, 2, function(x) mean(x, na.rm=T)))
  assay_median <- lapply(data, function(y) apply(y, 2, function(x) median(x, na.rm=T)))

  ## Set plotorder
  plotorder <- rep(NA, 2*length(data))
  plotorder[seq(1, length(plotorder), 2)] <- 1:length(data)
  plotorder[seq(2, length(plotorder), 2)] <- (length(data)+1):(2*length(data))

  ## PLOTS
  data_cv <- data.frame(CVs, CVs_r_m) ; data_melt <- melt(data_cv, id.vars=NULL)
  data_melt$variable <- factor(data_melt$variable, levels=levels(data_melt$variable)[plotorder])
  pdf(filename,
      width=12, height=12, useDingbats=F)
  par(mar=c(9,4,3,1))

  # CV boxplots
  boxplot(data_melt$value~data_melt$variable, data=data_melt, outcol=0, las=1, col=grey.colors(2),
          ylab="CVs [%]", xlab=NA, main="CVs between replicates \n one point = one antigen", names=NA)
  beeswarm(data_melt$value~data_melt$variable, data=data_melt, pch=16, cex=0.3, corral="gutter", add=T)
  legend("topleft",
         legend=c("True replicates",
                  paste0("False replicates (", iter," iterations)"),
                  "CV=10%"),
         fill=c(grey.colors(2), NA),
         border=c(rep("black", 2), 0),
         lty=c(rep(NA, 2), 2),
         col=c(rep(NA, 2), "red"))
  abline(h=10, col="red", lty=2)
  text(1:dim(data_cv)[2],par("usr")[3]-11, labels = levels(data_melt$variable),
       srt = 45, adj=c(1.1,1.1), xpd = TRUE, cex=0.9)
  mtext(paste0(rep(nrreplicates, each=2), " samples"), at=1:dim(data_cv)[2], side=1, line=0.5, cex=0.8)

  # Correlations
  par(mfrow=c(5,2), mar=c(4,4,3,1))
  for(l in 1:length(cor_samp_s)){
    dens_s <- density(cor_samp_s[[l]])
    dens_p <- density(cor_samp_p[[l]])
    dens_s_r <- density(cor_samp_s_r_m[[l]])
    dens_p_r <- density(cor_samp_p_r_m[[l]])

    plot(range(0, 1), range(dens_s$y, dens_s_r$y), type = "n", xlab = "Spearman's rho",
         ylab = "Density", main=paste0(names(cor_samp_s)[l],": ", nrreplicates[l], " samples"))
    lines(dens_s, col = "black")
    lines(dens_s_r, col = "cornflowerblue")
    legend("topleft", legend=c("True replicates", paste0("False replicates (", iter," iterations)")), lty=1, col=c("black","cornflowerblue"), cex=0.7)

    plot(range(0, 1), range(dens_p$y, dens_p_r$y), type = "n", xlab = bquote("Pearson's R"^"2"),
         ylab = "Density", main=paste0(names(cor_samp_s)[l],": ", nrreplicates[l], " samples"))
    lines(dens_p, col = "black")
    lines(dens_p_r, col = "cornflowerblue")
  }

  # Assay CV vs replicate CV
  par(mfrow=c(2,3))
  for(l in 1:length(assay_cv)){
    plot(CVs[,l], assay_cv[[l]], pch=16, cex=0.7, xlim=c(0,50),
         xlab="CVs of replicates, per antigen [%]", ylab="CVs of assay, per antigen [%]", main=colnames(CVs)[l])
    textxy(CVs[,l][which(CVs[,l] > 10)],
         assay_cv[[l]][which(CVs[,l] > 10)], cex=0.4, offset=0.55,
         labs=names(assay_cv[[l]])[which(CVs[,l] > 10)])
    abline(v=10, lty=2, col="grey")
    abline(h=10, lty=2, col="grey")
  }
  for(f in 1:(6-length(assay_cv))){
    frame()
  }

  # Assay CV vs assay max
  for(l in 1:length(assay_max)){
    plot(assay_max[[l]], assay_cv[[l]], pch=16, cex=0.7,
         xlab="Max signal intensity per antigen [MFI]", ylab="CVs of assay, per antigen [%]", main=names(assay_max)[l])
    abline(h=10, lty=2, col="grey")
  }
  for(f in 1:(6-length(assay_cv))){
    frame()
  }

  # Assay CV vs assay median
  for(l in 1:length(assay_median)){
    plot(assay_median[[l]], assay_cv[[l]], pch=16, cex=0.7,
         xlab="Median signal intensity per antigen [MFI]", ylab="CVs of assay, per antigen [%]", main=names(assay_median)[l])
    abline(h=10, lty=2, col="grey")
  }
  for(f in 1:(6-length(assay_cv))){
    frame()
  }

  # Assay CV vs assay mean
  for(l in 1:length(assay_mean)){
    plot(assay_mean[[l]], assay_cv[[l]], pch=16, cex=0.7,
         xlab="Mean signal intensity per antigen [MFI]", ylab="CVs of assay, per antigen [%]", main=names(assay_mean)[l])
    abline(h=10, lty=2, col="grey")
  }

  dev.off()
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
         width=4.5*length(perp), height=height, useDingbats=useDingbats)
}

#' Overview of signals in relation to neg control beads
#'
#' Plot overview of singals in relation to the neg control beads.
#' Based on output from Autoimmunity Profiling scoring function \code{\link[rappp:ap_scoring2]{ap_scoring2()}}.
#'
#' @param x List with at least four elements, see Deatils for naming and content.
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @details The x list needs to include at least the element
#'
#'     MFI = assay mfi, column names of negative control columns should include empty|bare|blank|his6abp|hisabp,
#'
#'     SCORE = scored data, column names of negative control columns should include empty|bare|blank|his6abp|hisabp,
#'
#'     CUTOFF_KEY = Cutoff key as data.frame with cutoff values, scores and colors.
#'
#'     BEADS = Beads info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any beads with no text (ie. "" or NA) or "NegControl" in such column will be included in the transformation.
#'
#' @export

ap_negbeads <- function(x,
                        filename="neg-control-beads.pdf", width=15, height=10, useDingbats=F){

  pdf(filename, width=width, height=height, useDingbats=useDingbats)

  layout(matrix(c(1,1,2,3), ncol=2, byrow=T))
  par(mar=c(4,4,4,5))

  if("Filtered" %in% colnames(x$BEADS)){
    plotdata <- x$MFI[, which(is.na(x$BEADS$Filtered) |
                                x$BEADS$Filtered == "" |
                                grepl("NegControl", x$BEADS$Filtered))]
    plotdata_score <- x$SCORE[, which(is.na(x$BEADS$Filtered) |
                                        x$BEADS$Filtered == "" |
                                        grepl("NegControl", x$BEADS$Filtered))]
  } else {
    plotdata <- x$MFI
    plotdata_score <- x$SCORE
  }

  plotcolor <- x$CUTOFF_KEY

  beeswarm(data.frame(t(plotdata)), log=T, corral="gutter", cex=0.5, las=2,
           pwcol=ifelse(grepl("empty|bare|blank", rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T),"magenta",
                        ifelse(grepl("his6abp|hisabp",rep(colnames(plotdata), dim(plotdata)[1]), ignore.case=T), "chartreuse1",
                               as.color(paste(plotcolor$color[t(plotdata_score)*10+1]), 0.4))),
           pwpch=rep(ifelse(grepl("empty|bare|blank|his6abp|hisabp", colnames(plotdata), ignore.case=T), 16, 1), dim(plotdata)[1]),
           ylab="Signal intensity [AU]", cex.axis=0.5)
  legend(par("usr")[2], 10^par("usr")[4], xpd=T, cex=0.7,
         legend=c("empty", "his6abp", rev(rownames(plotcolor))),
         col=c("magenta", "chartreuse1", paste(rev(plotcolor$color))),
         pch=c(rep(16, 2), rep(1, length(plotcolor$color))))

  par(pty="s")
  plot(apply(plotdata, 1, median, na.rm=T), plotdata[,grep("empty|bare|blank", colnames(plotdata), ignore.case=T)],
       las=1, xlab="Median signal per sample", ylab="Empty bead signal", main="Empty bead",
       col=paste(plotcolor$color[plotdata_score[,grep("empty|bare|blank", colnames(plotdata), ignore.case=T)]*10+1]))

  plot(apply(plotdata, 1, median, na.rm=T), plotdata[,grep("his6abp|hisabp", colnames(plotdata), ignore.case=T)],
       las=1, xlab="Median signal per sample", ylab="His6ABP bead signal", main="His6ABP bead",
       col=paste(plotcolor$color[plotdata_score[,grep("his6abp|hisabp", colnames(plotdata), ignore.case=T)]*10+1]))

  dev.off()
}

#' Calculate reactivity frequencies
#'
#' Calculate number of reactivities and the corresponding frequencies in Autoimmunity Profiling data.
#'
#' @param x List with at least two elements, see Details for naming and content.
#' @param samplegroups factor vector of groupings. Only samples with an assigned level are included in plots.
#'     If left as \code{NULL} (default), the all non-filtered, if filetring done otherwise all, will be assigned "Sample".
#' @param check.names logical, altered default from \code{\link[base:data.frame]{data.frame()}}.
#' @details
#'
#' The x list needs to include at least the elements:
#'
#'     SAMPLES = Sample info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any beads with no text (ie. "") in such column will be included.
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
                             check.names = FALSE) {

  if("BINARY_CO" %in%  names(x)){
    data_bin <- append(x$BINARY, list(Selected_co=x$BINARY_CO))
  } else {
    data_bin <- x$BINARY
  }

  if(is.null(samplegroups)){
    if("Filtered" %in% colnames(x$SAMPLES)){
      samplegroups <- factor(ifelse(is.na(x$SAMPLES$Filtered) | x$SAMPLES$Filtered == "", "Sample", NA))
    } else {
      samplegroups <- factor(rep("Sample", dim(x$SAMPLES)[1]))
    }
  }
  data_size <- table(samplegroups)
  n_ag <- lapply(data_bin, function(i) apply(i, 1, function(l) sum(!is.na(l))))

  # Calculate per antigen
  data_sum_ag <- lapply(data_bin, function(i) apply(i, 2, function(l) aggregate(l, by=list(samplegroups), FUN=sum)))
  names(data_sum_ag) <- names(data_bin)

  data_freq_ag <- lapply(1:length(data_sum_ag),
                         function(cutoff) lapply(data_sum_ag[[cutoff]],
                                                 function(antigen) round(antigen$x/data_size*100,1)))
  names(data_freq_ag) <- names(data_sum_ag)


  data_sum_ag <- lapply(data_sum_ag,
                        function(cutoff) data.frame(do.call(cbind,
                                                            lapply(cutoff,
                                                                   function(antigen) antigen$x)), check.names = check.names))
  data_sum_ag <- lapply(data_sum_ag, function(i) {
    rownames(i) <- levels(samplegroups) ;
    colnames(i) <- colnames(data_sum_ag[[length(data_sum_ag)]]) ; i } )

  data_freq_ag <- lapply(data_freq_ag, function(cutoff) data.frame(do.call(cbind, cutoff), check.names = check.names))
  data_freq_ag <- lapply(data_freq_ag, function(i) {
    colnames(i) <- colnames(data_freq_ag[[length(data_freq_ag)]]) ; i } )

  data_sum_ag <- do.call(rbind, data_sum_ag)
  data_freq_ag <- do.call(rbind, data_freq_ag)

  # Calculate per sample
  data_sum_samp <- lapply(data_bin,
                          function(i) apply(i, 1,
                                            function(l) sum(l, na.rm=T)))
  names(data_sum_samp) <- names(data_bin)

  data_freq_samp <- lapply(1:length(data_sum_samp), function(cutoff) round(data_sum_samp[[cutoff]]/n_ag[[cutoff]]*100,1))
  names(data_freq_samp) <- names(data_sum_samp)

  data_sum_samp <- do.call(cbind, data_sum_samp)
  data_freq_samp <- do.call(cbind, data_freq_samp)

  # Add to input
  output <- list(SAMPLEGROUPS=data.frame(Sample=x$SAMPLES$sample_name,
                                         Grouping=samplegroups),
                 REACTSUM_AG=data_sum_ag,
                 REACTFREQ_AG=data_freq_ag,
                 REACTSUM_SAMP=data_sum_samp,
                 REACTFREQ_SAMP=data_freq_samp)

  return(output)
}

#' Reactivity result plots
#'
#' Plot beeswarm, density and frequency plots for each antigen.
#' Based on output from Autoimmunity Profiling wrapper function \code{\link[rappp:ap_norm2]{ap_norm2()}}.
#'
#' @param x list with at least nine elements, see Deatils for naming and content.
#' @param samplegroups factor vector of groupings. Only samples with an assigned level are included in plots.
#'     If left as \code{NULL} (default), the all non-filtered if filtering has been done,
#'     otherwise all, will be assigned "Sample".
#'     Passed to \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} to calculate frequencies.
#' @param groupcolors colors for each group in samplegroups.
#' @param agtoplot indices for which antigens to plot, default is all.
#'     Character vector with column names of what to plot also ok.
#' @param filename string with filename and desired path, end with .pdf
#' @param height width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
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
#'     Any beads with no text (ie. "" or NA) in such column will be included.
#'
#'     BEADS = beads info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any beads with no text (ie. "" or NA) will be included in the transformation.
#'
#'     DENS = Density output used for cutoff selection,
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
                         groupcolors=1:6,
                         agtoplot=NULL,
                         filename="AntigenResults.pdf",
                         height=15,
                         useDingbats=F,
                         check.names=FALSE) {

  print("Calculating frequencies")
  react_summary <- ap_reactsummary2(x,
                                    samplegroups = samplegroups,
                                    check.names = check.names)
  samplegroups <- react_summary$SAMPLEGROUPS$Grouping
  data_size <- table(samplegroups)

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

    if(sum(grepl("Selected_co", rownames(react_summary$REACTSUM_AG))) > 0){
      data_sum <- react_summary$REACTSUM_AG[grep("Selected_co", rownames(react_summary$REACTSUM_AG)), ]
      data_freq <- react_summary$REACTFREQ_AG[grep("Selected_co", rownames(react_summary$REACTFREQ_AG)), ]
      data_freq_all <- react_summary$REACTFREQ_AG[-grep("Selected_co", rownames(react_summary$REACTFREQ_AG)), ]
    } else {
      data_sum <- react_summary$REACTSUM_AG
      data_freq <- react_summary$REACTFREQ_AG
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
      pdf(filename,
          width=ifelse(length(levels(samplegroups)) > 1, 20, 15), height=height, useDingbats=useDingbats)
      par(mfrow=c(4,3), mar=c(6,5,5,3), mgp=c(3,1,0))

      if(length(levels(samplegroups)) > 1){
        layout(matrix(1:16, nrow=4, byrow=T))
      } else {
        layout(rbind(c(1,2,2,3,3),
                     t(sapply(seq(3,9,3), function(x) c(1,2,2,3,3)+x))))
      }

      n=1
    for(a in agtoplot){
      tmp_ag <- colnames(data_cont)[a]
      print(paste("Plotting ag", n, "of", length(agtoplot),"(",tmp_ag,")"))
      n=n+1

      dens <- x$DENS[[tmp_ag]]
      tmp_which_co <- cutoffs$score[which(cutoffs$bead == tmp_ag)]*10+1
      tmp_cutoff <- cokey$xmad[tmp_which_co]

      # MADs Beeswarm, antigen score coloring
      plotdata <- data_cont[,tmp_ag]
      boxplot(plotdata~samplegroups, col="lightgrey", outcol=0, las=2,
              ylab="MADs [AU]", xaxt="n", xlab=NA,
              ylim=c(min(data_cont), ifelse(max(plotdata) > 50, max(plotdata), 50)))
      abline(h=tmp_cutoff, lty=2)
      beeswarm(plotdata~samplegroups, pch=16, corral="gutter", corralWidth=0.5, cex=0.8, add=T,
               pwcol=as.color(paste(cokey$color[data_score[,tmp_ag]*10+1]), 0.8))
      mtext(paste0("Above dashed line: "), adj=0.5,
            side=1, at=par("usr")[1], line=0, cex=0.7)
      mtext(paste0(data_sum[, grep(paste0("\\Q",tmp_ag,"\\E"), colnames(data_sum))], " of ", data_size,
                   " (", data_freq[, grep(paste0("\\Q",tmp_ag,"\\E"), colnames(data_freq))], "%)\n", levels(samplegroups)),
            side=1, at=1:length(levels(samplegroups)), line=1, cex=0.7)
      legend(par("usr")[2], par("usr")[4], legend=rev(c("<0",cokey$xmad[-1])),
             title=expression(bold("MADs cutoff")),
             pch=16, cex=0.6, bty="n", xjust=0.2, title.adj=4,
             col=rev(paste(cokey$color)), xpd=NA)
      mtext("Visualization of signals.", line=0.1, cex=0.65)

      # Histrogram & Density
      h <- hist(data_score[,tmp_ag], breaks=seq(min(cokey$score)-0.1,max(cokey$score)+0.1, 0.1), prob=T, right=F,
                main=NA, xlim=c(-0.1, max(cokey$score)+0.1), xlab="MADs cutoff\nDensity bandwidth = 0.1", xaxt="n")
      axis(1, labels=c("<0",cokey$xmad[-1]), at=h$breaks[-c(1, length(h$breaks))], cex.axis=0.8)
      abline(v=(tmp_which_co-1)/10, lty=2)
      mtext(tmp_ag, line=3, font=2)
      mtext("Distribution of binned values, \n algorithm assigned cutoff at dashed line.", line=0, cex=0.65)
      lines(dens,
            col="maroon")

      # Frequency
      plotdata <- data_freq_all[,grep(paste0("\\Q",tmp_ag,"\\E"), colnames(data_freq_all)), drop=F]
      if(length(levels(samplegroups)) > 1){
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
      abline(v=tmp_which_co, lty=2)
      matplot(plotdata, type="l", ylim=c(0,100), lty=1, lwd=2, add=T,
              col=groupcolors)

      if(length(levels(samplegroups)) > 1){
        mtext("Percentage of reactive samples per group at each exemplified cutoff.", line=1, cex=0.65)
      legend(par("usr")[1], par("usr")[4], horiz=T, yjust=0, xpd=NA, bty="n", cex=0.8,
             lty=1, lwd=2, col=groupcolors[seq_along(samplegroups)], legend=levels(samplegroups))
      } else {
        mtext("Percentage of reactive samples per group at each exemplified cutoff.", line=0.1, cex=0.65)
      }

      # Empty plot-spot until Fisher is added
      if(length(levels(samplegroups)) > 1){
        frame()
      }

      ##### NOT INCLUDED YET!
      # if(length(levels(samplegroups)) > 1){
      #   # Fisher line plot
      #   plotdata <- melt(fisher_p[[g]])
      #   plotdata <- plotdata[which(plotdata$Var2 == tmp_ag),]
      #   plotdata <- -log10(do.call(cbind,split(plotdata$value, plotdata$L1)))
      #
      #   matplot(plotdata, type="l", xaxt="n", las=1, lty=1, lwd=2,
      #           col=apply(comparisons[[g]], 2,
      #                     function(x) colorRampPalette(as.character(groupcolors))(3)[2]),
      #           ylim=c(ifelse(min(plotdata) < -log10(cofisher), min(plotdata), -log10(cofisher)),
      #                  ifelse(max(plotdata) > -log10(cofisher), max(plotdata), -log10(cofisher))),
      #           ylab="Fisher's exact test p-value (-log10)",
      #           xlab="Antigen specific MADs cutoff")
      #   axis(1, at=1:dim(cokey)[1], labels=cokey$xmad, cex.axis=0.8)
      #   abline(h=-log10(cofisher), lty=2)
      #   legend(par("usr")[1], par("usr")[4], horiz=T, legend=c(colnames(plotdata), paste0("-log10(",cofisher,")")),
      #          lty=c(rep(1, dim(plotdata)[2]), 2), lwd=1.7,
      #          col=c(apply(comparisons[[g]], 2,
      #                      function(x) colorRampPalette(as.character(groupcolors$Color[match(x, groupcolors$Subtype)]))(3)[2]),
      #                "black"),
      #          cex=0.55, xpd=NA, bty="n", yjust=0.1, seg.len=2.7)
      # }
    }
    dev.off()

    return(react_summary)
  }

#' Analysis summary
#'
#' Summarizes number of filtered/flagged beads/samples, amino acid lenghts and protein representation.
#'
#' @param x List with at least two elements, see Deatils for naming and content.
#' @details The x list needs to include at least the elements
#'
#'     SAMPLES = Sample info. Including column "sample_name" with LIMS-IDs.
#'
#'     BEADS = Beads info. Including columns:
#'
#'     "Type" with info about type of content on bead, at least including "PrEST" for PrESTs,
#'
#'     "PrEST.seq..aa." with amino acid sequences,
#'
#'     "Filtered" with filtering annotation, e.g. from other ap_-functions
#'
#'     "Flagged" with filtering annotation, e.g. from other ap_-functions
#'
#' @export

ap_summary <- function(x) {
  # Cohorts
  print("Samples per cohort")
  print(table(matrix(unlist(strsplit(as.character(x$SAMPLES$sample_name),"-")), ncol=2, byrow=T)[,1]))
  # Uniprot IDs
  print("Table of Uniprot IDs")
  print(sort(table(unlist(strsplit(as.character(x$BEADS$Uniprot[which(x$BEADS$Type == "PrEST")]), ";")))))
  print("Unique Uniprot IDs")
  print(length(table(unlist(strsplit(as.character(x$BEADS$Uniprot[which(x$BEADS$Type == "PrEST")]), ";")))))
  # Aminoacids
  print("Table of sequence lenghts")
  print(sort(apply(as.matrix(x$BEADS$PrEST.seq..aa.[grep("PrEST",x$BEADS$Type, ignore.case=T)],ncol=1), 1, nchar))) # Check aa-sequence length range
  print("min sequence lenghts")
  print(min(apply(as.matrix(x$BEADS$PrEST.seq..aa.[grep("PrEST",x$BEADS$Type, ignore.case=T)],ncol=1), 1, nchar))) # Check aa-sequence length range
  print("max sequence lenghts")
  print(max(apply(as.matrix(x$BEADS$PrEST.seq..aa.[grep("PrEST",x$BEADS$Type, ignore.case=T)],ncol=1), 1, nchar))) # Check aa-sequence length range
  print("median sequence lenghts")
  print(median(apply(as.matrix(x$BEADS$PrEST.seq..aa.[grep("PrEST",x$BEADS$Type, ignore.case=T)],ncol=1), 1, nchar))) # Check aa-sequence median length
  # Filtered & Flagged summary
  print("Filtered antigens")
  print(table(x$BEADS$Filtered))
  print("Flagged antigens")
  print(table(x$BEADS$Flagged))
  print("Filtered samples")
  print(table(x$SAMPLES$Filtered))
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
#' @param row.names if TRUE, the row names of the data frames are included in the Excel file worksheets.
#'     Deafult altered from \code{\link[WriteXLS:WriteXLS]{WriteXLS()}}.
#' @param filename string with filename and desired path, end with .xlsx.
#' @param ... arguments passed to \code{\link[WriteXLS:WriteXLS]{WriteXLS()}}.
#' @details The x list needs to include at least the elements specified under \code{elements}.
#'   It is recommended to append the output from \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} or
#'   \code{\link[rappp:ap_agresults]{ap_agresults()}} to the output from \code{\link[rappp:ap_norm2]{ap_norm2()}}
#'   and use the combined list as function input.
#'
#'   Exceptions for input element names:\cr
#'   - If an element is named BEADS in the input data, its name will be changed to ANTIGENS.
#'   Therefore, ANTIGENS is a default element to export while BEADS is not.\cr
#'   - Sums and frequencies at the cutoffs from \code{\link[rappp:ap_cutoff_selection2]{ap_cutoff_selection2()}}
#'   will be combined into a new element called REACTIVITIES
#'  if the output from \code{\link[rappp:ap_reactsummary2]{ap_reactsummary2()}} or
#'  \code{\link[rappp:ap_agresults]{ap_agresults()}} is included in the output
#'  Therefore, REACTIVITIES is a default element to export although it is not present in the input data.
#'
#'   If other list structures are used, it is most likely more convenient to just use
#'   \code{\link[WriteXLS:WriteXLS]{WriteXLS()}}, which this function is built on.
#'
#' @export

ap_excel <- function(x,
                     elements = c("MFI", "MADS", "SCORE", "BINARY_CO",
                                "REACTIVITIES", "ANTIGEN_CUTOFFS", "CUTOFF_KEY",
                                "ANTIGENS", "SAMPLES", "COUNT"),
                     filename = "DataOutput.xlsx",
                     row.names = TRUE,
                     ...) {
  excel <- x
  if("BEADS" %in% names(excel)){
  names(excel)[which(names(excel) == "BEADS")] <- "ANTIGENS"
  }

  if("REACTSUM_AG" %in% names(excel)){
  tmp_sum <- excel$REACTSUM_AG[grep("Selected", rownames(excel$REACTSUM_AG)),]
  rownames(tmp_sum) <- gsub("Selected_co","Sum", rownames(tmp_sum))
  }

  if("REACTSUM_AG" %in% names(excel)){
  tmp_freq <- excel$REACTFREQ_AG[grep("Selected", rownames(excel$REACTFREQ_AG)),]
  rownames(tmp_freq) <- gsub("Selected_co","Frequency", rownames(tmp_freq))
  }

  if(exists("tmp_sum") & exists("tmp_freq")){
  excel <- append(excel,
                  list(REACTIVITIES=data.frame(t(tmp_sum),t(tmp_freq))))
  }

  excel <- excel[match(elements,
                       names(excel))]

  WriteXLS(excel,
           ExcelFileName = filename,
           row.names = row.names, ...)
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
#' @return If MADlimits is provided a data.frame with three columns will be returned:
#'
#'    [,1] MADs cutoff value
#'
#'    [,2] Corresponding score value
#'
#'    [,3] Corresponding color using the Zissou1 palette in \code{\link[wesanderson]{wes_palette}}
#' @export

ap_cutoffs2image <- function(cutoffkey = NULL,
                             MADlimits = NULL) {

  if(!is.null(cutoffkey)){
    xmad_score <- cutoffkey
  } else if(!is.null(MADlimits)){
    xmad_score <- ap_cutoffs2(MADlimits = MADlimits)
  } else {
    warning("Either cutoffkey or MADlimits has to be provided.")
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
        if(start_aas[which.max(start_aas)-1] > (length(sequence)-aa_length+1)){
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


#' Align sequences
#'
#' Visualize how sequences align to a full length protein (exact matches).
#'
#' @param inputfile file with sequence information, see Details for needed formatting and columns.
#' @param ouputfile filename for the output file.
#' @param gene string with genename (or other name-identifier) for the alignment.
#' @param uniprot string with Uniprot ID for the alignment.
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
#' @export

align_sequences <- function(inputfile,
                            ouputfile,
                            gene = NULL,
                            uniprot = NULL,
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

  pdf(ouputfile,
      height=ifelse(length(all)*0.5 < 3, 5, length(all)*0.5),
      width=25, useDingbats=F)
  if(shouldtextplot){
    layout(matrix(c(1,1,2), ncol=1))
  } else {
    layout(matrix(1, ncol=1))
  }
  par(lend="butt", mar=c(2,2,2,25))
  y_at <- seq(0, length(all)*0.5-0.5, 0.5)
  ylabs <- ifelse(sequences$Type == "Peptide", paste0(sequences$Sequence1), "")
  aasteps <- ifelse(FLlength < 200, 10, ifelse(FLlength < 500, 50, 100))
  plot(0,0, col=0, xlim=c(1,FLlength), ylim=c(0,sum(sequences$Include)*0.5+0.5), xaxt="n",xlab=NA, las=1, yaxt="n", ylab=NA)
  axis(1, at=seq(1, FLlength, aasteps),
       labels=as.character(seq(1, FLlength, aasteps)),
       las=1)
  axis(4, at=y_at, labels=ylabs, las=2, cex=0.8)
  abline(v=seq(1, FLlength, aasteps), lty=2, col="lightgrey")
  abline(h=y_at, lty=2, col="lightgrey")
  n=0
  for(i in 1:length(all)){
    if(length(all[[i]]) != 2){
      x0=which(rollapply(all[[1]], length(all[[i]]), identical, all[[i]]))
      x1=(which(rollapply(all[[1]], length(all[[i]]), identical, all[[i]]))+length(all[[i]])-1)
      segments(x0=x0, y0=n,
               x1=x1, y1=n,
               lwd=5, col=as.character(sequences$Color[i]))
      text(x=ifelse(x0[1] < FLlength*0.75, x0[1], x1[1]), y=(n+0.1), adj=ifelse(x0[1] < FLlength*0.75, 0, 1),
           labels=paste0(names(all)[i]," (",x0,"-",x1,", ",length(all[[i]])," aa)"),
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
      text(x=ifelse(x0[1] < FLlength*0.75, x0[1], x1[1]), y=(n+0.1), adj=ifelse(x0[1] < FLlength*0.75, 0, 1),
           labels=paste0(names(all)[i]," (",x0_1,"-",x1_1,", ",x0_2,"-",x1_2,", ",
                         length(all[[i]][[1]]),"+",length(all[[i]][[2]])," aa)"),
           cex=0.9, offset=0, xpd=NA)
    }
    n=n+0.5
  }
  mtext(paste0(gene," (Uniprot ",uniprot,")"), side=3, font=2)

  if(shouldtextplot){
    ap_textplot(sequences[,textplotcolumns][-1,], show.rownames=F)
  }

  dev.off()
}
