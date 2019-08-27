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
#'     BEADS = Beads info (including Type-column with PrEST for PrESTs),
#'
#'     FILTERINFO = Vector with info on which filter steps has been done.
#' @return Updated input x with relevant filtering info and a pdf with plot (if \code{shouldplot=T}).
#' @export

ap_ct <- function(x, empty_bead, empty_co_multiple=3,
                  shouldplot=T, filename="coupling_efficiency.pdf", width=30, height=6, useDingbats=F, ...) {

    empty_co <- mean(x$CT[,empty_bead], ...) + empty_co_multiple*sd(x$CT[,empty_bead], ...)

    if(shouldplot){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)
      par(mar=c(12,4,2,8), cex.axis=0.8)
      layout(matrix(c(1,1,1,2), nrow=1))
      bs=beeswarm(x$CT, pch=16, las=2, corral="gutter", xaxt="n",
                  main="Copupling efficiency test", ylab="Signal intensity [MFI]",
                  pwcol=rep(ifelse(grepl("empty|bare|blank", colnames(x$CT),ignore.case=T), "orange",
                                   ifelse(grepl("his6abp", colnames(x$CT),ignore.case=T), "darkgreen",
                                          ifelse(grepl("higg", colnames(x$CT),ignore.case=T), "blue",
                                                 ifelse(grepl("ebna", colnames(x$CT),ignore.case=T), "purple",
                                                        ifelse(apply(x$CT, 2, mean) < empty_co, "red", "darkgrey"))))),
                            each=dim(x$CT)[1]))

      axis(1, at=bs$x[seq(1, dim(bs)[1], 3)], labels=unique(bs$x.orig), cex.axis=0.5, las=2, tick=F)
      abline(v=seq(min(bs$x)-0.5, max(bs$x)+0.5, 1), lty=2, col="lightgrey")
      legend(par("usr")[2], par("usr")[4],
             legend=c("Empty", "His6ABP", "ahIgG", "EBNA1", "Flagged"),
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
#'     "Sample" with sample names, preferably LIMS-IDs, where
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
                   shouldplot=T, filename="anti-humanIgX.pdf", width=10, height=6, useDingbats=F) {

    plotdata <- unlist(x$MFI[,IgX_bead])
    sampledata <- x$SAMPLES
    SamplesNames <- sampledata$Sample
    AssayNum <- sampledata$AssayNum

    cosIgG <- median(plotdata, na.rm=T)+cosfac*mad(plotdata, constant = 1, na.rm=T)
    tmp <- plotdata[grepl("empty|blank|buffer", SamplesNames, ignore.case=T)]
    cosIgG <- c(cosIgG, IgX_cutoff)

    which_lowIgG <- which(plotdata<cosIgG[3])

    if(shouldplot){
      pdf(filename, width=width, height=height, useDingbats=useDingbats)
      layout(matrix(c(1,1,1,6,7,
                      1,1,1,2,3,
                      1,1,1,4,5), nrow=3, byrow=T))
      par(mar=c(5,5,4,4))

      plot(1:length(plotdata), plotdata, cex=0.5, pch=c(16:18,6,8)[AssayNum],
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

      # Display samples with high but still outliers
      if(length(which(plotdata<cosIgG[2] & plotdata>cosIgG[3])) > 0) {
        plottext <- data.frame(AssayWell=sampledata$AssayWell,
                               InternalID=sampledata$Sample,
                               Subject=sampledata$tube_label,
                               MFI=plotdata)[which(plotdata<cosIgG[2] & plotdata>cosIgG[3]),]
        plottext <- plottext[order(plottext$MFI, decreasing=T),]

        if(dim(plottext)[1] > 20){
          ap_textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          ap_textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          mtext(paste0("anti-hIg", IgType, " MFI between ",cosIgG[3]," & ", cosIgG[2]), font=2, cex=0.5, xpd=NA, at=-0.5)
        } else {
          ap_textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          mtext(paste0("anti-hIg", IgType, " MFI between ",cosIgG[3]," & ", cosIgG[2]), font=2, cex=0.5)
          frame()
        }

      } else {
        frame()
        frame()
      }

      # Display and remove samples with low total IgG signal
      if(length(which_lowIgG) > 0) {
        plottext <- data.frame(AssayWell=sampledata$AssayWell,
                               InternalID=sampledata$Sample,
                               Subject=sampledata$tube_label,
                               MFI=plotdata)[which_lowIgG,]
        plottext <- plottext[order(plottext$MFI, decreasing=T),]

        if(dim(plottext)[1] > 20){
          ap_textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          ap_textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          mtext(paste0("anti-hIg", IgType, " MFI below ",cosIgG[3]), font=2, cex=0.5, xpd=NA, at=-0.5)
        } else {
          ap_textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          mtext(paste0("anti-hIg", IgType, " MFI below ",cosIgG[3]), font=2, cex=0.5)
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
#'     "Sample" with sample names, preferably LIMS-IDs, where
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
#'     "Sample" with sample names, preferably LIMS-IDs, where
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
    replicates <- data[[l]][grep("pool|rep|mix|commercial", samples[[l]]$Sample, ignore.case=T),]
    nrreplicates[l] <- dim(replicates)[1]
    CVs[,l] <- apply(replicates, 2, cv, digits=5, na.rm=T)

    cor_samp_s[[l]] <- cor(t(replicates), method="spearman", use="pairwise.complete.obs")
    cor_samp_s[[l]] <- cor_samp_s[[l]][upper.tri(cor_samp_s[[l]])]

    cor_samp_p[[l]] <- cor(t(log(replicates)), method="pearson", use="pairwise.complete.obs")^2
    cor_samp_p[[l]] <- cor_samp_p[[l]][upper.tri(cor_samp_p[[l]])]

    # Iterate over random sets of samples
    CVs_r <- matrix(NA, ncol=iter, nrow=dim(x$MFI)[2])
    tmp_data <- data[[l]][-grep("pool|rep|mix|commercial", samples[[l]]$Sample, ignore.case=T),]
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
#' @param z a list of data
#' @param perp a vector of values to test as perplexities, eg. c(2,5,10,50).
#'     sqrt(Nsamples) will always be included.
#' @param iterations number of iterations per tSNE
#' @param groups grouping for colors and lines
#' @param names text to display as point labels
#' @param legend TRUE or FALSE
#' @param legendname Title of legend
#' @param main Title of page
#' @param filename String with filename and desired path, end with .pdf
#' @param height height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.
#' @details The x list needs to include at least the element
#'     MFI = assay mfi,
#'
#' @export

tsne_perp <- function(z, perp=c(2,5,10,50), iterations=1000, groups, names,
                      legend=T, legendname=NULL, main=NULL,
                      filename="t-SNE_perplezities.pdf", height=16, useDingbats=F) {

  g <- rep(list(NULL), length(z)*(length(perp)+1))
  n=1
  for(l in 1:length(z)){
    perp <- sort(c(perp, round(sqrt(dim(z[[l]])[1]))))
    for(p in 1:length(perp)){
      tsne_tmp = Rtsne(z[[l]], check_duplicates=F, perplezity=perp[p], maz_iter=iterations)

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
#'     COKEY = Cutoff key as data.frame with cutoff values, scores and colors.
#'
#'     BEADS = Beads info, if any should be excluded then these should be annotated in a column called "Filtered".
#'     Any beads with no text (ie. "") or "NegControl" in such column will be included in the transformation.
#'
#' @export

ap_negbeads <- function(x,
                        filename="neg-control-beads.pdf", width=15, height=10, useDingbats=F){

  pdf(filename, width=width, height=height, useDingbats=useDingbats)

  layout(matrix(c(1,1,2,3), ncol=2, byrow=T))
  par(mar=c(4,4,4,5))

  if("Filtered" %in% colnames(x$BEADS)){
    plotdata <- x$MFI[, which(x$BEADS$Filtered == "" | grepl("NegControl", x$BEADS$Filtered))]
    plotdata_score <- x$SCORE[, which(x$BEADS$Filtered == "" | grepl("NegControl", x$BEADS$Filtered))]
  } else {
    plotdata <- x$MFI
    plotdata_score <- x$SCORE
  }

  plotcolor <- x$COKEY

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
