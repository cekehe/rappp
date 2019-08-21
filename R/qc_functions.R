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
#' @param width,height Width and height for pdf, see \link[grDevices]{pdf}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \link[grDevices]{pdf}.
#' @param ... Further arguments passed to \link[stats]{mean} and \link[stats]{sd}.
#' @details The x list needs to include at least the elements
#'
#'     CT = coupling test mfi,
#'
#'     BEADS = Beads info (including Type-column with PrEST for PrESTs),
#'
#'     FILTERINFO = Vector with info on whichfilter steps has been done.
#' @return
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
        textplot(tmp_text, mar=c(2,2,1,2),
                 show.rownames=F, show.colnames=F, hadj=0, valign="top", cex=0.8)
        mtext("Protein fragments with low coupling efficiency signal", font=2, cex=0.9, xpd=NA)
      } else {
        frame()
        mtext("No protein fragments displayed \n low coupling efficiency signal.", font=2, cex=0.7, line=-3)
      }

      dev.off()
    }

    x$BEADS <- data.frame(Flagged=ifelse(apply(x$CT, 2, mean) < empty_co &
                                           grepl("PrEST", x$BEADS$Type, ignore.case=T),
                                         "Coupling", ""),
                          x$BEADS)
    x$BEADS$Flagged <- paste(x$BEADS$Flagged)

    x$FILTERINFO <- c(x$FILTERINFO, "CouplingEfficiency")

    return(x)
  }


#' Loading control (FAR FROM DONE)
#'
#' Filter samples with low MFI for the anti-human IgX bead.
#'
#'
#' @param x List with at least four elements, see Deatils for naming and content.
#' @param IgX_bead Column index for empty bead.
#' @param IgType Which Imunoglobulin is measured, default is G.
#' @param IgX_cutoff MFI cutoff value for filtering.
#' @param cosfac Median absolute deviation multipliers in vector c(upper, lower),
#'     for drawing lines and detecting potential outliers.
#' @param shouldplot Logical, should a plot be made?
#' @param filename String with filename and desired path, end with .pdf
#' @param width,height Width and height for pdf, see \link[grDevices]{pdf}.
#' @param useDingbats Logical. Default is \code{FALSE}, compared to in default \link[grDevices]{pdf}.
#' @details The x list needs to include at least the elements:
#'
#'     MFI = assay mfi,
#'
#'     SAMPLES = Sample info. See below for required columns.
#'
#'     BEADS = Beads info (including Type-column with PrEST for PrESTs),
#'
#'     FILTERINFO = Vector with info on whichfilter steps has been done.
#'
#' The SAMPLES element needs at least the columns:
#'
#'     "Sample" with sample names, preferably LIMS-IDs, where
#'     replicates (named with one of pool|rep|mix|commercial)
#'     and blanks (named with one of empty|blank|buffer) are also stated,
#'
#'     "AssayNum" with assay number (vector with 1s if only one assay, support for up to 5 assys in one plot),
#'
#'     "Well384" with Well IDs, e.g A01, B01 etc.,
#'
#'     "tube_label" with alternative sample names, eg. from collaborator,
#'
#' @return
#' @export

ap_igx <- function(x, IgX_bead, IgType="G", IgX_cutoff=5000, cosfac=c(3, -3),
                   shouldplot=T, filename, width=10, height=6, useDingbats=F) {

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
        plottext <- data.frame(Well384=sampledata$Well384,
                               InternalID=sampledata$Sample,
                               Subject=sampledata$tube_label,
                               MFI=plotdata)[which(plotdata<cosIgG[2] & plotdata>cosIgG[3]),]
        plottext <- plottext[order(plottext$MFI, decreasing=T),]

        if(dim(plottext)[1] > 20){
          textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          mtext(paste("anti-hIg", IgType, " MFI between",cosIgG[3],"&", cosIgG[2]), font=2, cex=0.5, xpd=NA, at=-0.5)
        } else {
          textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          mtext(paste("anti-hIg", IgType, " MFI between",cosIgG[3],"&", cosIgG[2]), font=2, cex=0.5)
          frame()
        }

      } else {
        frame()
        frame()
      }

      # Display and remove samples with low total IgG signal
      if(length(which_lowIgG) > 0) {
        plottext <- data.frame(Well384=sampledata$Well384,
                               InternalID=sampledata$Sample,
                               Subject=sampledata$tube_label,
                               MFI=plotdata)[which_lowIgG,]
        plottext <- plottext[order(plottext$MFI, decreasing=T),]

        if(dim(plottext)[1] > 20){
          textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          mtext(paste("anti-hIg", IgType, " MFI below",cosIgG[3]), font=2, cex=0.5, xpd=NA, at=-0.5)
        } else {
          textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          mtext(paste("anti-hIg", IgType, " MFI below",cosIgG[3]), font=2, cex=0.5)
          frame()
        }

      } else {
        frame()
        frame()
      }
      dev.off()
    }

    # Sample list filtered
    tmp_remove <- rownames(sampledata)[which_lowIgG]
    tmp_remove <- tmp_remove[-grep("empty", tmp_remove, ignore.case=T)]
    x$SAMPLES$Filtered <- ifelse(rownames(x$SAMPLES) %in% tmp_remove,
                                         paste0(x$SAMPLES$Filtered, ", hIg", IgType, ""),
                                         x$SAMPLES$Filtered)
    x$SAMPLES$Filtered <- gsub("^, ", "", x$SAMPLES$Filtered)

  x$FILTERINFO <- c(x$FILTERINFO, "CouplingEfficiency")

  return(x)
}
