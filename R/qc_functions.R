#' Coupling efficency test
#'
#' Flags beads with signal similar to empty bead in coupling test, produces plot if wanted.
#'
#'
#' @param x List with at least three elements:
#'     CT = coupling test mfi,
#'     BEADS = Beads info (including Type-column with PrEST for PrESTs),
#'     FILTERINFO = Vector with info on whichfilter steps has been done.
#' @param empty_bead Column index for empty bead.
#' @param empty_co_multiple Number of sd above empty for cutoff.
#' @param shouldplot Logical, should a plot be made?
#' @param ... Further arguments passed do \link[stats]{mean} and \link[stats]{sd}
#' @return
#' @export

ap_ct <- function(x, empty_bead, empty_co_multiple=3, shouldplot=T, ...) {

    empty_co <- mean(x$CT[,empty_bead], ...) + empty_co_multiple*sd(x$CT[,empty_bead], ...)

    if(shouldplot){
      pdf(paste0(analysis,"Coupling-efficiency_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
          width=30, height=6, useDingbats=F)
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
#' @param x List with at least three elements:
#'     MFI = assay mfi,
#'     SAMPLES = Sample info
#'     BEADS = Beads info (including Type-column with PrEST for PrESTs),
#'     FILTERINFO = Vector with info on whichfilter steps has been done.
#' @param IgX_bead Column index for empty bead.
#' @param IgG_cutoff MFI cutoff value for filtering.
#' @param shouldplot Logical, should a plot be made?
#' @param ... Further arguments passed do \link[stats]{mean} and \link[stats]{sd}
#' @return
#' @export

ap_igx <- function(x, IgX_bead, IgG_cutoff=5000, shouldplot=T, ...) {
  cosfac <- c(3, -3)
  which_lowIgG <- rep(list(NULL), length(assay_list))
  which_hIgGremove <- rep(list(NULL), length(assay_list))
  for(l in 1:length(assay_list)){
    plotdata <- unlist(assay_list[[l]][,grep("hIgG", colnames(assay_list[[l]]), ignore.case=T)])
    sampledata <- samples_list[[l]]
    SamplesNames <- rownames(sampledata) # INPUT NEEDED: Change to column with blank/pool/replicate information
    AssayNum <- samples_list[[l]]$AssayNum # INPUT NEEDED: Change to column with assay run information information, or create a vector here
    cosIgG <- median(plotdata, na.rm=T)+cosfac*mad(plotdata, constant = 1, na.rm=T)
    tmp <- plotdata[grepl("empty|blank|buffer", SamplesNames, ignore.case=T)]
    cosIgG <- c(cosIgG, 5000) #max(tmp)+3*mad(tmp, constant = 1))

    if(shouldplot){
      pdf(paste0(analysis,names(count_list)[l],"_ahIgG_",format(Sys.time(),"%Y-%m-%d_%H%M%S"),".pdf"),
          width=10, height=6, useDingbats=F)
      layout(matrix(c(1,1,1,6,7,
                      1,1,1,2,3,
                      1,1,1,4,5), nrow=3, byrow=T))
      par(mar=c(5,5,4,4))

      plot(1:length(plotdata), plotdata, cex=0.5, pch=c(16:18)[AssayNum],
           xlab="Samples in analysis order",ylab="Signal intensity (MFI)",main="Total hIgG",
           col=ifelse(grepl("empty|blank|buffer", SamplesNames, ignore.case=T),2,
                      ifelse(grepl("pool|rep|mix|commercial", SamplesNames, ignore.case=T),5, 4)))
      if(length(unique(AssayNum)) > 1){
        legend(par("usr")[1], par("usr")[4], horiz=T, yjust=0.1, bty="n",
               legend=c("Sample","Replicate","Buffer", paste("Assay", unique(AssayNum))),
               cex=0.8,
               pch=c(rep(NA, 3), c(16:18)[unique(AssayNum)]),
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
    }

    # Display samples with high but still outliers
    if(shouldplot){
      if(length(which(plotdata<cosIgG[2] & plotdata>cosIgG[3])) > 0) {
        # INPUT NEEDED: Change to match which columns you want to display, number of columns is optional
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
          mtext(paste("anti-hIgG MFI between",cosIgG[3],"&", cosIgG[2]), font=2, cex=0.5, xpd=NA, at=-0.5)
        } else {
          textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top", cex=0.6)
          mtext(paste("anti-hIgG MFI between",cosIgG[3],"&", cosIgG[2]), font=2, cex=0.5)
          frame()
        }

      } else {
        frame()
        frame()
      }
    }

    # Display and remove samples with low total IgG signal
    if(length(which(plotdata<cosIgG[3])) > 0) {
      which_lowIgG[[l]] <- which(plotdata<cosIgG[3])
      # INPUT NEEDED: Change to match which columns you want to display, number of columns is optional
      plottext <- data.frame(Well384=sampledata$Well384,
                             InternalID=sampledata$Sample,
                             Subject=sampledata$tube_label,
                             MFI=plotdata)[which(plotdata<cosIgG[3]),]
      plottext <- plottext[order(plottext$MFI, decreasing=T),]

      if(shouldplot){
        if(dim(plottext)[1] > 20){
          textplot(plottext[1:20,],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          textplot(plottext[21:dim(plottext)[1],],
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          mtext(paste("anti-hIgG MFI below",cosIgG[3]), font=2, cex=0.5, xpd=NA, at=-0.5)
        } else {
          textplot(plottext,
                   halign="left", show.rownames=F, hadj=0, cmar=1, valign="top")
          mtext(paste("anti-hIgG MFI below",cosIgG[3]), font=2, cex=0.5)
          frame()
        }
      }

    } else {
      if(shouldplot){
        frame()
        frame()
      }
    }
    if(shouldplot){ dev.off() }

    # Sample list filtered
    tmp_remove <- rownames(sampledata)[which_lowIgG[[l]]] ; tmp_remove <- tmp_remove[-grep("empty", tmp_remove, ignore.case=T)]
    samples_list[[l]]$Filtered <- ifelse(rownames(samples_list[[l]]) %in% tmp_remove,
                                         paste0(samples_list[[l]]$Filtered, ", hIgG"),
                                         samples_list[[l]]$Filtered)
    samples_list[[l]]$Filtered <- gsub("^, ", "", samples_list[[l]]$Filtered)
  }

  x$FILTERINFO <- c(x$FILTERINFO, "CouplingEfficiency")

  return(x)
}
