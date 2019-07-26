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

