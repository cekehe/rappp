#' Display text information in a graphics plot.
#'
#' This function displays text output in a graphics window.
#' It is the equivalent of 'print' except that the output is displayed as a plot.
#'
#' @details Based on \code{\link[gplots:textplot]{textplot()}} with slightly altered code.
#' Altered code found onlne (https://gist.github.com/johncolby/1482973).
#' Additional alteration done by CH, including making working methods, see code.
#'
#' The input should be one of classes \code{matrix}, \code{data.frame}, \code{vectors} longer than 1,
#' single \code{character} string, single \code{integer} or single \code{numeric} value.
#'
#' @export

setGeneric("ap_textplot",
           function(object, halign="center", valign="center", cex, cmar=1.5, ... ) {
             if(is.matrix(object) || is.data.frame(object) || is.vector(object) && length(object)>1) {

             } else if(is.character(object) || is.numeric(object) || is.integer(object)) {

             } else {
               stop("Object needs to be a matrix, data.frame, vector longer than 1 or
                    single character string, integer or numeric value.") }
           }
)

# Function to be used for classes matrix, data.frame and vectors longer than 1. (setMethod individually below)
.tables <- function(object,
                    halign=c("center","left","right"),
                    valign=c("center","top","bottom"),
                    cex, cmar=cmar, rmar=0.5,
                    show.rownames=TRUE, show.colnames=TRUE,
                    hadj=1,
                    vadj=1,
                    mar= c(1,1,4,1)+0.1,
                    col.data=par("col"),
                    col.rownames=par("col"),
                    col.colnames=par("col"),
                    ... ){

  if(is.vector(object) && length(object)>1)
    object <- t(as.matrix(object))
  else
    object <- as.matrix(object)


  # check dimensions of col.data, col.rownames, col.colnames
  if(length(col.data)==1)
    col.data <- matrix(col.data, nrow=nrow(object), ncol=ncol(object))
  else
    if( nrow(col.data)!=nrow(object) || ncol(col.data)!=ncol(object) )
      stop("Dimensions of 'col.data' do not match dimensions of 'object'.")

  if(length(col.rownames)==1)
    col.rownames <- rep(col.rownames, nrow(object))

  if(length(col.colnames)==1)
    if(show.rownames)
      col.colnames <- rep(col.colnames, ncol(object)+1)
    else
      col.colnames <- rep(col.colnames, ncol(object))

    halign=match.arg(halign)
    valign=match.arg(valign)

    opar <- par()[c("mar","xpd","cex")]
    on.exit( par(opar) )
    par(mar=mar, xpd=FALSE )

    # setup plot area
    plot.new()
    plot.window(xlim=c(0,1),ylim=c(0,1), log = "", asp=NA)



    # add 'r-style' row and column labels if not present
    if( is.null(colnames(object) ) )
      colnames(object) <- paste( "[,", 1:ncol(object), "]", sep="" )
    if( is.null(rownames(object)) )
      rownames(object) <- paste( "[", 1:nrow(object), ",]", sep="")


    # extend the matrix to include row and column labels
    if( show.rownames )
    {
      object <- cbind( rownames(object), object )
      col.data <- cbind( col.rownames, col.data )

    }
    if( show.colnames )
    {
      object <- rbind( colnames(object), object )
      col.data <- rbind( col.colnames, col.data )
    }

    # set the character size
    if( missing(cex) )
    {
      cex <- 1.0
      lastloop <- FALSE
    }
    else
    {
      lastloop <- TRUE
    }

    iterations <- 10000
    for (i in 1:iterations) # max 20 iteration in original, increased in ap_textplot to better fine tune char size. /CH
    {
      if(!lastloop){ #added to let xpos part be run even if cex is predefined. /CH
        # oldcex <- cex # Not necessary with changed loop break check. /CH

        width  <- sum(
          apply( object, 2,
                 function(x) max(strwidth(x,cex=cex) ) )
        ) +
          strwidth('M', cex=cex) * cmar * ncol(object)

        height <- strheight('M', cex=cex) * nrow(object) * (1 + rmar)

        # if(lastloop) break # Moved to below xpos etc /CH

        cex <- cex / max(width,height)
      }

      # Changed check for loop, now based on total column width below /CH
      #   if (abs(oldcex - cex) < 0.001)
      #   {
      #     lastloop <- TRUE
      #   }
      # }

      # compute the individual row and column heights
      rowheight<-strheight("M",cex=cex) * (1 + rmar)  # Changed from W to M /CH
      colwidth<- apply( object, 2, function(XX) max(strwidth(XX, cex=cex)) ) +
        strwidth("M")*cmar # Changed from W to M /CH

      width  <- sum(colwidth)
      height <- rowheight*nrow(object)

      # setup x alignment
      if(halign=="left")
        xpos <- 0
      else if(halign=="center")
        xpos <- 0 + (1-width)/2
      else #if(halign=="right")
        xpos <- 0 + (1-width)

      # setup y alignment
      if(valign=="top")
        ypos <- 1
      else if (valign=="center")
        ypos <- 1 - (1-height)/2
      else #if (valign=="bottom")
        ypos <- 0 + height

      # Added instead of (abs(oldcex - cex) < 0.001) check /CH
      if ((xpos + width) <= 1 & (ypos - height) >= 0)
      {
        lastloop <- TRUE
      }

      if(lastloop) break # moved from above /CH
    }

    x <- xpos
    y <- ypos

    # iterate across elements, plotting them
    xpos<-x
    for(i in 1:ncol(object)) {
      xpos <- xpos + hadj*colwidth[i]
      for(j in 1:nrow(object)) {
        ypos<-y-(j-1)*rowheight
        if( (show.rownames && i==1) || (show.colnames && j==1) )
          text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=2,
               col=col.data[j,i], ... )
        else
          text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=1,
               col=col.data[j,i], ... )
      }
      xpos <- xpos + (1-hadj)*colwidth[i]
    }

    par(opar)
}

setMethod("ap_textplot", "matrix", .tables)
setMethod("ap_textplot", "data.frame", .tables)
setMethod("ap_textplot", "vector", .tables)

## Function to be used for classes character, numeric and integer, ie. vectors of length 1. (setMethod individually below)
.singles <- function (object,
                      halign = c("center", "left", "right"),
                      valign = c("center", "top", "bottom"),
                      cex, fixed.width=TRUE,
                      cspace=1,
                      lspace=1,
                      mar=c(0,0,3,0)+0.1,
                      tab.width=8,
                      ...)
{
  object <- paste(object,collapse="\n",sep="")
  object <- replaceTabs(object, width=tab.width)

  halign = match.arg(halign)
  valign = match.arg(valign)
  plot.new()

  opar <- par()[c("mar","xpd","cex","family")]
  on.exit( par(opar) )

  par(mar=mar,xpd=FALSE )
  if(fixed.width)
    par(family="mono")

  plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "", asp = NA)

  slist   <- unlist(lapply(object, function(x) strsplit(x,'\n')))
  slist   <- lapply(slist, function(x) unlist(strsplit(x,'')))

  slen    <- sapply(slist, length)
  slines  <- length(slist)

  if (missing(cex))
  {
    lastloop <- FALSE
    cex <- 1
  }
  else
    lastloop <- TRUE


  for (i in 1:20)
  {
    oldcex <- cex
    #cat("cex=",cex,"\n")
    #cat("i=",i,"\n")
    #cat("calculating width...")
    cwidth  <- max(sapply(unlist(slist), strwidth,  cex=cex)) * cspace
    #cat("done.\n")
    #cat("calculating height...")
    cheight <- max(sapply(unlist(slist), strheight, cex=cex)) * ( lspace + 0.5 )
    #cat("done.\n")

    width <- strwidth(object, cex=cex)
    height <- strheight(object, cex=cex)

    if(lastloop) break

    cex <- cex  / max(width, height)

    if (abs(oldcex - cex) < 0.001)
    {
      lastloop <- TRUE
    }

  }

  if (halign == "left")
    xpos <- 0
  else if (halign == "center")
    xpos <- 0 + (1 - width)/2
  else xpos <- 0 + (1 - width)

  if (valign == "top")
    ypos <- 1
  else if (valign == "center")
    ypos <- 1 - (1 - height)/2
  else ypos <- 1 - (1 - height)

  text(x=xpos, y=ypos, labels=object, adj=c(0,1),
       cex=cex, ...)

  par(opar)
  invisible(cex)
}

setMethod("ap_textplot", "character", .singles)
setMethod("ap_textplot", "numeric", .singles)
setMethod("ap_textplot", "integer", .singles)


## Old code which doesn't work as methods properly:
# ap_textplot <- function(object, halign="center", valign="center", cex, ... ){
#   UseMethod('ap_textplot')
# }

#
# ap_textplot.default <- function(object,
#                              halign=c("center","left","right"),
#                              valign=c("center","top","bottom"),
#                              cex, ... ){
#
#   if (is.matrix(object) || (is.vector(object) && length(object)>1) )
#     return(ap_textplot.matrix(object, halign, valign, cex, ... ))
#
#   halign <- match.arg(halign)
#   valign <- match.arg(valign)
#
#   ap_textplot.character(object, halign,  valign, cex, ...)
# }

# ap_textplot.data.frame <- function(object,
#                                 halign=c("center","left","right"),
#                                 valign=c("center","top","bottom"),
#                                 cex, ... ) {
#   ap_textplot.matrix(object, halign, valign, cex, ... )
# }


# ap_textplot.matrix <- function(object,
#                             halign=c("center","left","right"),
#                             valign=c("center","top","bottom"),
#                             cex, cmar=2, rmar=0.5,
#                             show.rownames=TRUE, show.colnames=TRUE,
#                             hadj=1,
#                             vadj=1,
#                             mar= c(1,1,4,1)+0.1,
#                             col.data=par("col"),
#                             col.rownames=par("col"),
#                             col.colnames=par("col"),
#                             ... ){
#
#   if(is.vector(object))
#     object <- t(as.matrix(object))
#   else
#     object <- as.matrix(object)
#
#   # check dimensions of col.data, col.rownames, col.colnames
#   if(length(col.data)==1)
#     col.data <- matrix(col.data, nrow=nrow(object), ncol=ncol(object))
#   else
#     if( nrow(col.data)!=nrow(object) || ncol(col.data)!=ncol(object) )
#       stop("Dimensions of 'col.data' do not match dimensions of 'object'.")
#
#   if(length(col.rownames)==1)
#     col.rownames <- rep(col.rownames, nrow(object))
#
#   if(length(col.colnames)==1)
#     if(show.rownames)
#       col.colnames <- rep(col.colnames, ncol(object)+1)
#   else
#     col.colnames <- rep(col.colnames, ncol(object))
#
#   halign=match.arg(halign)
#   valign=match.arg(valign)
#
#   opar <- par()[c("mar","xpd","cex")]
#   on.exit( par(opar) )
#   par(mar=mar, xpd=FALSE )
#
#   # setup plot area
#   plot.new()
#   plot.window(xlim=c(0,1),ylim=c(0,1), log = "", asp=NA)
#
#
#
#   # add 'r-style' row and column labels if not present
#   if( is.null(colnames(object) ) )
#     colnames(object) <- paste( "[,", 1:ncol(object), "]", sep="" )
#   if( is.null(rownames(object)) )
#     rownames(object) <- paste( "[", 1:nrow(object), ",]", sep="")
#
#
#   # extend the matrix to include row and column labels
#   if( show.rownames )
#   {
#     object <- cbind( rownames(object), object )
#     col.data <- cbind( col.rownames, col.data )
#
#   }
#   if( show.colnames )
#   {
#     object <- rbind( colnames(object), object )
#     col.data <- rbind( col.colnames, col.data )
#   }
#
#   # set the character size
#   if( missing(cex) )
#   {
#     cex <- 1.0
#     lastloop <- FALSE
#   }
#   else
#   {
#     lastloop <- TRUE
#   }
#
#   iterations <- 10000
#   for (i in 1:iterations) # max 20 iteration in original, increased in ap_textplot to better fine tune char size. /CH
#   {
#     if(!lastloop){ #added to let xpos part be run even if cex is predefined. /CH
#      # oldcex <- cex # Not necessary with changed loop break check. /CH
#
#     width  <- sum(
#       apply( object, 2,
#              function(x) max(strwidth(x,cex=cex) ) )
#     ) +
#       strwidth('M', cex=cex) * cmar * ncol(object)
#
#     height <- strheight('M', cex=cex) * nrow(object) * (1 + rmar)
#
#     # if(lastloop) break # Moved to below xpos etc /CH
#
#     cex <- cex / max(width,height)
#     }
#
#     # Changed check for loop, now based on total column width below /CH
#     #   if (abs(oldcex - cex) < 0.001)
#     #   {
#     #     lastloop <- TRUE
#     #   }
#     # }
#
#     # compute the individual row and column heights
#     rowheight<-strheight("M",cex=cex) * (1 + rmar)  # Changed from W to M /CH
#     colwidth<- apply( object, 2, function(XX) max(strwidth(XX, cex=cex)) ) +
#       strwidth("M")*cmar # Changed from W to M /CH
#
#     width  <- sum(colwidth)
#     height <- rowheight*nrow(object)
#
#     # setup x alignment
#     if(halign=="left")
#       xpos <- 0
#     else if(halign=="center")
#       xpos <- 0 + (1-width)/2
#     else #if(halign=="right")
#       xpos <- 0 + (1-width)
#
#     # setup y alignment
#     if(valign=="top")
#       ypos <- 1
#     else if (valign=="center")
#       ypos <- 1 - (1-height)/2
#     else #if (valign=="bottom")
#       ypos <- 0 + height
#
#     # Added instead of (abs(oldcex - cex) < 0.001) check /CH
#     if ((xpos + width) <= 1 & (ypos - height) >= 0)
#     {
#       lastloop <- TRUE
#     }
#
#     if(lastloop) break # moved from above /CH
#   }
#
#   x <- xpos
#   y <- ypos
#
#   # iterate across elements, plotting them
#   xpos<-x
#   for(i in 1:ncol(object)) {
#     xpos <- xpos + hadj*colwidth[i]
#     for(j in 1:nrow(object)) {
#       ypos<-y-(j-1)*rowheight
#       if( (show.rownames && i==1) || (show.colnames && j==1) )
#         text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=2,
#              col=col.data[j,i], ... )
#       else
#         text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=1,
#              col=col.data[j,i], ... )
#     }
#     xpos <- xpos + (1-hadj)*colwidth[i]
#   }
#
#   par(opar)
# }

# ap_textplot.character <- function (object,
#                                 halign = c("center", "left", "right"),
#                                 valign = c("center", "top", "bottom"),
#                                 cex, fixed.width=TRUE,
#                                 cspace=1,
#                                 lspace=1,
#                                 mar=c(0,0,3,0)+0.1,
#                                 tab.width=8,
#                                 ...)
# {
#   object <- paste(object,collapse="\n",sep="")
#   object <- replaceTabs(object, width=tab.width)
#
#   halign = match.arg(halign)
#   valign = match.arg(valign)
#   plot.new()
#
#   opar <- par()[c("mar","xpd","cex","family")]
#   on.exit( par(opar) )
#
#   par(mar=mar,xpd=FALSE )
#   if(fixed.width)
#     par(family="mono")
#
#   plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "", asp = NA)
#
#   slist   <- unlist(lapply(object, function(x) strsplit(x,'\n')))
#   slist   <- lapply(slist, function(x) unlist(strsplit(x,'')))
#
#   slen    <- sapply(slist, length)
#   slines  <- length(slist)
#
#   if (missing(cex))
#   {
#     lastloop <- FALSE
#     cex <- 1
#   }
#   else
#     lastloop <- TRUE
#
#
#   for (i in 1:20)
#   {
#     oldcex <- cex
#     #cat("cex=",cex,"\n")
#     #cat("i=",i,"\n")
#     #cat("calculating width...")
#     cwidth  <- max(sapply(unlist(slist), strwidth,  cex=cex)) * cspace
#     #cat("done.\n")
#     #cat("calculating height...")
#     cheight <- max(sapply(unlist(slist), strheight, cex=cex)) * ( lspace + 0.5 )
#     #cat("done.\n")
#
#     width <- strwidth(object, cex=cex)
#     height <- strheight(object, cex=cex)
#
#     if(lastloop) break
#
#     cex <- cex  / max(width, height)
#
#     if (abs(oldcex - cex) < 0.001)
#     {
#       lastloop <- TRUE
#     }
#
#   }
#
#   if (halign == "left")
#     xpos <- 0
#   else if (halign == "center")
#     xpos <- 0 + (1 - width)/2
#   else xpos <- 0 + (1 - width)
#
#   if (valign == "top")
#     ypos <- 1
#   else if (valign == "center")
#     ypos <- 1 - (1 - height)/2
#   else ypos <- 1 - (1 - height)
#
#   text(x=xpos, y=ypos, labels=object, adj=c(0,1),
#        cex=cex, ...)
#
#   par(opar)
#   invisible(cex)
# }
