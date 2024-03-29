#' Coefficient of Variation
#'
#' This function computes the Coefficient of Variation (CV) of the values in x.
#'
#' @details The argument na.rm can be included. Default is FALSE in both
#'     \code{\link[base:mean]{mean()}} and \code{\link[stats:sd]{sd()}}.
#'     If set to TRUE, then missing values are removed before computation proceeds.
#'
#' @param x Vector or matrix with values.
#' @param format If the output should be as "percent" (default) or "decimal".
#' @param digits Integer indicating the number of decimal places (see \code{\link[base:round]{round()}})
#' @param ... Further arguments passed do \code{\link[base:mean]{mean()}} and \code{\link[stats:sd]{sd()}}
#' @return The CV of the input values.
#' @export

cv <- function(x, format="percent", digits=2, ...) {
  format <- match.arg(format, c("percent", "decimal"))

  output <- sd(x, ...)/mean(x, ...)

  if(format == "percent") { output <- output*100 }

  return(round(output, digits=digits))
}

#' Print a 96-well plate layout
#'
#' Create a 96-well plate layout from vector format.
#'
#' @param x values to print in cells/wells
#' @param n indicies for where the values should be printed, default is all wells.
#' @export

layout96 <- function(x, n=1:96) {
  wells96 <- matrix(NA, nrow=8, ncol=12,
                    dimnames=list(LETTERS[1:8], 1:12))
  wells96[n] <- x
  View(wells96)
  return(wells96)
}

#' Print a 384-well plate layout
#'
#' Create a 384-well plate layout from vector format.
#'
#' @param x values to print in cells/wells
#' @param n indicies for where the values should be printed, default is all wells.
#' @export

layout384 <- function(x, n=1:384) {
  wells384 <- matrix(NA, nrow=16, ncol=24,
                    dimnames=list(LETTERS[1:16], 1:24))
  wells384[n] <- x
  View(wells384)
  return(wells384)
}

#' Scatterplot Matrices
#'
#' A matrix of scatterplots is produced.
#' Alternative version of \code{\link[graphics:pairs]{pairs()}} (default S3 method) where all axes are on the bottom and left sides.
#' Argument information copied from \code{\link[graphics:pairs]{pairs()}}.
#' NB! Better functionality for layout alterations is in the works.
#'
#' @details Please see \code{\link[graphics:pairs]{pairs()}}
#'
#' @param x the coordinates of points given as numeric columns of a matrix or data frame.
#'    Logical and factor columns are converted to numeric in the same way that
#'    \code{\link[base:data.matrix]{data.matrix()}} does.
#' @param labels the names of the variables.
#' @param panel \code{function(x, y, ...)} which is used to plot the contents of each panel of the display.
#' @param ... arguments to be passed to or from methods.
#'     Also, graphical parameters (\code{\link[graphics:par]{par()}}) can be given as can arguments to
#'     \code{plot} such as \code{main}.
#'     \code{par("oma")} will be set appropriately unless specified.
#' @param horInd,verInd The (numerical) indices of the variables to be plotted on the horizontal and vertical axes respectively.
#' @param lower.panel,upper.panel separate panel functions (or \code{NULL}) to be used below and above the diagonal respectively.
#' @param diag.panel optional \code{function(x, ...)} to be applied on the diagonals.
#' @param text.panel optional \code{function(x, y, labels, cex, font, ...)} to be applied on the diagonals.
#' @param label.pos y position of labels in the text panel.
#' @param line.main if \code{main} is specified, \code{line.main} gives the \code{line argument} to
#'     \code{\link[graphics:mtext]{mtext()}}
#'     which draws the title. You may want to specify \code{oma} when changing \code{line.main}.
#' @param cex.labels,font.labels graphics parameters for the text panel.
#' @param row1attop logical. Should the layout be matrix-like with row 1 at the top, or graph-like with row 1 at the bottom?
#'     The latter (non default) leads to a basically symmetric scatterplot matrix.
#' @param gap distance between subplots, in margin lines.
#' @param log a character string indicating if logarithmic axes are to be used, see
#'     \code{\link[graphics:plot.default]{plot.default()}}
#'     or a numeric vector of indices specifying the indices of those
#'     variables where logarithmic axes should be used for both x and y.
#'     \code{log = "xy"} specifies logarithmic axes for all variables.
#' @export

pairs2 <- function (x, labels, panel = points,
                    # lower.panel.args = NULL,
                    # upper.panel.args = NULL,
                    # diag.panel.args = NULL,
                    ..., horInd = 1:nc, verInd = 1:nc,
                    lower.panel = panel, upper.panel = panel, diag.panel = NULL,
                    text.panel = textPanel, label.pos = 0.5 + has.diag/3, line.main = 3,
                    cex.labels = NULL, font.labels = 1, row1attop = TRUE, gap = 1,
                    log = "")
{
  # if(is.null(args)){
  #   args <- list(lower.panel = NULL,
  #                upper.panel = NULL,
  #                diag.panel = NULL)
  # } else if (!(lower.panel %in% names(args))){
  #   args <- append(args, list(lower.panel = NULL))
  # } else if (!(upper.panel %in% names(args))){
  #   args <- append(args, list(upper.panel = NULL))
  # }  else if (!(diag.panel %in% names(args))){
  #   args <- append(args, list(diag.panel = NULL))
  # }

  # if(is.null(args)){
  #   lower.panel.args = NULL
  #   upper.panel.args = NULL
  #   diag.panel.args = NULL
  # } else {
  #   if(lower.panel %in% names(args)){
  #     lower.panel.args = args$lower.panel
  #   }
  #   if(upper.panel %in% names(args)){
  #     upper.panel.args = args$upper.panel
  #   }
  #   if(diag.panel %in% names(args)){
  #     diag.panel.args = args$diag.panel
  #   }
  # }

  if (doText <- missing(text.panel) || is.function(text.panel))
    textPanel <- function(x = 0.5, y = 0.5, txt, cex, font) text(x,
                                                                 y, txt, cex = cex, font = font)
  localAxis <- function(side, x, y, xpd, bg, col.axis = NULL, main,
                        oma, ...) {
    xpd <- NA
    if (side%%2L == 1L && xl[j])
      xpd <- FALSE
    if (side%%2L == 0L && yl[i])
      xpd <- FALSE
    if (side%%2L == 1L)
      Axis(x, side = side, xpd = xpd, ...)
    else Axis(y, side = side, xpd = xpd, ...)
  }
  localPlot <- function(..., main, oma, font.main, cex.main) plot(...)
  localLowerPanel <- function(..., main, oma, font.main, cex.main) lower.panel(...)
  localUpperPanel <- function(..., main, oma, font.main, cex.main) upper.panel(...)
  localDiagPanel <- function(..., main, oma, font.main, cex.main) diag.panel(...)

  # localLowerPanel <- function(..., main, oma, font.main, cex.main) {
  #   if(is.null(lower.panel.args)) {
  #     lower.panel(...)
  #   } else {
  #     lower.panel(..., lower.panel.args)
  #   }
  # }
  #
  # if(is.null(upper.panel.args)) {
  #   localUpperPanel <- function(..., main, oma, font.main, cex.main) { upper.panel(...) }
  # } else {
  #   localUpperPanel <- function(..., main, oma, font.main, cex.main, upper.panel.args=upper.panel.args) {
  #     do.call(upper.panel, args=c(..., upper.panel.args)) }
  # }
  #
  # localDiagPanel <- function(..., main, oma, font.main, cex.main)  {
  #   if(is.null(diag.panel.args)) {
  #     diag.panel(...)
  #   } else {
  #     diag.panel(..., diag.panel.args)
  #   }
  # }

  dots <- list(...)
  nmdots <- names(dots)
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]]))
        x[[i]] <- as.numeric(x[[i]])
      if (!is.numeric(unclass(x[[i]])))
        stop("non-numeric argument to 'pairs'")
    }
  }
  else if (!is.numeric(x))
    stop("non-numeric argument to 'pairs'")
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel))
    lower.panel <- match.fun(lower.panel)
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel))
    upper.panel <- match.fun(upper.panel)
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel))
    diag.panel <- match.fun(diag.panel)
  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  nc <- ncol(x)
  if (nc < 2L)
    stop("only one column in the argument to 'pairs'")
  if (!all(horInd >= 1L && horInd <= nc))
    stop("invalid argument 'horInd'")
  if (!all(verInd >= 1L && verInd <= nc))
    stop("invalid argument 'verInd'")
  if (doText) {
    if (missing(labels)) {
      labels <- colnames(x)
      if (is.null(labels))
        labels <- paste("var", 1L:nc)
    }
    else if (is.null(labels))
      doText <- FALSE
  }
  oma <- if ("oma" %in% nmdots)
    dots$oma
  main <- if ("main" %in% nmdots)
    dots$main
  if (is.null(oma))
    oma <- c(4, 4, if (!is.null(main)) 6 else 4, 4)
  opar <- par(mfcol = c(length(horInd), length(verInd)), mar = rep.int(gap/2,
                                                                       4), oma = oma)
  on.exit(par(opar))
  dev.hold()
  on.exit(dev.flush(), add = TRUE)
  xl <- yl <- logical(nc)
  if (is.numeric(log))
    xl[log] <- yl[log] <- TRUE
  else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }
  ni <- length(iSet <- if (row1attop) horInd else rev(horInd))
  nj <- length(jSet <- verInd)
  for (j in jSet) for (i in iSet) {
    l <- paste0(if (xl[j])
      "x"
      else "", if (yl[i])
        "y"
      else "")
    localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE,
              type = "n", ..., log = l)
    if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
      box()
      j.odd <- (match(j, jSet) + !row1attop)%%2L
      i.odd <- (match(i, iSet) + !row1attop)%%2L

      ## Added code
      # draw x-axis
      if (i == iSet[ni] & j != jSet[nj])
        localAxis(1L, x[, j], x[, i], ...)
      # draw y-axis
      if (j == jSet[1L] & i != iSet[1L])
        localAxis(2L, x[, j], x[, i], ...)

      ## In original code
      # if (i == iSet[1L] && (!j.odd || !has.upper || !has.lower))
      #   localAxis(3L, x[, j], x[, i], ...)
      # if (i == iSet[ni] && (j.odd || !has.upper || !has.lower))
      #   localAxis(1L, x[, j], x[, i], ...)
      # if (j == jSet[1L] && (!i.odd || !has.upper || !has.lower))
      #   localAxis(2L, x[, j], x[, i], ...)
      # if (j == jSet[nj] && (i.odd || !has.upper || !has.lower))
      #   localAxis(4L, x[, j], x[, i], ...)

      mfg <- par("mfg")
      if (i == j) {
        if (has.diag)
          localDiagPanel(as.vector(x[, i]), ...)
        if (doText) {
          par(usr = c(0, 1, 0, 1))
          if (is.null(cex.labels)) {
            l.wid <- strwidth(labels, "user")
            cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
          }
          xlp <- if (xl[i])
            10^0.5
          else 0.5
          ylp <- if (yl[j])
            10^label.pos
          else label.pos
          text.panel(xlp, ylp, labels[i], cex = cex.labels,
                     font = font.labels)
        }
      }
      else if (i < j)
        localLowerPanel(as.vector(x[, j]), as.vector(x[,
                                                       i]), ...)
      else localUpperPanel(as.vector(x[, j]), as.vector(x[,
                                                          i]), ...)
      if (any(par("mfg") != mfg))
        stop("the 'panel' function made a new plot")
    }
    else par(new = FALSE)
  }
  if (!is.null(main)) {
    font.main <- if ("font.main" %in% nmdots)
      dots$font.main
    else par("font.main")
    cex.main <- if ("cex.main" %in% nmdots)
      dots$cex.main
    else par("cex.main")
    mtext(main, 3, line.main, outer = TRUE, at = 0.5, cex = cex.main,
          font = font.main)
  }
  invisible(NULL)
}


#' Panel function for pairs-functions
#'
#' Writes spearman and pearson correlation values,
#' and the respective p-values, on chosen panel.
#' Values are color coded (see Deatils for key).\cr\cr
#' For linear input! See \code{\link[rappp:panel.corSPpLog]{panel.corSPpLog()}} for log-transformed input.
#'
#' @param x,y values to be correlated, automatically fed to function
#'     \code{\link[graphics:pairs]{pairs()}} or \code{\link[rappp:pairs2]{pairs2()}}
#' @param digits how many significant digits are to be used,
#'     passed \code{\link[base:format]{format()}} with altered default.
#' @param cex.cor text size, automatically calculated if left empty.
#' @details If used with \code{\link[rappp:pairs2]{pairs2()}},
#'     use for upper panel:\cr
#'     \code{pairs2(MATRIX, upper.panel=panel.corSPp)}
#'
#'     Color code:\cr
#'     \code{plot(rep(1, 11),10*seq(0,1,0.1)+1 ,col=rainbow(20)[10:20], pch=16, cex=3,
#'     xaxt="n", yaxt="n", ylab="Absolute correlation value", xlab=NA)\cr
#'     axis(2, at=1:11, seq(0,1,0.1), las=1)}
#' @export

panel.corSPp <- function(x, y, digits = 2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  rs <- cor.test(x, y, method="spearman", use="pairwise.complete.obs")
  rp <- cor.test(x, y, method="pearson", use="pairwise.complete.obs")

  txt_s1 <- format(c(rs$estimate, 0.123456789), digits = digits)[1]
  txt_s2 <- format(c(rs$p.value, 0.123456789), digits = digits, scientific=T)[1]
  txt_s <- bquote(bold(rho~" = "~.(txt_s1)~", p = "~.(txt_s2)))

  txt_p1 <- format(c(rp$estimate, 0.123456789), digits = digits)[1]
  txt_p2 <- format(c(rp$p.value, 0.123456789), digits = digits, scientific=T)[1]
  txt_p <- bquote(bold("r"~" = "~.(txt_p1)~", p = "~.(txt_p2)))
  if(missing(cex.cor)) { cex.cor <- log(1.5)/(dim(x)[2]/3) }
  colrs <- rainbow(20)[10:20]
  rs[which(is.na(rs))] <- 0
  rp[which(is.na(rp))] <- 0
  coltxt_s <- colrs[floor(10*abs(as.numeric(txt_s1))+1)]
  coltxt_p <- colrs[floor(10*abs(as.numeric(txt_p1))+1)]
  text(0.5, 0.6, txt_s, cex = cex.cor, col=coltxt_s, font=2)
  text(0.5, 0.4, txt_p, cex = cex.cor, col=coltxt_p, font=2)
}

#' Panel function for pairs-functions
#'
#' Writes spearman and pearson correlation values,
#' and the respective p-values, on chosen panel.
#' Values are color coded (see Deatils for key).\cr\cr
#' For log-transformed (log="xy" in pairs-function) input!
#' See \code{\link[rappp:panel.corSPp]{panel.corSPp()}} for linear input.
#'
#' @param x,y values to be correlated, automatically fed to function
#'     \code{\link[graphics:pairs]{pairs()}} or \code{\link[rappp:pairs2]{pairs2()}}
#' @param digits how many significant digits are to be used,
#'     passed \code{\link[base:format]{format()}} with altered default.
#' @param cex.cor text size, automatically calculated if left empty.
#' @details If used with \code{\link[rappp:pairs2]{pairs2()}},
#'     use for upper panel:\cr
#'     \code{pairs2(MATRIX, upper.panel=panel.corSPpLog, log="xy")}
#'
#'     Color code:\cr
#'     \code{plot(rep(1, 11),10*seq(0,1,0.1)+1 ,col=rainbow(20)[10:20], pch=16, cex=3,
#'     xaxt="n", yaxt="n", ylab="Absolute correlation value", xlab=NA)\cr
#'     axis(2, at=1:11, seq(0,1,0.1), las=1)}
#' @export

panel.corSPpLog <- function(x, y, digits = 2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  rs <- cor.test(x, y, method="spearman", use="pairwise.complete.obs")
  rp <- cor.test(log(x), log(y), method="pearson", use="pairwise.complete.obs")

  txt_s1 <- format(c(rs$estimate, 0.123456789), digits = digits)[1]
  txt_s2 <- format(c(rs$p.value, 0.123456789), digits = digits, scientific=T)[1]
  txt_s <- bquote(bold(rho~" = "~.(txt_s1)~", p = "~.(txt_s2)))

  txt_p1 <- format(c(rp$estimate, 0.123456789), digits = digits)[1]
  txt_p2 <- format(c(rp$p.value, 0.123456789), digits = digits, scientific=T)[1]
  txt_p <- bquote(bold("r"~" = "~.(txt_p1)~", p = "~.(txt_p2)))
  if(missing(cex.cor)) { cex.cor <- log(1.5)/(dim(x)[2]/3) }
  colrs <- rainbow(20)[10:20]
  rs[which(is.na(rs))] <- 0
  rp[which(is.na(rp))] <- 0
  coltxt_s <- colrs[floor(10*abs(as.numeric(txt_s1))+1)]
  coltxt_p <- colrs[floor(10*abs(as.numeric(txt_p1))+1)]
  text(exp(1)+0.3, (exp(1)+0.5), txt_s, cex = cex.cor, col=coltxt_s, font=2)
  text(exp(1)+0.3, (exp(1)-0.5), txt_p, cex = cex.cor, col=coltxt_p, font=2)
}

