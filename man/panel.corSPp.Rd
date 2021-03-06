% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/general_small_functions.R
\name{panel.corSPp}
\alias{panel.corSPp}
\title{Panel function for pairs-functions}
\usage{
panel.corSPp(x, y, digits = 2, cex.cor)
}
\arguments{
\item{x, y}{values to be correlated, automatically fed to function
\code{\link[graphics:pairs]{pairs()}} or \code{\link[rappp:pairs2]{pairs2()}}}

\item{digits}{how many significant digits are to be used,
passed \code{\link[base:format]{format()}} with altered default.}

\item{cex.cor}{text size, automatically calculated if left empty.}
}
\description{
Writes spearman and pearson correlation values,
and the respective p-values, on chosen panel.
Values are color coded (see Deatils for key).\cr\cr
For linear input! See \code{\link[rappp:panel.corSPpLog]{panel.corSPpLog()}} for log-transformed input.
}
\details{
If used with \code{\link[rappp:pairs2]{pairs2()}},
    use for upper panel:\cr
    \code{pairs2(MATRIX, upper.panel=panel.corSPp)}

    Color code:\cr
    \code{plot(rep(1, 11),10*seq(0,1,0.1)+1 ,col=rainbow(20)[10:20], pch=16, cex=3,
    xaxt="n", yaxt="n", ylab="Absolute correlation value", xlab=NA)\cr
    axis(2, at=1:11, seq(0,1,0.1), las=1)}
}
