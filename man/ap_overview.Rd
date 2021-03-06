% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ch_functions.R
\name{ap_overview}
\alias{ap_overview}
\title{Signal overview}
\usage{
ap_overview(x, shouldpdf = TRUE, filename = "Signal_overview.pdf",
  width = 25, height = 15, useDingbats = FALSE, ...)
}
\arguments{
\item{x}{List with at least one elements, see Deatils for naming and content.}

\item{shouldpdf}{Logical, should it plot to pdf?}

\item{filename}{String with filename and desired path, end with .pdf}

\item{width, height}{Width and height for pdf, see \code{\link[grDevices:pdf]{pdf()}}.}

\item{useDingbats}{Logical. Default is \code{FALSE}, compared to in default \code{\link[grDevices:pdf]{pdf()}}.}

\item{...}{Further arguments passed to \code{\link[graphics:boxplot]{boxplot()}}.}
}
\description{
Boxplots of signals per antigen and sample in different orders.
}
\details{
The x list needs to include at least the element
    MFI = assay mfi.

    Data points with the value NA or 0 will be set to 1 for the plotting to allow for
    logarithmic scale without filtering any beads or samples.

Note: The function plots to a layout containing three areas.
}
