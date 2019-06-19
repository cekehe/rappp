#' Beeswarm plots for SBA data using ggplot
#'
#' Makes beeswarm (jitter) plots for all analytes (beads) in the input object \code{data}.
#' Grouping is performed according to the \code{grouping_var} column in the df \code{sample_info}.
#'
#' @param data A df with data from an SBA experiment. May be transformed or not.
#' @param filename A string with the desired name of the output \code{.pdf}.
#' @param sample_info A df with the sample info.
#' @param bead_info A df with the bead info.
#' @param grouping_var A string with the name of the column in \code{sample_info} containing the desired grouping variable.
#' @param sample_name_var A string with the name of the column in \code{sample_info} containing the desired sample identifier.
#' @param xlab A character vector containing the desired grouping labels to be printed.
#' @export


beeswarm_ggplot <- function(data, filename, sample_info, bead_info, grouping_var, sample_name_var, xlab) {

  data_long <- data %>%
    cbind(sample_info[[sample_name_var]], sample_info[[grouping_var]]) %>%
    `colnames<-`(c(colnames(.)[1:(ncol(.)-2)], sample_name_var, grouping_var)) %>%
    gather(key = "Bead_name", value = "MAD", -sample_name_var, -grouping_var)

  beeswarm_plotting <- function(i) {
    ggplot(data_long) +
      facet_wrap_paginate(. ~ Bead_name, ncol = 3, nrow = 4, page = i) +
      geom_jitter(aes(x = DiagICD10c6, y = MAD), width = 0.3, size = 0.7) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_x_discrete(labels = xlab)
  }

  pdf(paste0(filename ,".pdf"), paper = "a4", height = 12, width = 12, useDingbats = FALSE)

  print({
    p <- map(1:ceiling(ncol(data)/12), beeswarm_plotting)
    p
  })

  dev.off()
}
