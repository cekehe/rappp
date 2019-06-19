#' beeswarm
#'
#' Beeswarm plots for SBA data.
#'
#' @details
#'
#' @param x
#' @return
#' @examples
#' @export


beeswarm <- function(data, filename, sample_info, bead_info, grouping_var, sample_name_var, xlab) {

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

  map(1:ceiling(ncol(data)/12), beeswarm_plotting)

  dev.off()
}
