#' xMAD beeswarm
#'
#' Beeswarm plots for xMAD data.
#'
#' @details
#'
#' @param x
#' @return
#' @examples
#' @export


beeswarm <- function(data, sample_info, bead_info) {
  pdf("test.pdf", paper = "a4", height = 12, width = 12, useDingbats = FALSE)

  ggplot(data) +
    geom_jitter(aes(x = sample_info$DiagICD10c6, y = B.001_34300_HSPA14), width = 0.3)

  dev.off()
}
