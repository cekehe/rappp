#' --------------------------------------------------------------------------------------------------
#' Bead count with base
#' --------------------------------------------------------------------------------------------------
#' Investigate if you have enough beads in each well of each beadID
#'
#'
#'
#'
#' @section Warnings:
#' This is a warning
#' @param bead_count_df A dataframe with bead count. Each row is a well, each column is a bead.
#' @param bead_count_limit How many beads do you need? Can be vector if you want several lines.
#' @param data_in_file_name Yes or No.
#' @return Plot with bead count per sample and per bead
#' @examples
#' bead_count(df)
#' @export


bead_count_base <- function(bead_count_df, bead_count_limit = 32, path_for_output = "./", date_in_file_name = T, ...) {

  file_name_and_path <- paste0(path_for_output, paste0(Sys.Date() ," Bead count with base.pdf"))

  #par(oma= c(0.5, 0.5, 0.5, 0.5), mar=c(8.1, 8.1, 4.1, 4.1))
  pdf(file = file_name_and_path, width= 16, height= 8)


  # --------------------------------------------------------------------------------------------------
  # Per antibody
  # --------------------------------------------------------------------------------------------------
  boxplot(bead_count_df, main = "Bead count", cex.main = 2, names = F, xaxt = "n",
          pch = "", cex = 0.5,
          #log = "y",
          las=2
  )

  beeswarm(bead_count_df, add = TRUE, cex = 0.5, pch = 16, col = "grey62", corral = "gutter",  ...)

   # To get the boxes above the beeswarm
  boxplot(bead_count_df, main = "", cex.main = 2.5, names = F, xaxt = "n",
          pch = "", cex = 0.5,
          las=2, add = T
  )

  abline(h= bead_count_limit, col= "skyblue2")
  mtext("Number of beads", side = 2, line = 2.8, las = 3, cex = 1.2, font = 2)
  mtext(paste0("BeadID (n = ", length(colnames(bead_count_df)), ")"), side = 1, line = 2, las = 1, cex = 1.2, font = 2)

  # --------------------------------------------------------------------------------------------------
  # Per sample
  # --------------------------------------------------------------------------------------------------
  boxplot(t(bead_count_df), main = "Bead count", cex.main = 2, names = F, xaxt = "n",
          pch = "", cex = 0.5,
          #log = "y",
          las=2
  )

  # beeswarm(t(bead_count_df), add = TRUE, cex = 0.5, pch = 16, col = "grey62", corral = "gutter",  ...)
  #
  # # To get the boxes above the beeswarm
  # boxplot(t(bead_count_df), main = "", cex.main = 2.5, names = F, xaxt = "n",
  #         pch = "", cex = 0.5,
  #         log = "y",las=2, add = T
  # )

  abline(h= bead_count_limit, col= "skyblue2")
  mtext("Number of beads", side = 2, line = 2.8, las = 3, cex = 1.2, font = 2)
  mtext(paste0("Samples (n = ", length(rownames(bead_count_df)), ")"), side = 1, line = 2, las = 1, cex = 1.2, font = 2)

dev.off()
}




