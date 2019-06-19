#'
#'
#' Bead count with ggplot
#'
#' Investigate if you have enough beads in each well of each beadID
#'
#'
#'
#'
#' @section Warnings:
#' This is a warning
#' @param bead_count_df A dataframe with bead count. Each row is a well, each column is a bead.
#' @param bead_count_limit How many beads do you need? Can be vector if you want several lines.
#' @param path_for_output Path for output
#' @param date_in_file_name Yes or No
#' @return Plot with bead count per sample and per bead
#' @examples
#' bead_count_base(df)
#' @export


bead_count_base <- function(bead_count_df, bead_count_limit = 32, path_for_output = "./", date_in_file_name = T, ...) {


  file_name_and_path <- ifelse(date_in_file_name == T, paste0(path_for_output, paste0(Sys.Date() ," Bead count with base.pdf")),
                               paste0(path_for_output, paste0(" Bead count with base.pdf")))

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

  beeswarm::beeswarm(bead_count_df, add = TRUE, cex = 0.5, pch = 16, col = "grey62", corral = "gutter",  ...)

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

  beeswarm::beeswarm(data.frame(t(bead_count_df)), add = TRUE, cex = 0.5, pch = 16, col = "grey62", corral = "gutter",  ...)

  # To get the boxes above the beeswarm
  boxplot(t(bead_count_df), main = "", cex.main = 2.5, names = F, xaxt = "n",
          pch = "", cex = 0.5,
          log = "y",las=2, add = T
  )

  abline(h= bead_count_limit, col= "skyblue2")
  mtext("Number of beads", side = 2, line = 2.8, las = 3, cex = 1.2, font = 2)
  mtext(paste0("Samples (n = ", length(rownames(bead_count_df)), ")"), side = 1, line = 2, las = 1, cex = 1.2, font = 2)

dev.off()
}







#'
#'
#' Bead count with ggplot
#'
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
#' @param path_for_output Path for output
#' @param date_in_file_name Yes or No
#' @return Plot with bead count per sample and per bead
#' @examples
#' bead_count_ggplot(df)
#' @export



bead_count_ggplot <- function(bead_count_df, bead_count_limit = 32, path_for_output = "./", date_in_file_name = T, ...) {
  new_theme <- theme(
    axis.title=element_text(color="grey62",face="bold", size = 20),
    plot.title=element_text(color="black",face="bold", size = 30),
    axis.text.x = element_text(color="white"),
    axis.text.y = element_text(color="black", size = 14),
    axis.ticks.length = unit(0, "cm"),
    plot.subtitle=element_text(color="Pink"),
    panel.grid=element_blank()
  )

  file_name_and_path <- ifelse(date_in_file_name == T, paste0(path_for_output, paste0(Sys.Date() ," Bead count with ggplot.pdf")),
                               paste0(path_for_output, paste0(" Bead count with ggplot.pdf")))

  df_bead_count_samples_in_col <- cbind(bead_count_df, rownames(bead_count_df)); colnames(df_bead_count_samples_in_col) <- c(colnames(bead_count_df), "Samples")

  df_bead_count_long <- gather(df_bead_count_samples_in_col, key = "BeadID", value = "Number_of_beads", colnames(df_bead_count_samples_in_col), na.rm = FALSE, convert = FALSE, -"Samples")
  df_bead_count_long$BeadID <- as.factor(df_bead_count_long$BeadID)

  pdf(file = file_name_and_path, width= 16, height= 8)

  first_plot <- ggplot(data = df_bead_count_long) +
    labs(title="Bead count")+
    geom_jitter(mapping = aes(x=BeadID, y=Number_of_beads), width = 0.2, cex = 0.5, col = "grey62")+
    geom_boxplot(mapping = aes(x=BeadID, y=Number_of_beads))+
    geom_hline(yintercept = bead_count_limit, show.legend = NA, col = "skyblue2")+
    theme_bw()+
    new_theme

  print(first_plot)


  second_plot <- ggplot(data = df_bead_count_long) +
    labs(title="Bead count")+
    geom_jitter(mapping = aes(x=Samples, y=Number_of_beads), width = 0.2, cex = 0.5, col = "grey62")+
    geom_boxplot(mapping = aes(x=Samples, y=Number_of_beads))+
    geom_hline(yintercept = bead_count_limit, show.legend = NA, col = "skyblue2")+
    theme_bw()+
    new_theme

  print(second_plot)

  dev.off()
}

