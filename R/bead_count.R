#' Bead count
#'
#' Bead count
#'
#' Investigate if you have enough beads in each well of each beadID
#'
#'
#' @section Warnings:
#' This is a warning
#' @param bead_count_file A dataframe with bead count. Each row is one well, each column is a bead.
#' @return Plot with bead count per sample and per bead
#' @examples
#' bead_count(df)
#' @export

bead_count <- function(x) {
  for(i in 1:x){
    cat(rep( letters[i], times = i),"\n")
  }
  # Check if second object can be called by our function
  cat(second_object)
}


