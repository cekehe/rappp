#' Coefficient of Variation
#'
#' This function computes the Coefficient of Variation (CV) of the values in x.
#'
#' @details The argument na.rm can be included. Default is FALSE in both \link[stats]{mean} and \link[stats]{sd}.
#' If set to TRUE, then missing values are removed before computation proceeds.
#'
#' @param x Vector or matrix with values.
#' @param format If the output should be as "percent" (default), or "decimal".
#' @param digits Integer indicating the number of decimal places (see \link[stats]{round})
#' @param ... Further arguments passed do \link[stats]{mean} and \link[stats]{sd}
#' @return List of MADs, with same structure as input list.
#' @export

cv <- function(x, format="percent", digits=2, ...) {
  format <- match.arg(format, c("percent", "decimal"))

  output <- sd(x, ...)/mean(x, ...)

  if(format == "percent") { output <- output*100 }

  return(round(output, digits=digits))
}
