#' @title Fit polynomial growth curve
#'
#' @param gcdata A data set from gc.data()
#'
#' @details
#' This function fits polynomial growth curve and outputs a model object that can be used for
#' parameter estimation and plot.
#'
#' @return A model object
#' @export
#'
#' @examples
#' lm1 <- gc.poly(df.test)
#'
#'
#'
gc.poly <- function(gcdata) {
  model.curve <- stats::lm(logOD ~ time + time2 + time3 + time4, data = gcdata)                     # 4th order polynomial regression (assumes independent obs)
  return(model.curve)                                                                               # return result
}
