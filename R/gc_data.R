#' @title Compile data set for growth curve fitting
#'
#' @description
#' This function convert raw data from microplate readers to a long format time-structured data set
#' that can be used to fit polynomial growth curve.
#'
#' @param data Raw data from microplate readers or data in similar format
#' @param var.OD Variable name of OD measure(s). If multiple columns exist, use c() for input.
#' @param min The starting point of bacterial growth measurement
#' @param max The end-time of bacterial growth measurement
#' @param interval Time interval of data collections. Make sure the units are the same in all three
#' time inputs.
#'
#'
#' @details
#' This function only takes time-structured data, i.e. OD needs to be measured at the same time
#' interval. It will provide a long format data set with 4 time variables: t, t^2, t^3, and t^4.
#' Missing data is allowed.
#'
#' @import dplyr tidyr
#' @return
#' A data set with 4 time variables: t, t^2, t^3, and t^4.
#' @export
#'
#' @examples
#' df.test <- gc.data(df, "C1", min = 20/60, max = 1100/60, interval = 20/60)
#' df.test <- gc.data(df, c("C1", "C2", "C3"), min = 20/60, max = 1100/60, interval = 20/60)
#'
#'
#'
gc.data <- function(data, var.OD, min = NULL, max = NULL, interval = NULL) {
  stopifnot("Input must be a data frame" = is.data.frame(data))
  if(is.null(min)) {min = 0}

  gcdata <- data %>%
    dplyr::mutate(time = seq(min, max, interval),
                  time2 = time^2,
                  time3 = time^3,
                  time4 = time^4,
                  dplyr::across(var.OD, ~ log(.x), .names = "{.col}.log")) %>%
    tidyr::pivot_longer(., var.OD, names_to = "replicate", values_to = "logOD") %>%
    dplyr::select(replicate, time, logOD, time2, time3, time4)

  return(gcdata)
}
