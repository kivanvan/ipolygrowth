#' @title Plot polynomial growth curve
#'
#' @description
#' This function plot the observed bacterial growth data with a line of fitted polynomial growth curve.
#'
#' @param gcdata A data set from gc.data()
#' @param model.curve A polynomial regression model
#' @param color Color of fitted polynomial growth curve. The default color is blue.
#' @param title Title of the plot. If not specified, there will be no title.
#'
#' @import dplyr
#' @return A plot
#' @export
#'
#' @examples
#' gc.plot(df.test, lm1)
#' gc.plot(df.test, lm1, color = "red", title = "My curve")
#'
#'
#'
gc.plot <- function(gcdata, model.curve, color = "blue", title = NULL) {
  gcdata <- gcdata %>%
    dplyr::mutate(., fit = stats::predict(model.curve))

  gcdata2 <- gcdata %>%
    dplyr::select(time, fit) %>%
    dplyr::distinct()

  p <-
    ggplot2::ggplot()+
    ggplot2::geom_line(data= gcdata2, ggplot2::aes(x=time, y=fit), color=color)+
    ggplot2::geom_point(data= gcdata, ggplot2::aes(x=time, y=logOD))+
    ggplot2::labs(x = "time", y= "log OD", title = title)+
    ggplot2::theme_classic()
  print(p)
}

