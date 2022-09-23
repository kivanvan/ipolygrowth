#' Calculate growth curve parameters
#'
#' @param model.curve A polynomial regression model
#' @param df.input A data set with observed OD values
#' @param min.time Starting time of bacterial growth
#' @param max.time Ending time of bacterial growth
#'
#' @return A data set with growth curve parameters
#' @export
#'
#' @examples
#'df.test.result <- gc.parameters(lm1, df.test, 20/60, 1100/60)
#'
gc.parameter <- function(model.curve, df.input, min.time, max.time) {
  df.input <- df.input %>%
    dplyr::mutate(fit = stats::predict(model.curve))

  df.result <- data.frame(
    #id = NA_integer_, visit = NA_character_,
    peak.growth.rate = NA_real_, peak.growth.time = NA_real_, start.peak.growth = NA_real_,
    max.od = NA_real_, max.od.time = NA_real_, fit.max.od = NA_real_, fit.max.od.time = NA_real_,
    doubling.time = NA_real_,
    observedOD.begining = NA_real_, observedOD.end = NA_real_,
    predictedOD.begining = NA_real_, predictedOD.end = NA_real_,
    B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)

  df.result$B0 <- A0 <- as.numeric(coef(model.curve)[1])                                            # extract betas
  df.result$B1 <- A1 <- as.numeric(coef(model.curve)[2])
  df.result$B2 <- A2 <- as.numeric(coef(model.curve)[3])
  df.result$B3 <- A3 <- as.numeric(coef(model.curve)[4])
  df.result$B4 <- A4 <- as.numeric(coef(model.curve)[5])


  ## calculate peak growth rate ##
  # first derivative of 4th degree polynomial is A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3)
  # second derivative of 4th degree polynomial is A2 * 2 + A3 * (3 * (2 * x)) + A4 * (4 * (3 * x^2))
  x.roots <- polyroot(c(A2 * 2, A3 * (3 * (2)), A4 * (4 * (3))))                                    # find roots of second derivative; peak growth times
  x <- c(x.roots,0)                                                                                 # add 0 to calculate the slope at 0

  pgr <- Re(exp(eval(A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3))))                         # calculate growth rate (slope) at peak growth times (real)
  pgr.index <- which(pgr == max(pgr[pgr>0 & Re(x) >= min.time & Re(x) <= max.time]))[1]             # choose growth time with maximum growth rate and peak growth time between 0 and 18 hrs-[[change 0 and 18 to min and max time]]

  df.result$peak.growth.time <- pgt <- Re(x[pgr.index])                                             # peak growth time
  df.result$peak.growth.rate <- pgr[pgr.index]                                                      # peak growth rate at peak growth time
  df.result$doubling.time <- Re(log(2)/eval(A1 + A2*(2*x[pgr.index]) + A3*(3*x[pgr.index]^2) + A4*(4*x[pgr.index]^3)))## doubling time = ln(2)/growth rate
  df.result$observedOD.begining <- exp(df.input$logOD[1])
  df.result$observedOD.end <- exp(df.input$logOD[length(df.input$logOD)])
  df.result$predictedOD.begining <- exp(df.input$fit[1])
  df.result$predictedOD.end <- exp(df.input$fit[length(df.input$fit)])


  ## calculate start of peak growth ##
  # calculate the time when the line of intercept from the polynomial model intersects with the linear line of peak growth at peak growth time
  y0 <- log(df.result$predictedOD.begining)                                                         # line at predicted logOD at 0 hrs
  y1 <- eval(expression(A0 + A1 * pgt + A2 * pgt^2 + A3 * pgt^3 + A4 * pgt^4))                      # calculate predicted log OD at peak growth time
  df.result$start.peak.growth <- solve(pgr[pgr.index], y0-y1+pgr[pgr.index]*pgt)


  ## assign max OD to result data set ##
  df.result$max.od <- max(df.input$logOD, na.rm = TRUE)
  df.result$max.od.time <- df.input$time[which(df.input$logOD == df.result$max.od)][1]              # when there are multiple times, select the first time

  ## calculate the max OD ##
  if (any((df.input$diffOD<0.0035) & (df.input$diffOD>0), na.rm=TRUE)) {
    if (any(df.input$time[(df.input$diffOD<0.0035) & (df.input$diffOD>0)] > df.result$peak.growth.time, na.rm=TRUE)) {
      # find times that satisfy the requirement: difference in OD between two time points is less than 0.0035 and larger than 0
      temp.time <- df.input$time[(df.input$diffOD<0.0035) & (df.input$diffOD>0)]
      # select the first time that is later than the peak growth time
      temp <- temp.time[temp.time > df.result$peak.growth.time][1]
      # return the corresponding fitted OD value
      df.result$fit.max.od <- df.input$fit[df.input$time == temp]
    } else {df.result$fit.max.od <- max(df.input$fit)}
  } else {df.result$fit.max.od <- max(df.input$fit)}

  df.result$fit.max.od.time <- df.input$time[df.input$fit == df.result$fit.max.od][1]

  return(df.result)
}
