#' @title Growth curve parameters
#'
#' @description
#' This function takes a polynomial regression and calculates a series of bacterial growth curve
#' parameters.
#'
#'
#' @param model A polynomial regression model
#' @param gcdata A data set from gc.data()
#' @param min The starting point of bacterial growth measurement [Do I need this?]
#' @param max The end-time of bacterial growth measurement
#' @param threshold.pct The percentage of the log OD range that determines the maximum difference
#' in log OD for max OD calculation.
#' The default is 0.2%.
#'
#'
#' @details
#' The peak growth time is calculated as the solution of the second derivative of the polynomial
#' model that is between the start and end time point. The peak growth rate is the exponentiation of
#' the slope of polynomial curve at the peak growth time.
#'
#'
#' @import dplyr
#' @return
#' A data set with growth curve parameters:
#' peak growth rate, time to peak growth rate, start time of peak growth (if not at the starting
#' time), log of max OD (observed and predicted), time to max OD (observed and predicted), doubling
#' time, log OD in the beginning (observed and predicted), log OD at the end (observed and predicted),
#' and beta coefficients of time from polynomial regression.
#'
#' If replicates are imported, the observed log of max OD, observed beginning log OD and end log OD
#' will be calculated as the mean across replicates.
#'
#' Note that peak growth rate is back calculated to the original OD unit.
#' @export
#'
#'
#' @examples
#'  df.test.result <- gc.parameter(lm1, df.test, 20/60, 1100/60)
#'
#'
gc.parameter <- function(model, gcdata, min, max, threshold.pct = 0.2/100) {
  gcdata2 <- gcdata %>%
    dplyr::mutate(fit = stats::predict(model)) %>%
    dplyr::select(time, fit) %>%
    dplyr::distinct() %>%                                                                           # create a data frame for predicted log OD and time
    dplyr::mutate(diffOD = c(0, diff(fit, lag = 1)))                                                # difference in fitted log OD

  gcdata3 <- gcdata %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(log.OD = mean(logOD, na.rm = T))                                               # take avg of observed OD across replicates

  gcdata3 <- dplyr::left_join(gcdata3, gcdata2, by = "time")                                        # merge observed data with predicted data

  df.result <- data.frame(
    #id = NA_integer_, visit = NA_character_,
    peak.growth.rate = NA_real_, peak.growth.time = NA_real_, start.peak.growth = NA_real_,
    max.od = NA_real_, max.od.time = NA_real_, fit.max.od = NA_real_, fit.max.od.time = NA_real_,
    doubling.time = NA_real_,
    observedOD.begining = NA_real_, observedOD.end = NA_real_,
    predictedOD.begining = NA_real_, predictedOD.end = NA_real_,
    B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)

  df.result$B0 <- A0 <- as.numeric(stats::coef(model)[1])                                           # extract betas
  df.result$B1 <- A1 <- as.numeric(stats::coef(model)[2])
  df.result$B2 <- A2 <- as.numeric(stats::coef(model)[3])
  df.result$B3 <- A3 <- as.numeric(stats::coef(model)[4])
  df.result$B4 <- A4 <- as.numeric(stats::coef(model)[5])


  ## calculate peak growth rate ##
  # first derivative of 4th degree polynomial is A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3)
  # second derivative of 4th degree polynomial is A2 * 2 + A3 * (3 * (2 * x)) + A4 * (4 * (3 * x^2))
  x.roots <- polyroot(c(A2 * 2, A3 * (3 * (2)), A4 * (4 * (3))))                                    # find roots of second derivative; peak growth times
  x <- c(x.roots,0)                                                                                 # add 0 to calculate the slope at 0

  pgr <- Re(exp(eval(A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3))))                         # calculate growth rate (slope) at peak growth times (real)
  pgr.index <- which(pgr == max(pgr[pgr>0 & Re(x) >= 0 & Re(x) <= max]))[1]                         # choose growth time with maximum growth rate and peak growth time between the start and end point

  df.result$peak.growth.time <- pgt <- Re(x[pgr.index])                                             # peak growth time
  df.result$peak.growth.rate <- pgr[pgr.index]                                                      # peak growth rate at peak growth time
  df.result$doubling.time <- Re(log(2)/eval(A1 + A2*(2*x[pgr.index]) + A3*(3*x[pgr.index]^2) + A4*(4*x[pgr.index]^3)))## doubling time = ln(2)/growth rate
  df.result$observedOD.begining <- exp(gcdata3$log.OD[1])
  df.result$observedOD.end <- exp(gcdata3$log.OD[length(gcdata3$log.OD)])
  df.result$predictedOD.begining <- exp(gcdata3$fit[1])
  df.result$predictedOD.end <- exp(gcdata3$fit[length(gcdata3$fit)])


  ## calculate start of peak growth ##
  # calculate the time when the line of intercept from the polynomial model intersects with the linear line of peak growth at peak growth time
  if (df.result$peak.growth.time == 0) {
    df.result$start.peak.growth <- 0
  } else {
    y0 <- log(df.result$predictedOD.begining)                                                       # line at predicted logOD at 0 hrs
    y1 <- eval(expression(A0 + A1 * pgt + A2 * pgt^2 + A3 * pgt^3 + A4 * pgt^4))                    # calculate predicted log OD at peak growth time
    df.result$start.peak.growth <- solve(log(pgr[pgr.index]), log(pgr[pgr.index])*pgt-(y1-y0))
    }


  ## assign max OD to result data set ##
  df.result$max.od <- max(gcdata3$log.OD, na.rm = TRUE)
  df.result$max.od.time <- gcdata3$time[which(gcdata3$log.OD == df.result$max.od)][1]               # when there are multiple times, select the first time

  ## calculate the max OD ##
  threshold <- (range(gcdata$logOD)[2] - range(gcdata$logOD)[1]) * threshold.pct                    # calculate cutoff based on range

  if (any((gcdata3$diffOD < threshold) & (gcdata3$diffOD > 0), na.rm=TRUE)) {
    if (any(gcdata3$time[(gcdata3$diffOD < threshold) & (gcdata3$diffOD > 0)] > df.result$peak.growth.time, na.rm=TRUE)) {

      temp.time <- gcdata3$time[(gcdata3$diffOD < threshold) & (gcdata3$diffOD > 0)]                # find times that satisfy the requirement: difference in OD between two time points is less than threshold and larger than 0

      temp <- temp.time[temp.time > df.result$peak.growth.time][1]                                  # select the first time that is later than the peak growth time

      df.result$fit.max.od <- gcdata3$fit[gcdata3$time == temp]                                     # return the corresponding fitted OD value

    } else {df.result$fit.max.od <- max(gcdata3$fit)}
  } else {df.result$fit.max.od <- max(gcdata3$fit)}

  df.result$fit.max.od.time <- gcdata3$time[gcdata3$fit == df.result$fit.max.od][1]

  return(df.result)
}

