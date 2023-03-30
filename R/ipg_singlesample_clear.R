#' @Title Growth parameter estimates for single sample
#'
#' @description
#' Add descriptions
#'
#' @param df.input Data contains time and measurements from a single sample. Data needs to be in long format.
#' @param var.time Name of the time variable
#' @param y Name of the measure variable
#' @param max.od.select A threshold for maximum detection in exponential phase. Default is 0.2%. See explanation below.
#'
#' @details
#' Explain the max.od.select
#'
#' @return A list that contains a table of estimates, the polynomial model, a table of beta coefficients, and a table of fitted values.
#' @import dplyr stats
#'
#' @examples
#' Add examples
#'
ipg_singlesample <- function(df.input, var.time, y, max.od.select = 0.2/100) {
  gcdata <- df.input %>%
    dplyr::rename(time = .data[[var.time]]) %>%
    dplyr::mutate(time_sq = time^2,
                  time_cb = time^3,
                  time_qt = time^4) %>%
    dplyr::arrange(time)

  curve <- stats::lm(as.formula(paste(y, "~ time + time_sq + time_cb + time_qt")), data = gcdata)   # 4th order polynomial regression (assumes independent obs)

  gcdata2 <- gcdata %>%
    dplyr::mutate(fit = stats::predict(curve)) %>%
    dplyr::select(time, fit) %>%
    dplyr::distinct() %>%                                                                           # create a data frame for predicted log OD and time
    dplyr::mutate(diffOD = c(0, diff(fit, lag = 1)))                                                # difference in fitted log OD

  gcdata3 <- gcdata2                                                                                # [if gcdata3 isn't needed, change the downstream code]

  tb.beta <- data.frame(B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)
  tb.beta$B0 <- A0 <- as.numeric(stats::coef(curve)[1])                                             # extract betas
  tb.beta$B1 <- A1 <- as.numeric(stats::coef(curve)[2])
  tb.beta$B2 <- A2 <- as.numeric(stats::coef(curve)[3])
  tb.beta$B3 <- A3 <- as.numeric(stats::coef(curve)[4])
  tb.beta$B4 <- A4 <- as.numeric(stats::coef(curve)[5])

  tb.result <- data.frame(
    peak.growth.rate = NA_real_, peak.growth.time = NA_real_, start.peak.growth = NA_real_,
    fit.max.od = NA_real_, fit.max.od.time = NA_real_,doubling.time = NA_real_)


  ## calculate peak growth rate ##
  # first derivative of 4th degree polynomial is A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3)
  # second derivative of 4th degree polynomial is A2 * 2 + A3 * (3 * (2 * x)) + A4 * (4 * (3 * x^2))
  x.roots <- polyroot(c(A2 * 2, A3 * (3 * (2)), A4 * (4 * (3))))                                    # find roots of second derivative; peak growth times
  x <- c(x.roots,0)                                                                                 # add 0 to calculate the slope at 0

  pgr <- Re(eval(A1 + A2 * (2 * x) + A3 * (3 * x^2) + A4 * (4 * x^3)))                              # calculate growth rate (slope) at peak growth times (real)
  max <- max(gcdata$time, na.rm = T)
  pgr.index <- which(pgr == max(pgr[pgr>0 & Re(x) >= 0 & Re(x) <= max]))[1]                         # choose growth time with maximum growth rate and peak growth time between the start and end point

  tb.result$peak.growth.time <- pgt <- Re(x[pgr.index])                                             # peak growth time
  tb.result$peak.growth.rate <- pgr[pgr.index]                                                      # peak growth rate at peak growth time
  tb.result$doubling.time <- Re(log(2)/eval(A1 + A2*(2*x[pgr.index]) + A3*(3*x[pgr.index]^2) + A4*(4*x[pgr.index]^3)))## doubling time = ln(2)/growth rate


  ## calculate start of peak growth ##
  # calculate the time when the line of intercept from the polynomial model intersects with the linear line of peak growth at peak growth time
  if (tb.result$peak.growth.time == 0) {
    tb.result$start.peak.growth <- 0
  } else {
    y0 <- gcdata2$fit[1]                                                                            # line at predicted logOD at 0 hrs
    y1 <- eval(expression(A0 + A1 * pgt + A2 * pgt^2 + A3 * pgt^3 + A4 * pgt^4))                    # calculate predicted log OD at peak growth time
    tb.result$start.peak.growth <- solve(pgr[pgr.index], pgr[pgr.index]*pgt-(y1-y0))
  }


  ## calculate the max OD ##
  r <- gcdata %>% dplyr::select(any_of(y)) %>% range()
  threshold <- (r[2] - r[1]) * max.od.select                    # calculate cutoff based on range

  if (any((gcdata3$diffOD < threshold) & (gcdata3$diffOD > 0), na.rm=TRUE)) {
    if (any(gcdata3$time[(gcdata3$diffOD < threshold) & (gcdata3$diffOD > 0)] > tb.result$peak.growth.time, na.rm=TRUE)) {

      temp.time <- gcdata3$time[(gcdata3$diffOD < threshold) & (gcdata3$diffOD > 0)]                # find times that satisfy the requirement: difference in OD between two time points is less than threshold and larger than 0

      temp <- temp.time[temp.time > tb.result$peak.growth.time][1]                                  # select the first time that is later than the peak growth time

      tb.result$fit.max.od <- gcdata3$fit[gcdata3$time == temp]                                     # return the corresponding fitted OD value

    } else {tb.result$fit.max.od <- max(gcdata3$fit)}
  } else {tb.result$fit.max.od <- max(gcdata3$fit)}

  tb.result$fit.max.od.time <- gcdata3$time[gcdata3$fit == tb.result$fit.max.od][1]

  colnames(tb.result) <- c("peak growth rate", "peak growth time", "start of peak growth", "max od (fit)", "max od time (fit)", "doubling time")

  output <- list(estimates = tb.result,
                 model = curve,
                 betas = tb.beta,
                 fitted = gcdata2 %>% dplyr::select(time, fit))

  return(output)
}
