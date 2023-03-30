#' @title Growth parameter estimates for multiple samples
#'
#' @description
#' Add descriptions
#'
#' @param df.input Data contains time, measurements, and sample identifier from multiple samples. Data needs to be in long format.
#' @param factors Sample identifier
#' @param time Name of the time variable
#' @param y Name of the measurement variable
#' @param max.od.select A vector of threshold for maximum detection in exponential phase. See explanation below.
#'
#' @details
#' Explain the max.od.select
#'
#' @return A list that contains a table of estimates, polynomial models, a table of beta coefficients, and a table of fitted values.
#' @import dplyr stats
#'
#' @examples
#' Add examples
#'
ipg_multisample <- function(df.input, factors, time, y, max.od.select) {
  ls.df <- split(df.input, df.input[, factors], sep = "^")                                          # split datasets; it's separated by .

  df.factor <- data.frame(order = names(ls.df)) %>%
    tidyr::separate_wider_delim(., "order", "^", names = factors)                                   # unique combination of stratified factors

  duptimes <- df.input[, c(factors, time)] %>%
   # duptimes <- df.input %>% dplyr::select(tidyselect::all_of(vf), time) %>%
    dplyr::distinct() %>% dplyr::group_by(dplyr::across(tidyselect::all_of(factors))) %>%
    # dplyr::distinct() %>% dplyr::group_by(.data[[factors]]) %>%
    dplyr::summarise(n = n())
  duptimes <- duptimes$n                                                                            # get duplicated times for each data subset [edit]
  df.factor2 <- df.factor %>% dplyr::slice(rep(1:nrow(.), duptimes))                                # duplicate factors for fitted value

  ls.model <- ls.fit <- vector(mode = 'list', length = nrow(df.factor))                             # create lists to store model and fitted values
  names(ls.model) <- do.call(paste, df.factor[names(df.factor)])                                    # match component name with stratified factors

  output.estimates <- data.frame(
#    df.factor,
    peak.growth.rate = NA_real_, peak.growth.time = NA_real_, start.peak.growth = NA_real_,
    fit.max.od = NA_real_, fit.max.od.time = NA_real_,doubling.time = NA_real_)                     # create a dataset to store param estimates together

  output.betas <- data.frame(
#    df.factor,
    B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)                      # create a dataset to store beta cof together

  output.fitted <- data.frame(
    time = NULL, fit = NULL
  )


  for (i in 1:length(ls.df)) {
    temp <- ipg_singlesample(ls.df[[i]], time, y, max.od.select[i]) # requires a vector of threshold
    ls.model[[i]] <- temp$model
    # print(ncol(output.estimates))
    # print(temp$estimates)
    output.estimates[i, ] <- temp$estimates
    # print(ncol(output.estimates))
    # print(temp$estimates)
    # print(output.estimates)
    output.betas[i, ] <- temp$betas
    output.fitted <- rbind(output.fitted, temp$fitted)
  }

  colnames(output.estimates) <- c("peak growth rate", "peak growth time", "start of peak growth", "max od (fit)", "max od time (fit)", "doubling time")
  output.estimates <- cbind(df.factor, output.estimates)
  output.betas <- cbind(df.factor, output.betas)
  output.fitted <- cbind(df.factor2, output.fitted)

  ls.output <- list(
    model = ls.model,
    estimates = output.estimates,
    betas = output.betas,
    fitted = output.fitted
  )
  return(ls.output)
}
