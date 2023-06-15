#' @title Growth parameter estimates for multiple samples
#'
#' @description
#' This function estimates growth parameters for a dataset with multiple samples. Technical replicates are allowed in the dataset.
#'
#' @param data Input data containing the time and dependent variable from a single sample. Data needs to be in long format.
#' Remove rows with missing in the dependent variable.
#' @param id Sample ID to indicate the biological replicates. When exists more than one ID, put a vector of variable names of the IDs
#' @param time.name Name of the time variable
#' @param y.name Name of the dependent variable
#' @param epsilon A threshold for maximum detection of dependent variable in exponential phase. Default is 0.2%.
#' This number will be multiplied by the range of the dependent variable to establish a threshold for maximum detection.
#'
#' @details
#' The function uses the same approach to estimate growth parameters as in `ipg_singlesample()`.
#'
#' @return A list that contains a table of estimates, polynomial models, a table of beta coefficients, and a table of fitted values, all by sample identifier(s).
#' Growth parameters include peak growth rate, peak growth time, lag time, max y (fit), max y time (fit), and doubling time at the peak growth.
#'
#' @importFrom rlang .data
#'
#' @examples
#' library(dplyr)
#' data <- growthrates::bactgrowth
#' out.multisample <- ipg_multisample(data, c("strain", "conc"), "time", "value", 0.2/100)
#' @export
#'
#'
ipg_multisample <- function(data, id, time.name, y.name, epsilon) {
  n <- . <- NULL
  X <- 0

  ls.df <- split(data, data[, id], sep = "^")                                                       # split datasets; name is separated by ^

  df.factor <- data.frame(order = names(ls.df)) %>%
    tidyr::separate_wider_delim("order", delim = "^", names = id)                                   # unique combination of stratified id by row;

  ls.model <- vector(mode = 'list', length = nrow(df.factor))                                       # create lists to store model and fitted values
  names(ls.model) <- do.call(paste, df.factor[names(df.factor)])                                    # match component name with stratified id

  output.estimates <- data.frame(
    peak.growth.rate = NA_real_, peak.growth.time = NA_real_,doubling.time = NA_real_,
    start.peak.growth = NA_real_, fit.max.y = NA_real_, fit.max.y.time = NA_real_)                  # create a dataset to store param estimates together

  output.betas <- data.frame(
    B0 = NA_real_, B1 = NA_real_, B2 = NA_real_, B3 = NA_real_, B4 = NA_real_)                      # create a dataset to store beta coefficient

  output.fitted <- data.frame(time = NULL, fit = NULL)                                              # create a dataset to store fitted values

  duptimes <- data[, c(id, time.name)] %>%                                                          # There will be issue when missing in y
    dplyr::distinct() %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(id))) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop_last")
    duptimes <- duptimes$n                                                                          # get duplicated times for each data subset
    df.factor2 <- df.factor %>% dplyr::slice(rep(1:nrow(.), duptimes))                              # duplicate id for fitted value

  for (i in 1:length(ls.df)) {
    tryCatch(
      ipg_singlesample(ls.df[[i]], time.name, y.name, epsilon[i]),
      message = function(cnd) {
        X <<- X+1
      },
      finally = {
        temp <- suppressMessages(ipg_singlesample(ls.df[[i]], time.name, y.name, epsilon[i]))  # requires a vector of threshold
        ls.model[[i]] <- temp$model
        output.estimates[i, ] <- temp$estimates
        output.betas[i, ] <- temp$betas
        output.fitted <- rbind(output.fitted, temp$fitted)
      }
    )
  }

  if (X>0) {
    message(paste0("max y time is equal to the largest value of \"", time.name, "\" in ", X, " samples."))
  }


  ## clean output ##----
  colnames(output.estimates) <- c("peak growth rate", "peak growth time", "doubling time", "lag time", "max y", "max y time")
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
