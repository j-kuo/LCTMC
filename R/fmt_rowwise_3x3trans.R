#' @title Formats data for model fitting
#'
#' @description Takes a data frame object in long format and convert into a row-wise transition format.
#' Input data should meet certain requirements indicated in the **Note** section
#'
#' @param data a data frame object. Data should be stored in long format. See example for data structure
#'
#' @return a data.frame object which contains the data in row-wise transitions format.
#'
#' @note Four conditions of the input data frame are checked, if any fails, then error will be returned
#' \enumerate{
#'   \item there must be a column named "id" in the input data
#'   \item there must be columns named "state_at_obsTime" and "obsTime" in the input data
#'   \item input data cannot have columns named "x0" and "w0" as these are reserved for the function to create a intercept term variable
#'   \item each person indicated by the "id" column should have at least 2 observation minimum
#' }
#'
#' @export
#'
#' @example inst/examples/ex_fmt_rowwise_3x3trans.R

fmt_rowwise_3x3trans = function(data = data.frame()) {
  # check ... column names
  if (!"id" %in% colnames(data)) {
    stop("`id` should be a column name in `data` serving as the individual level identifier")
  }
  if (!all(c("state_at_obsTime", "obsTime") %in% colnames(data))) {
    stop("`state_at_obsTime` and `obsTime` should be column names in `data`")
  }
  if (any(c('x0', 'w0') %in% colnames(data))) {
    stop("`x0` and `w0` cannot be column names in `data`. These are reserved for the intercept terms")
  }

  # check ... all people must have a least 1 transition
  temp = do.call(`c`, Map(`nrow`, split(x = data, f = data$id)))
  temp_obs_count = unname(temp)
  temp_id = names(temp)

  # messages
  if (any(temp_obs_count == 1)) {
    x = as.character(temp_id[temp_obs_count == 1])

    if (length(x) > 1) {
      err_msg = paste("ID=(", x[1], ", ", x[2], ", and possibly others) only have 1 observation in the input data frame.", sep = "")
    } else {
      err_msg = paste("ID=(", x, ") only have 1 observation in the input data frame.", sep = "")
    }
    stop(paste(err_msg, "\n  ",  "All `data$id` should have at least two observations", sep = ""))
  }

  # `arrange` by ID
  data = data[order(data$id), ]

  # mutate outside of `group_by`
  data$x0 = 1
  data$w0 = 1
  data$y1 = data$state_at_obsTime
  data$t1 = data$obsTime

  # `mutate` within `group_by`
  lead.group_by = function(x) c(x[2:length(x)], NA)
  data$y2 = do.call(`c`, Map(`lead.group_by`, split(x = data$state_at_obsTime, f = data$id)))
  data$t2 = do.call(`c`, Map(`lead.group_by`, split(x = data$obsTime, f = data$id)))

  # `mutate` outside of `group_by`
  data$dt = data$t2 - data$t1

  # `filter`
  data = data[!is.na(data$y2), ]

  # `mutate` outside of `group_by`
  data$trans.1_1 = (data$y1 == 1) * (data$y2 == 1) # 1-1
  data$trans.1_2 = (data$y1 == 1) * (data$y2 == 2) # 1-2
  data$trans.1_3 = (data$y1 == 1) * (data$y2 == 3) # 1-3
  data$trans.2_1 = (data$y1 == 2) * (data$y2 == 1) # 2-1
  data$trans.2_2 = (data$y1 == 2) * (data$y2 == 2) # 2-2
  data$trans.2_3 = (data$y1 == 2) * (data$y2 == 3) # 2-3

  # `select` minus
  data = data[, !(colnames(data) %in% c("state_at_obsTime", "y1", "y2", "t1", "t2"))]

  # `select`
  everything_else = colnames(data)[!colnames(data) %in% c("id", "obsTime")]
  data = data[, c("id", "obsTime", everything_else)]

  # return
  rownames(data) = NULL
  data
}
