#' @title Formats data for model fitting
#'
#' @description Takes a data frame object in long format and convert into a row-wise transition format.
#' Input data should meet certain requirements indicated in the **Note** section
#'
#' @param data a data frame object. Data should be stored in long format. See example for data structure
#' @param type a character scalar. Equals either "2x2" or "3x3" to indicate the number of category in the outcome variable.
#' @param X_names a named character vector. Hosting the names of covariates for the CTMC model. It should be a column in `data`. \cr
#' Best set to x0, x1, x2 to avoid errors
#' @param W_names a named character vector. Hosting names of covariates for the latent class component. It should be a column in `data`. \cr
#' Best set to w0, w1, w2 to avoid errors
#' @param scaling a named numeric vector indicating which covariate are scaled by how much. Set to 1 for no scaling. \cr
#' For example, `scaling = c(x0 = 1, x1 = 0.01, x2 = 1, w0 = 1, w1 = 1, w2 = 1, dt = 0.002)`
#' @param ... the following are optional parameters\
#' \describe{
#'   \item{trace}{a logical scalar, if TRUE, function will print which parameters were scaled. Default is FALSE which does not print messages}
#' }
#'
#' @return a list object containing the 4 major data objects needed to perform model fitting. The objects are the following:
#' \describe{
#'   \item{Xmat}{a numeric matrix of the "X" covariates which are variables that affect the CTMC process}
#'   \item{Wmat}{a numeric matrix of the "W" covariates which are variables that affect the latent class probabilities}
#'   \item{dt}{a numeric vector of the time interval between observations}
#'   \item{df_trans}{a data.frame object containing the binary transition indicator variables}
#' }
#'
#' @note This is the first step out of six of fitting a latent class CTMC model (i.e., data preparation). \cr\cr
#' Four conditions of the input data frame, `data`, are checked. If any fails, then error will be returned
#' \enumerate{
#'   \item there must be a column named "id" in the input data
#'   \item there must be columns named "state_at_obsTime" and "obsTime" in the input data
#'   \item input data cannot have columns named "x0" and "w0" as these are reserved for the function to create a intercept term variable
#'   \item each person indicated by the "id" column should have at least 2 observation minimum
#' }
#'
#' @seealso [lctmc_2x2()], [lctmc_3x3()]
#'
#' @example inst/examples/ex_fmt_rowwise_trans.R

fmt_rowwise_trans = function(data = data.frame(),
                             type = c("2x2", "3x3"),
                             X_names = c("x"),
                             W_names = c("w"),
                             scaling = c(x = 1, w = 1),
                             ...) {
  ### optional args
  opt_args = list(...)
  trace = ifelse(is.null(opt_args$trace), FALSE, opt_args$trace)

  ### check (1)
  if (!"id" %in% colnames(data)) {
    stop("`id` should be a column name in `data` serving as the individual level identifier")
  }
  if (!all(c("state_at_obsTime", "obsTime") %in% colnames(data))) {
    stop("`state_at_obsTime` and `obsTime` should be column names in `data`")
  }
  if (any(c('x0', 'w0') %in% colnames(data))) {
    stop("`x0` and `w0` cannot be column names in `data`. These are reserved for the intercept terms")
  }
  if (length(type) != 1 || !(type %in% c("2x2", "3x3"))) {
    stop("`type` should be a length 1 character vector. Equals either '2x2' or '3x3'")
  }

  ### check (2) ... compute N per id
  temp = do.call(`c`, Map(`nrow`, split(x = data, f = data$id)))
  temp_obs_count = unname(temp)
  temp_id = names(temp)

  ### check (2) ... error if anyone has less than 2 obs
  if (any(temp_obs_count == 1)) {
    x = as.character(temp_id[temp_obs_count == 1])

    if (length(x) > 1) {
      err_msg = paste("ID=(", x[1], ", ", x[2], ", and possibly others) only have 1 observation in the input data frame.", sep = "")
    } else {
      err_msg = paste("ID=(", x, ") only have 1 observation in the input data frame.", sep = "")
    }
    stop(paste(err_msg, "\n  ",  "All `data$id` should have at least two observations", sep = ""))
  }

  ### `arrange` by ID
  data = data[order(data$id), ]

  ### mutate outside of `group_by`
  data$x0 = 1
  data$w0 = 1
  data$y1 = data$state_at_obsTime
  data$t1 = data$obsTime

  ### `mutate` within `group_by`
  lead.group_by = function(x) c(x[2:length(x)], NA)
  data$y2 = do.call(`c`, Map(`lead.group_by`, split(x = data$state_at_obsTime, f = data$id)))
  data$t2 = do.call(`c`, Map(`lead.group_by`, split(x = data$obsTime, f = data$id)))

  ### `mutate` outside of `group_by`
  data$dt = data$t2 - data$t1

  ### `filter`
  data = data[!is.na(data$y2), ]

  ### `mutate` outside of `group_by` ~ common for both 2x2 and 3x3
  data$trans.1_1 = (data$y1 == 1) * (data$y2 == 1) # 1-1
  data$trans.1_2 = (data$y1 == 1) * (data$y2 == 2) # 1-2
  data$trans.2_1 = (data$y1 == 2) * (data$y2 == 1) # 2-1
  data$trans.2_2 = (data$y1 == 2) * (data$y2 == 2) # 2-2

  ### `mutate` outside of `group_by` ~ only 3x3
  if (type == "3x3") {
    data$trans.1_3 = (data$y1 == 1) * (data$y2 == 3) # 1-3
    data$trans.2_3 = (data$y1 == 2) * (data$y2 == 3) # 2-3
  }

  ### `select` minus
  data = data[, !(colnames(data) %in% c("state_at_obsTime", "y1", "y2", "t1", "t2"))]

  ### `select`
  everything_else = colnames(data)[!colnames(data) %in% c("id", "obsTime")]
  data = data[, c("id", "obsTime", everything_else)]

  ### no row names
  rownames(data) = NULL

  ### covariate scaling
  if (any(scaling != 1)) {
    ## loop through each variable that need to be scaled
    for (sc in seq_along(scaling)) {
      if (scaling[sc] != 1) {
        # msg
        if (trace) {
          sc_name = names(scaling)[sc]
          sc_factor = sprintf("%.3f", scaling[sc])
          cat(" * `", sc_name, "` being scaled by a factor of ", sc_factor, "\n", sep = "")
        }

        # scaling
        data[[names(scaling)[sc]]] = data[[names(scaling)[sc]]] * scaling[sc]
      }
    }
  } else {
    if (trace) {
      cat(" * no scaling were applied", "\n", sep = "")
    }
  }

  ### create X mat
  Xmat = as.matrix(data[, X_names])

  ### create W mat
  Wmat = unique(data[, c("id", W_names)])
  Wmat = as.matrix(Wmat[!colnames(Wmat) %in% "id"])

  ### dt
  dt = data[["dt"]]

  ### data frame with transition indicators
  data = data[!colnames(data) %in% c(X_names, W_names, "dt")]

  ### return
  list(Xmat = Xmat, Wmat = Wmat, dt = dt, df_trans = data)
}
