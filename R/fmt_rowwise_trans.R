#' @title Formats data for model fitting
#'
#' @description Takes a data frame object in long format and convert into a row-wise transition format.
#'
#' @param data a data frame object. Data should be stored in long format. See example for data structure
#' @param type a character scalar. Equals either "2x2" or "3x3" to indicate the number of category in the outcome variable.
#' @param X_names a named character vector. Hosting the names of covariates for the CTMC model. It should be a column in `data`. \cr
#' Best set to x0, x1, x2 to avoid errors
#' @param W_names a named character vector. Hosting names of covariates for the latent class component. It should be a column in `data`. \cr
#' Best set to w0, w1, w2 to avoid errors
#' @param scaling a named numeric vector indicating which covariates are scaled by how much. Set to 1 for no scaling. \cr
#' For example, `scaling = c(x0 = 1, x1 = 0.01, x2 = 1, w0 = 1, w1 = 1, w2 = 1, dt = 0.002)`
#' @param trace a logical scalar, if TRUE, function will print which parameters were scaled.
#' Default is FALSE which does not print messages
#'
#' @return a list object containing the 4 major data objects needed to perform model fitting. The objects are the following:
#' \describe{
#'   \item{Xmat}{a numeric matrix of the "X" covariates which are variables that affect the CTMC process}
#'   \item{Wmat}{a numeric matrix of the "W" covariates which are variables that affect the latent class probabilities}
#'   \item{dt}{a numeric vector of the time interval between observations}
#'   \item{df_trans}{a data.frame object containing the binary transition indicator variables}
#' }
#'
#' @note This is the first step out of six of fitting a latent class CTMC model (i.e., data preparation).
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]; [create_controls()]
#'
#' @example inst/examples/ex_fmt_rowwise_trans.R

fmt_rowwise_trans = function(data = data.frame(),
                             X_names = c("x0"),
                             W_names = c("x0"),
                             type = c("2x2", "3x3"),
                             scaling = c("x0" = 1),
                             trace = FALSE) {
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
