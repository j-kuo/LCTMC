#' @title Checks Model Specification of LCTMC
#'
#' @description To be used only within `lctmc_2x2()` or `lctmc_3x3()`. Takes the exact same arguments as
#' those functions and performs several checks to ensure the model is specified correctly.
#'
#' @param data reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param X_names reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param W_names reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param K reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param par_constraint reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param controls reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param parallel_optim reference [lctmc_2x2()] or [lctmc_3x3()]
#' @param MyModelName reference [lctmc_2x2()] or [lctmc_3x3()]
#'
#' @return NULL
#'
#' @note This function checks the following conditions:
#' \enumerate{
#'   \item `controls` is created from the `LCTMC::create_controls()` function
#'   \item `parallel_optim` is a list object containing the correct element names
#'   \item `MyModelName` is a character scalar, also checks for lengthy names
#'   \item `scaling` parameter in `controls$fmt_data` matches `X_names` & `W_names`
#'   \item `par_constraint` is either NULL or only zero constrain
#'   \item `K` is an integer greater than 1
#'   \item `controls$init02$which_step1` has the correct options
#'   \item `data` contains a column called `id` as the individual identifying numbers
#'   \item `data` contains `obsTime` & `state_at_obsTime`
#'   \item `data` does not contain `x0` & `w0` as those are reserved for intercepts
#'   \item `type` argument within `controls` has the correct options
#'   \item `id` within `data` appears more than once (i.e., at least 2 obs or equivalently 1 transition)
#' }
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]

check_lctmc = function(data,
                       X_names,
                       W_names,
                       K,
                       par_constraint,
                       controls,
                       parallel_optim,
                       MyModelName) {
  ### checks (1)
  if (!"lctmc_control" %in% class(controls)) {
    stop("`controls` should be a 'lctmc_control' object obtained from `LCTMC::create_controls()`")
  }
  if (!is.list(parallel_optim) || !all(c("run", "cl") %in% names(parallel_optim))) {
    stop("`parallel_optim` should be a list with elements 'run' and 'cl'")
  }
  if (!is.character(MyModelName)) {
    stop("`MyModelName` should be a character scalar.")
  }
  if (nchar(MyModelName) > 85) {
    stop("`MyModelName` is too long! Pick a shorter name (max length is 85 characters)")
  }
  if (!all(names(controls$fmt_data$scaling) %in% c("dt", X_names, W_names))) {
    stop("names of `scaling` must match with variable names specified in `X_names` and `W_names`")
  }
  if (!is.null(par_constraint) && any(par_constraint != 0)) {
    stop("Currently the 'LCTMC' package only supports constraints equal to 0")
  }
  if (K <= 1 || !is.integer(K)) {
    stop("`K` should be greater than 1, and it should be of an 'integer' class object (e.g., K = 3L)")
  }
  if ((length(controls$init02$which_step1) != 1) || !(controls$init02$which_step1 %in% c("all", "best"))) {
    stop("`which_step1` should be a length 1, either '100%' or 'best'")
  }

  if (!"id" %in% colnames(data)) {
    stop("`id` should be a column name in `data` serving as the individual level identifier")
  }
  if (!all(c("state_at_obsTime", "obsTime") %in% colnames(data))) {
    stop("`state_at_obsTime` and `obsTime` should be column names in `data`")
  }
  if (any(c('x0', 'w0') %in% colnames(data))) {
    stop("`x0` and `w0` cannot be column names in `data`. These are reserved for the intercept terms")
  }
  if (length(controls$type) != 1 || !(controls$type %in% c("2x2", "3x3"))) {
    stop("`type` should be a length 1 character vector. Equals either '2x2' or '3x3'")
  }

  ### checks (2)
  temp = do.call(`c`, Map(`nrow`, split(x = data, f = data$id)))
  temp_obs_count = unname(temp)
  temp_id = names(temp)
  if (any(temp_obs_count == 1)) {
    err_msg = paste("Some IDs only appear once in the input data frame, `data`.", sep = "")
    err_msg = paste(err_msg, " All `data$id` should have at least two observations.", sep = "")
    stop(err_msg)
  }
}
