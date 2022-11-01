#' @title Creates the control object for fitting LCTMC
#'
#' @description Create a custom object that acts like a list which serve as object that holds all control
#' parameters for fitting the LCTMC model.
#'
#' @param type a character scalar. Equals either "2x2" or "3x3" to indicate the number of category in the outcome variable.
#' @param ... can be the following values:
#' \describe{
#'   \item{scaling}{controls what factor covariate should be scaled to ensure algorithm stability. \cr
#'         See [fmt_rowwise_trans()] for more info}
#'   \item{trace}{controls whether data formatting function should print function details. \cr
#'         See [fmt_rowwise_trans()] for more info}
#'   \item{N_sub}{controls how many subjects to use for step 1 of initial value generation. \cr
#'         See [gen_inits01_lctmc_2x2()] for more info}
#'   \item{pct_keep}{controls how much of left/right tail should be truncated for the k-means algorithm. \cr
#'         See [gen_inits01_lctmc_2x2()] for more info}
#'   \item{parallelize.init01}{controls whether step 1 of initial value generation should be parallelized. \cr
#'         See [gen_inits01_lctmc_2x2()] for more info}
#'   \item{arg6}{...}
#'   \item{arg7}{...}
#'   \item{arg8}{...}
#'   \item{arg9}{...}
#'   \item{arg9}{...}
#'   \item{arg9}{...}
#' }
#'
#' @return a custom class object that acts like a list. The output contains the following:
#' \describe{
#'   \item{fmt_data}{contains elements: `scaling` and `trace`}
#'   \item{init01}{...}
#'   \item{init02}{...}
#'   \item{EM}{...}
#'   \item{SE}{...}
#'   \item{rescale}{...}
#'   \item{type}{controls whether contorls is for a 2x2 model or a 3x3 model.}
#' }
#'
#' @note XXXXXXXXXXXXXXXX
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]
#'
#' @export
#'
#' @example inst/examples/ex_fmt_rowwise_trans.R

create_controls = function(type, ...) {
  ### unpack ...
  control_args = list(...)

  ### initialize control object
  ctrl = vector(mode = "list", length = 6)
  names(ctrl) = c("fmt_data", "init01", "init02", "EM", "SE", "rescale")

  ### check `type`
  if (!type %in% c("2x2", "3x3")) {
    stop("`type` should be either '2x2' or '3x3'")
  }

  ### `type`
  ctrl$type = type

  ### format data controls
  control_args$trace = ifelse(is.null(control_args$trace), FALSE, control_args$trace)
  ctrl[["fmt_data"]] = list(
    scaling = control_args$scaling,
    trace = control_args$trace
  )

  ### init01 controls
  control_args$N_sub = ifelse(is.null(control_args$N_sub), Inf, control_args$N_sub)
  ctrl[["init01"]] = list(
    N_sub = control_args$N_sub,
    pct_keep = control_args$pct_keep,
    paralleize = control_args$parallelize.init01
  )

  ### init 02 controls
  control_args$which_step1 = ifelse(is.null(control_args$which_step1), "best", control_args$which_step1)
  ctrl[["init02"]] = list(
    which_step1 = control_args$which_step1
  )

  ### EM controls
  ctrl[["EM"]] = list(

  )

  ### SE approx controls
  ctrl[["SE"]] = list(

  )

  ### SE approx controls
  ctrl[["rescale"]] = list(

  )



  ### return
  class(ctrl) = append("lctmc_control", class(ctrl))
  return(ctrl)
}
