#' @title Creates the control object for fitting LCTMC
#'
#' @description Create a custom object that acts like a list which serve as object that holds all control
#' parameters for fitting the LCTMC model.
#'
#' @param type a character scalar. Equals either "2x2" or "3x3" to indicate the number of category in the outcome variable.
#' @param ... can be the following values:
#' \describe{
#'   \item{data}{certain parameters will require an input data frame to compute the default values. \cr
#'         Use `data` to specify this data frame.}
#'   \item{scaling}{controls what factor covariate should be scaled to ensure algorithm stability. \cr
#'         See [fmt_rowwise_trans()] for more info. Default value depends on `data`}
#'   \item{trace}{controls whether data formatting function should print function details. \cr
#'         See [fmt_rowwise_trans()] for more info}
#'   \item{N_sub}{controls how many subjects to use for step 1 of initial value generation. \cr
#'         See [gen_inits01_lctmc_2x2()] for more info}
#'   \item{pct_keep}{controls how much of left/right tail should be truncated for the k-means algorithm. \cr
#'         See [gen_inits01_lctmc_2x2()] for more info}
#'   \item{parallelize.init01}{controls whether step 1 of initial value generation should be parallelized. \cr
#'         See [gen_inits01_lctmc_2x2()] for more info}
#'   \item{which_step1}{a character scalar. Controls which kmeans result to use for step 2 of initial value generation. \cr
#'         Default is "best" which uses the parameter vector that yields the highest \eqn{\log(P(Y))}.
#'         The alternative is "all" which uses all of the kmeans estimate without truncating either tail.}
#'   \item{EM.maxit}{controls how many EM iterations the algorithm will perform. \cr
#'         See [EM_lctmc_2x2()] for more info}
#'   \item{EM.ELL_tol}{controls the convergence tolerance on the expected conditional log-likelihood value. \cr
#'         See [EM_lctmc_2x2()] for more info}
#'   \item{EM.LPY_tol}{controls the convergence tolerance on the observed log-likelihood value. \cr
#'         See [EM_lctmc_2x2()] for more info}
#'   \item{EM.par_tol}{controls the convergence tolerance on the magnitude change in parameter values. \cr
#'         See [EM_lctmc_2x2()] for more info}
#'   \item{LBFGSB.fnscale}{controls the `fnscale` argument for `optim()`. \cr
#'         See [EM_lctmc_2x2()] & [optim()] for more info. Default value depends on `data`.}
#'   \item{LBFGSB.maxit}{controls the `maxit` argument for `optim()`. \cr
#'         See [EM_lctmc_2x2()] & [optim()] for more info}
#'   \item{LBFGSB.factr}{controls the `factr` argument for `optim()`. \cr
#'         See [EM_lctmc_2x2()] & [optim()] for more info}
#'   \item{solve_tol}{controls the tolerance value for detecting linear dependency when `solve()` is called. \cr
#'         See [get_SE_lctmc_2x2()] & [solve()] for more info}
#'   \item{symmetric_tol}{controls the tolerance value when checking if the hessian matrix is symmetric. \cr
#'         See [get_SE_lctmc_2x2()] for more info}
#'   \item{eigen0_tol}{controls the tolerance value for treating eigen values as essentially 0. \cr
#'         See [get_SE_lctmc_2x2()] & [eigen()] for more info}
#' }
#'
#' @return a custom class object that acts like a list. The output contains the following:
#' \describe{
#'   \item{fmt_data}{contains elements: `scaling` and `trace`.}
#'   \item{init01}{contains elements: `N_sub`, `pct_keep`, and `parallelize.init01`}
#'   \item{init02}{contains elements: `which_step1`.}
#'   \item{EM}{contains elements: `EM.maxit`, `EM.ELL_tol`, `EM.LPY_tol`, `EM.par_tol`,
#'         `LBFGS.fnscale`, `LBFGS.maxit`, and `LBFGSB.factr`.}
#'   \item{SE}{contains elements: `solve_tol`, `symmetric_tol`, and `eigen0_tol`.}
#'   \item{rescale}{currently there are no control options for rescaling parameters.}
#'   \item{type}{controls whether contorls is for a '2x2' model or a '3x3' model.}
#' }
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]
#'
#' @export
#'
#' @example inst/examples/ex_create_controls.R

create_controls = function(type, ...) {
  ### unpack ...
  control_args = list(...)

  ### data arg
  if (is.null(control_args$scaling)) {
    if (is.null(control_args$data)) {
      scaling.default = 1
    } else {
      temp = Map(f = function(t) t[2:length(t)]-t[1:(length(t)-1)],
                 split(x = control_args$data$obsTime, f = control_args$data$id))
      temp = do.call(`c`, temp)
      scaling.default = c(
        dt = 1/max(temp),
        sapply(control_args$data[c("x1", "x2", "w1", "w2")], function(x) 1/max(x))
      )
      scaling.default = abs(scaling.default)
    }

    control_args$scaling = scaling.default
  }

  if (is.null(control_args$LBFGSB.fnscale)) {
    if (is.null(control_args$data)) {
      LBFGSB.fnscale.default = -100
    } else {
      LBFGSB.fnscale.default = -nrow(control_args$data)
    }

    control_args$LBFGSB.fnscale = LBFGSB.fnscale.default
  }


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
  control_args$trace = ifelse(is.null(control_args$trace), TRUE, control_args$trace)
  ctrl[["fmt_data"]] = list(
    scaling = control_args$scaling,
    trace = control_args$trace
  )

  ### init01 controls
  control_args$N_sub = ifelse(is.null(control_args$N_sub), Inf, control_args$N_sub)
  if (is.null(control_args$pct_keep)) {
    control_args$pct_keep = seq(0.4, 1.0, 0.001)
  }
  control_args$parallelize.init01 = ifelse(is.null(control_args$parallelize.init01), TRUE, control_args$parallelize.init01)
  ctrl[["init01"]] = list(
    N_sub = control_args$N_sub,
    pct_keep = control_args$pct_keep,
    parallelize = control_args$parallelize.init01
  )

  ### init 02 controls
  control_args$which_step1 = ifelse(is.null(control_args$which_step1), "best", control_args$which_step1)
  ctrl[["init02"]] = list(
    which_step1 = control_args$which_step1
  )

  ### EM controls
  control_args$EM.maxit = ifelse(is.null(control_args$EM.maxit), 10, control_args$EM.maxit)
  control_args$EM.ELL_tol = ifelse(is.null(control_args$EM.ELL_tol), 1e-2, control_args$EM.ELL_tol)
  control_args$EM.LPY_tol = ifelse(is.null(control_args$EM.LPY_tol), 1e-4, control_args$EM.LPY_tol)
  control_args$EM.par_tol = ifelse(is.null(control_args$EM.par_tol), 1e-4, control_args$EM.par_tol)
  control_args$LBFGSB.maxit = ifelse(is.null(control_args$LBFGSB.maxit), 1e4, control_args$LBFGSB.maxit)
  control_args$LBFGSB.factr = ifelse(is.null(control_args$LBFGSB.factr), 1e-4, control_args$LBFGSB.factr)
  ctrl[["EM"]] = list(
    EM.maxit = control_args$EM.maxit,
    EM.ELL_tol = control_args$EM.ELL_tol,
    EM.LPY_tol = control_args$EM.LPY_tol,
    EM.par_tol = control_args$EM.par_tol,
    LBFGSB.fnscale = control_args$LBFGSB.fnscale,
    LBFGSB.maxit = control_args$LBFGSB.maxit,
    LBFGSB.factr = control_args$LBFGSB.factr
  )

  ### SE approx controls
  control_args$solve_tol = ifelse(is.null(control_args$solve_tol), (.Machine$double.eps)^2, control_args$solve_tol)
  control_args$symmetric_tol = ifelse(is.null(control_args$symmetric_tol), 5e-11, control_args$symmetric_tol)
  control_args$eigen0_tol = ifelse(is.null(control_args$eigen0_tol), 1e-10, control_args$eigen0_tol)
  ctrl[["SE"]] = list(
    solve_tol = control_args$solve_tol,
    symmetric_tol = control_args$symmetric_tol,
    eigen0_tol = control_args$eigen0_tol
  )

  ### SE approx controls
  ctrl[["rescale"]] = NULL

  ### return
  class(ctrl) = append("lctmc_control", class(ctrl))
  ctrl
}
