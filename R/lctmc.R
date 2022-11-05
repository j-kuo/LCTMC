#' @title Fit a Latent Class CTMC model
#'
#' @description Fits a latent class CTMC where we assume the data generating process is a CTMC when conditioned on a un-observed latent variable.
#'
#' @name lctmc
#'
#' @param data a data frame object with data stored in long-format
#' @param X_names a named character vector. Hosting the names of covariates for the CTMC model. It should be a column in `data`. \cr
#' Best set to x0, x1, x2 to avoid errors
#' @param W_names a named character vector. Hosting names of covariates for the latent class component. It should be a column in `data`. \cr
#' Best set to w0, w1, w2 to avoid errors
#' @param K an integer scalar. Use this variable to tell the function how many latent classes there should be. \cr
#' Note that the number of latent classes will affect the number of parameters in the model.
#' @param par_constraint a named numeric vector to indicate which parameter is constrained. Set equal to NULL for unconstrained model. \cr
#' For example, `c(alpha1.1 = 0)` constraints the parameter 'alpha1.1' to be a constant 0. **NOTE:** Current version of the code will *only* work with constrains equal to 0.
#' @param controls a "lctmc_control" object obtained from [create_controls()]. Through this function user will be able to control various
#' aspect of the algorithm. See the documentation [create_controls()] for details on what can be controlled.
#' @param parallel_optim a list object telling the function whether parallel process should be used. \cr
#' The list should contain **two** elements: \cr
#' (1) `run` a logical scalar, if TRUE then this function will use parallel processing. If FALSE, then the `cl` argument is ignored. \cr
#' (2) `cl` is an object obtained from the `parallel` package, for example \cr `cl = parallel::makeCluster(spec = 2)`
#' @param MyModelName a character scalar. Gives the current model fitting process a name. This name will be used when the function is logging the algorithm progress.
#'
#' @return A list object containing the following elements:
#' \itemize{
#'   \item **init01** a list object obtained from the `gen_inits01_lctmc` functions
#'   \item **init02** a named numeric vector obtained from the `gen_inits02_lctmc` functions
#'   \item **EM** a list object obtained from the `EM_lctmc` functions
#'   \item **SE** a list object obtained from the `get_SE_lctmc` functions
#'   \item **K** an integer scalar that is identical to the input argument, `K`.
#'   \item **n_pars** an integer scalar that indicates the number of parameters estimated
#'         (total number of parameter minus number of constrained parameters).
#'   \item **X_names** a character vector equivalent to the input argument `X_names`
#'   \item **W_names** a character vector equivalent to the input argument `W_names`
#'   \item **run_time** a "difftime" object. Indicating the total algorithm run time.
#' }
#'
#' @note The model fitting process essentially breaks down into the following steps:
#' \enumerate{
#'   \item data processing, formatting, scaling
#'   \item generate initial values step 1, individualized fitting
#'   \item generate initial values step 2, optimization on \eqn{log(P(Y))}
#'   \item EM algorithm to obtain tighter estimate
#'   \item hessian approximation for SE
#'   \item re-scale parameters back to original scale
#' }
#'
#' @seealso [fmt_rowwise_trans()]; [gen_inits01_lctmc_2x2()]; [gen_inits02_lctmc_2x2()];
#' [EM_lctmc_2x2()]; [get_SE_lctmc_2x2()]; [rescale_theta()]
#'
#' @example inst/examples/ex_running_lctmc.R
NULL

#' @rdname lctmc
#' @export
lctmc_2x2 = function(data = data.frame(),
                     X_names = c(),
                     W_names = c(),
                     K = integer(),
                     par_constraint,
                     controls = create_controls(),
                     parallel_optim = list(run = FALSE, cl = NA),
                     MyModelName = "lctmc_2x2") {
  ### checking
  if (!all(names(controls$fmt_data$scaling) %in% c("dt", X_names, W_names))) {
    stop("names of `scaling` must match with variable names specified in `X_names` and `W_names`")
  }
  if (any(par_constraint != 0)) {
    stop("Currently the 'LCTMC' package only supports constraints equal to 0")
  }
  if (K <= 1 || !is.integer(K)) {
    stop("`K` should be greater than 1, and it should be of an 'integer' class object (e.g., K = 3L)")
  }
  if ((length(controls$init02$which_step1) != 1) || !(controls$init02$which_step1 %in% c("all", "best"))) {
    stop("`which_step1` should be a length 1, either '100%' or 'best'")
  }

  ### starting time
  t0 = Sys.time()
  cat("RUN DATE: ", as.character(Sys.Date()), "\n", sep = "")
  cat("--------------------\n\n", sep = "")

  ### process data + scaling if any
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "format", MyModelName = MyModelName)

  my_df = fmt_rowwise_trans(
    data = data,
    type = controls$type,
    X_names = X_names,
    W_names = W_names,
    scaling = controls$fmt_data$scaling,
    trace = controls$fmt_data$trace
  )
  my_df.Xmat = my_df[["Xmat"]]
  my_df.Wmat = my_df[["Wmat"]]
  my_df.dt = my_df[["dt"]]
  my_df = my_df[["df_trans"]]

  trace_lctmc_progress(section = "tail1", type = "format", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### generate initial values ~ step 1 (k-mean)
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "init1", MyModelName = MyModelName)

  my_model.init01 = gen_inits01_lctmc_2x2(
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    K = K,
    par_constraint = par_constraint,
    parallel_optim = parallel_optim,
    N_sub = controls$init01$N_sub,
    pct_keep = controls$init01$pct_keep,
    parallelize = controls$init01$parallelize
  )

  trace_lctmc_progress(section = "tail1", type = "init1", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### generate initial values ~ step 2 (direct optim)
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "init2", MyModelName = MyModelName)
  cat(" * using `which_step1` = '", controls$init02$which_step1, "' as the intial values for Step 2 \n", sep = "")

  if (controls$init02$which_step1 == "all") {
    step2_inits = my_model.init01[["step1_full"]]$theta
  } else {
    step2_inits = my_model.init01[["step1_best"]]
  }

  my_model.init02 = gen_inits02_lctmc_2x2(
    step2_inits = step2_inits,
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    K = K,
    par_constraint = par_constraint,
    parallel_optim = parallel_optim
  )

  trace_lctmc_progress(section = "tail1", type = "init2", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### fit EM
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "em", MyModelName = MyModelName)

  my_model.EM = EM_lctmc_2x2(
    # theta
    EM_inits = my_model.init02,
    par_constraint = par_constraint,
    K = K,
    # data
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    # controls (EM parts)
    EM.maxit = controls$EM$EM.maxit,
    ELL_diff.tol = controls$EM$EM.ELL_tol,
    LPY_diff.tol = controls$EM$EM.LPY_tol,
    par_diff.tol = controls$EM$EM.par_tol,
    # controls (L-BFGS-B parts)
    LBFGSB.maxit = controls$EM$LBFGSB.maxit,
    fnscale = controls$EM$LBFGSB.fnscale,
    factr = controls$EM$LBFGSB.factr,
    parallel_optim = parallel_optim
  )

  trace_lctmc_progress(section = "tail1", type = "em", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### compute information matrix
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "se", MyModelName = MyModelName)

  my_model.SE = get_SE_lctmc_2x2(
    em = my_model.EM,
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    K = K,
    par_constraint = par_constraint,
    solve.tol = controls$SE$solve_tol,
    symmetric.tol = controls$SE$symmetric_tol,
    eigen0.tol = controls$SE$eigen0_tol
  )

  trace_lctmc_progress(section = "tail1", type = "se", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### re-scaling estimates
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "rescale", MyModelName = MyModelName)

  my_model.rescale = rescale_theta(
    df_theta = my_model.SE$SE,
    v_theta = my_model.init01$step1_full$theta,
    scaling = controls$fmt_data$scaling,
    mult_vars = c("mle_theta", "SE", "L_CI", "U_CI"),
    add_vars = c("mle_theta", "L_CI", "U_CI")
  )
  my_model.SE$SE = my_model.rescale$df_theta
  my_model.init01$step1_full$theta = my_model.rescale$v_theta

  trace_lctmc_progress(section = "tail1", type = "rescale", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### close parallel connection
  if (parallel_optim$run) {
    parallel::stopCluster(cl = parallel_optim$cl)
  }

  ### misc. output for convenience
  data = data.frame(data)
  n_pars = nrow(my_model.SE$Covariance)
  tf = Sys.time()

  ### output
  out = list(
    init01 = my_model.init01,
    init02 = my_model.init02,
    EM = my_model.EM,
    SE = my_model.SE,
    K = K,
    n_pars = n_pars,
    X_names = X_names,
    W_names = W_names,
    run_time = difftime(tf, t0, units = 'hour')
  )
  class(out) = c("lctmc", "lctmc_2x2", "list")
  return(out)
}

#' @rdname lctmc
#' @export
lctmc_3x3 = function(data = data.frame(),
                     X_names = c(),
                     W_names = c(),
                     K = integer(),
                     par_constraint,
                     controls = create_controls(),
                     parallel_optim = list(run = FALSE, cl = NA),
                     MyModelName = "lctmc_3x3") {
  ### checking
  if (!all(names(controls$fmt_data$scaling) %in% c("dt", X_names, W_names))) {
    stop("names of `scaling` must match with variable names specified in `X_names` and `W_names`")
  }
  if (any(par_constraint != 0)) {
    stop("Currently the 'LCTMC' package only supports constraints equal to 0")
  }
  if (K <= 1 || !is.integer(K)) {
    stop("`K` should be greater than 1, and it should be of an 'integer' class object (e.g., K = 3L)")
  }
  if ((length(controls$init02$which_step1) != 1) || !(controls$init02$which_step1 %in% c("all", "best"))) {
    stop("`which_step1` should be a length 1, either '100%' or 'best'")
  }

  ### starting time
  t0 = Sys.time()
  cat("RUN DATE: ", as.character(Sys.Date()), "\n", sep = "")
  cat("-------------------- \n\n", sep = "")

  ### process data + scaling if any
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "format", MyModelName = MyModelName)

  my_df = fmt_rowwise_trans(
    data = data,
    type = controls$type,
    X_names = X_names,
    W_names = W_names,
    scaling = controls$fmt_data$scaling,
    trace = controls$fmt_data$trace
  )
  my_df.Xmat = my_df[["Xmat"]]
  my_df.Wmat = my_df[["Wmat"]]
  my_df.dt = my_df[["dt"]]
  my_df = my_df[["df_trans"]]

  trace_lctmc_progress(section = "tail1", type = "format", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### generate initial values ~ step 1 (k-mean)
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "init1", MyModelName = MyModelName)

  my_model.init01 = gen_inits01_lctmc_3x3(
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    K = K,
    par_constraint = par_constraint,
    parallel_optim = parallel_optim,
    N_sub = controls$init01$N_sub,
    pct_keep = controls$init01$pct_keep,
    parallelize = controls$init01$parallelize
  )

  trace_lctmc_progress(section = "tail1", type = "init1", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### generate initial values ~ step 2 (direct optim)
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "init2", MyModelName = MyModelName)
  cat(" * using `which_step1` = '", controls$init02$which_step1, "' as the intial values for Step 2 \n", sep = "")

  if (controls$init02$which_step1 == "all") {
    step2_inits = my_model.init01[["step1_full"]]$theta
  } else {
    step2_inits = my_model.init01[["step1_best"]]
  }

  my_model.init02 = gen_inits02_lctmc_3x3(
    step2_inits = step2_inits,
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    K = K,
    par_constraint = par_constraint,
    parallel_optim = parallel_optim
  )

  trace_lctmc_progress(section = "tail1", type = "init2", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### fit EM
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "em", MyModelName = MyModelName)

  my_model.EM = EM_lctmc_3x3(
    # theta
    EM_inits = my_model.init02,
    par_constraint = par_constraint,
    K = K,
    # data
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    # controls (EM parts)
    EM.maxit = controls$EM$EM.maxit,
    ELL_diff.tol = controls$EM$EM.ELL_tol,
    LPY_diff.tol = controls$EM$EM.LPY_tol,
    par_diff.tol = controls$EM$EM.par_tol,
    # controls (L-BFGS-B parts)
    LBFGSB.maxit = controls$EM$LBFGSB.maxit,
    fnscale = controls$EM$LBFGSB.fnscale,
    factr = controls$EM$LBFGSB.factr,
    parallel_optim = parallel_optim
  )

  trace_lctmc_progress(section = "tail1", type = "em", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### compute information matrix
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "se", MyModelName = MyModelName)

  my_model.SE = get_SE_lctmc_3x3(
    em = my_model.EM,
    df = my_df,
    df_Xmat = my_df.Xmat,
    df_Wmat = my_df.Wmat,
    df_dt = my_df.dt,
    K = K,
    par_constraint = par_constraint,
    solve.tol = controls$SE$solve_tol,
    symmetric.tol = controls$SE$symmetric_tol,
    eigen0.tol = controls$SE$eigen0_tol
  )

  trace_lctmc_progress(section = "tail1", type = "se", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### re-scaling estimates
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "rescale", MyModelName = MyModelName)

  my_model.rescale = rescale_theta(
    df_theta = my_model.SE$SE,
    v_theta = my_model.init01$step1_full$theta,
    scaling = controls$fmt_data$scaling,
    mult_vars = c("mle_theta", "SE", "L_CI", "U_CI"),
    add_vars = c("mle_theta", "L_CI", "U_CI")
  )
  my_model.SE$SE = my_model.rescale$df_theta
  my_model.init01$step1_full$theta = my_model.rescale$v_theta

  trace_lctmc_progress(section = "tail1", type = "rescale", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### close parallel connection
  if (parallel_optim$run) {
    parallel::stopCluster(cl = parallel_optim$cl)
  }

  ### misc. output for convenience
  data = data.frame(data)
  n_pars = nrow(my_model.SE$Covariance)
  tf = Sys.time()

  ### output
  out = list(
    init01 = my_model.init01,
    init02 = my_model.init02,
    EM = my_model.EM,
    SE = my_model.SE,
    K = K,
    n_pars = n_pars,
    X_names = X_names,
    W_names = W_names,
    run_time = difftime(tf, t0, units = 'hour')
  )
  class(out) = c("lctmc", "lctmc_3x3", "list")
  return(out)
}
