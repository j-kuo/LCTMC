#' @title Fit a Latent Class CTMC model (2x2)
#'
#' @description A latent class CTMC (general 2x2 case) where we assume the data generating process is a CTMC when conditioned on a un-observed latent variable.
#'
#' @param data a data frame object with data stored in long-format
#' @param dt_scale a named numeric scalar. It is a scaling parameter for observation time intervals. \cr
#' For example, `dt_scale = c(dt = 1/10)` would converts a time unit from 10 to 1.
#' @param x_scale a named numeric vector. It is a scaling parameters for covariates that affect the CTMC process. \cr
#' For example, `x_scale = c(x0 = 1, x1 = 1/5, x2 = 1)` would convert x1 from 10 to 2.
#' @param w_scale a named numeric vector. It is a scaling parameters for covariates that affect the latent class component of the model  \cr
#' For example, `w_scale = c(w0 = 1, x1 = 1/3, x2 = 1)` would convert w1 from 12 to 4.
#' @param K an integer scalar. Use this variable to tell the function how many latent classes there should be. \cr
#' Note that the number of latent classes will affect the number of parameters in the model.
#' @param X_names a named character vector. Hosting the names of covariates for the CTMC model. It should be a column in `data`. Best set to x0, x1, x2 to avoid errors
#' @param W_names a named character vector. Hosting names of covariates for the latent class component. It should be a column in `data`. Best set to w0, w1, w2 to avoid errors
#' @param par_constraint a named numeric vector to indicate which parameter is constrained. Set equal to NULL for unconstrained model. \cr
#' For example, `c(alpha1.1 = 0)` constraints the parameter 'alpha1.1' to be a constant 0. **NOTE:** Current version of the code will *only* work with constrains equal to 0.
#' @param N_sub a numeric scalar. This is used for step 1 of initial value generation where the algorithm fits a traditional CTMC model for each individual. \cr
#' Fitting the model for *all* individuals might have long run time without improvement in the accuracy of the estimation.
#' Thus, using this argument to set a maximum number of individuals to use for the initial value generation could save some computation time.
#' @param pct_keep a numeric vector where each element of the vector ranges from 0 to 1 (ideally, 0.50 to 0.95).
#' This argument controls what percentage of the individual effect should be used for the K-means algorithm for initial value generation. The algorithm will consider all percentages specified in this vector. \cr
#' For example, for `pct_keep = c(0.5)`, after individuals effects are estimated, only the 25th to 75th percentile are used for the K-Means algorithm to obtain cluster level estimates. \cr
#' Additionally, note that the threshold "1" is always appended to `pct_keep`, so it will always at least consider the 100% case.
#' @param parallelize a logical scalar. Set to TRUE if we want the for-loop for the individual-wise CTMC to be parallelized
#' @param which_step1 a character scalar. Either "100%" or "best". The former will use 100% of the individual CTCM effects to perform the K-means algorithm. \cr
#' The latter will compute the \eqn{log(P(Y))} value for all possible K-means result and use the thetas that yields the largest \eqn{log(P(Y))}.
#' @param theta.names a list of parameter names, e.g., `list(c("alpha1.1", "alpha2.1"), c("beta0.21_1", "beta1.21_1", "beta2.21_1"))`. \cr
#' This argument is used to determine the number of ECM sub-steps. Parameters within the same element of the list will be optimized together as a group. \cr
#' User can use the function `gen_theta_names()` to generate this with ease.
#' @param EM_controls a list object holding the arguments to tune the EM algorithm.
#' The following elements are necessary for running the algorithm: \cr
#' (1) `maxit` a numeric scalar. This specifies the max number of EM iterations \cr
#' (2) `ELL_tol` a numeric scalar. This is the tolerance of the conditional expected log likelihood \cr
#' (3) `LPY_tol` a numeric scalar. This is the tolerance of observed log likelihood, \eqn{log(P(Y))} \cr
#' (4) `par_tol` a numeric scalar. This is the tolerance of parameter changes per EM iteration \cr \cr
#' Note that for all 3 tolerance parameters, the smaller the value, the more accurate the estimate will be. But it will also lead to longer computation time.
#' @param optim_controls a list object holding the arguments to control the L-BFGS optimization at each ECM step.
#' The following elements are necessary for running the algorithm: \cr
#' (1) `fnscale` a numeric scalar. This value scales the objective function. Additionally, its sign determines whether the algorithm is performing a maximization or minimization task. \cr
#' (2) `maxit` a numeric scalar. This value specifies the max number of L-BFGS iterations. \cr
#' (3) `factr` a numeric scalar. This value controls the tolerance of L-BFGS optimization. The smaller in magnitude this argument is the more precise the optimization algorithm will be.
#' @param test_if_global_optim a list containing two elements: `test` and `true_params`. \cr
#' (1) `test` is a logical scalar. It indicates whether the function should check whether the MLE is actually a global optimal point.
#' This is done by comparing the observed log likelihood at MLE vs. the log likelihood at the true parameter value. Hence, the second element of the list is `true_params`.
#' For this reason, it only makes sense to check when fitting simulated datasets where we would know the true parameter values. \cr
#' (2) `true_params` should be generated using the "LCTMC.simulate" package via the `gen_true_param()` function.
#' @param parallel_optim a list object telling the function whether parallel process should be used for the Step 2 of the initial value generation. \cr
#' The list should contain **two** elements: \cr
#' (1) `run` a logical scalar, if TRUE then this function will use parallel processing. If FALSE, then the `cl` argument is ignored. \cr
#' (2) `cl` is an object obtained from the `parallel` package, for example \cr `cl = parallel::makeCluster(spec = 2)`
#' @param MyModelName a character scalar. Gives the current model fitting process a name. This name will be used when the function is logging the algorithm progress.
#'
#' @return A list object containing the 4 elements:
#' \itemize{
#'   \item **inits** a list object obtained from the `gen_inits_lctmc_2x2()` function
#'   \item **EM** a list object obtained from the `EM_lctmc_2x2()` function
#'   \item **SE** a list object obtained from the `get_SE_lctmc_2x2()` function
#'   \item **global_optim** a numeric scalar that indicates whether the MLE has converged to the global optimal point. See notes on the input argument for `test_if_global_optim`.
#'   \item **data** a data.frame object that is identical to the input argument, `data`.
#'   \item **K** an integer scalar that is identical to the input argument, `K`.
#'   \item **n_pars** an integer scalar that indicates the number of parameters estimated (total number of parameter minus number of constrained parameters).
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
#'   \item test for global optimal if specified to be tested
#' }
#'
#' @export
#'
#' @example inst/examples/ex_lctmc_2x2.R

lctmc_2x2 = function(data = data.frame(),
                     dt_scale = c(),
                     x_scale = c(),
                     w_scale = c(),
                     K = integer(),
                     X_names = c(),
                     W_names = c(),
                     par_constraint,
                     N_sub,
                     pct_keep = c(),
                     parallelize = FALSE,
                     which_step1 = c("100%", "best"),
                     theta.names = list(),
                     EM_controls = list(),
                     optim_controls = list(),
                     test_if_global_optim = list(test = FALSE, true_params = NA),
                     parallel_optim = list(run = FALSE, cl = NA),
                     MyModelName = "lctmc_2x2") {
  ### checking
  if (any(names(x_scale) != X_names)) {
    stop("names of `x_scale` must match with variable names specified in `X_names`")
  }
  if (any(names(w_scale) != W_names)) {
    stop("names of `w_scale` must match with variable names specified in `W_names`")
  }
  if (!all(names(par_constraint) %in% unlist(theta.names))) {
    stop("some constrained parameters are not specified in the `theta.names` argument")
  }
  if (K <= 1 || !is.integer(K)) {
    stop("`K` should be greater than 1, and it should be of an 'integer' class object (e.g., K = 3L)")
  }

  ### starting time
  t0 = Sys.time()
  cat("RUN DATE: ", as.character(Sys.Date()), "\n", sep = "")
  cat("--------------------\n\n", sep = "")

  ### data (transition indicators), X matrix, W matrix, and dt
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "format", MyModelName = MyModelName)
  my_df = fmt_rowwise_2x2trans(data = data)
  scaling = c(dt_scale, x_scale, w_scale)

  ### if any scaling factor is not 1
  if (any(scaling != 1)) {
    ## loop through each variable that need to be scaled
    for (sc in seq_along(scaling)) {
      if (scaling[sc] != 1) {
        # msg
        cat(" * `", names(scaling)[sc], "` being scaled by a factor of ", sprintf("%.3f", scaling[sc]), "\n", sep = "")

        # scaling
        my_df[[names(scaling)[sc]]] = my_df[[names(scaling)[sc]]] * scaling[sc]
      }
    }

    ## obtain X mat, an observation level matrix
    my_df.Xmat = as.matrix(my_df[, X_names])

    ## obtain W mat, an individual level matrix
    my_df.Wmat = unique(my_df[, c("id", W_names)])
    my_df.Wmat = as.matrix(my_df.Wmat[!colnames(my_df.Wmat) %in% "id"])

    ## obtain dt, an observation level vector
    my_df.dt = my_df$dt

    ## obtain data frame with row-wise transition indicators
    my_df = my_df[!colnames(my_df) %in% c(X_names, W_names, "dt")]
  } else {
    cat(" * no scaling were applied", "\n", sep = "")
  }
  trace_lctmc_progress(section = "tail1", type = "format", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### generate initial values (step1 + step2)
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "init1", MyModelName = MyModelName)
  my_model.inits = gen_inits_lctmc_2x2(
    df = my_df,
    Xmat = my_df.Xmat,
    Wmat = my_df.Wmat,
    dt = my_df.dt,
    N_sub = N_sub,
    pct_keep = pct_keep,
    par_constraint = par_constraint,
    K = K,
    returns = 2L,
    parallelize = parallelize,
    which_step1 = which_step1,
    parallel_optim = parallel_optim,
    MyModelName = MyModelName
  )
  trace_lctmc_progress(section = "tail1", type = "init1", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### fit EM
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "em", MyModelName = MyModelName)
  my_model.EM = EM_lctmc_2x2(
    # theta
    theta.init = my_model.inits$step2,
    theta.names = theta.names,
    par_constraint = par_constraint,
    K = K,
    # data
    df = my_df,
    df.Xmat = my_df.Xmat,
    df.Wmat = my_df.Wmat,
    df.dt = my_df.dt,
    # controls (EM parts)
    EM.maxit = EM_controls$maxit,
    ELL_diff.tol = EM_controls$ELL_tol,
    LPY_diff.tol = EM_controls$LPY_tol,
    par_diff.tol = EM_controls$par_tol,
    # controls (L-BFGS-B parts)
    L_BFGS_B.ctrl = optim_controls,
    parallel_optim = parallel_optim,
    MyModelName = MyModelName
  )
  trace_lctmc_progress(section = "tail1", type = "em", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### compute information matrix
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "se", MyModelName = MyModelName)
  my_model.SE = get_SE_lctmc_2x2(
    em = my_model.EM,
    df = my_df,
    df.Xmat = my_df.Xmat,
    df.Wmat = my_df.Wmat,
    df.dt = my_df.dt,
    par_constraint = par_constraint,
    K = K,
    MyModelName = MyModelName
  )
  trace_lctmc_progress(section = "tail1", type = "se", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### re-scaling estimates
  trace_lctmc_progress(section = "header1", type = "format", MyModelName = MyModelName)
  trace_lctmc_progress(section = "header2", type = "rescale", MyModelName = MyModelName)
    ## constants
    MULTS_pars = c("mle_theta", "SE", "L_CI", "U_CI")
    ADD_pars = c("mle_theta", "L_CI", "U_CI")
    coef_df = my_model.SE$SE
    inits = my_model.inits$step1_full

    ## re-scale alphas
    # df of MLE
    coef_df[grepl(pattern = "alpha1", coef_df$names), MULTS_pars] = coef_df[grepl(pattern = "alpha1", coef_df$names), MULTS_pars] * scaling['w1']
    coef_df[grepl(pattern = "alpha2", coef_df$names), MULTS_pars] = coef_df[grepl(pattern = "alpha2", coef_df$names), MULTS_pars] * scaling['w2']
    # vector of step1 full
    inits[grepl(pattern = "alpha1", names(inits))] = inits[grepl(pattern = "alpha1", names(inits))] * scaling['w1']
    inits[grepl(pattern = "alpha2", names(inits))] = inits[grepl(pattern = "alpha2", names(inits))] * scaling['w2']
    cat(" * alphas (1 & 2) re-scaled \n", sep = "")

    ## re-scale slope of q_rs
    # df of MLE
    coef_df[grepl(pattern = "beta1", coef_df$names), MULTS_pars] = coef_df[grepl(pattern = "beta1", coef_df$names), MULTS_pars] * scaling['x1']
    coef_df[grepl(pattern = "beta2", coef_df$names), MULTS_pars] = coef_df[grepl(pattern = "beta2", coef_df$names), MULTS_pars] * scaling['x2']
    # vector of step1 full
    inits[grepl(pattern = "beta1", names(inits))] = inits[grepl(pattern = "beta1", names(inits))] * scaling['x1']
    inits[grepl(pattern = "beta2", names(inits))] = inits[grepl(pattern = "beta2", names(inits))] * scaling['x2']
    cat(" * betas (1 & 2) re-scaled \n", sep = "")

    ## re-scale intercept of q_rs
    # df of MLE
    coef_df[grepl(pattern = "beta0", coef_df$names), ADD_pars] = coef_df[grepl(pattern = "beta0", coef_df$names), ADD_pars] + log(scaling['dt'])
    # vector of step1 full
    inits[grepl(pattern = "beta0", names(inits))] = inits[grepl(pattern = "beta0", names(inits))] + log(scaling['dt'])
    cat(" * beta (0) re-scaled \n", sep = "")

    ## update
    my_model.SE$SE = coef_df
    my_model.inits$step1_full = inits
  trace_lctmc_progress(section = "tail1", type = "rescale", ref_t = t0, MyModelName = MyModelName)
  trace_lctmc_progress(section = "tail2", type = "format", ref_t = t0, MyModelName = MyModelName)


  ### test if global optimum has been reached by comparing with true parameter values (i.e., only used when evaluating simulation result)
  global_optim = NULL
  if (test_if_global_optim$test && K == 3) {
    ## generate theta names for the bik() function
    theta.names.bik = gen_theta_names(K = K, type = "2x2", purpose = "bik")

    ## various data objects (non-scaled) ~ this is because true parameter values are not scaled
    my_df.test = fmt_rowwise_2x2trans(data = data)
    my_df.Xmat.test = as.matrix(my_df.test[, X_names])
    my_df.Wmat.test = unique(my_df.test[, c("id", W_names)])
    my_df.Wmat.test = as.matrix(my_df.Wmat.test[!colnames(my_df.Wmat.test) %in% c("id")])
    my_df.dt.test = my_df.test$dt
    my_df.test = my_df.test[!colnames(my_df.test) %in% c(X_names, W_names, "dt")]

    ## get true parameter as a named vector
    mle_tht = my_model.SE$SE$mle_theta
    names(mle_tht) = my_model.SE$SE$names
    df_true_tht = align_MLE_2x2(true = test_if_global_optim$true_params, mle = mle_tht, K = K)
    true_tht = df_true_tht$true_theta
    names(true_tht) = df_true_tht$names

    ## compute LPY at True params
    bik_at_true = bik_all_2x2(
      theta = true_tht,
      data = my_df.test,
      Xmat = my_df.Xmat.test,
      Wmat = my_df.Wmat.test,
      dt = my_df.dt.test,
      K = K,
      P.rs = FALSE,
      theta.names = theta.names.bik
    )
    LPY_at_true = sum(log(Reduce(`+`, bik_at_true)))

    ## compute LPY at MLE params
    bik_at_mle = bik_all_2x2(
      theta = mle_tht,
      data = my_df.test,
      Xmat = my_df.Xmat.test,
      Wmat = my_df.Wmat.test,
      dt = my_df.dt.test,
      K = K,
      P.rs = FALSE,
      theta.names = theta.names.bik
    )
    LPY_at_mle = sum(log(Reduce(`+`, bik_at_mle)))

    ## result
    if ((LPY_at_true - LPY_at_mle) >= 1) {
      # when "mle" is NOT global
      global_optim = 0
    } else {
      # when "mle" is global
      global_optim = 1
    }
  }


  ### close parallel connection
  if (parallel_optim$run) {
    parallel::stopCluster(cl = parallel_optim$cl)
  }

  ### misc. output for convenience
  data = data.frame(data)
  n_pars = length(unlist(theta.names)) - length(par_constraint)
  tf = Sys.time()

  ### output
  out = list(
    inits = my_model.inits,
    EM = my_model.EM,
    SE = my_model.SE,
    global_optim = global_optim,
    data = data,
    K = K,
    n_pars = n_pars,
    X_names = X_names,
    W_names = W_names,
    run_time = difftime(tf, t0, units = 'hour')
  )
  class(out) = c("lctmc", "lctmc_2x2", "list")
  return(out)
}
