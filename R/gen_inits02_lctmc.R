#' @title Step 2/2 for generating initial value for model fitting
#'
#' @description Performs step 2 of the initial value generation. This is a direct optimization on the observed data likelihood function.
#' We condition on the number of latent clusters and sum over joint probability between the observed data and the unobserved data. \cr
#' That is,
#' \deqn{
#'   P(Y_{i}) = \sum_{k=1}^{K} P(Y_{i} \cap Z_{i}=k)
#' }
#' See **Note** section for more info.
#'
#' @name gen_inits02_lctmc
#'
#' @param step2_inits a named numeric vector. This should be a vector of model parameters. \cr
#' The names of the vector should be parameters names and the values are the initial value for direct likelihood optimization
#' @param df a data frame object containing row-wise transition data as binary variables. \cr
#' For example, if `trans.2_1` equals 1 then it means the observation was a transition from stage 2 to stage 1 within `df_dt` amount of time.
#' @param df_Xmat a matrix object housing the covariates that affect the CTMC portion of the model. \cr
#' This matrix should have the same number of rows as the data frame object, `df`
#' @param df_Wmat a matrix object housing the covariates that affect the latent classification part of the model. \cr
#' This matrix should have number of rows equal to unique number of individuals in the data frame object, `df`
#' @param df_dt a numeric vector housing the length of time interval between observations.
#' This vector's length should be equal to number of rows in the data frame object, `df`
#' @param K an integer scalar. Use this variable to tell the function how many latent classes there should be. \cr
#' @param par_constraint See documentation in [lctmc_2x2()] or [lctmc_3x3()]
#' @param parallel_optim See documentation in [lctmc_2x2()] or [lctmc_3x3()]
#' @param LBFGSB.maxit a numeric scalar. This is the maximum number of iterations set for the L-BFGS-B algorithm.
#' See the `maxit` argument under [optim()]'s **Details** section for more info.
#' @param LBFGSB.factr a numeric scalar. This is the convergence tolerance value for the L-BFGS-B algorithm.
#' See the `factr` argument under [optim()]'s **Details** section for more info.
#'
#' @return a named numeric vector. This is the vector of model parameters that maximizes the observed data log-likelihood function (see **Note** section).
#'
#' @note This is the third step out of six of fitting a latent class CTMC model (i.e., generate initial values via direct optimization). \cr\cr
#' This step is equivalent to obtaining the MLE. It performs numerical optimization on the observed data log likelihood function:
#' \deqn{
#'   log(P(Y)) = \sum_{i} log(P(Y_{i}))
#' }
#' where
#' \deqn{
#'   P(Y_{i}) = \sum_{k} P(Y_{i}=y_{i} \cap Z_{i}=k)
#' }
#' and the summand could be simplified as:
#' \deqn{
#'   P(Y_{i} \cap Z_{i}) = P(Z_{i}) \cdot P(Y_{i} | Z_{i})
#' }
#'
#' @seealso [lctmc_2x2()]; [gen_inits01_lctmc_2x2()]; [EM_lctmc_2x2()]
#'
#' @importFrom optimParallel optimParallel
#'
#' @example inst/examples/ex_running_lctmc.R
NULL

#' @rdname gen_inits02_lctmc
gen_inits02_lctmc_2x2 = function(step2_inits,
                                 df,
                                 df_Xmat,
                                 df_Wmat,
                                 df_dt,
                                 K,
                                 par_constraint,
                                 parallel_optim,
                                 LBFGSB.maxit,
                                 LBFGSB.factr) {
  ### run optim in parallel (?)
  optim2 = ifelse(parallel_optim$run, optimParallel::optimParallel, stats::optim)

  ### constants needed for optimization
  theta.names.bik = gen_theta_names(K = K, type = "2x2", purpose = "bik")
  constraint_index = names(par_constraint)[names(par_constraint) %in% names(step2_inits)]

  ### STEP 2  ~~>  optimization
  if (parallel_optim$run) {
    step2_out = optim2(
      par = step2_inits,
      function(par) {
        ## constrain parameters
        par[constraint_index] = par_constraint
        ## compute log(P(Y))
        y = bik_all_2x2(
          theta = par,
          data = df,
          Xmat = df_Xmat,
          Wmat = df_Wmat,
          dt = df_dt,
          K = K,
          P.rs = FALSE,
          theta.names = theta.names.bik
        )
        y = impute_bik(x = y)
        ## `bi = y$bi1 + y$bi2 + ... + y$biK`
        bi = Reduce(`+`, y)
        ## return
        sum(log(bi))
      },
      method = "L-BFGS-B",
      control = list(fnscale = -length(df_dt), trace = 1, maxit = LBFGSB.maxit, factr = LBFGSB.factr),
      parallel = list(cl = parallel_optim$cl, forward = FALSE, loginfo = FALSE)
    )
  }
  if (!parallel_optim$run) {
    step2_out = optim2(
      par = step2_inits,
      function(par) {
        ## constrain parameters
        par[constraint_index] = par_constraint
        ## compute log(P(Y))
        y = bik_all_2x2(
          theta = par,
          data = df,
          Xmat = df_Xmat,
          Wmat = df_Wmat,
          dt = df_dt,
          K = K,
          P.rs = FALSE,
          theta.names = theta.names.bik
        )
        y = impute_bik(x = y)
        ## `bi = y$bi1 + y$bi2 + ... + y$biK`
        bi = Reduce(`+`, y)
        ## return
        sum(log(bi))
      },
      method = "L-BFGS-B",
      control = list(fnscale = -length(df_dt), trace = 1, maxit = LBFGSB.maxit, factr = LBFGSB.factr)
    )
  }
  step2_out = step2_out$par

  ### STEP 2  ~~> exits
  out = step2_out
  class(out) = append("lctmc_2x2.inits02", class(out))
  return(out)
}

#' @rdname gen_inits02_lctmc
gen_inits02_lctmc_3x3 = function(step2_inits,
                                 df,
                                 df_Xmat,
                                 df_Wmat,
                                 df_dt,
                                 K,
                                 par_constraint,
                                 parallel_optim,
                                 LBFGSB.maxit,
                                 LBFGSB.factr) {
  ### run optim in parallel (?),

  optim2 = ifelse(parallel_optim$run, optimParallel::optimParallel, stats::optim)

  ### constants needed for optimization
  theta.names.bik = gen_theta_names(K = K, type = "3x3", purpose = "bik")
  constraint_index = names(par_constraint)[names(par_constraint) %in% names(step2_inits)]

  ### STEP 2  ~~>  optimization
  if (parallel_optim$run) {
    step2_out = optim2(
      par = step2_inits,
      function(par) {
        ## constrain parameters
        par[constraint_index] = par_constraint
        ## compute log(P(Y))
        y = bik_all_3x3(
          theta = par,
          data = df,
          Xmat = df_Xmat,
          Wmat = df_Wmat,
          dt = df_dt,
          K = K,
          P.rs = FALSE,
          theta.names = theta.names.bik
        )
        y = impute_bik(x = y)
        ## `bi = y$bi1 + y$bi2 + ... + y$biK`
        bi = Reduce(`+`, y)
        ## return
        sum(log(bi))
      },
      method = "L-BFGS-B",
      control = list(fnscale = -length(df_dt), trace = 1, maxit = LBFGSB.maxit, factr = LBFGSB.factr),
      parallel = list(cl = parallel_optim$cl, forward = FALSE, loginfo = FALSE)
    )
  }
  if (!parallel_optim$run) {
    step2_out = optim2(
      par = step2_inits,
      function(par) {
        ## constrain parameters
        par[constraint_index] = par_constraint
        ## compute log(P(Y))
        y = bik_all_3x3(
          theta = par,
          data = df,
          Xmat = df_Xmat,
          Wmat = df_Wmat,
          dt = df_dt,
          K = K,
          P.rs = FALSE,
          theta.names = theta.names.bik
        )
        y = impute_bik(x = y)
        ## `bi = y$bi1 + y$bi2 + ... + y$biK`
        bi = Reduce(`+`, y)
        ## return
        sum(log(bi))
      },
      method = "L-BFGS-B",
      control = list(fnscale = -length(df_dt), trace = 1, maxit = LBFGSB.maxit, factr = LBFGSB.factr)
    )
  }
  step2_out = step2_out$par

  ### STEP 2  ~~> exits
  out = step2_out
  class(out) = append("lctmc_3x3.inits02", class(out))
  return(out)
}
