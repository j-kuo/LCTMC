#' @title Compare log likelihood between MLE and true parameters
#'
#' @description **SIMULATION STUDY ONLY** \cr
#' Input a fitted LCTMC model and provide some true parameter object (obtained from the 'LCTMC.simulate' package).
#' This function will then compute the log observed data log-likelihood evaluated at the MLE and the true parameter value.
#' If the MLE is truly a maximum likelihood estimate, then it should yield a larger log-likelihood value.
#'
#' @name test_global_optim
#'
#' @param m a 'lctmc_2x2' or 'lctmc_3x3' object obtained from the `lctmc` functions
#' @param true_param a list object housing the true parameter values. This object can be obtained from
#' the 'LCTMC.simulate' package via the `gen_true_param()` function.
#' @param tol a numeric scalar. Tolerance value when comparing the log-likelihood values.
#' The smaller this value, the more precise the comparison.
#' @param data a data.frame object. This should be the data set that the LCTMC model fitted on.
#'
#' @return A list object of length 3. The elements of this list are:
#' \describe{
#'   \item{global_optim}{a logical scalar. TRUE if the \eqn{\log(L_{mle}) - \log(L_{true}) \geq \epsilon},
#'   where \eqn{\epsilon} is the `tol` input argument.}
#'   \item{L_mle}{a numeric scalar. The observed log-likelihood value when evaluated at the MLE.}
#'   \item{L_true}{a numeric scalar. The observed log-likelihood value when evaluated at the true parameter values.}
#' }
#'
#' @note This function is specifically meant to be used during a simulation study. \cr
#' Currently, this function only supports when `K = 3` as the 'LCTMC.simulate' package only
#' simulates the case of 3 latent classes.
#'
#' @export
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]
test_global_optim = function(m, true_param, tol, data) {
  UseMethod("test_global_optim")
}

#' @rdname test_global_optim
#' @export
test_global_optim.lctmc_2x2 = function(m, true_param, tol, data) {
  ## generate theta names for the bik() function
  theta.names.bik = gen_theta_names(K = 3L, type = "2x2", purpose = "bik")

  ## various data objects (non-scaled) ~ this is because true parameter values are not scaled
  my_df.test = fmt_rowwise_trans(
    data = data,
    type = "2x2",
    X_names = m$X_names,
    W_names = m$W_names,
    scaling = 1,
    trace = FALSE
  )
  my_df.Xmat.test = my_df.test$Xmat
  my_df.Wmat.test = my_df.test$Wmat
  my_df.dt.test = my_df.test$dt
  my_df.test = my_df.test$df_trans

  ## get MLE as a named vector
  mle_tht = m$SE$SE$mle_theta
  names(mle_tht) = m$SE$SE$names

  ## get true parameter as a named vector
  df_true_tht = align_MLE_2x2(true = true_param, mle = mle_tht, K = 3)
  true_tht = df_true_tht$true_theta
  names(true_tht) = df_true_tht$names

  ## compute LPY at True params
  bik_at_true = bik_all_2x2(
    theta = true_tht,
    data = my_df.test,
    Xmat = my_df.Xmat.test,
    Wmat = my_df.Wmat.test,
    dt = my_df.dt.test,
    K = 3L,
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
    K = 3L,
    P.rs = FALSE,
    theta.names = theta.names.bik
  )
  LPY_at_mle = sum(log(Reduce(`+`, bik_at_mle)))

  ## result
  global_optim = (LPY_at_mle - LPY_at_true) >= tol

  ## output
  out = list(global_optim = global_optim, L_mle = LPY_at_mle, L_true = LPY_at_true)
  return(out)
}

#' @rdname test_global_optim
#' @export
test_global_optim.lctmc_3x3 = function(m, true_param, tol, data) {
  ## generate theta names for the bik() function
  theta.names.bik = gen_theta_names(K = 33L, type = "3x3", purpose = "bik")

  ## various data objects (non-scaled) ~ this is because true parameter values are not scaled
  my_df.test = fmt_rowwise_trans(
    data = data,
    type = "3x3",
    X_names = m$X_names,
    W_names = m$W_names,
    scaling = 1,
    trace = FALSE
  )
  my_df.Xmat.test = my_df.test$Xmat
  my_df.Wmat.test = my_df.test$Wmat
  my_df.dt.test = my_df.test$dt
  my_df.test = my_df.test$df_trans

  ## get MLE as a named vector
  mle_tht = m$SE$SE$mle_theta
  names(mle_tht) = m$SE$SE$names

  ## get true parameter as a named vector
  df_true_tht = align_MLE_3x3(true = true_param, mle = mle_tht, K = 3)
  true_tht = df_true_tht$true_theta
  names(true_tht) = df_true_tht$names

  ## compute LPY at True params
  bik_at_true = bik_all_3x3(
    theta = true_tht,
    data = my_df.test,
    Xmat = my_df.Xmat.test,
    Wmat = my_df.Wmat.test,
    dt = my_df.dt.test,
    K = 3L,
    P.rs = FALSE,
    theta.names = theta.names.bik
  )
  LPY_at_true = sum(log(Reduce(`+`, bik_at_true)))

  ## compute LPY at MLE params
  bik_at_mle = bik_all_3x3(
    theta = mle_tht,
    data = my_df.test,
    Xmat = my_df.Xmat.test,
    Wmat = my_df.Wmat.test,
    dt = my_df.dt.test,
    K = 3L,
    P.rs = FALSE,
    theta.names = theta.names.bik
  )
  LPY_at_mle = sum(log(Reduce(`+`, bik_at_mle)))

  ## result
  global_optim = (LPY_at_mle - LPY_at_true) >= tol

  ## output
  out = list(global_optim = global_optim, L_mle = LPY_at_mle, L_true = LPY_at_true)
  return(out)
}
