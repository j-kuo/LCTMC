#' @title Computes the expected sojourn time
#'
#' @description Uses the parameters estimated from the latent class CTMC model (`lctmc_2x2()` or `lctmc_3x3()`),
#' and computes the expected sojourn time for some user specified covariate values
#'
#' @name compute_sojourn
#'
#' @param m a 'lctmc_2x2' or 'lctmc_3x3' object obtained from the `lctmc` functions
#' @param x1 a numeric scalar. It is the covariate corresponding to the parameters `beta1.12_k`, `beta1.21_k`, and `beta1.23_k`
#' @param x2 a numeric scalar. It is the covariate corresponding to the parameters `beta2.12_k`, `beta2.21_k`, and `beta2.23_k`
#'
#' @return A list object of length equal to the number of latent classes specified when fitting the LCTMC model.
#' (See the documentation for `K` at [lctmc_2x2()]). \cr
#' Each element of the list is a named numeric vector of the expected sojourn time for the respective latent class
#'
#' @note The sojourn time is waiting time before a transition happens.
#' In CTMC models the sojourn is assumed to follow an exponential distribution. \cr
#' In the latent class CTMC model, the distribution is only exponential when we condition on the latent variable, \eqn{Z}
#'
#' @export
#'
#' @seealso [compute_Q()]
#'
#' @example inst/examples/ex_compute_sojourn.R
compute_sojourn = function(m, x1, x2) {
  UseMethod("compute_sojourn")
}

#' @rdname compute_sojourn
#' @export
compute_sojourn.lctmc_2x2 = function(m, x1, x2) {
  ## checks
  if (!all(is.numeric(c(x1, x2)))) {
    stop("`x1`, `x2` should both be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1)) {
    stop("`x1`, `x2` should both be length 1")
  }

  ## calls the compute Q function to generate transition rates
  fn_call = match.call()
  fn_call[[1]] = as.name("compute_Q")
  Qlist = eval(fn_call)

  ## loop through each latent class
  sojoun_list = list()
  for (k in seq_along(Qlist)) {
    # get current class's transition rates
    Qmat.k = Qlist[[k]]
    sojourn.k = 1 / (-1 * diag(Qmat.k))

    # append current class
    sojoun_list[[names(Qlist)[k]]] = sojourn.k
  }

  ## return
  return(sojoun_list)
}

#' @rdname compute_sojourn
#' @export
compute_sojourn.lctmc_3x3 = function(m, x1, x2) {
  ## checks
  if (!all(is.numeric(c(x1, x2)))) {
    stop("`x1`, `x2` should all both numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1)) {
    stop("`x1`, `x2` should all both length 1")
  }

  ## calls the compute Q function to generate transition rates
  fn_call = match.call()
  fn_call[[1]] = as.name("compute_Q")
  Qlist = eval(fn_call)

  ## loop through each latent class
  sojoun_list = list()
  for (k in seq_along(Qlist)) {
    # get current class's transition rates
    Qmat.k = Qlist[[k]]
    sojourn.k = 1 / (-1 * diag(Qmat.k))

    # append current class
    sojoun_list[[names(Qlist)[k]]] = sojourn.k
  }

  ## return
  return(sojoun_list)
}
