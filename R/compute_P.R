#' @title Computes the transition probability matrix, \eqn{P}
#'
#' @description Uses the parameters estimated from the latent class CTMC model (`lctmc_2x2()` or `lctmc_3x3()`),
#' and computes the transition **probability** matrices for some user specified covariate values
#'
#' @name compute_P
#'
#' @param m a 'lctmc_2x2' or 'lctmc_3x3' object obtained from the `lctmc` functions
#' @param x1 a numeric scalar. It is the covariate corresponding to the parameters `beta1.12_k`, `beta1.21_k`, and `beta1.23_k`
#' @param x2 a numeric scalar. It is the covariate corresponding to the parameters `beta2.12_k`, `beta2.21_k`, and `beta2.23_k`
#' @param dt a numeric scalar. This value indicates the time difference between observations.
#'
#' @return A list object of length equal to the number of latent classes specified when fitting the LCTMC model.
#' (See the documentation for `K` at [lctmc_2x2()]). \cr
#' Each element of the list is the transition probability matrix for the corresponding latent class.
#'
#' @note The transition **probability** matrix is directly related to the transition **rate** matrix, \eqn{Q}. \cr
#' The relation is the following:
#' \deqn{
#'   P = \exp(Q \cdot t)
#' }
#' where \eqn{Q} is obtained from the `compute_Q()` function, and \eqn{t} is equivalent to the input argument `dt`.
#' The exponential of a square matrix is defined as the Taylor expansion,
#' \deqn{
#'   \exp(Qt) = \sum_{n=0}^{\infty}\frac{1}{n!}(Qt)^{n}
#' }
#'
#' @export
#'
#' @seealso [compute_Q()]; [get_P_2x2()]; [lctmc_2x2()];
#'
#' @example inst/examples/ex_compute_statistics.R
compute_P = function(m, x1, x2, dt) {
  UseMethod("compute_P")
}

#' @rdname compute_P
#' @export
compute_P.lctmc_2x2 = function(m, x1, x2, dt) {
  ## checks
  if (!all(is.numeric(c(x1, x2, dt)))) {
    stop("`x1`, `x2`, `dt` should all be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1) | length(dt) != 1) {
    stop("`x1`, `x2`, `dt` should all be length 1")
  }
  if (!(dt >= 0)) {
    stop("`dt` cannot be a negative number!")
  }

  ## calls the compute Q function to generate transition rates
  fn_call = match.call()
  fn_call[[1]] = as.name("compute_Q")
  fn_call$dt = NULL
  Qlist = eval(fn_call)

  ## loop through each latent class
  Plist = list()
  for (k in seq_along(Qlist)) {
    # get current class's transition rates
    Qmat.k = Qlist[[k]]
    q_12 = Qmat.k[1, 2]
    q_21 = Qmat.k[2, 1]

    # get current class's probability matrix
    Pmat.k = get_P_2x2(q12 = q_12, q21 = q_21, dt = dt)
    Pmat.k = matrix(c(Pmat.k$P11, Pmat.k$P12, Pmat.k$P21, Pmat.k$P22),
                    byrow = TRUE, nrow = 2, ncol = 2)
    colnames(Pmat.k) = rownames(Pmat.k) = paste("state_", 1:2, sep = "")

    # append current class
    Plist[[names(Qlist)[k]]] = Pmat.k
  }

  ## return
  return(Plist)
}

#' @rdname compute_P
#' @export
compute_P.lctmc_3x3 = function(m, x1, x2, dt) {
  ## checks
  if (!all(is.numeric(c(x1, x2, dt)))) {
    stop("`x1`, `x2`, `dt` should all be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1) | length(dt) != 1) {
    stop("`x1`, `x2`, `dt` should all be length 1")
  }
  if (!(dt >= 0)) {
    stop("`dt` cannot be a negative number!")
  }

  ## calls the compute Q function to generate transition rates
  fn_call = match.call()
  fn_call[[1]] = as.name("compute_Q")
  fn_call$dt = NULL
  Qlist = eval(fn_call)

  ## loop through each latent class
  Plist = list()
  for (k in seq_along(Qlist)) {
    # get current class's transition rates
    Qmat.k = Qlist[[k]]
    q_12 = Qmat.k[1, 2]
    q_21 = Qmat.k[2, 1]
    q_23 = Qmat.k[2, 3]

    # get current class's probability matrix
    Pmat.k = get_P_3x3(q12 = q_12, q21 = q_21, q23 = q_23, dt = dt)
    Pmat.k = matrix(c(Pmat.k$P11, Pmat.k$P12, Pmat.k$P13,
                      Pmat.k$P21, Pmat.k$P22, Pmat.k$P23,
                      0, 0, 1),
                    byrow = TRUE, nrow = 3, ncol = 3)
    colnames(Pmat.k) = rownames(Pmat.k) = paste("state_", 1:3, sep = "")

    # append current class
    Plist[[names(Qlist)[k]]] = Pmat.k
  }

  ## return
  return(Plist)
}
