#' @title Computes the transition rate matrix, \eqn{Q}
#'
#' @description Uses the parameters estimated from the latent class CTMC model (`lctmc_2x2()` or `lctmc_3x3()`),
#' and computes the transition **rate** matrices for some user specified covariate values
#'
#' @name compute_Q
#'
#' @param m a 'lctmc_2x2' or 'lctmc_3x3' object obtained from the `lctmc` functions
#' @param x1 a numeric scalar. It is the covariate corresponding to the parameters `beta1.12_k`, `beta1.21_k`, and `beta1.23_k`
#' @param x2 a numeric scalar. It is the covariate corresponding to the parameters `beta2.12_k`, `beta2.21_k`, and `beta2.23_k`
#'
#' @return A list object of length equal to the number of latent classes specified when fitting the LCTMC model.
#' (See the documentation for `K` at [lctmc_2x2()]). \cr
#' Each element of the list is the transition rate matrix for the corresponding latent class
#'
#' @note The transition rate matrix is also commonly referred to as the "generator" matrix for CTMC models.
#' It satisfies some mathematical properties:
#' \itemize{
#'   \item \eqn{q_{ij} \geq 0}
#'   \item \eqn{q_{ii} = -\sum_{j}^{M}{q_{ij}}}
#' }
#'
#' @export
#'
#' @seealso [lctmc_2x2()]; [lctmc_3x3()]
#'
#' @example inst/examples/ex_compute_Q.R
compute_Q = function(m, x1, x2) {
  UseMethod("compute_Q")
}


#' @rdname compute_Q
#' @export
compute_Q.lctmc_2x2 = function(m, x1, x2) {
  ## checks
  if (!all(is.numeric(c(x1, x2)))) {
    stop("`x1` and `x2` should both be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1)) {
    stop("`x1` and `x2` should both be length 1")
  }

  ## constants
  x0 = 1
  df_theta = m$SE$SE

  ## loop through each latent class
  Qlist = list()
  for (k in 1:m$K) {
    # subset to current class
    df_theta.k = df_theta[grepl(pattern = paste("_", k, "$", sep = ""), df_theta$names), ]

    # compute tranistion rates for current class
    q_12 = as.numeric(exp(c(x0, x1, x2) %*% df_theta.k$mle_theta[grepl(pattern = "12", df_theta.k$names)]))
    q_21 = as.numeric(exp(c(x0, x1, x2) %*% df_theta.k$mle_theta[grepl(pattern = "21", df_theta.k$names)]))

    # build matrix for current class
    Qmat.k = matrix(c(0, q_12, q_21, 0), byrow = TRUE, nrow = 2, ncol = 2)
    diag(Qmat.k) = -rowSums(Qmat.k)
    colnames(Qmat.k) = rownames(Qmat.k) = paste("state_", 1:2, sep = "")

    # append current class
    Qlist[[paste("k", k, sep = "")]] = Qmat.k
  }

  ## return
  return(Qlist)
}

#' @rdname compute_Q
#' @export
compute_Q.lctmc_3x3 = function(m, x1, x2) {
  ## checks
  if (!all(is.numeric(c(x1, x2)))) {
    stop("`x1` and `x2` should both be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1)) {
    stop("`x1` and `x2` should both be length 1")
  }

  ## constants
  x0 = 1
  df_theta = m$SE$SE

  ## loop through each latent class
  Qlist = list()
  for (k in 1:m$K) {
    # subset to current class
    df_theta.k = df_theta[grepl(pattern = paste("_", k, "$", sep = ""), df_theta$names), ]

    # compute tranistion rates for current class
    q_12 = as.numeric(exp(c(x0, x1, x2) %*% df_theta.k$mle_theta[grepl(pattern = "12", df_theta.k$names)]))
    q_21 = as.numeric(exp(c(x0, x1, x2) %*% df_theta.k$mle_theta[grepl(pattern = "21", df_theta.k$names)]))
    q_23 = as.numeric(exp(c(x0, x1, x2) %*% df_theta.k$mle_theta[grepl(pattern = "23", df_theta.k$names)]))

    # build matrix for current class
    Qmat.k = matrix(c(0, q_12, 0,
                      q_21, 0, q_23,
                      0, 0, 0), byrow = TRUE, nrow = 3, ncol = 3)
    diag(Qmat.k) = -rowSums(Qmat.k)
    colnames(Qmat.k) = rownames(Qmat.k) = paste("state_", 1:3, sep = "")

    # append current class
    Qlist[[paste("k", k, sep = "")]] = Qmat.k
  }

  ## return
  return(Qlist)
}
