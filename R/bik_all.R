#' @title Computes the quantity \eqn{b_{ik}}
#'
#' @description Computes \eqn{b_{ik}} for \eqn{k \geq 2}. See **Note** section for more info
#'
#' @name bik_all
#'
#' @param theta a named vector of numeric values for the model parameters
#' @param data a data frame object containing the binary transition indicator variables, i.e., `trans.1_1`, `trans.1_2`, etc.
#' @param Xmat a matrix object containing covariates "X" which are the covariates for the CTMC model
#' @param Wmat a matrix object containing covariates "W" which are the covariates for the latent class component of the model
#' @param dt a vector of numeric values for the time interval between observations
#' @param K a integer scalar. Use this variable to tell the function how many latent classes there should be, The number of latent classes will affect the number of parameter in the model.
#' @param P.rs a logical scalar. If TRUE, then \eqn{P_{rs(k)}} is returned as a list object where each element of the list varies by `k`.
#' If FALSE, then \eqn{b_{ik}} is returned as a list object where each element of the list varies by `k`.
#' @param theta.names a list parameter names. It is nested by the value of `K`. This list object can be generated using `gen_theta_names(..., purpose = "bik")`. \cr
#' Because this function is often using within the `optim()` function, it is inefficient to generate the names every iteration via `paste()`, hence having this pre-generated will help with speed.
#'
#' @return A list object. Depending on input value for `P.rs`. \cr\cr
#' If `P.rs` is specified to be TRUE, then the transition probability from state `r` to state `s` within `dt` amount of time is return.
#' As such, each element of `P.rs` is the same length as the input argument `dt`. \cr
#' If `P.rs` is set to FALSE however, then the quantity \eqn{b_{ik}} is returned (see note section for more info)
#'
#' @note This function computes the transition rates internally using the input arguments (`theta`, `Xmat`, `Wmat`).
#' Furthermore, this functions combines the `get_P` and the `Li` functions to obtain the quantities \eqn{P_{rs(k)}} or \eqn{b_{ik}}. \cr\cr
#' \eqn{P_{rs(k)}} is the transition probability from state `r` to state `s` within `dt` amount of time, condition on the latent variable **Z**. If \eqn{\delta t = t_{2}-t_{1}} then,
#' \deqn{
#'   P_{rs(k)}(\delta t) = P(Y(t_2)=s | Y(t_1)=r \cap Z=k)
#' }
#' Similarly, \eqn{b_{ik}} is the following
#' \deqn{
#'   b_{ik} = \pi_{k} \cdot \prod_{j}{P_{rs(k)}(\delta t_{ij})}
#' }
#' where \eqn{r=y_{i(j-1)}} and \eqn{s=y_{ij}}
#'
#' @seealso [get_P_2x2()]; [Li_2x2()]; [fmt_rowwise_trans()], [gen_theta_names()]; [impute_bik()]
#'
#' @example inst/examples/ex_bik_all.R
NULL

#' @rdname bik_all
bik_all_2x2 = function(theta = c(),
                       data = data.frame(),
                       Xmat = matrix(),
                       Wmat = matrix(),
                       dt = c(),
                       K,
                       P.rs = FALSE,
                       theta.names) {
  ## allocate some variables spaces
  e.pi.list = list()
  P.rs.list = list()
  bik.list = list()
  pi_K = 0

  ## loop for each class
  for (k in 1:K) {
    # alphas (only up to K-1)
    if (k < K) {
      e.pi_k = exp(as.numeric(Wmat %*% theta[theta.names[[k]]$alpha]))
      e.pi.list[[k]] = e.pi_k
      pi_K = pi_K + e.pi_k
    } else {
      e.pi.list[[k]] = 1
    }

    # betas for transition 1 --> 2
    q.12_k = exp(as.numeric(Xmat %*% theta[theta.names[[k]]$q12]))

    # betas for transition 2 --> 1
    q.21_k = exp(as.numeric(Xmat %*% theta[theta.names[[k]]$q21]))

    # P.rs
    P.rs_k = get_P_2x2(dt = dt, q12 = q.12_k, q21 = q.21_k)
    P.rs.list[[k]] = P.rs_k
  }

  ## if only P.rs is needed then return it here
  if (P.rs) {
    # return
    names(P.rs.list) = paste("P.rs_", 1:K, sep = "")
    return(P.rs.list)
  } else {
    # get last class's probability mass
    pi_K = 1 / (1 + pi_K)

    # loop through each class to obtain bik
    for (k in 1:K) {
      pi = pi_K * e.pi.list[[k]]
      bik.list[[k]] = pi * Li_2x2(P = P.rs.list[[k]], data = data)
    }

    # return
    bik.list
  }
}

#' @rdname bik_all
bik_all_3x3 = function(theta = c(),
                       data = data.frame(),
                       Xmat = matrix(),
                       Wmat = matrix(),
                       dt = c(),
                       K,
                       P.rs = FALSE,
                       theta.names) {
  ## allocate some variables spaces
  e.pi.list = list()
  P.rs.list = list()
  bik.list = list()
  pi_K = 0

  ## loop for each class
  for (k in 1:K) {
    # alphas (only up to K-1)
    if (k < K) {
      e.pi_k = exp(as.numeric(Wmat %*% theta[theta.names[[k]]$alpha]))
      e.pi.list[[k]] = e.pi_k
      pi_K = pi_K + e.pi_k
    } else {
      e.pi.list[[k]] = 1
    }

    # betas for transition 1 --> 2
    q.12_k = exp(as.numeric(Xmat %*% theta[theta.names[[k]]$q12]))

    # betas for transition 2 --> 1
    q.21_k = exp(as.numeric(Xmat %*% theta[theta.names[[k]]$q21]))

    # betas for transition 2 --> 3
    q.23_k = exp(as.numeric(Xmat %*% theta[theta.names[[k]]$q23]))

    # P.rs
    P.rs_k = get_P_3x3(dt = dt, q12 = q.12_k, q21 = q.21_k, q23 = q.23_k)
    P.rs.list[[k]] = P.rs_k
  }

  ## if only P.rs is needed then return it here
  if (P.rs) {
    # return
    names(P.rs.list) = paste("P.rs_", 1:K, sep = "")
    return(P.rs.list)
  } else {
    # get last class's probability mass
    pi_K = 1 / (1 + pi_K)

    # loop through each class to obtain bik
    for (k in 1:K) {
      pi = pi_K * e.pi.list[[k]]
      bik.list[[k]] = pi * Li_3x3(P = P.rs.list[[k]], data = data)
    }

    # return
    bik.list
  }
}
