#' @title Computes transition probabilities
#'
#' @description Computes the transition probabilities for the respective CTMC model (2x2 or 3x3). \cr
#' The only transition rates in this model are \eqn{q_{12}, q_{21} > 0}
#'
#' @name get_P
#' @param q12 a numeric vector for the transition rate from stage 1 to 2
#' @param q21 a numeric vector for the transition rate from stage 2 to 1
#' @param q23 a numeric vector for the transition rate from stage 2 to 3 (only applicable to 3x3 models)
#' @param log_q12 a numeric vector for the **log** transition rate from stage 1 to 2
#' @param log_q21 a numeric vector for the **log** transition rate from stage 2 to 1
#' @param log_q23 a numeric vector for the **log** transition rate from stage 2 to 3 (only applicable to 3x3 models)
#' @param dt a vector of numeric values for the time interval between observations
#'
#' @return a data.frame object with each column being one of the transition probability in the respective model. \cr\cr
#' For the 2x2 case, the possible transitions are `1-1, 1-2, 2-1, 2-2`,
#' and the corresponding transition probabilities are: `P11, P12, P21, P22`. \cr
#' For the 3x3 case, the possible transitions are: `1-1, 1-2, 1-3, 2-1, 2-2, 2-3`
#' and the corresponding transition probabilities are: `P11, P12, P13, P21, P22, P23`.
#' In the case of exactly-observed data on **state 3**, the transition densities are `P13_exact, P23_exact`.
#'
#' @note The probabilities will be undefined when \eqn{q_{21}=0} and \eqn{q_{12} = q_{23} \geq 0}.
#' Thus, a naive solution is implemented via the `impute_bik()` function which defaults NA values to some small values as penalty. \cr\cr
#' As far as function specification goes, the input argument should be of the same length to ensure regular behavior.
#' Additionally, this function does not need to depend on the number of latent classes because what we are computing is
#' the transition probability conditioned on the given latent class. \cr\cr
#' At the moment, the log rates are not being used as part of the computation.
#'
#' @seealso [Li_2x2()]; [bik_all_2x2()]; [impute_bik()]
#'
#' @example inst/examples/ex_get_P.R
NULL

#' @rdname get_P
get_P_2x2 = function(q12 = numeric(),
                     q21 = numeric(),
                     q23 = NA,
                     log_q12 = numeric(),
                     log_q21 = numeric(),
                     log_q23 = NA,
                     dt = numeric()) {
  # constants for both
  K = q12 + q21
  K_exp = (1-exp(-K*dt)) / K

  # compute P11 and P12
  p12 = q12 * K_exp
  p11 = 1 - p12

  # compute P21 and P22
  p21 = q21 * K_exp
  p22 = 1 - p21


  # output as df
  data.frame(
    # from 1 --> {1,2}
    P11 = p11,
    P12 = p12,
    # from 2 --> {1,2}
    P21 = p21,
    P22 = p22
  )
}

#' @rdname get_P
get_P_3x3 = function(q12 = numeric(),
                     q21 = numeric(),
                     q23 = numeric(),
                     log_q12 = numeric(),
                     log_q21 = numeric(),
                     log_q23 = numeric(),
                     dt = numeric()) {
  # constants for both
  K = q12 + q21 + q23
  discr = sqrt(K^2 - 4*q12*q23)
  r1 = (-K + discr) / 2
  r2 = (-K - discr) / 2

  # transition from state 1
  A = r2 / discr
  B = -1 - A
  A_exp = A*exp(r1*dt)
  B_exp = B*exp(r2*dt)
  p12 = (r1*A_exp + r2*B_exp) / q23
  p13 = A_exp + B_exp + 1
  p11 = 1 - p12 - p13

  # transitions from state 2
  C = (q23 + r2) / discr
  D = -1 - C
  C_exp = C*exp(r1*dt)
  D_exp = D*exp(r2*dt)
  CD_exp = C_exp + D_exp
  p22 = (r1*C_exp + r2*D_exp) / q23
  p23 = CD_exp + 1
  p21 = -p22 - CD_exp


  # output as df
  data.frame(
    # from 1 --> {1,2,3}
    P11 = p11,
    P12 = p12,
    P13 = p13,
    P13_exact = p12*q23,
    # from 2 --> {1,2,3}
    P21 = p21,
    P22 = p22,
    P23 = p23,
    P23_exact = p22*q23
  )
}
