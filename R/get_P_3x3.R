#' @title Computes transition probabilities (3x3)
#'
#' @description Computes the transition probabilities for a special 3x3 CTMC model. \cr
#' Specifically, the case with \eqn{q_{13} = q_{31} = q_{32} = q_{33} = 0}
#'
#' @param q12 a numeric vector for the transition rate from stage 1 to 2
#' @param q21 a numeric vector for the transition rate from stage 2 to 1
#' @param q23 a numeric vector for the transition rate from stage 2 to 3
#' @param dt a numeric vector for the time difference between observations
#'
#' @return a data.frame object with 8 (6+2) elements where each element is a transition probability/density. \cr
#' For this special 3x3 case, the possible transitions are: `1-1, 1-2, 1-3, 2-1, 2-2, 2-3`
#' and the corresponding transition probabilities are: `P11, P12, P13, P21, P22, P23`. \cr
#' In the case of exactly-observed data on **state 3**, the transition densities are `P13_exact, P23_exact`.
#'
#' @note The input argument should be of the same length, or their lengths should be multiples of one another. \cr
#' In addition, this function does not need to depend on the number of latent classes because what we are computing is
#' the transition probability conditioned on the given latent class. \cr
#' Additionally, when the absorbing state (Y=3) is exactly observed, such as in the case of deaths, the likelihood should use transition densities instead of probabilities.
#'
#' @seealso [Li_3x3()], [bik_all_3x3()], [get_P_2x2()]
#'
#' @example inst/examples/ex_get_P_3x3.R

get_P_3x3 = function(q12 = c(), q21 = c(), q23 = c(), dt = c()) {
  # constants for both
  K = q12 + q21 + q23
  discr = sqrt(K^2 - 4*q12*q23)
  r1 = (-K + discr) / 2
  r2 = (-K - discr) / 2

  # transition from state 1
  A = r2 / discr
  B = -1 - A
  A_exp = A * exp(r1*dt)
  B_exp = B * exp(r2*dt)
  p12 = (r1*A_exp + r2*B_exp) / q23
  p13 = A_exp + B_exp + 1
  p11 = 1 - p12 - p13

  # transitions from state 2
  C = (q23 + r2) / discr
  D = -1 - C
  C_exp = C * exp(r1*dt)
  D_exp = D * exp(r2*dt)
  p21 = -(r1*C_exp + r2*D_exp) / q23 - C_exp - D_exp
  p23 = C_exp + D_exp + 1
  p22 = 1 - p21 - p23

  # output as df
  data.frame(
    # P 1 --> {1,2,3}
    P11 = p11,
    P12 = p12,
    P13 = p13,
    P13_exact = p12*q23,
    # P 2 --> {1,2,3}
    P21 = p21,
    P22 = p22,
    P23 = p23,
    P23_exact = p22*q23
  )
}
