#' @title Computes transition probabilities (2x2)
#'
#' @description Computes the transition probabilities for the general 2x2 CTMC model. \cr
#' The only transition rates in this model are \eqn{q_{12}, q_{21} > 0}
#'
#' @param q12 a numeric vector for the transition rate from stage 1 to 2
#' @param q21 a numeric vector for the transition rate from stage 2 to 1
#' @param dt a numeric vector for the time difference between observations
#'
#' @return a data.frame object with 4 elements where each element is a resulting transition probability. \cr
#' For the 2x2 case, the possible transitions are `1-1, 1-2, 2-1, 2-2` \cr
#' and the corresponding transition probabilities are: `P11, P12, P21, P22`.
#'
#' @note The input argument should be of the same length, or their lengths should be multiples of one another. \cr
#' In addition, this function does not need to depend on the number of latent classes because what we are computing is
#' the transition probability conditioned on the given latent class.
#'
#' @example inst/examples/ex_get_P_2x2.R

get_P_2x2 = function(q12 = c(), q21 = c(), dt = c()) {
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
    # P 1 --> {1,2}
    P11 = p11,
    P12 = p12,
    # P 2 --> {1,2}
    P21 = p21,
    P22 = p22
  )
}
