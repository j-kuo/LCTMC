#' @title Individual likelihood of the CTMC process
#'
#' @description Computes the individual-level likelihood conditioned on the latent variable, \eqn{P(Y_{i} | Z_{i})}
#'
#' @name Li
#'
#' @param P A data frame object containing the transition probabilities obtained from the `get_P` functions.
#' @param data A data frame object containing the binary transition indicator variables named as `trans.1_1`, `trans.1_2`, etc.
#' This data frame object can be obtained from `fmt_rowwise_trans()`
#'
#' @return A vector of likelihood/probability (i.e., numeric) values. Each element is the likelihood of each individual indicated by the `id` variable in `data`.
#' Length of output should be equal to number of unique IDs.
#'
#' @note This function does not need to depend on the number of latent classes because what we are computing is the likelihood conditioned on the given latent class. \cr
#'
#' \deqn{
#'   L_{i} = \prod_{j}^{n_{i}} L_{ij}
#' }
#' where \eqn{L_{ij}} is the likelihood of the \eqn{j^{th}} transition of the person \eqn{i^{th}} person.
#'
#' @importFrom collapse fprod
#'
#' @seealso [get_P_2x2()]; [bik_all_2x2()], [fmt_rowwise_trans()]
#'
#' @example inst/examples/ex_Li.R
NULL

#' @rdname Li
Li_2x2 = function(P = data.frame(), data = data.frame()) {
  ### Lij given z
  Lij = data$trans.1_1*P$P11 + data$trans.1_2*P$P12 + data$trans.2_1*P$P21 + data$trans.2_2*P$P22

  ### Li given z
  collapse::fprod(Lij, data$id, use.g.names = FALSE)
}

#' @rdname Li
Li_3x3 = function(P = data.frame(), data = data.frame()) {
  ### Lij given z
  Lij = data$trans.1_1*P$P11 + data$trans.1_2*P$P12 + data$trans.1_3*P$P13_exact +
    data$trans.2_1*P$P21 + data$trans.2_2*P$P22 + data$trans.2_3*P$P23_exact

  ### Li given z
  collapse::fprod(Lij, data$id, use.g.names = FALSE)
}
