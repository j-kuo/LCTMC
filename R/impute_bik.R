#' @title Imputes \eqn{b_{ik}}
#'
#' @description This function imputes the quantity \eqn{b_{ik}} when calculating the likelihood. It works hand-in-hand with the `bik_all_2x2()` and `bik_all_3x3()` functions. \cr
#' During numerical optimization, sometimes the quantity \eqn{b_{ik}} could be evaluated to be a non-positive value (even though it should be greater than or equal to 0).
#' This function is usually called whenever \eqn{b_{ik}} is to be imputed. \cr
#' It performs the imputation by finding the minimum value greater than 0, then multiplying it by penalizing factor, \eqn{\epsilon}.
#'
#' @param x a list object obtained from `bik_all_2x2()` or `bik_all_3x3()`
#' @param eps a numeric scalar which is used to determine the imputation value to replace negative or zero values. \cr
#' For example, with `eps = 0.001` and `bi1 = c(0.25, -1.0, 0.5)` will replace the `-1` value to `0.25*0.001 = 0.00025` because `0.25` is the minimum value excluding the `-1`
#' @param EPS a numeric scalar which is used to determine the imputation value in the case when the entire vector is negative or zero.
#'
#' @return a list that has the same structure as the input argument `x`, except any negative or zero values are imputed
#'
#' @note This function should not be called externally in most use cases. it is typically called within the EM algorithm and the initial value generation process.
#'
#' @example inst/examples/ex_impute_bik.R

impute_bik = function(x = list(), eps = 1e-3, EPS = 1e-24) {
  ### loop through each element of `x`
  for (i in seq_along(x)) {
    ## current vector in list `x`
    v = x[[i]]
    v[is.na(v)] = -1 # NA should be imputed

    ## impute
    if (sum(v > 0, na.rm = T) > 0) {
      # if at least some positive values are in `v` then impute with `min(bik)*eps`
      v[v <= 0] = min(v[v > 0]) * eps
      # if any still 0 impute with `min(bik)`. This happens if `min(bik)` is extremely small
      v[v == 0] = min(v[v > 0])
    } else {
      # if the entire vector is non-positive impute with `EPS`
      v = rep(EPS, length(v))
    }

    ## update
    x[[i]] = v
  }

  ### return
  x
}
