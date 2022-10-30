#' @title Prints the fitted coefficients of a "lctmc" object
#'
#' @description Prints estimated coefficient for a latent class CTMC model
#'
#' @param x a "lctmc" object obtained from either the `lctmc_2x2()` or the `lctmc_3x3()` function
#' @param ... NULL
#'
#' @return NULL
#'
#' @exportS3Method
#'
#' @example inst/examples/ex_print.R

print.lctmc = function(x, ...){
  print(x$SE$SE)

  dot_length = nchar(nrow(x$SE$SE)) + 5*(1 + nchar("beta2.23_3")) + 7
  cat(rep(".", times = dot_length), "\n", sep = "")
  cat(rep(".", times = dot_length), "\n", sep = "")

  print(
    data.frame(
      type = ifelse(any(grepl(pattern = "2x2", class(x))), "2x2", "3x3"),
      K = x$K,
      n_pars = x$n_pars,
      hess_code = x$SE$hess_code,
      covariance_code = x$SE$covariance_code,
      run_time = paste(sprintf("%.2f", x$run_time), units(x$run_time), sep = " ")
    )
  )
}
