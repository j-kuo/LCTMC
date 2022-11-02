#' @title Prints the fitted coefficients of a "lctmc" object
#'
#' @description Prints estimated coefficient for a latent class CTMC model
#'
#' @param x a "lctmc" object obtained from either the [lctmc_2x2()] or the [lctmc_3x3()] function
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


#' @title Prints the LCTMC control object
#'
#' @description Prints the control elements of a "lctmc_control" object
#'
#' @param x a "lctmc_control" object obtained from either the [create_controls()]
#' @param ... NULL
#'
#' @return NULL
#'
#' @exportS3Method
#'
#' @example inst/examples/ex_print.R
print.lctmc_control = function(x, ...){
  categories = c("fmt_data", "init01", "init02", "EM", "SE", "rescale")

  for (ctg in categories) {
    if (!is.null(x[[ctg]])) {
      cat("$", ctg, "\n", sep = "")
      longest = max(nchar(names(x[[ctg]])))

      for (i in seq_along(x[[ctg]])) {
        arg = x[[ctg]][[i]]
        arg_name = names(x[[ctg]])[i]

        white_space = rep(" ", times = longest-nchar(arg_name)+1)
        if (length(arg) > 1) {
          min_arg = sprintf("%.3f", min(arg))
          max_arg = sprintf("%.3f", max(arg))
          arg = paste("[", min_arg, "-", max_arg, "] length=", length(arg), sep = "")
        }
        cat("  *", white_space, arg_name, ": ", paste(arg, collapse = ";"), "\n", sep = "")
      }
    } else {
      cat("$", ctg, "\n", sep = "")
      cat("  NULL", "\n", sep = "")
    }

    cat(rep("=-=", 16), "\n", sep = "")
  }

  cat(rep(" ", 13), "Model Type ~ ", x[["type"]], "\n", sep = "")
}
