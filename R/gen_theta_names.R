#' @title Generates model parameter names
#'
#' @description This function provides a convenient way of generating the parameter names of the latent class CTMC model. \cr
#' This is useful during the EM algorithm where we perform ECM steps by elements of this function's output list object. \cr
#' Or it is useful whenever \eqn{b_{ik}} should be computed, as those functions reference the parameter names internally,
#' and it would be inefficient to create the names via `paste()` every time
#'
#' @param K an integer scalar. Used to determine the number of latent classes the model is fitting.
#' @param type a character scalar which is used to determine the type of latent class CTMC model, either `type = "2x2"` or `type = "3x3"`
#' @param purpose a character scalar to indicate whether the output list object should be formatted for the EM algorithm or for computing \eqn{b_{ik}}
#'
#' @return a list object containing the name of parameters
#'
#' @export
#'
#' @example inst/examples/ex_gen_theta_names.R

gen_theta_names = function(K = integer(), type = c("2x2", "3x3"), purpose = c("em", "bik")) {
  ### checks
  if (!is.integer(K) || K <= 1) {
    stop("`K` should be an integer value greater than 1")
  }
  if (length(type) != 1 || !(type %in% c("2x2", "3x3"))) {
    stop("`type` should be length 1 and it should be either '2x2' or '3x3'")
  }
  if (length(purpose) != 1 || !(purpose %in% c("em", "bik"))) {
    stop("`purpose` should be length 1 and it should be either 'em' or 'bik'")
  }

  ### 2x2 model
  if (type == "2x2") {
    ## type ('em')
    i = 1
    names.a = names.b = list()
    for (k in 1:K) {
      # alpha names
      if (k < K) {
        names.a[[k]] = paste(c('alpha0', 'alpha1', 'alpha2'), ".", k, sep = "")
      }
      # beta names
      for (t in c("12", "21")) {
        names.b[[i]] = paste(c("beta0", "beta1", "beta2"), ".", t, "_", k, sep = "")
        i = i + 1
      }
    }
    par_names = c(names.a, names.b)


    ## type ('bik)
    theta_names = list()
    for (k in 1:K) {
      if (k < K) {
        theta_names[[k]] = list(
          alpha = paste(c('alpha0', "alpha1", "alpha2"), ".", k, sep = ""),
          q12 = paste(c('beta0.12', "beta1.12", "beta2.12"), "_", k, sep = ""),
          q21 = paste(c('beta0.21', "beta1.21", "beta2.21"), "_", k, sep = "")
        )
      } else {
        theta_names[[k]] = list(
          q12 = paste(c('beta0.12', "beta1.12", "beta2.12"), "_", k, sep = ""),
          q21 = paste(c('beta0.21', "beta1.21", "beta2.21"), "_", k, sep = "")
        )
      }
    }
  }

  ### 3x3 model
  if (type == "3x3") {
    ## type (1)
    i = 1
    names.a = names.b = list()
    for (k in 1:K) {
      # alpha names
      if (k < K) {
        names.a[[k]] = paste(c('alpha0', 'alpha1', 'alpha2'), ".", k, sep = "")
      }
      # beta names
      for (t in c("12", "21", "23")) {
        names.b[[i]] = paste(c("beta0", "beta1", "beta2"), ".", t, "_", k, sep = "")
        i = i + 1
      }
    }
    par_names = c(names.a, names.b)

    ## type ('bik)
    theta_names = list()
    for (k in 1:K) {
      if (k < K) {
        theta_names[[k]] = list(
          alpha = paste(c('alpha0', "alpha1", "alpha2"), ".", k, sep = ""),
          q12 = paste(c('beta0.12', "beta1.12", "beta2.12"), "_", k, sep = ""),
          q21 = paste(c('beta0.21', "beta1.21", "beta2.21"), "_", k, sep = ""),
          q23 = paste(c('beta0.23', "beta1.23", "beta2.23"), "_", k, sep = "")
        )
      } else {
        theta_names[[k]] = list(
          q12 = paste(c('beta0.12', "beta1.12", "beta2.12"), "_", k, sep = ""),
          q21 = paste(c('beta0.21', "beta1.21", "beta2.21"), "_", k, sep = ""),
          q23 = paste(c('beta0.23', "beta1.23", "beta2.23"), "_", k, sep = "")
        )
      }
    }
  }

  ### output
  if (purpose == "em") {
    return(par_names)
  } else if (purpose == "bik") {
    return(theta_names)
  } else {
    return(NA)
  }
}
