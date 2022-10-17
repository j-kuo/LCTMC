#' @title Re-order the MLE of estimated parameters
#'
#' @description This function re-orders the estimated MLE from fitting the latent CTMC model. The latent class model does not have a unique solution
#' because ordering the latent classes are interchangeable as long as the respective transition rate matrices are consistent.
#' To see this take a simple example of a three-category random variable with probability mass:
#' \deqn{
#'   P(X=1) = 1/2\\
#'   P(X=2) = 1/3\\
#'   P(X=3) = 1/6
#' }
#' The solutions could be parametrized as either (1/2, 1/3), (1/2, 1/6), or (1/3, 1/6)
#' as long as the category stays consistent (i.e., X=1 is attached to probability mass of 1/2, X=2 with probability 1/3, etc.). \cr
#' Therefore, this function's purpose is to re-order the classes such that each classes are sorted by the `alpha0` parameter in descending order
#' (i.e., `alpha0.1 > alpha0.2 > ...`)
#'
#' @name align_MLE
#'
#' @param true a list object obtained from `gen_true_param()` of the 'LCTMC.simulate' package.
#' When this argument is specified we can compare the true parameter values vs. the estimated MLE. (only used for simulation study purposes). \cr
#' Note user can create this object without the 'LCTMC.simulate', as long as the list is structured correctly.
#' @param mle a named numeric vector for the MLE of model parameters.
#' @param K an integer scalar. Use this variable to tell the function how many latent classes there should be. \cr
#' Note that the number of latent classes will affect the number of parameters in the model, thus the argument `theta.name` should be in sync with `K`
#'
#' @return A data frame object with the following columns:
#' \itemize{
#'   \item **names** - the name of the corresponding parameters.
#'   \item **true_theta** the true parameter value for the underlying data generation process.
#'   \item **mle_theta** - the MLE of the corresponding parameters.
#' }
#'
#' @note Note that the argument `true` can be left unspecified so that we are only re-arranging the vector of MLE.
#'
#' @seealso [lctmc_2x2()], [[lctmc_3x3()]]
#'
#' @example inst/examples/ex_align_MLE.R
NULL

#' @rdname align_MLE
#' @export
align_MLE_2x2 = function(true, mle = c(), K = 3L) {
  ## check MLE
  if (missing(mle)) {
    stop("must specify a named numeric vector for the `mle` argument")
  }
  if (!is.numeric(mle) || is.null(names(mle))) {
    stop("`mle` must be a named numeric vector where the names are model parameter names")
  }

  ## MLEs' pi
  mle.alpha = c(K, 0,0,0)
  for (k in (K-1):1) {
    mle.alpha_k = c(mle[grepl(pattern = paste("^alpha0\\.", k, "$", sep = ""), names(mle))], # alpha0 ~ intercept
                    mle[grepl(pattern = paste("^alpha1\\.", k, "$", sep = ""), names(mle))], # alpha1 ~ 1st covariate
                    mle[grepl(pattern = paste("^alpha2\\.", k, "$", sep = ""), names(mle))]) # alpha2 ~ 2nd covariate
    mle.alpha = rbind(c(k, mle.alpha_k), mle.alpha, deparse.level = 0)
  }
  colnames(mle.alpha) = c("class", "alpha0", "alpha1", "alpha2")

  ## MLEs' beta
  mle.beta = mle[!grepl(pattern = "^alpha\\d\\.\\d$", names(mle))]
  mle.beta_names = names(mle.beta)
  mle.beta_names_index = list()
  for (k in 1:K) {
    beta.12_k.index = which(grepl(pattern = paste("12_", k, "$", sep = ""), mle.beta_names))
    beta.21_k.index = which(grepl(pattern = paste("21_", k, "$", sep = ""), mle.beta_names))
    mle.beta_names_index[[as.character(k)]] = c(beta.12_k.index, beta.21_k.index)
  }
  mle.beta_names = mle.beta_names[unlist(mle.beta_names_index)]

  ## re-arrange alphas ... ACCORDING to magnitude of alpha0
  mle.alpha.ref_row = mle.alpha[which.min(mle.alpha[, 'alpha0']), 2:4]
  for (i in seq_along(mle.alpha[, 1])) {
    mle.alpha[i, 2:4] = mle.alpha[i, 2:4] - mle.alpha.ref_row
  }
  mle.alpha = mle.alpha[order(-mle.alpha[, 'alpha0']), ]
  mle.orders = as.character(mle.alpha[, 'class'])

  ## change names accordingly
  mle.beta_new = c()
  for (i in seq_along(mle.orders)) {
    class = mle.orders[i]
    index = mle.beta_names_index[[class]]
    mle.beta_new = c(mle.beta_new, mle.beta[index])
  }

  ## re-parametrize alphas (create parameter names)
  mle.alpha_names = expand.grid(c("alpha0", "alpha1", "alpha2"), 1:(K-1))
  mle.alpha_names = paste(mle.alpha_names$Var1, mle.alpha_names$Var2, sep = ".")

  ## re-parametrize alphas (re-ordering)
  mle.alpha_new = c()
  for (k in 1:(K-1)) {
    mle.alpha_new = append(mle.alpha_new, mle.alpha[k, 2:4])
  }
  mle.alpha_new = as.numeric(mle.alpha_new)
  names(mle.alpha_new) = mle.alpha_names

  ## vector of parameter names
  names.mle = c(mle.alpha_names, mle.beta_names)

  ## everything True parameter, if it's specified
  if (!missing("true")) {
    if (is.list(true)) {
      # True parameters' pi
      true.alpha = matrix(c(1, true$pi$pi.Z1$alpha0, true$pi$pi.Z1$alpha1, true$pi$pi.Z1$alpha2,
                            2, true$pi$pi.Z2$alpha0, true$pi$pi.Z2$alpha1, true$pi$pi.Z2$alpha2,
                            3, 0, 0, 0),
                          byrow = TRUE, nrow = 3)
      colnames(true.alpha) = c("class", "alpha0", "alpha1", "alpha2")

      # true parameters' beta
      true.beta = c(beta0.12_1 = log(true$r0$r0.Z1$q12), beta1.12_1 = true$beta$beta.Z1$q12[1], beta2.12_1 = true$beta$beta.Z1$q12[2],
                    beta0.21_1 = log(true$r0$r0.Z1$q21), beta1.21_1 = true$beta$beta.Z1$q21[1], beta2.21_1 = true$beta$beta.Z1$q21[2],
                    beta0.12_2 = log(true$r0$r0.Z2$q12), beta1.12_2 = true$beta$beta.Z2$q12[1], beta2.12_2 = true$beta$beta.Z2$q12[2],
                    beta0.21_2 = log(true$r0$r0.Z2$q21), beta1.21_2 = true$beta$beta.Z2$q21[1], beta2.21_2 = true$beta$beta.Z2$q21[2],
                    beta0.12_3 = log(true$r0$r0.Z3$q12), beta1.12_3 = true$beta$beta.Z3$q12[1], beta2.12_3 = true$beta$beta.Z3$q12[2],
                    beta0.21_3 = log(true$r0$r0.Z3$q21), beta1.21_3 = true$beta$beta.Z3$q21[1], beta2.21_3 = true$beta$beta.Z3$q21[2])
      true.beta_names = names(true.beta)
      true.beta_names_index = list()
      for (k in 1:3) {
        beta.12_k.index = which(grepl(pattern = paste("12_", k, "$", sep = ""), true.beta_names))
        beta.21_k.index = which(grepl(pattern = paste("21_", k, "$", sep = ""), true.beta_names))
        true.beta_names_index[[as.character(k)]] = c(beta.12_k.index, beta.21_k.index)
      }
      true.beta_names = true.beta_names[unlist(true.beta_names_index)]

      # checks
      if (!all(dim(true.alpha) == dim(mle.alpha))) {
        stop("Mis-matching parameters: `true.alpha` and `mle.alpha` have different dimensions")
      }
      if (length(true.beta) != length(mle.beta)) {
        stop("Mis-matching parameters: `true.beta` and `mle.beta` have mis-matching length")
      }

      # re-arrange alphas ... ACCORDING to magnitude of alpha0
      true.alpha.ref_row = true.alpha[which.min(true.alpha[, 'alpha0']), 2:4]
      for (i in seq_along(true.alpha[, 1])) {
        true.alpha[i, 2:4] = true.alpha[i, 2:4] - true.alpha.ref_row
      }
      true.alpha = true.alpha[order(-true.alpha[, 'alpha0']), ]
      true.orders = as.character(true.alpha[, 'class'])

      # change names accordingly
      true.beta_new = c()
      for (i in seq_along(true.orders)) {
        class = true.orders[i]
        index = true.beta_names_index[[class]]
        true.beta_new = c(true.beta_new, true.beta[index])
      }

      # re-parametrize pi (create param names)
      true.alpha_names = expand.grid(c("alpha0", "alpha1", "alpha2"), 1:(K-1))
      true.alpha_names = paste(true.alpha_names$Var1, true.alpha_names$Var2, sep = ".")

      # re-parametrize pi (re-ordering)
      true.alpha_new = c()
      for (k in 1:(K-1)) {
        true.alpha_new = append(true.alpha_new, true.alpha[k, 2:4])
      }
      true.alpha_new = as.numeric(true.alpha_new)
      names(true.alpha_new) = true.alpha_names

      # final checks ---> matching parameter names
      names.true = c(true.alpha_names, true.beta_names)
      if (!all(names.true == names.mle)) {
        stop("`true` argument was specified, but somehow mis-matching names were produced")
      }
    } else {
      stop("`true` must be a list object as created by `LCTMC.simulate::gen_true_param()`")
    }
  }else{
    true.alpha_names = rep(NA, length(mle.alpha_names))
    true.beta_names = rep(NA, length(mle.beta_names))
    true.alpha_new = rep(NA, length(mle.alpha_new))
    true.beta_new = rep(NA, length(mle.beta_new))
  }

  ## return
  data.frame(
    names = names.mle,
    true_theta = as.numeric(c(true.alpha_new, true.beta_new)),
    mle_theta = as.numeric(c(mle.alpha_new, mle.beta_new))
  )
}

#' @rdname align_MLE
#' @export
align_MLE_3x3 = function(true, mle = numeric(), K = 3L) {
  ## check MLE
  if (missing(mle)) {
    stop("must specify a named numeric vector for the `mle` argument")
  }
  if (!is.numeric(mle) || is.null(names(mle))) {
    stop("`mle` must be a named numeric vector where the names are model parameter names")
  }

  ## MLEs' pi
  mle.alpha = c(K, 0, 0, 0)
  for (k in (K-1):1) {
    mle.alpha_k = c(mle[grepl(pattern = paste("^alpha0\\.", k, "$", sep = ""), names(mle))], # alpha0 ~ intercept
                    mle[grepl(pattern = paste("^alpha1\\.", k, "$", sep = ""), names(mle))], # alpha1 ~ 1st covariate
                    mle[grepl(pattern = paste("^alpha2\\.", k, "$", sep = ""), names(mle))]) # alpha2 ~ 2nd covariate
    mle.alpha = rbind(c(k, mle.alpha_k), mle.alpha, deparse.level = 0)
  }
  colnames(mle.alpha) = c("class", "alpha0", "alpha1", "alpha2")

  ## MLEs' beta
  mle.beta = mle[!grepl(pattern = "^alpha\\d\\.\\d$", names(mle))]
  mle.beta_names = names(mle.beta)
  mle.beta_names_index = list()
  for (k in 1:K) {
    beta.12_k.index = which(grepl(pattern = paste("12_", k, "$", sep = ""), mle.beta_names))
    beta.21_k.index = which(grepl(pattern = paste("21_", k, "$", sep = ""), mle.beta_names))
    beta.23_k.index = which(grepl(pattern = paste("23_", k, "$", sep = ""), mle.beta_names))
    mle.beta_names_index[[as.character(k)]] = c(beta.12_k.index, beta.21_k.index, beta.23_k.index)
  }
  mle.beta_names = mle.beta_names[unlist(mle.beta_names_index)]

  ## re-arrange alphas ... ACCORDING to magnitude of alpha0
  mle.alpha.ref_row = mle.alpha[which.min(mle.alpha[, 'alpha0']), 2:4]
  for (i in seq_along(mle.alpha[, 1])) {
    mle.alpha[i, 2:4] = mle.alpha[i, 2:4] - mle.alpha.ref_row
  }
  mle.alpha = mle.alpha[order(-mle.alpha[, 'alpha0']), ]
  mle.orders = as.character(mle.alpha[, 'class'])

  ## change names accordingly
  mle.beta_new = c()
  for (i in seq_along(mle.orders)) {
    class = mle.orders[i]
    index = mle.beta_names_index[[class]]
    mle.beta_new = c(mle.beta_new, mle.beta[index])
  }

  ## re-parametrize alphas (create parameter names)
  mle.alpha_names = expand.grid(c("alpha0", "alpha1", "alpha2"), 1:(K-1))
  mle.alpha_names = paste(mle.alpha_names$Var1, mle.alpha_names$Var2, sep = ".")

  ## re-parametrize alphas (re-ordering)
  mle.alpha_new = c()
  for (k in 1:(K-1)) {
    mle.alpha_new = append(mle.alpha_new, mle.alpha[k, 2:4])
  }
  mle.alpha_new = as.numeric(mle.alpha_new)
  names(mle.alpha_new) = mle.alpha_names

  ## vector of parameter names
  names.mle = c(mle.alpha_names, mle.beta_names)

  ## everything True parameter, if it's specified
  if (!missing("true")) {
    if (is.list(true)) {
      # True parameters' pi
      true.alpha = matrix(c(1, true$pi$pi.Z1$alpha0, true$pi$pi.Z1$alpha1, true$pi$pi.Z1$alpha2,
                            2, true$pi$pi.Z2$alpha0, true$pi$pi.Z2$alpha1, true$pi$pi.Z2$alpha2,
                            3, 0, 0, 0),
                          byrow = TRUE, nrow = 3)
      colnames(true.alpha) = c("class", "alpha0", "alpha1", "alpha2")

      # True parameters' beta
      true.beta = c(beta0.12_1 = log(true$r0$r0.Z1$q12), beta1.12_1 = true$beta$beta.Z1$q12[1], beta2.12_1 = true$beta$beta.Z1$q12[2],
                    beta0.21_1 = log(true$r0$r0.Z1$q21), beta1.21_1 = true$beta$beta.Z1$q21[1], beta2.21_1 = true$beta$beta.Z1$q21[2],
                    beta0.23_1 = log(true$r0$r0.Z1$q23), beta1.23_1 = true$beta$beta.Z1$q23[1], beta2.23_1 = true$beta$beta.Z1$q23[2],
                    beta0.12_2 = log(true$r0$r0.Z2$q12), beta1.12_2 = true$beta$beta.Z2$q12[1], beta2.12_2 = true$beta$beta.Z2$q12[2],
                    beta0.21_2 = log(true$r0$r0.Z2$q21), beta1.21_2 = true$beta$beta.Z2$q21[1], beta2.21_2 = true$beta$beta.Z2$q21[2],
                    beta0.23_2 = log(true$r0$r0.Z2$q23), beta1.23_2 = true$beta$beta.Z2$q23[1], beta2.23_2 = true$beta$beta.Z2$q23[2],
                    beta0.12_3 = log(true$r0$r0.Z3$q12), beta1.12_3 = true$beta$beta.Z3$q12[1], beta2.12_3 = true$beta$beta.Z3$q12[2],
                    beta0.21_3 = log(true$r0$r0.Z3$q21), beta1.21_3 = true$beta$beta.Z3$q21[1], beta2.21_3 = true$beta$beta.Z3$q21[2],
                    beta0.23_3 = log(true$r0$r0.Z3$q23), beta1.23_3 = true$beta$beta.Z3$q23[1], beta2.23_3 = true$beta$beta.Z3$q23[2])
      true.beta_names = names(true.beta)
      true.beta_names_index = list()
      for (k in 1:3) {
        beta.12_k.index = which(grepl(pattern = paste("12_", k, "$", sep = ""), true.beta_names))
        beta.21_k.index = which(grepl(pattern = paste("21_", k, "$", sep = ""), true.beta_names))
        beta.23_k.index = which(grepl(pattern = paste("23_", k, "$", sep = ""), true.beta_names))
        true.beta_names_index[[as.character(k)]] = c(beta.12_k.index, beta.21_k.index, beta.23_k.index)
      }
      true.beta_names = true.beta_names[unlist(true.beta_names_index)]

      # checks
      if (!all(dim(true.alpha) == dim(mle.alpha))) {
        stop("Mis-matching parameters: `true.alpha` and `mle.alpha` have different dimensions")
      }
      if (length(true.beta) != length(mle.beta)) {
        stop("Mis-matching parameters: `true.beta` and `mle.beta` have mis-matching length")
      }

      # re-arrange alphas ... ACCORDING to magnitude of alpha0
      true.alpha.ref_row = true.alpha[which.min(true.alpha[, 'alpha0']), 2:4]
      for (i in seq_along(true.alpha[, 1])) {
        true.alpha[i, 2:4] = true.alpha[i, 2:4] - true.alpha.ref_row
      }
      true.alpha = true.alpha[order(-true.alpha[, 'alpha0']), ]
      true.orders = as.character(true.alpha[, 'class'])

      # change names accordingly
      true.beta_new = c()
      for (i in seq_along(true.orders)) {
        class = true.orders[i]
        index = true.beta_names_index[[class]]
        true.beta_new = c(true.beta_new, true.beta[index])
      }

      # re-parametrize pi (create param names)
      true.alpha_names = expand.grid(c("alpha0", "alpha1", "alpha2"), 1:(K-1))
      true.alpha_names = paste(true.alpha_names$Var1, true.alpha_names$Var2, sep = ".")

      # re-parametrize pi (re-ordering)
      true.alpha_new = c()
      for (k in 1:(K-1)) {
        true.alpha_new = append(true.alpha_new, true.alpha[k, 2:4])
      }
      true.alpha_new = as.numeric(true.alpha_new)
      names(true.alpha_new) = true.alpha_names

      # final checks ---> matching parameter names
      names.true = c(true.alpha_names, true.beta_names)
      if (!all(names.true == names.mle)) {
        stop("`true` argument was specified, but somehow mis-matching names were produced")
      }
    } else {
      stop("`true` must be a list object as created by `LCTMC.simulate::gen_true_param()`")
    }
  } else {
    true.alpha_names = rep(NA, length(mle.alpha_names))
    true.beta_names = rep(NA, length(mle.beta_names))
    true.alpha_new = rep(NA, length(mle.alpha_new))
    true.beta_new = rep(NA, length(mle.beta_new))
  }

  ## return
  data.frame(
    names = names.mle,
    true_theta = as.numeric(c(true.alpha_new, true.beta_new)),
    mle_theta = as.numeric(c(mle.alpha_new, mle.beta_new))
  )
}
