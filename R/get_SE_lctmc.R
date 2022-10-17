#' @title SE Approximation for the Latent CTMC model
#'
#' @description Uses the large sample approximation for to estimate the standard error of the MLE (see **Note** for more info). \cr
#' The result of the approximation is the variance-covariance matrix for the vector of the estimated MLE,
#' where the diagonal elements are the variance and off-diagonals are the covariance.
#'
#' @name get_SE_lctmc
#'
#' @param em a list object of with the custom class 'lctmc_2x2.mle' or 'lctmc_3x3.mle'. This type of object is the output of the `EM_lctmc` functions
#' @param df a data frame object containing the binary row-wise transition indicator variables
#' @param df.Xmat a matrix object with same number of rows as `df`. This matrix object should contain the covariates which affect the CTMC part of the model
#' @param df.Wmat a matrix object with number of rows equal to the unique number of individuals in `df`. \cr
#' This matrix object should contain the covariates which affect the latent class probability part of the model.
#' @param df.dt a numeric vector with length equal to the number of rows as `df`. This vector contains the time difference between observations.
#' @param par_constraint a named numeric vector to indicate which parameter is constrained. Set equal to NULL for unconstrained model. \cr
#' For example, `c(alpha1.1 = 0)` constraints the parameter 'alpha1.1' to be a constant 0. **NOTE:** Current version of the code will *only* work with constrains equal to 0.
#' @param K the number of categories the latent class variable has.
#' @param solve.tol a numeric scalar, typically a small decimal value. It is the tolerance for detecting linear dependencies in the hessian matrix. \cr
#' Defaults to `(.Machine$double.eps)^2` if not specified.
#' @param symmetric.tol a numeric scalar. Tolerane value for checking symmetric matrix. \cr Default is 5e-11
#' @param eigen0.tol a numeric scalar. Tolerance value for eigenvalues, any values smaller than this will be treated as 0. \cr Default is 1e-10
#'
#' @return a list object containing 3 elements:
#' \itemize{
#'   \item `SE` is a data frame object containing columns for: the MLE, the approximated SE, and the 95% confidence interval for the MLE
#'   \item `covariance_code` a numeric scalar that can take be one of three values. \cr
#'          0: when the covariance matrix is neither positive definite or semi positive definite \cr
#'          1: when the covariance matrix is semi positive definite \cr
#'          2: when the covariance matrix is positive definite
#'   \item `hess_code` similar to `covariance_code` \cr
#'          0: when the hessian is neither negative definite or semi negative definite (saddle point) \cr
#'          1: when the covariance matrix is semi positive definite (inconclusive result) \cr
#'          2: when the covariance matrix is positive definite (a local optimal point)
#'   \item `Covariance` a matrix object which is the estimated variance covariance matrix for the estimated parameters. \cr
#'          Note that the estimated SE is simply square root of the diagonal elements.
#' }
#'
#' @note This is step five out of six for fitting a latent class CTMC model (i.e., SE approximation via the hessian matrix). \cr\cr
#' This method relies on the large sample theory that
#' \deqn{
#'   \hat\theta \sim N(\theta, -H^{-1}(\hat\theta))
#' }
#' where, \eqn{H(\hat\theta)} is the hessian of the observed log likelihood function evaluated at the MLE. And \eqn{\hat\theta} is the estimated MLE.
#' Since this SE approximation relies on \eqn{N \rarr \infty}, user should be cautious when considering sample size and number of model parameters.
#'
#' @importFrom numDeriv hessian
#'
#' @example inst/examples/ex_get_SE_lctmc.R
NULL

#' @rdname get_SE_lctmc
get_SE_lctmc_2x2 = function(em,
                            df,
                            df.Xmat,
                            df.Wmat,
                            df.dt,
                            par_constraint,
                            K,
                            solve.tol = (.Machine$double.eps)^2,
                            symmetric.tol = 5e-11,
                            eigen0.tol = 1e-10) {
  ### perform checks
  if (!("lctmc_2x2.mle" %in% class(em))) {
    stop("`em` should be a custom class list object 'lctmc_2x2.mle' obtained from `EM_lctmc_2x2()`")
  }

  ### extract the EM run with best log(P(Y)) value
  best_index = sapply(em, function(x) x$LPY_value)
  best_index = max(which(best_index == max(best_index))) # in case of ties, get the last index

  ### align results with true parameter values
  df.theta = align_MLE_2x2(mle = em[[best_index]]$pars_value, K = K)
  df.theta = df.theta[!colnames(df.theta) %in% c("true_theta")]

  ### msg
  cat(" * best EM run occurred at step ", best_index, "/", length(em), "\n", sep = "")

  ### get vector of the MLE (without constrained elements)
  mle = df.theta$mle_theta[!(df.theta$names %in% names(par_constraint))]
  mle.names = df.theta$names[!(df.theta$names %in% names(par_constraint))]

  ### numerical derivative to get Hessian matrix
  theta.names.bik = gen_theta_names(K = K, type = "2x2", purpose = "bik")
  hess = numDeriv::hessian(
    func = function(x) {
      ## theta vector, adding back the constrained elements
      names(x) = mle.names
      x = append(x, par_constraint)
      ## compute log(P(Y))
      bik_all.mle = bik_all_2x2(
        theta = x,
        data = df,
        Xmat = df.Xmat,
        Wmat = df.Wmat,
        dt = df.dt,
        K = K,
        theta.names = theta.names.bik
      )
      bik_all.mle = impute_bik(bik_all.mle, eps = 1e-3, EPS = 1e-24)
      ## sum over all class
      bi = Reduce(`+`, bik_all.mle)
      ## return
      sum(log(bi))
    },
    x = mle,
    method.args = list(r = 4) # use 6 for more accuracy but longer run time
  )

  ### Hessian and covariance matrix
  colnames(hess) = rownames(hess) = mle.names
  cov_mat = -1 * solve(a = hess, tol = solve.tol)

  ### check for covariance matrix
  covariance_code = -1
  if (isSymmetric(cov_mat, tol = symmetric.tol)) {
    ## eigenvalues
    covariance.eigen = eigen(cov_mat, only.values = TRUE)$values
    covariance.eigen[abs(covariance.eigen) < eigen0.tol] = 0

    ## check definiteness
    if (all(covariance.eigen > 0)) {
      covariance_code = 2
      cat(" * the covariance matrix is positive definite ~ at least a local maxima is reached \n")
    } else if (all(covariance.eigen >= 0)) {
      covariance_code = 1
      cat(" * the covariance matrix is positive semi-definite ~ this is inconclusive \n")
    } else {
      covariance_code = 0
      cat(" * the covariance matrix is not positive (semi-)definite ~ this is a saddle point \n")
    }
  }

  ### check for hessian matrix
  hess_code = -1
  if (isSymmetric(hess, tol = symmetric.tol)) {
    ## eigenvalues
    hess.eigen = eigen(hess, only.values = TRUE)$values
    hess.eigen[abs(hess.eigen) < eigen0.tol] = 0

    ## check definiteness
    if (all(hess.eigen < 0)) {
      hess_code = 2
    } else if (all(hess.eigen <= 0)) {
      hess_code = 1
    } else {
      hess_code = 0
    }

    ## does the hessian matrix agree with covariance matrix (?)
    if (hess_code == covariance_code) {
      cat(" * the hessian matrix is in agreement with the covariance matrix \n")
    } else {
      cat(" * the hessian matrix is NOT in agreement with the covariance matrix \n")
    }
  }

  ### a data frame with parameter names & respective SEs
  df.se = data.frame(names = colnames(cov_mat), SE = sqrt(diag(cov_mat)))

  ### alpha=0.05 critical value
  z_crit = stats::qnorm(p = 1 - 0.05/2, lower.tail = TRUE, log.p = FALSE)

  ### compute SE: SQRT( diag of -H^(-1) )
  df.theta = merge(df.theta, df.se, by = 'names', all.x = TRUE, sort = FALSE)
  df.theta$L_CI = df.theta$mle_theta - z_crit * df.theta$SE
  df.theta$U_CI = df.theta$mle_theta + z_crit * df.theta$SE

  ### output
  out = list(SE = df.theta, covariance_code = covariance_code, hess_code = hess_code, Covariance = cov_mat)
  class(out) = c("lctmc_2x2.se", "list")
  return(out)
}

#' @rdname get_SE_lctmc
get_SE_lctmc_3x3 = function(em,
                            df,
                            df.Xmat,
                            df.Wmat,
                            df.dt,
                            par_constraint,
                            K,
                            solve.tol = (.Machine$double.eps)^2,
                            symmetric.tol = 5e-11,
                            eigen0.tol = 1e-10) {
  ### perform checks
  if (!("lctmc_3x3.mle" %in% class(em))) {
    stop("`em` should be a custom class list object 'lctmc_3x3.mle' obtained from `EM_lctmc_3x3()`")
  }

  ### extract the EM run with best log(P(Y)) value
  best_index = sapply(em, function(x) x$LPY_value)
  best_index = max(which(best_index == max(best_index))) # in case of ties, get the last index

  ### align results with true parameter values
  df.theta = align_MLE_3x3(mle = em[[best_index]]$pars_value, K = K)
  df.theta = df.theta[!colnames(df.theta) %in% c("true_theta")]

  ### msg
  cat(" * best EM run occurred at step ", best_index, "/", length(em), "\n", sep = "")

  ### get vector of the MLE (without constrained elements)
  mle = df.theta$mle_theta[!(df.theta$names %in% names(par_constraint))]
  mle.names = df.theta$names[!(df.theta$names %in% names(par_constraint))]

  ### numerical derivative to get Hessian matrix
  theta.names.bik = gen_theta_names(K = K, type = "3x3", purpose = "bik")
  hess = numDeriv::hessian(
    func = function(x) {
      ## theta vector, adding back the constrained elements
      names(x) = mle.names
      x = append(x, par_constraint)
      ## compute log(P(Y))
      bik_all.mle = bik_all_3x3(
        theta = x,
        data = df,
        Xmat = df.Xmat,
        Wmat = df.Wmat,
        dt = df.dt,
        K = K,
        theta.names = theta.names.bik
      )
      bik_all.mle = impute_bik(bik_all.mle, eps = 1e-3, EPS = 1e-24)
      ## sum over all class
      bi = Reduce(`+`, bik_all.mle)
      ## return
      sum(log(bi))
    },
    x = mle,
    method.args = list(r = 4) # use 6 for more accuracy but longer run time
  )

  ### Hessian and covariance matrix
  colnames(hess) = rownames(hess) = mle.names
  cov_mat = -1 * solve(a = hess, tol = solve.tol)

  ### check for covariance matrix
  covariance_code = -1
  if (isSymmetric(cov_mat, tol = symmetric.tol)) {
    ## eigenvalues
    covariance.eigen = eigen(cov_mat, only.values = TRUE)$values
    covariance.eigen[abs(covariance.eigen) < eigen0.tol] = 0

    ## check definiteness
    if (all(covariance.eigen > 0)) {
      covariance_code = 2
      cat(" * the covariance matrix is positive definite ~ at least a local maxima is reached \n")
    } else if (all(covariance.eigen >= 0)) {
      covariance_code = 1
      cat(" * the covariance matrix is positive semi-definite ~ this is inconclusive \n")
    } else {
      covariance_code = 0
      cat(" * the covariance matrix is not positive (semi-)definite ~ this is a saddle point \n")
    }
  }

  ### check for hessian matrix
  hess_code = -1
  if (isSymmetric(hess, tol = symmetric.tol)) {
    ## eigenvalues
    hess.eigen = eigen(hess, only.values = TRUE)$values
    hess.eigen[abs(hess.eigen) < eigen0.tol] = 0

    ## check definiteness
    if (all(hess.eigen < 0)) {
      hess_code = 2
    } else if (all(hess.eigen <= 0)) {
      hess_code = 1
    } else {
      hess_code = 0
    }

    ## does the hessian matrix agree with covariance matrix (?)
    if (hess_code == covariance_code) {
      cat(" * the hessian matrix is in agreement with the covariance matrix \n")
    } else {
      cat(" * the hessian matrix is NOT in agreement with the covariance matrix \n")
    }
  }

  ### a data frame with parameter names & respective SEs
  df.se = data.frame(names = colnames(cov_mat), SE = sqrt(diag(cov_mat)))

  ### alpha=0.05 critical value
  z_crit = stats::qnorm(p = 1 - 0.05/2, lower.tail = TRUE, log.p = FALSE)

  ### compute SE: SQRT( diag of -H^(-1) )
  df.theta = merge(df.theta, df.se, by = 'names', all.x = TRUE, sort = FALSE)
  df.theta$L_CI = df.theta$mle_theta - z_crit * df.theta$SE
  df.theta$U_CI = df.theta$mle_theta + z_crit * df.theta$SE

  ### output
  out = list(SE = df.theta, covariance_code = covariance_code, hess_code = hess_code, Covariance = cov_mat)
  class(out) = c("lctmc_3x3.se", "list")
  return(out)
}

