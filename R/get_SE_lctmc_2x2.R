#' @title SE Approximation for the Latent CTMC model (2x2)
#'
#' @description Uses the large sample approximation for to estimate the standard error of the MLE (see **Note** for more info). \cr
#' The result of the approximation is the variance-covariance matrix for the vector of the estimated MLE, where the diagonal elements are the variance and off-diagonals are the covariance.
#'
#' @param em a list object of with custom class 'lctmc_2x2.mle'. This type of object is the output of the function `EM_lctmc_2x2()`
#' @param df a data frame object containing the binary row-wise transition indicator variables (usually obtained from `fmt_rowwise_2x2trans()`)
#' @param df.Xmat a matrix object with same number of rows as `df`. This matrix object should contain the covariates which affect the CTMC part of the model
#' @param df.Wmat a matrix object with number of rows equal to the unique number of individuals in `df`.
#' This matrix object should contain the covariates which affect the latent class probability part of the model.
#' @param df.dt a numeric vector with length equal to the number of rows as `df`. This vector contains the time difference between observations.
#' @param hessian.round a numeric scalar. This variable is used to specified how many decimal point should the hessian approximation be rounded to.
#' This is necessary because by definition the hessian matrix should be symmetric. However, due to it being obtained from numerical approximation, it will not always be symmetric (rounding errors).
#' So, this variable is used to round to a reasonable number of decimal places to ensure the hessian is symmetric. \cr
#' If this number is too small (i.e., large margin of error), then the hessian will be inaccurate. Contrary, if this number is too large, the approximated hessian will not be symmetric and an error will be returned.
#' @param par_constraint a named numeric vector to indicate which parameter is constrained. Set equal to NULL for unconstrained model. \cr
#' For example, `c(alpha1.1 = 0)` constraints the parameter 'alpha1.1' to be a constant 0. **NOTE:** Current version of the code will *only* work with constrains equal to 0.
#' @param K the number of categories the latent class variable has.
#' @param solve.tol a numeric scalar, typically a small decimal value. It is the tolerance for detecting linear dependencies in the hessian matrix. Defaults to `(.Machine$double.eps)^2` if not specified.
#' @param MyModelName a character scalar. Gives the current model fitting process a name. This name will be used when the function is logging the algorithm progress.
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
#' @note This function uses the large sample theory
#' \deqn{
#'   \hat\theta \sim N(\theta, -H^{-1}(\hat\theta))
#' }
#' where, \eqn{H(\hat\theta)} is the hessian of the observed log likelihood function evaluated at the MLE. And \eqn{\hat\theta} is the estimated MLE.
#' Since this SE approximation relies on large sample, user should be cautious when considering sample size and number of model parameters.
#'
#' @importFrom numDeriv hessian
#' @importFrom matrixcalc is.positive.definite
#' @importFrom matrixcalc is.positive.semi.definite
#' @importFrom matrixcalc is.negative.definite
#' @importFrom matrixcalc is.negative.semi.definite
#'
#' @export
#'
#' @example inst/examples/ex_get_SE_lctmc_2x2.R

get_SE_lctmc_2x2 = function(em,
                            df,
                            df.Xmat,
                            df.Wmat,
                            df.dt,
                            hessian.round,
                            par_constraint,
                            K,
                            solve.tol = (.Machine$double.eps)^2,
                            MyModelName) {
  ### perform checks
  if (!("lctmc_2x2.mle" %in% class(em))) {
    stop("`em` should be a custom class list object 'lctmc_2x2.mle' obtained from `EM_lctmc_2x2()`")
  }

  ### align results with true parameter values
  best_index = which.max(sapply(em, function(x) x$LPY_value))
  df.theta = align_MLE_2x2(mle = em[[best_index]]$pars_value, K = K)
  df.theta = df.theta[!colnames(df.theta) %in% c("true_theta")]

  ### msg
  cat("####~{", MyModelName, "}~", paste(rep("#", times = 100 - 8 - nchar(MyModelName)), sep = ""),"\n", sep = "")
  cat(paste(rep("#", times = 27), sep = ""), " Computing Hessian Matrix for SE approx. ", paste(rep("#", times = 32), sep = ""), "\n", sep = "")
  cat(" - best EM run occurred at step ", best_index, "/", length(em), "\n", sep = "")

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
  cov_mat = round(x = cov_mat, digits = hessian.round)

  ### check for positive definite-ness
  covariance_code = -1
  hess_code = -1
  tryCatch(
    ## note: 2, 1, 0 = local extrmum, inconclusive, saddle point
    expr = {
      # perform check 1: Covariance matrix
      if (matrixcalc::is.positive.definite(cov_mat)) {
        covariance_code = 2
      } else if (matrixcalc::is.positive.semi.definite(cov_mat)) {
        covariance_code = 1
      } else {
        covariance_code = 0
      }

      # perform check 2: Hessian matrix
      if (matrixcalc::is.negative.definite(hess)) {
        hess_code = 2
      } else if (matrixcalc::is.negative.semi.definite(hess)) {
        hess_code = 1
      } else {
        hess_code = 0
      }
    },
    error = function(e) cat("CAUTION: hessian/covariance matrix might not be symmetric \n\n", sep = "")
  )

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
