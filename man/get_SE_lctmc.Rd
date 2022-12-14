% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_SE_lctmc.R
\name{get_SE_lctmc}
\alias{get_SE_lctmc}
\alias{get_SE_lctmc_2x2}
\alias{get_SE_lctmc_3x3}
\title{SE Approximation for the Latent CTMC model}
\usage{
get_SE_lctmc_2x2(
  em,
  df,
  df_Xmat,
  df_Wmat,
  df_dt,
  K,
  par_constraint,
  solve.tol = (.Machine$double.eps)^2,
  symmetric.tol = 5e-11,
  eigen0.tol = 1e-10
)

get_SE_lctmc_3x3(
  em,
  df,
  df_Xmat,
  df_Wmat,
  df_dt,
  K,
  par_constraint,
  solve.tol = (.Machine$double.eps)^2,
  symmetric.tol = 5e-11,
  eigen0.tol = 1e-10
)
}
\arguments{
\item{em}{a list object of with the custom class 'lctmc_2x2.EM' or 'lctmc_3x3.EM'. This type of object is the output of the \code{EM_lctmc} functions}

\item{df}{a data frame object containing row-wise transition data as binary variables. \cr
For example, if \code{trans.2_1} equals 1 then it means the observation was a transition from stage 2 to stage 1 within \code{df_dt} amount of time.}

\item{df_Xmat}{a matrix object housing the covariates that affect the CTMC portion of the model. \cr
This matrix should have the same number of rows as the data frame object, \code{df}}

\item{df_Wmat}{a matrix object housing the covariates that affect the latent classification part of the model. \cr
This matrix should have number of rows equal to unique number of individuals in the data frame object, \code{df}}

\item{df_dt}{a numeric vector housing the length of time interval between observations.
This vector's length should be equal to number of rows in the data frame object, \code{df}}

\item{K}{an integer scalar. Use this variable to tell the function how many latent classes there should be. \cr}

\item{par_constraint}{See documentation in \code{\link[=lctmc_2x2]{lctmc_2x2()}} or \code{\link[=lctmc_3x3]{lctmc_3x3()}}}

\item{solve.tol}{a numeric scalar, typically a small decimal value. It is the tolerance for detecting linear dependencies in the hessian matrix. \cr
Defaults to \code{(.Machine$double.eps)^2} if not specified.}

\item{symmetric.tol}{a numeric scalar. Tolerance value for checking symmetric matrix. \cr
Default is 5e-11}

\item{eigen0.tol}{a numeric scalar. Tolerance value for eigenvalues, any values smaller than this will be treated as 0. \cr
Default is 1e-10}
}
\value{
a list object containing 4 elements:
\describe{
\item{SE}{a data frame object containing columns for:
\eqn{\hat{\theta}_{mle}}, \eqn{SE(\hat{\theta}_{mle})}, the 95\% CI, and the Wald's test \eqn{p} value
for testing \eqn{H_{o}: \hat{\theta}_{mle} = 0}}
\item{covariance_code}{a numeric scalar that can take be one of three values: \cr
0: when the covariance matrix is neither positive definite or semi positive definite \cr
1: when the covariance matrix is semi positive definite \cr
2: when the covariance matrix is positive definite}
\item{hess_code}{similar to \code{covariance_code}: \cr
0: when the hessian is neither negative definite or semi negative definite (saddle point) \cr
1: when the covariance matrix is semi positive definite (inconclusive result) \cr
2: when the covariance matrix is positive definite (a local optimal point)}
\item{Covariance}{a matrix object which is the estimated variance covariance matrix for the estimated parameters. \cr
Note that \eqn{SE(\hat{\theta}_{mle})} is simply square root of the diagonal elements.}
}
}
\description{
Uses the large sample approximation for to estimate the standard error of the MLE (see \strong{Note} for more info). \cr
The result of the approximation is the variance-covariance matrix for the vector of the estimated MLE,
where the diagonal elements are the variance and off-diagonals are the covariance.
}
\note{
This is step five out of six for fitting a latent class CTMC model (i.e., SE approximation via the hessian matrix). \cr\cr
This method relies on the large sample theory that
\deqn{
  \hat\theta \sim N(\theta, -H^{-1}(\hat\theta))
}
where, \eqn{H(\hat\theta)} is the hessian of the observed log likelihood function evaluated at the MLE. And \eqn{\hat\theta} is the estimated MLE.
Since this SE approximation relies on \eqn{N \rarr \infty}, user should be cautious when considering sample size and number of model parameters.
}
\examples{
# The example below demonstrates how to specify the 'lctmc' functions correctly
#   to perform the latent class CTMC model.
# -------------------------------------------------------------------------------------
# Both the 2x2 and 3x3 examples are listed below, however it should be noted
#   that the models may take a while to run (approx. 10-30min on a high spec laptop)
# -------------------------------------------------------------------------------------
# In the 2x2 model, we constrained the model so that the covariate effects on
#   the latent class probability are equal to 0.
# -------------------------------------------------------------------------------------
# In the 3x3 model, we left the model un-constrained.
# -------------------------------------------------------------------------------------
# In this code, we also use a for-loop to fit models with 2, 3, and 4 latent classes
#   the actual number of latent classes within these data sets is 2.
# -------------------------------------------------------------------------------------
# The model with the smallest BIC is the "best" model.
# -------------------------------------------------------------------------------------

\dontrun{
  ## this is a 2x2 example, with 3 latent classes
  data("example_df2x2", package = "LCTMC")
  ctrl_2x2 = LCTMC::create_controls(type = "2x2", data = example_df2x2)

  m2x2_list = list()
  for (k in as.integer(2:4)) {
    m2x2 = LCTMC::lctmc_2x2(
      # data
      data = example_df2x2,
      # general model specification
      K = k,
      X_names = c('x0', 'x1', 'x2'),
      W_names = c('w0', 'w1', 'w2'),
      par_constraint = c(alpha1.1 = 0, alpha2.1 = 0),
      # misc.
      controls = ctrl_2x2,
      parallel_optim = list(
        run = TRUE, cl = parallel::makeCluster(spec = parallel::detectCores()-1)
      ),
      MyModelName = paste("My 2x2 (K=", k, ") model", sep = "")
    )

    m2x2_list[[paste("k", k, sep = "")]] = m2x2
  }

  BIC = sapply(
    X = m2x2_list[-1],
    FUN = function(x) {
      log_like = LCTMC::test_global_optim(m = x, data = example_df2x2)$L_mle
      k = x$n_pars
      n = x$n_trans
      -2*log_like + k*log(n)
    }
  )
  BIC
  which.min(BIC)


  # - # - # - # - # - # - # - # - #


  ## this is a 3x3 example, with 3 latent classes
  data("example_df3x3", package = "LCTMC")
  ctrl_3x3 = LCTMC::create_controls(type = "3x3", data = example_df3x3)

  m3x3_list = list()
  for (k in as.integer(2:4)) {
    m3x3 = LCTMC::lctmc_3x3(
      # data
      data = example_df3x3,
      # general model specification
      K = k,
      X_names = c('x0', 'x1', 'x2'),
      W_names = c('w0', 'w1', 'w2'),
      par_constraint = NULL,
      # misc.
      controls = ctrl_3x3,
      parallel_optim = list(
        run = T,
        cl = parallel::makeCluster(spec = parallel::detectCores()-1)
      ),
      MyModelName = paste("My 3x3 (K=", k, ") model", sep = "")
    )

    m3x3_list[[paste("k", k, sep = "")]] = m3x3
  }

  BIC = sapply(
    X = m3x3_list[-1],
    FUN = function(x) {
      log_like = LCTMC::test_global_optim(m = x, data = example_df3x3)$L_mle
      k = x$n_pars
      n = x$n_trans
      -2*log_like + k*log(n)
    }
  )
  BIC
  which.min(BIC)
}
}
\seealso{
\code{\link[=lctmc_2x2]{lctmc_2x2()}}; \code{\link[=EM_lctmc_2x2]{EM_lctmc_2x2()}}; \code{\link[=rescale_theta]{rescale_theta()}}
}
