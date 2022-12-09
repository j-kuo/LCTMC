#' @title Compute transition rate ratios (RR)
#'
#' @description Computes the within-class and between-class transition rate ratios from the \eqn{Q} matrix. \cr
#' The within-class RR compares the transition rate *within* latent class (i.e., holding the latent class constant),
#' when covariates are increased by 1-unit. \cr
#' The between-class rate ratios compares the transition *across* latent class,
#' while holding covariates constant at the user specified values.
#'
#' @name compute_RR
#'
#' @param m a 'lctmc_2x2' or 'lctmc_3x3' object obtained from the `lctmc` functions
#' @param x1 a numeric scalar. It is the covariate corresponding to the parameters `beta1.12_k`, `beta1.21_k`, and `beta1.23_k`
#' @param x2 a numeric scalar. It is the covariate corresponding to the parameters `beta2.12_k`, `beta2.21_k`, and `beta2.23_k`
#' @param ci_method a character scalar. Either the "delta" method or "endpoint" method.
#'
#' @return A list object of length two.
#' The first element contains the RR for within-class transition rates, for 1-unit increase in covariate values.
#' The second element contains the RR for between-class comparison. Holding all covariates constant at the user specified values.
#'
#' @note For within-class comparisons, the increment in covariate is always 1-unit.
#' The RR for 1-unit increase in \eqn{x_{1}}, for the transition from stage \eqn{r} to \eqn{s} of class \eqn{k} is the following:
#' \deqn{
#'   RR_{rs}^{(k)(W)} = \exp(\beta_{1rs}^{(k)})
#' }
#' where \eqn{\beta_{1rs}^{(k)}} is the coefficient associated with \eqn{x_{1}}. \cr
#' For between-class comparisons, the last class is always used as the referent group (denoted as \eqn{K'}).
#' The RR for the transition from stage \eqn{r} to \eqn{s} of class \eqn{k} vs. the referent class is:
#' \deqn{
#'   RR_{rs}^{(k)(B)} = \exp((\beta_{rs}^{(k)} - \beta_{rs}^{(K')}) \cdot X)
#' }
#' where \eqn{X} is the vector of covariates specified by the input arguments `x1` and `x2`.
#'
#' @export
#'
#' @seealso [compute_Q()]
#'
#' @example inst/examples/ex_compute_statistics.R
compute_RR = function(m, x1, x2, ci_method) {
  UseMethod("compute_RR")
}

#' @rdname compute_RR
#' @export
compute_RR.lctmc_2x2 = function(m, x1, x2, ci_method) {
  ## checks
  if (!all(is.numeric(c(x1, x2)))) {
    stop("`x1`, `x2` should both be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1)) {
    stop("`x1`, `x2` should both be length 1")
  }

  ## data frame of coefficients
  df_theta = m$SE$SE

  ## calls the compute Q function to generate transition rates
  fn_call = match.call()
  fn_call[[1]] = as.name("compute_Q")
  fn_call[["ci_method"]] = NULL
  Qlist = eval(fn_call)

  ## within-clsas RR
  z_crit = stats::qnorm(p = 1 - 0.05/2, mean = 0, sd = 1)
  RR_within_list = list()
  RR_L_CI_within_list = list()
  RR_U_CI_within_list = list()
  for (k in 1:m$K) {
    df_theta.k = df_theta[grepl(pattern = paste("_", k, "$", sep = ""), df_theta$names), ]

    RR = L_CI = U_CI = c()
    for (rs in c("12", "21")) {
      # get correct betas as data frame
      betas_rs = df_theta.k[grepl(pattern = rs, df_theta.k$names), ]
      betas_rs_x1 = betas_rs[grepl(pattern = "beta1", betas_rs$names), ]
      betas_rs_x2 = betas_rs[grepl(pattern = "beta2", betas_rs$names), ]

      # rate ratios
      RR_rs_x1 = exp(betas_rs_x1$mle_theta)
      RR_rs_x2 = exp(betas_rs_x2$mle_theta)

      # conf. interval
      if (ci_method == "delta") {
        RR_L_CI_rs_x1 = RR_rs_x1 - z_crit * (betas_rs_x1$SE*RR_rs_x1) # x1
        RR_U_CI_rs_x1 = RR_rs_x1 + z_crit * (betas_rs_x1$SE*RR_rs_x1)
        RR_L_CI_rs_x2 = RR_rs_x2 - z_crit * (betas_rs_x2$SE*RR_rs_x2) # x2
        RR_U_CI_rs_x2 = RR_rs_x2 + z_crit * (betas_rs_x2$SE*RR_rs_x2)
      } else if (ci_method == "endpoint") {
        RR_L_CI_rs_x1 = exp(betas_rs_x1$L_CI) # x1
        RR_U_CI_rs_x1 = exp(betas_rs_x1$U_CI)
        RR_L_CI_rs_x2 = exp(betas_rs_x2$L_CI) # x2
        RR_U_CI_rs_x2 = exp(betas_rs_x2$U_CI)
      }

      # combine x1 & x2
      rr = c(RR_rs_x1, RR_rs_x2)
      l_ci = c(RR_L_CI_rs_x1, RR_L_CI_rs_x2)
      u_ci = c(RR_U_CI_rs_x1, RR_U_CI_rs_x2)

      # make names
      r_to_s = sub("(\\d{1})(\\d{1})", "\\1to\\2", rs)
      names(rr) = paste("RR_", r_to_s, "_", c("x1", "x2"), sep = "")
      names(l_ci) = paste("L_CI_", r_to_s, "_", c("x1", "x2"), sep = "")
      names(u_ci) = paste("U_CI_", r_to_s, "_", c("x1", "x2"), sep = "")

      # combine all combo of x1 & x2
      RR = append(RR, rr)
      L_CI = append(L_CI, l_ci)
      U_CI = append(U_CI, u_ci)
    }

    RR_within_list[[paste("k", k, sep = "")]] = RR
    RR_L_CI_within_list[[paste("k", k, sep = "")]] = L_CI
    RR_U_CI_within_list[[paste("k", k, sep = "")]] = U_CI
  }

  ## between-class RR
  RR_between_list = vector("list", length = m$K)
  for (k in m$K:1) {
    ## current class's Q matrix
    Qmat.k = Qlist[[k]]

    ## rates
    rates = c(RR_12 = Qmat.k[1, 2],
              RR_21 = Qmat.k[2, 1])

    ## if k is referent, create referent class rates
    if (k == m$K) {
      rates_ref = rates
    }

    ## append rate ratios
    RR_between_list[[k]] = rates / rates_ref
  }
  names(RR_between_list) = paste("k", 1:m$K, sep = "")

  ## return
  out = list(within.RR = RR_within_list,
             within.RR.L_CI = RR_L_CI_within_list,
             within.RR.U_CI = RR_U_CI_within_list,
             between.RR = RR_between_list)
  return(out)
}

#' @rdname compute_RR
#' @export
compute_RR.lctmc_3x3 = function(m, x1, x2, ci_method) {
  ## checks
  if (!all(is.numeric(c(x1, x2)))) {
    stop("`x1`, `x2` should both be numeric values")
  }
  if ((length(x1) != 1) | (length(x2) != 1)) {
    stop("`x1`, `x2` should both be length 1")
  }

  ## data frame of coefficients
  df_theta = m$SE$SE

  ## calls the compute Q function to generate transition rates
  fn_call = match.call()
  fn_call[[1]] = as.name("compute_Q")
  fn_call[["ci_method"]] = NULL
  Qlist = eval(fn_call)

  ## within-clsas RR
  z_crit = stats::qnorm(p = 1 - 0.05/2, mean = 0, sd = 1)
  RR_within_list = list()
  RR_L_CI_within_list = list()
  RR_U_CI_within_list = list()
  for (k in 1:m$K) {
    df_theta.k = df_theta[grepl(pattern = paste("_", k, "$", sep = ""), df_theta$names), ]

    RR = L_CI = U_CI = c()
    for (rs in c("12", "21", "23")) {
      # get correct betas as data frame
      betas_rs = df_theta.k[grepl(pattern = rs, df_theta.k$names), ]
      betas_rs_x1 = betas_rs[grepl(pattern = "beta1", betas_rs$names), ]
      betas_rs_x2 = betas_rs[grepl(pattern = "beta2", betas_rs$names), ]

      # rate ratios
      RR_rs_x1 = exp(betas_rs_x1$mle_theta)
      RR_rs_x2 = exp(betas_rs_x2$mle_theta)

      # conf. interval
      if (ci_method == "delta") {
        RR_L_CI_rs_x1 = RR_rs_x1 - z_crit * (betas_rs_x1$SE*RR_rs_x1) # x1
        RR_U_CI_rs_x1 = RR_rs_x1 + z_crit * (betas_rs_x1$SE*RR_rs_x1)
        RR_L_CI_rs_x2 = RR_rs_x2 - z_crit * (betas_rs_x2$SE*RR_rs_x2) # x2
        RR_U_CI_rs_x2 = RR_rs_x2 + z_crit * (betas_rs_x2$SE*RR_rs_x2)
      } else if (ci_method == "endpoint") {
        RR_L_CI_rs_x1 = exp(betas_rs_x1$L_CI) # x1
        RR_U_CI_rs_x1 = exp(betas_rs_x1$U_CI)
        RR_L_CI_rs_x2 = exp(betas_rs_x2$L_CI) # x2
        RR_U_CI_rs_x2 = exp(betas_rs_x2$U_CI)
      }

      # combine x1 & x2
      rr = c(RR_rs_x1, RR_rs_x2)
      l_ci = c(RR_L_CI_rs_x1, RR_L_CI_rs_x2)
      u_ci = c(RR_U_CI_rs_x1, RR_U_CI_rs_x2)

      # make names
      r_to_s = sub("(\\d{1})(\\d{1})", "\\1to\\2", rs)
      names(rr) = paste("RR_", r_to_s, "_", c("x1", "x2"), sep = "")
      names(l_ci) = paste("L_CI_", r_to_s, "_", c("x1", "x2"), sep = "")
      names(u_ci) = paste("U_CI_", r_to_s, "_", c("x1", "x2"), sep = "")

      # combine all combo of x1 & x2
      RR = append(RR, rr)
      L_CI = append(L_CI, l_ci)
      U_CI = append(U_CI, u_ci)
    }

    RR_within_list[[paste("k", k, sep = "")]] = RR
    RR_L_CI_within_list[[paste("k", k, sep = "")]] = L_CI
    RR_U_CI_within_list[[paste("k", k, sep = "")]] = U_CI
  }

  ## between-class RR
  RR_between_list = vector("list", length = m$K)
  for (k in m$K:1) {
    ## current class's Q matrix
    Qmat.k = Qlist[[k]]

    ## rates
    rates = c(RR_12 = Qmat.k[1, 2],
              RR_21 = Qmat.k[2, 1],
              RR_23 = Qmat.k[2, 3])

    ## if k is referent, create referent class rates
    if (k == m$K) {
      rates_ref = rates
    }

    ## append rate ratios
    RR_between_list[[k]] = rates / rates_ref
  }
  names(RR_between_list) = paste("k", 1:m$K, sep = "")

  ## return
  out = list(within.RR = RR_within_list,
             within.RR.L_CI = RR_L_CI_within_list,
             within.RR.U_CI = RR_U_CI_within_list,
             between.RR = RR_between_list)
  return(out)
}
