#' @title Re-scales model parameters
#'
#' @description This function scales the model parameter by some user specified scaling factors.
#' It can take either parameters in data.frame format or in vector form.
#'
#' @param df_theta a data frame object with data stored in long-format
#' @param v_theta a named numeric scalar. It is a scaling parameter for observation time intervals. \cr
#' @param scaling a named numeric vector indicating the factor of scaling. \cr
#' For example (1), if `scaling = c(x1 = 0.5)`, then the coefficient associated with 'x1' are multiplied by a factor of 0.5. \cr
#' For example (2), if `scaling = c(dt = 0.5)`, then the coefficient associated with 'dt' are added by a factor of log(0.5).
#' @param mult_vars a character vector, indicating which columns in `df_theta` needs multiplicative re-scaling. \cr
#' Default is `c("mle_theta", "SE", "L_CI", "U_CI")`
#' @param add_vars a character vector, indicating which columns in `df_theta` needs additive re-scaling. \cr
#' Default is `c("mle_theta", "L_CI", "U_CI")`.
#'
#' @return A list object containing the 2 elements:
#' \itemize{
#'   \item **df_theta** a data.frame object with its columns re-scaled accordingly. Its dimension is equal to the input `df_theta`. \cr
#'   If the input `df_theta` is NULL then the output will also be NULL
#'   \item **v_theta** a named numeric vector with its elements re-scaled accordingly. Its length is equal to the input `v_theta`. \cr
#'   If the input `v_theta` is NULL then the output will also be NULL
#' }
#'
#' @note This is the final step of fitting a latent class CTMC model (i.e., scaling estimated parameters back to the original data's units).
#'
#' @seealso [lctmc_2x2()] [get_SE_lctmc_2x2()]
#'
#' @example inst/examples/ex_running_lctmc.R

rescale_theta = function(df_theta = NULL,
                         v_theta = NULL,
                         scaling = numeric(),
                         mult_vars = c("mle_theta", "SE", "L_CI", "U_CI"),
                         add_vars = c("mle_theta", "L_CI", "U_CI")) {
  ### checks
  if (!is.data.frame(df_theta)) {
    stop("`df_theta` should be a data.frame object")
  }
  if (!all(c("names", mult_vars, add_vars) %in% colnames(df_theta))) {
    err = paste("all of the following variables should be a column in `df_theta` \n",
                paste(c("names", mult_vars, add_vars), collapse = ""), sep = "")
    stop(err)
  }
  if (!is.numeric(v_theta) || is.null(names(v_theta))) {
    stop("`v_theta` should be a named numeric vector")
  }

  ### df_theta
  if (!is.null(df_theta)) {
    ## all the indices
    a1_index = grepl(pattern = "alpha1", df_theta$names)
    a2_index = grepl(pattern = "alpha2", df_theta$names)
    b1_index = grepl(pattern = "beta1", df_theta$names)
    b2_index = grepl(pattern = "beta2", df_theta$names)
    b0_index = grepl(pattern = "beta0", df_theta$names)

    ## re-scale
    df_theta[a1_index, mult_vars] = df_theta[a1_index, mult_vars] * scaling['w1']
    df_theta[a2_index, mult_vars] = df_theta[a2_index, mult_vars] * scaling['w2']
    cat(" * `df_theta` ~ alphas (1 & 2) re-scaled \n", sep = "")

    df_theta[b1_index, mult_vars] = df_theta[b1_index, mult_vars] * scaling['x1']
    df_theta[b2_index, mult_vars] = df_theta[b2_index, mult_vars] * scaling['x2']
    cat(" * `df_theta` ~ betas (1 & 2) re-scaled \n", sep = "")

    df_theta[b0_index, add_vars] = df_theta[b0_index, add_vars] + log(scaling['dt'])
    cat(" * `df_theta` ~ log(r0) re-scaled \n", sep = "")
  }

  ### v_theta
  if (!is.null(v_theta)) {
    ## all the indices
    a1_index = grepl(pattern = "alpha1", names(v_theta))
    a2_index = grepl(pattern = "alpha2", names(v_theta))
    b1_index = grepl(pattern = "beta1", names(v_theta))
    b2_index = grepl(pattern = "beta2", names(v_theta))
    b0_index = grepl(pattern = "beta0", names(v_theta))

    ## re-scale
    v_theta[a1_index] = v_theta[a1_index] * scaling['w1']
    v_theta[a2_index] = v_theta[a2_index] * scaling['w2']
    cat(" * `v_theta` ~ alphas (1 & 2) re-scaled \n", sep = "")

    v_theta[b1_index] = v_theta[b1_index] * scaling['x1']
    v_theta[b2_index] = v_theta[b2_index] * scaling['x2']
    cat(" * `v_theta` ~ betas (1 & 2) re-scaled \n", sep = "")

    v_theta[b0_index] = v_theta[b0_index] + log(scaling['dt'])
    cat(" * `v_theta` ~ log(r0) re-scaled \n", sep = "")
  }

  ### return
  return(list(df_theta = df_theta, v_theta = v_theta))
}
