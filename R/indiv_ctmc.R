#' @title Estimate the CTMC parameters for each individual
#'
#' @description This function is part of the initial value generation process, where we first fit a CTMC model to each person in the input dataset.
#' Then the estimated effects are pooled and the K-means algorithm will be performed later to split the estimate into `K` clusters.
#'
#' @name indiv_ctmc
#'
#' @param subject_list a vector of ID numbers (charcter or numeric) matching with the `df$id`. The IDs in this vector should be unique
#' @param df a data frame object transformed via `fmt_rowwise_trans()`
#' @param Xmat a matrix of covariates that affects the CTMC process
#' @param dt a numeric vector holding the time interval between observations
#' @param trace a logical scalar. If TRUE then function will print progress update every 500 individuals fitted. If FALSE, no messages are printed.
#'
#' @return a matrix object housing the individual estimated CTMC parameter. \cr
#' The number of rows is equal to the number of individuals fitted. \cr
#' The number of columns is equal to 7, where the first six columns are the estimated coefficients and the last columns indicate whether convergence criteria is met for each person.
#'
#' @seealso [gen_inits01_lctmc_2x2()]; [gen_inits01_lctmc_3x3()]
#'
#' @example inst/examples/ex_indiv_ctmc.R
NULL

#' @rdname indiv_ctmc
indiv_ctmc_2x2 = function(subject_list, df, Xmat, dt, trace) {
  ### allocate space for each person's q_rs estimates
  person_logQ = matrix(NA, nrow = length(subject_list), ncol = 7)

  ### note that number of params does NOT change with `K`
  colnames(person_logQ) = c("beta0.12", "beta1.12", "beta2.12", # 1 --> 2
                            "beta0.21", "beta1.21", "beta2.21", # 2 --> 1
                            "convergence")

  ### for-loop for individualized fitting
  for (i in seq_along(subject_list)) {
    ## get i^th person data
    i_indices = df$id == subject_list[i]
    df.i = df[i_indices, ]
    dt.i = dt[i_indices]
    Xmat.i = Xmat[i_indices, ]

    ## uniform random initial values for each person
    par_init = stats::runif(n = 6, min = -0.15, 0.15)
    names(par_init) = c("beta0.12", "beta1.12", "beta2.12",
                        "beta0.21", "beta1.21", "beta2.21")

    ## get estimate of person i
    opt = stats::optim(
      par = par_init,
      fn = function(par) {
        # compute rates & P mat
        log_q12 = as.numeric(Xmat.i %*% par[1:3])
        log_q21 = as.numeric(Xmat.i %*% par[4:6])
        q12 = exp(log_q12)
        q21 = exp(log_q21)
        P = get_P_2x2(
          q12 = q12,
          q21 = q21,
          dt = dt.i
        )

        # log-likelihood
        L.ij = df.i$trans.1_1 * P$P11 + df.i$trans.1_2 * P$P12 +
          df.i$trans.2_1 * P$P21 + df.i$trans.2_2 * P$P22

        # impute
        L.ij[is.na(L.ij)] = -1
        if (any(L.ij <= 0)) {
          if (all(L.ij <= 0)) {
            L.ij = rep(1e-24, length(L.ij))
          } else {
            L.ij[L.ij <= 0] = min(L.ij[L.ij > 0]) * 1e-3
            L.ij[L.ij <= 0] = min(L.ij[L.ij > 0])
          }
        }

        # output
        sum(log(L.ij))
      },
      method = "Nelder-Mead",
      control = list(fnscale = -50, maxit = 500)
    )

    ## append result
    person_logQ[i, ] = c(opt$par, opt$convergence)

    ## msg
    mod_500 = (i %% 500) == 0
    if (trace && mod_500) {
      step1_counter = paste(i, "/", length(subject_list), sep = "")
      cat("   - finished person# ", step1_counter, "\n", sep = "")
    }
  }

  ### return
  rownames(person_logQ) = subject_list
  return(person_logQ)
}

#' @rdname indiv_ctmc
indiv_ctmc_3x3 = function(subject_list, df, Xmat, dt, trace) {
  ### allocate space for each person's q_rs estimates
  person_logQ = matrix(NA, nrow = length(subject_list), ncol = 10)

  ### note: number of params does NOT change with `K`
  colnames(person_logQ) = c("beta0.12", "beta1.12", "beta2.12", # 1 --> 2
                            "beta0.21", "beta1.21", "beta2.21", # 2 --> 1
                            "beta0.23", "beta1.23", "beta2.23", # 2 --> 3
                            "convergence")

  ### for-loop for individualized fitting
  for (i in seq_along(subject_list)) {
    ## get i^th person data
    i_indices = df$id == subject_list[i]
    df.i = df[i_indices, ]
    dt.i = dt[i_indices]
    Xmat.i = Xmat[i_indices, ]

    ## uniform random initial values for each person
    par_init = stats::runif(n = 9, min = -0.15, 0.15)
    names(par_init) = c("beta0.12", "beta1.12", "beta2.12",
                        "beta0.21", "beta1.21", "beta2.21",
                        "beta0.23", "beta1.23", "beta2.23")

    ## get estimate of person i
    opt = stats::optim(
      par = par_init,
      fn = function(par) {
        # compute rates & P mat
        log_q12 = as.numeric(Xmat.i %*% par[1:3])
        log_q21 = as.numeric(Xmat.i %*% par[4:6])
        log_q23 = as.numeric(Xmat.i %*% par[7:9])
        q12 = exp(log_q12)
        q21 = exp(log_q21)
        q23 = exp(log_q23)
        P = get_P_3x3(
          q12 = q12,
          q21 = q21,
          q23 = q23,
          dt = dt.i
        )

        # log-likelihood
        L.ij = df.i$trans.1_1*P$P11 + df.i$trans.1_2*P$P12 + df.i$trans.1_3*P$P13_exact +
          df.i$trans.2_1*P$P21 + df.i$trans.2_2*P$P22 + df.i$trans.2_3*P$P23_exact

        # impute
        L.ij[is.na(L.ij)] = -1
        if (any(L.ij <= 0)) {
          if (all(L.ij <= 0)) {
            L.ij = rep(1e-24, length(L.ij))
          } else {
            L.ij[L.ij <= 0] = min(L.ij[L.ij > 0]) * 1e-3
            L.ij[L.ij <= 0] = min(L.ij[L.ij > 0])
          }
        }

        # output
        sum(log(L.ij))
      },
      method = "Nelder-Mead",
      control = list(fnscale = -50, maxit = 500)
    )

    ## append result
    person_logQ[i, ] = c(opt$par, opt$convergence)

    ## msg
    mod_500 = (i %% 500) == 0
    if (trace && mod_500) {
      step1_counter = paste(i, "/", length(subject_list), sep = "")
      cat("   - finished person# ", step1_counter, "\n", sep = "")
    }
  }

  ### return
  rownames(person_logQ) = subject_list
  return(person_logQ)
}
