#' @title EM-algorithm for a latent CTMC model
#'
#' @description The EM algorithm attempts to maximize the likelihood by iteratively optimize the conditional expected log-likelihood.
#' This optimization approach is used to obtained the MLE of the latent class CTMC model.
#' Convergence for the EM algorithm is met when the improvement in the conditional expected log-likelihood does not improve further
#'
#' @param theta.init a named numeric vector, where the names are the model parameter names. Use this argument to specified the initial values for the EM algorithm
#' @param theta.names a list of character vectors, where each element of this list is a character vector specifying the names of model parameters. \cr
#' Note: model parameters grouped within the same element will be simultaneously optimized during the ECM step.
#' @param par_constraint a named numeric vector to indicate which parameter is constrained. Set equal to NULL for unconstrained model. \cr
#' For example, `c(alpha1.1 = 0)` constraints the parameter 'alpha1.1' to be a constant 0. **NOTE:** Current version of the code will *only* work with constrains equal to 0.
#' @param K an integer scalar. Used to determine the number of latent classes the model is fitting.
#' @param df a data frame object holding the binary row-wise transition indicator variables
#' @param df.Xmat a matrix object housing the covariates that affect the CTMC portion of the model. This matrix should have the same number of rows as the data frame object, `df`
#' @param df.Wmat a matrix object housing the covariates that affect the latent classification part of the model. This matrix should have number of rows equal to unique number of individuals in the data frame object, `df`
#' @param df.dt a numeric vector housing the length of time interval between observations. This vector's length should be equal to number of rows in the data frame object, `df`
#' @param EM.maxit a numeric scalar. Use this argument to specify the maximum iteration the EM algorithm.
#' Setting this variable too large will cause longer run time if optimal point is harder to reach. However, too low will lead to non-optimal points.
#' Typically, somewhere between 100-300 iterations is adequate depending on the problem.
#' @param ELL_diff.tol a numeric scalar that is greater than 0. This is the tolerance value for the expected conditional log-likelihood. \cr
#' The smaller this value, the more precise the estimate, however the longer run time will be needed.
#' @param LPY_diff.tol a numeric scalar that is greater than 0. This is the tolerance value for the log observed likelihood, \eqn{log(P(Y))}. \cr
#' The smaller this value, the more precise the estimate, however the longer run time will be needed.
#' @param par_diff.tol a numeric scalar that is greater than 0. This is the tolerance value for the changes in model parameter by iterations. \cr
#' The smaller this value, the more precise the estimate, however the longer run time will be needed.
#' @param L_BFGS_B.ctrl a list object holding the control parameter for the L-BFGS optimization algorithm. \cr
#' Note that the L-BFGS is utilized during each ECM step, where we attempt to maximize the conditional expected log-likelihood function with respect to a subset of the model parameters. \cr
#' This list object should contain **three** elements (see `?optim` for more info): \cr
#' (1) `fnscale` is a negative real number, this is used to scale the value of the objective function. Setting it to negative implies that a maximization is being performed. \cr
#' (2) `maxit` is the maximum number of iterations the algorithm will run before terminating \cr
#' (3) `factr` is tolerance parameter for algorithm convergence, the smaller the value, the better accracy but also longer run time
#' @param parallel_optim a list object telling the function whether parallel process should be used for the Step 2 of the initial value generation. \cr
#' The list should contain **two** elements: \cr
#' (1) `run` a logical scalar, if TRUE then this function will use parallel processing. If FALSE, then the `cl` argument is ignored. \cr
#' (2) `cl` is an object obtained from the `parallel` package, for example \cr `cl = parallel::makeCluster(spec = 2)`
#' @param MyModelName a character scalar. Gives the current model fitting process a name. This name will be used when the function is logging the algorithm progress.
#'
#' @return A list object containing as many elements as the number of EM iteration that were performed
#' Each element is a sub list object with the following elements:
#' \itemize{
#'   \item **iter** is a numeric scalar, indicating the count of current iteration
#'   \item **pars_value** is a named numeric vector for current iteration's estimated parameter values
#'   \item **ELL_value** is a numeric scalar for the current iteration's expected conditional log-likelihood evaluated at `pars_value`
#'   \item **LPY_value** is a numeric scalar for the current iteration's \eqn{log(P(Y))} evaluated at `pars_value`
#'   \item **par_diff.max** is a numeric scalar comparing previous iteration's `pars_value` and current iteration's `pars_value`. This value is maximum difference between the iterations. \cr
#'         This is one of the values to evaluate EM algorithm's convergence. When `par_diff.max` < `par_diff.tol`, one of the convergence criteria is met.
#'   \item **ELL_diff** is a numeric scalar comparing previous iteration and current iteration's expected conditional log-likelihood. \cr
#'         This is one of the values to evaluate EM algorithm's convergence. When `ELL_diff` < `ELL_diff.tol`, one of the convergence criteria is met.
#'   \item **LPY_diff** is a numeric scalar comparing previous iteration and current iteration's \eqn{log(P(Y))}. \cr
#'         This is one of the values to evaluate EM algorithm's convergence. When `LPY_diff` < `LPY_diff.tol`, one of the convergence criteria is met.
#'   \item **time_at_finish** is a 'difftime' object telling the time from start of algorithm to current iteration
#' }
#'
#' @note This function is actually the Expected Conditional Maximization (ECM) algorithm.
#' within each EM step, subsets of parameters are optimize one set at a time, the number of ECM step is determined by `theta.names`. \cr
#' Note that step 2 of `gen_inits_lctmc_2x2()` will likely produce a good approximation of the MLE. so when the EM algorithm requires large number of iterations
#' it is likely that we are getting a non-global optimal point.
#'
#' @importFrom optimParallel optimParallel
#'
#' @export
#'
#' @example inst/examples/ex_EM_lctmc_2x2.R

EM_lctmc_2x2 = function(theta.init,
                        theta.names,
                        par_constraint,
                        K,
                        df,
                        df.Xmat,
                        df.Wmat,
                        df.dt,
                        EM.maxit,
                        ELL_diff.tol,
                        LPY_diff.tol,
                        par_diff.tol,
                        L_BFGS_B.ctrl,
                        parallel_optim,
                        MyModelName) {
  ### checks
  if ((nrow(df) != nrow(df.Xmat)) || (nrow(df.Xmat) != length(df.dt))) {
    stop("Mis-matching dimensions in either `df`, `df.Xmat`, or `df.dt`")
  }
  if (length(unique(df$id)) != nrow(df.Wmat)) {
    stop("Number of unique ID in `df` does not matches with number of individuals in `df.Wmat`")
  }
  if (!is.numeric(df.dt)) {
    stop("`dt` must be a numeric vector (indicating time intervals)")
  }
  if (length(theta.init) != length(unlist(theta.names))) {
    stop("Mis-matching dimensions in `theta.init` and `theta.names`")
  }

  ### load `optimParallel`
  if (parallel_optim$run) {
    optim2 = optimParallel::optimParallel
  } else {
    optim2 = stats::optim
  }

  ### initialize
  t0 = Sys.time()
  theta.names.bik = gen_theta_names(K = K, type = "2x2", purpose = "bik")

  ### apply constraints to initial values
  constraint_index = names(par_constraint)[names(par_constraint) %in% names(theta.init)]
  theta.init[constraint_index] = par_constraint

  ### step 1: function evaluated at OLD theta
  bik_all.old = bik_all_2x2(theta = theta.init, data = df, Xmat = df.Xmat, Wmat = df.Wmat, dt = df.dt, K = K, P.rs = FALSE, theta.names = theta.names.bik)
  bik_all.old = impute_bik(x = bik_all.old, eps = 1e-3, EPS = 1e-24)

  ### step 1: compute log(P(Y)) and ELL
  denom.0 = Reduce(`+`, bik_all.old)
  numer.0 = Reduce(`+`, Map(f = function(x) x*log(x), bik_all.old))
  opt = list(par = theta.init, value = sum(numer.0/denom.0))
  LPY.0 = sum(log(denom.0))

  ### initialize output object
  EM_output = list()
  EM_output[[1]] = list(
    iter = 1,
    # values
    pars_value = opt$par, ELL_value = opt$value, LPY_value = LPY.0,
    # diffs
    par_diff.max = NA, ELL_diff = NA, LPY_diff = NA,
    # record time elapsed
    time_at_finish = difftime(Sys.time(), t0, units = "mins")
  )

  ### step 2 to step `EM.maxit`
  EM.i = 1
  EM.condition = T
  while (EM.condition) {
    # increment while-loop counter + msg
    EM.i = EM.i + 1
    cat(" * ~~~~~ Starting EM Step:", rep(" ", 4-nchar(EM.i)), EM.i, " - - - (with ", length(theta.names), " sub-steps) ",
        rep("~", 44),
        " \n", rep(" ", 3), sep = "")

    # function evaluated at OLD theta
    ELL_prev = EM_output[[EM.i-1]]$ELL_value
    LPY_prev = EM_output[[EM.i-1]]$LPY_value
    par_prev = EM_output[[EM.i-1]]$pars_value
    if (EM.i > 2) {
      # computes the "new" `bik_all.old` using theta from last round
      par_prev = par_curr
      bik_all.old = bik_all_2x2(theta = par_prev, data = df, Xmat = df.Xmat, Wmat = df.Wmat, dt = df.dt, K = K, P.rs = FALSE, theta.names = theta.names.bik)
      bik_all.old = impute_bik(x = bik_all.old, eps = 1e-3, EPS = 1e-24)
      # get denom
      denom.old = Reduce(`+`, bik_all.old)
    } else {
      # when `EM.i == 2` ... seems redundant but this part is necessary
      par_curr = par_prev
      denom.old = denom.0
    }

    # optim
    for (p.i in 1:length(theta.names)) {
      # print message for current ECM step -- within step `EM.i`
      cat(p.i, ".. ", sep = "")
      p = theta.names[[p.i]]
      # if running parallel, need to export variaable
      if (parallel_optim$run) {
        opt = optim2(
          par = par_curr[p],
          fn = function(par) {
            # my theta
            my_theta = c(par, par_curr[!names(par_curr) %in% p])
            # apply constraints
            my_theta[constraint_index] = par_constraint

            # function evaluated at theta
            bik_all.theta = bik_all_2x2(
              theta = my_theta,
              data = df,
              Xmat = df.Xmat,
              Wmat = df.Wmat,
              dt = df.dt,
              K = K,
              P.rs = FALSE,
              theta.names = theta.names.bik
            )
            bik_all.theta = impute_bik(x = bik_all.theta, eps = 1e-3, EPS = 1e-24)

            # compute expected log likelihood given obs. data
            numer = Reduce(`+`, Map(f = function(x, y) x*log(y), bik_all.old, bik_all.theta))
            out = sum(numer/denom.old)

            # return
            out
          },
          method = "L-BFGS-B",
          control = L_BFGS_B.ctrl,
          parallel = list(cl = parallel_optim$cl, forward = TRUE, loginfo = FALSE)
        )
      }
      if (!parallel_optim$run) {
        opt = optim2(
          par = par_curr[p],
          fn = function(par) {
            # my theta
            my_theta = c(par, par_curr[!names(par_curr) %in% p])
            # apply constraints
            my_theta[constraint_index] = par_constraint

            # function evaluated at theta
            bik_all.theta = bik_all_2x2(
              theta = my_theta,
              data = df,
              Xmat = df.Xmat,
              Wmat = df.Wmat,
              dt = df.dt,
              K = K,
              P.rs = FALSE,
              theta.names = theta.names.bik
            )
            bik_all.theta = impute_bik(x = bik_all.theta, eps = 1e-3, EPS = 1e-24)

            # compute expected loglikelihood given obs. data
            numer = Reduce(`+`, Map(f = function(x, y) x*log(y), bik_all.old, bik_all.theta))
            out = sum(numer/denom.old)

            # return
            out
          },
          method = "L-BFGS-B",
          control = L_BFGS_B.ctrl
        )
      }

      # update within-step for ECM
      par_curr[p] = opt$par
    }

    # compute P(Y) based on OLD parameter estimates
    LPY_prev = sum(log(Reduce(`+`, bik_all.old)))

    # compute P(Y) based on NEW parameter estimates
    b = bik_all_2x2(
      theta = par_curr,
      data = df,
      Xmat = df.Xmat,
      Wmat = df.Wmat,
      dt = df.dt,
      K = K,
      P.rs = FALSE,
      theta.names = theta.names.bik
    )
    b = impute_bik(x = b, eps = 1e-3, EPS = 1e-24)
    LPY_curr = sum(log(Reduce(`+`, b)))

    # print message for ending `EM.i` step
    cat("\n   ", paste(rep("~", 97), collapse = ""),
        "\n   Maximum Theta Difference:", sprintf("%.5f", max(abs(par_curr - par_prev))), "(", names(par_curr)[which.max(abs(par_curr - par_prev))], ")",
        "\n             ELL Difference:", sprintf("%.5f", (opt$value - ELL_prev)), "(", opt$value, ")",
        "\n             LPY Difference:", sprintf("%.5f", (LPY_curr - LPY_prev)), "(", LPY_curr, ")",
        "\n\n", sep = "")

    # update
    EM_output[[EM.i]] = list(
      iter = EM.i,
      # store values
      pars_value = par_curr, ELL_value = opt$value, LPY_value = LPY_curr,
      # compute convergence criteria
      par_diff.max = max(abs(par_curr - par_prev)), ELL_diff = (opt$value - ELL_prev), LPY_diff = (LPY_curr - LPY_prev),
      # misc
      time_at_finish = difftime(Sys.time(), t0, units = "mins")
    )

    # termination condition (1)
    # continue if (LPY_diff > tol or LPY_diff between [0, tol] while ELL_diff > tol)
    check_1 = EM_output[[EM.i]]$LPY_diff > LPY_diff.tol
    check_2.a = 0 <= EM_output[[EM.i]]$LPY_diff
    check_2.b = EM_output[[EM.i]]$LPY_diff <= LPY_diff.tol
    check_2.c = EM_output[[EM.i]]$ELL_diff > ELL_diff.tol
    if (check_1) {
      EM.c1 = TRUE
    } else if (check_2.a && check_2.b && check_2.c) {
      EM.c1 = TRUE
    } else {
      EM.c1 = FALSE
    }

    # termination condition (2)
    # continue if (condition 1 decide to terminate, and LPY_diff between [0, tol] while par_diff > tol)
    check_3.a = !(EM.c1)
    check_3.b = 0 <= EM_output[[EM.i]]$LPY_diff
    check_3.c = EM_output[[EM.i]]$LPY_diff <= LPY_diff.tol
    check_3.d = EM_output[[EM.i]]$par_diff.max > par_diff.tol
    if (check_3.a && check_3.b && check_3.c && check_3.d) {
      EM.c2 = TRUE
    } else {
      EM.c2 = FALSE
    }

    # termination condition (3)
    # if `maxit` is reached terminate regardless of other conditions
    EM.c3 = (EM.i < EM.maxit)

    # finalize condition
    EM.condition = (EM.c1 || EM.c2) && EM.c3
  }

  ### output as custom class
  class(EM_output) = c("lctmc_2x2.mle", "list")
  EM_output
}
