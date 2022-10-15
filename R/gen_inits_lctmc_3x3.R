#' @title Generates initial value for model fitting
#'
#' @description This function generates initial value for the Latent Class CTMC model (special 3x3 case). This is a two-step process, where we first use the K-means algorithm on the individuals' effect to create cluster level estimates.
#' Then, in the second step, a direct optimization on the observed data likelihood using Step 1's estimate is performed to obtain a set of improved initial values.
#' The output of this function should then be fed to the EM algorithm to obtain a tighter Maximum Likelihood Estimate.
#'
#' @param df a data frame object containing row-wise transition data as binary variables. \cr
#' For example, if `trans.2_1` equals 1 then it means the observation was a transition from stage 2 to stage 1 within `dt` amount of time.
#' @param Xmat a matrix object housing the covariates that affect the CTMC portion of the model. \cr
#' This matrix should have the same number of rows as the data frame object, `data`
#' @param Wmat a matrix object housing the covariates that affect the latent classification part of the model. \cr
#' This matrix should have number of rows equal to unique number of individuals in the data frame object, `data`
#' @param dt a numeric vector housing the length of time interval between observations. This vector's length should be equal to number of rows in the data frame object, `data`
#' @param N_sub a numeric scalar. This is used for Step 1 of initial value generation where the algorithm fits the traditional CTMC model for each individual. \cr
#' Fitting the model for *all* individuals might have long run time without improvement in the accuracy of the estimation. Thus, setting a maximum number of individuals to use for the initial value generation could shorten run without sacrificing in estimation accuracy.
#' @param pct_keep a numeric vector where each element of the vector ranges from 0 to 1 (ideally at minimum 0.50). \cr
#' This argument controls what percentage of the individual level estimated effects to be used for the K-means algorithm for initial value generation.
#' For example, for pct_keep = c(0.8), after individuals effects are estimated, only the 10th to 90th percentile are used for the K-Means algorithm to obtain cluster level estimates.
#' @param par_constraint a named numeric vector to indicate which parameter is constrained. Set equal to NULL for unconstrained model. \cr
#' For example, `c(alpha1.1 = 0)` constraints the parameter 'alpha1.1' to be a constant 0. **NOTE:** Current version of the code will *only* work with constrains equal to 0.
#' @param K the number of categories the latent class variable has
#' @param returns an integer scalar, either 1L or 2L. Other values will cause function to return NULL values. \cr
#' This argument determines if the function should return only the Step 1 estimates (`returns = 1L`) or both steps (`returns = 2L`)
#' @param parallelize a logical scalar. Set to TRUE if we want the for-loop for the individual-wise CTMC to be parallelized
#' @param which_step1 a character scalar. Either "100%" or "best". The former will use 100% of the individual CTCM effects to perform the K-means algorithm. \cr
#' The latter will compute the \eqn{log(P(Y))} value for all possible K-means result and use the thetas that yields the largest \eqn{log(P(Y))}.
#' @param parallel_optim a list object telling the function whether parallel process should be used for the Step 2 of the initial value generation. \cr
#' The list should contain **two** elements: \cr
#' (1) `run` a logical scalar, if TRUE then this function will use parallel processing. If FALSE, then the `cl` argument is ignored. \cr
#' (2) `cl` is an object obtained from the `parallel` package, for example \cr `cl = parallel::makeCluster(spec = 2)`
#' @param MyModelName a character scalar. Gives the current model fitting process a name. This name will be used when the function is logging the algorithm progress.
#'
#' @return A list containing the estimates obtained from step1 (K-means) or both steps (K-Means + direct optimization on observed data likelihood). \cr\cr
#' Step 1 of the algorithm outputs 3 elements:
#' \itemize{
#'   \item **step1** - a list of possible initial values, each element is the result of the K-means algorithm by trimming the individual estimate by the `pct_keep` argument
#'   \item **step1_full** - a single named numeric vector from the K-means using 100% of individual estimates (no trimming outliers, when `pct_keep = 1`)
#'   \item **step1_best** - a single named numeric vector that yielded the best observed log-likelihood value by testing trying each vector from `step1`
#' }
#' Step 2 of the algorithm outputs a named numeric vector obtained by maximizing the observed log-likelihood function (see **Note** section). This step is only performed is the argument `returns` is set to `2L`
#'
#' @note Step 1 is an individual level fitting algorithm and then uses K-means algorithm to perform clustering based on the individuals' effects
#' Step2 is equivalent to obtaining the MLE. It performs numerical optimization on the observed log-likelihood function:
#' \deqn{
#'   P(Y) = \sum_{i} log(P(Y_{i}))
#' }
#' where
#' \deqn{
#'   P(Y_{i}) = \sum_{k} P(Y_{i}=y_{i} \cap Z_{i}=k)
#' }
#' and the summand could be simplified via Bayes' Theorem:
#' \deqn{
#'   P(Y_{i} \cap Z_{i}) = P(Z_{i}) \cdot P(Y_{i} | Z_{i})
#' }
#'
#' @importFrom nnet multinom
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom foreach registerDoSEQ
#' @importFrom doParallel registerDoParallel
#'
#' @export
#'
#' @example inst/examples/ex_gen_inits_lctmc_3x3.R

gen_inits_lctmc_3x3 = function(df,
                               Xmat,
                               Wmat,
                               dt,
                               K,
                               N_sub,
                               pct_keep,
                               par_constraint,
                               returns = 2L,
                               parallelize,
                               which_step1,
                               parallel_optim,
                               MyModelName) {
  ### perform checks
  if (!all(pct_keep >= 0 & pct_keep <= 1)) {
    stop("`pct_keep` must be between 0 and 1, inclusive")
  }
  if (length(returns) != 1 || !is.integer(returns)) {
    stop("`returns` must be an integer scalar")
  }

  ### if "1" is not in `pct_keep`, append it
  pct_keep = append(pct_keep, 1)
  pct_keep = unique(pct_keep)

  ### load 'optimParallel'
  if (parallel_optim$run) {
    optim2 = optimParallel::optimParallel
  } else {
    optim2 = stats::optim
  }

  ### list of unique subject
  subject_list = unique(df$id)

  ### for Step 1, uses at max `N_sub` subjects
  set.seed(123)
  step1_N = min(N_sub, length(subject_list))
  randID_select = sample(seq_along(subject_list), size = step1_N, replace = FALSE)
  subject_list = subject_list[randID_select]

  ### STEP 1  ~~>  start
  cat("####~{", MyModelName, "}~", paste(rep("#", times = 100 - 8 - nchar(MyModelName)), sep = ""),"\n", sep = "")
  cat(paste(rep("#", times = 27), sep = ""), " Generating Initial Values -- Step 1 ", paste(rep("#", times = 36), sep = ""), "\n", sep = "")
  cat(" * using `N_sub` = ", step1_N, "/", length(unique(df$id)), " subjects for Step 1\n", sep = "")

  ### STEP 1  ~~>  estimate individual-wise CTMC
  if (!parallel_optim$run || !parallelize) {
    ## msg
    cat(" * starting non-parallelized individual fitting \n", sep = "")

    ## fit
    person_logQ.orig = indiv_ctmc_3x3(
      subject_list = subject_list,
      df = df,
      dt = dt,
      Xmat = Xmat,
      trace = TRUE
    )

  } else if (parallel_optim$run && parallelize) {
    ## msg
    cat(" * starting parallelized individual fitting \n", sep = "")

    ## register so `%dopar%` runs in parallel
    doParallel::registerDoParallel(cl = parallel_optim$cl)

    ## split `subject_id` into 300 groups
    n_groups = 300
    indiv_per_group = ceiling(length(subject_list)/n_groups)
    groups = rep(1:n_groups, each = indiv_per_group)
    if (length(groups) > length(subject_list)) {
      groups = groups[-(indiv_per_group * (1:(length(groups)-length(subject_list))))]
    }
    subject_list.split = split(x = subject_list, f = groups)

    ## run the 300 groups in parallel
    g = NULL
    person_logQ.orig = foreach::foreach(g = seq_along(subject_list.split), .inorder = FALSE, .combine = 'rbind') %dopar% {
      # ID list of current group
      subject_list.sub = subject_list.split[[g]]

      # fit
      logQ = indiv_ctmc_3x3(
        subject_list = subject_list.sub,
        df = df,
        dt = dt,
        Xmat = Xmat,
        trace = FALSE
      )

      # return
      logQ
    }

    ## reset foreach
    foreach::registerDoSEQ()

  } else {
    stop("Something went wrong at step 1 of initial value generation")
  }

  ### STEP 1  ~~>  loop through each specified `pct_keep`
  step1_out = list()
  for (i in seq_along(pct_keep)) {
    ## resets
    person_logQ = person_logQ.orig

    ## get min and max limits
    pk = pct_keep[i]
    limit.min = (1-pk) / 2
    limit.max = 1 - limit.min

    ## cut outlying estimates -- determined by `pct_keep`
    keep_obs = as.numeric(person_logQ[, 'convergence'] %in% c(0:1)) # 0 for converged, 1 for `maxit` reached
    for (c in 1:9) {
      # compute lower & upper quantile at current `pct_keep`
      lower_qnt = stats::quantile(person_logQ[, c], limit.min)
      upper_qnt = stats::quantile(person_logQ[, c], limit.max)

      # logical: which obs are within the limits
      within_qnts = (lower_qnt <= person_logQ[, c]) & (person_logQ[, c] <= upper_qnt)

      # update vector which is used later to truncate observations
      keep_obs = keep_obs & within_qnts
    }
    person_logQ = person_logQ[keep_obs, ]

    ## K-Means with `K` clusters
    set.seed(1)
    km = stats::kmeans(person_logQ[, 1:9], centers = K)

    ## re-arrange output such that most popular cluster always "1"
    km.cl = km$cluster
    km.order = order(-table(km.cl))
    for (k in 1:K) {
      km$cluster[km.cl == km.order[k]] = k
    }
    km$centers = km$centers[km.order, ]
    rownames(km$centers) = 1:K

    ## `Wmat.df` is used for estimating alpha parameters
    Wmat_sub = Wmat[randID_select, ]
    Wmat.df = data.frame(cbind(Wmat_sub[keep_obs, ], cl = km$cluster))

    ## if `par_constraint` is not NULL, then apply constrain, otherwise skip this step
    if (!is.null(par_constraint)) {
      # if alpha1.1 is constrained, then alpha1 for all classes will be constrained
      pc.par_names = names(par_constraint)[grepl(pattern = "\\.1$", names(par_constraint))]

      # loop through each ~ the sub(...) command extracts the covariate number e.g., alpha2.1 --> "2"
      for (pc in pc.par_names) {
        pc.data_name = paste("w", sub("(^alpha)(\\d)(\\.)(.+$)", "\\2", pc), sep = "")
        Wmat.df[, pc.data_name] = par_constraint[pc] * Wmat.df[, pc.data_name]
      }
    }

    ## multinomial model to determine alphas (`trace = FALSE` to suppress messages)
    multi_logistic = nnet::multinom(cl ~ w1 + w2, data = Wmat.df, trace = FALSE)
    multi_logistic.coef = stats::coef(multi_logistic)

    ## exception to handle `K=2L` case where stats::coef(multinom) output a vector instead of matrix
    if (K == 2) {
      multi_logistic.coef = as.matrix(multi_logistic.coef)
      multi_logistic.coef = t(multi_logistic.coef)
    }

    ## formatting before re-ordering
    colnames(multi_logistic.coef) = rownames(multi_logistic.coef) = NULL
    multi_logistic.coef = rbind(rep(0, ncol(multi_logistic.coef)), multi_logistic.coef)

    ## forces the K^th class to be referent in parametrization
    pi.alphas = c()
    for (k in 1:(K-1)) {
      pi.alphas = append(pi.alphas, multi_logistic.coef[k, ] - multi_logistic.coef[K, ])
    }

    ## create named vector for alpha parameters
    pi.alphas.names = expand.grid(c("alpha0", "alpha1", "alpha2"), 1:(K-1))
    pi.alphas.names = paste(pi.alphas.names$Var1, ".", pi.alphas.names$Var2, sep = "")

    ## put theta vector together for output: alpha coefficients
    names(pi.alphas) = pi.alphas.names
    theta = list(pi.alphas)

    ## put theta vector together for output: beta coefficients
    for (k in 1:K) {
      v = km$centers[k, ]
      names(v) = paste(names(v), "_", k, sep = "")
      theta[[k+1]] = v
    }

    ## append
    out_name = paste("pct_keep=", sprintf("%.3f", pk), sep = "")
    step1_out[[out_name]] = list(theta = unlist(theta), person_logQ = person_logQ)
  }
  step1_out[["person_logQ.orig"]] = person_logQ.orig

  ### STEP 1  ~~>  extract theta for the iteration where `pct_keep == 1`
  step1_out.full = step1_out[["pct_keep=1.000"]]$theta

  ### STEP 1  ~~>  determine the best option from Step 1
  constraint_index = names(par_constraint)[names(par_constraint) %in% names(step1_out[[1]]$theta)]
  theta.names.bik = gen_theta_names(K = K, type = "3x3", purpose = "bik")
  step1_out.log_PY = sapply(
    X = step1_out[-length(step1_out)],
    FUN = function(x) {
      ## apply constraint (NULLs will work)
      x.theta = x$theta
      x.theta[constraint_index] = par_constraint
      ## compute log(P(Y))
      y = bik_all_3x3(theta = x.theta, data = df, Xmat = Xmat, Wmat = Wmat, dt = dt, K = K, theta.names = theta.names.bik)
      ## `bi = y$bi1 + y$bi2 + ... y$biK`
      bi = Reduce(`+`, y)
      sum(log(bi))
    }
  )
  step1_out.best_index = which.max(step1_out.log_PY)
  step1_out.best = step1_out[[step1_out.best_index]]$theta

  ### STEP 1 ~~> constrain step 1 parameters
  step1_out.best[constraint_index] = par_constraint
  step1_out.full[constraint_index] = par_constraint

  ### STEP 1 ~~> print concluding messages
  if (!is.null(par_constraint)) {
    cat(" * ", "constrained parameters based on `par_constraint` \n", sep = "")
  } else {
    cat(" * ", "no parameter constraints provided, `par_constraint` is NULL \n", sep = "")
  }
  cat(" * best set of theta is when `pct_keep` = ", pct_keep[step1_out.best_index], "\n", sep = "")
  cat(" * log(P(Y)) = ", step1_out.log_PY[step1_out.best_index], " when evaluated at step1's best theta (`pct_keep`=", pct_keep[step1_out.best_index], ") \n\n", sep = "")

  ### STEP 1 ~~> exits. If specified to return only STEP 1, then exit here
  if (returns == 1L) {
    out = list(step1 = step1_out, step1_full = step1_out.full, step1_best = step1_out.best)
    class(out) = c("lctmc_3x3.inits", "list")
    return(out)
  }


  ### STEP 2  ~~>  initial value, use 'best' or '100%' (?)
  if (which_step1 == "100%" && length(which_step1) == 1) {
    step2_inits = step1_out.full
  } else if (which_step1 == "best" && length(which_step1) == 1) {
    step2_inits = step1_out.best
  } else {
    stop("ERROR: `which_step1` should be either '100%' or 'best'")
  }

  ### STEP 2  ~~>  start msg
  cat("####~{", MyModelName, "}~", paste(rep("#", times = 100 - 8 - nchar(MyModelName)), sep = ""),"\n", sep = "")
  cat(paste(rep("#", times = 27), sep = ""), " Generating Initial Values -- Step 2 ", paste(rep("#", times = 36), sep = ""), "\n", sep = "")
  cat(" * using `which_step1` = '", which_step1, "' as the intial values for Step 2 \n", sep = "")

  ### STEP 2  ~~>  optimization
  if (parallel_optim$run) {
    step2_out = optim2(
      par = step2_inits,
      function(par) {
        ## constrain parameters
        par[constraint_index] = par_constraint
        ## compute log(P(Y))
        y = bik_all_3x3(theta = par, data = df, Xmat = Xmat, Wmat = Wmat, dt = dt, K = K, theta.names = theta.names.bik)
        y = impute_bik(y, eps = 1e-3, EPS = 1e-24)
        ## `bi = y$bi1 + y$bi2 + ... + y$biK`
        bi = Reduce(`+`, y)
        ## return
        sum(log(bi))
      },
      method = "L-BFGS-B",
      control = list(fnscale = -length(dt), trace = 1, maxit = 1000, factr = 1e-5),
      parallel = list(cl = parallel_optim$cl, forward = FALSE, loginfo = FALSE)
    )
  }
  if (!parallel_optim$run) {
    step2_out = optim2(
      par = step2_inits,
      function(par) {
        ## constrain parameters
        par[constraint_index] = par_constraint
        ## compute log(P(Y))
        y = bik_all_3x3(theta = par, data = df, Xmat = Xmat, Wmat = Wmat, dt = dt, K = K, theta.names = theta.names.bik)
        y = impute_bik(y, eps = 1e-3, EPS = 1e-24)
        ## `bi = y$bi1 + y$bi2 + ... + y$biK`
        bi = Reduce(`+`, y)
        ## return
        sum(log(bi))
      },
      method = "L-BFGS-B",
      control = list(fnscale = -length(dt), trace = 1, maxit = 1000, factr = 1e-5)
    )
  }
  step2_out = step2_out$par
  cat("\n")

  ### STEP 2  ~~> exits
  if (returns == 2L) {
    out = list(
      step1 = step1_out,
      step1_full = step1_out.full,
      step1_best = step1_out.best,
      step2 = step2_out
    )
    class(out) = c("lctmc_3x3.inits", "list")
    return(out)
  }
}
