#' @title Step 1/2 for generating initial value for model fitting
#'
#' @description Performs step 1 of the initial value generation. This is a two-stage k-means procedure, where stage 1 involves fitting a
#' traditional CTMC model to each individual in the `df` data frame. Then the individual effects are pooled and fed into a k-means algorithm to split
#' the estimate into `K` clusters.
#'
#' @name gen_inits01_lctmc
#'
#' @param df a data frame object containing row-wise transition data as binary variables. \cr
#' For example, if `trans.2_1` equals 1 then it means the observation was a transition from stage 2 to stage 1 within `df_dt` amount of time.
#' @param df_Xmat a matrix object housing the covariates that affect the CTMC portion of the model. \cr
#' This matrix should have the same number of rows as the data frame object, `df`
#' @param df_Wmat a matrix object housing the covariates that affect the latent classification part of the model. \cr
#' This matrix should have number of rows equal to unique number of individuals in the data frame object, `df`
#' @param df_dt a numeric vector housing the length of time interval between observations.
#' This vector's length should be equal to number of rows in the data frame object, `df`
#' @param K an integer scalar. Use this variable to tell the function how many latent classes there should be. \cr
#' @param par_constraint See documentation in [lctmc_2x2()] or [lctmc_3x3()]
#' @param parallel_optim See documentation in [lctmc_2x2()] or [lctmc_3x3()]
#' @param N_sub a numeric scalar. This is used for step 1 of initial value generation where the algorithm fits a traditional CTMC model for each individual. \cr
#' Fitting the model for *all* individuals might have long run time without improvement in the accuracy of the estimation.
#' Thus, using this argument to set a maximum number of individuals to use for the initial value generation could save some computation time.
#' @param pct_keep a numeric vector where each element of the vector ranges from 0 to 1.
#' This argument controls what percentage of the individual effect should be used for the K-means algorithm for initial value generation.
#' The algorithm will consider all percentages specified in this vector. \cr
#' For example, for `pct_keep = c(0.5)`, after individuals effects are estimated, only the 25th to 75th percentile are fed into the K-Means algorithm to obtain cluster level estimates. \cr
#' Additionally, note that the threshold "1" is always appended to `pct_keep`, so it will always consider the 100% case.
#' @param parallelize a logical scalar. Set to TRUE if we want the for-loop for the individual-wise CTMC to be parallelized
#'
#' @return A list containing the estimates obtained from the K-means algorithm.
#' This step outputs 3 elements:
#' \itemize{
#'   \item **step1_full** - a list object containing the estimation from using 100% of individual effects (no trimming outliers, when `pct_keep = 1`). \cr
#'         First is the the theta estimates, and second is a 'kmeans' object from `kmeans()`.
#'   \item **step1_best** - a single named numeric vector that yielded the best observed log-likelihood value by testing trying each vector from `step1`
#' }
#'
#' @note This is the second step out of six of fitting a latent class CTMC model (i.e., generate initial values via a two stage k-means algorithm). \cr\cr
#' This is also part 1 of initial value generation,
#' which is an individual level fitting algorithm and then uses K-means method to perform clustering based on the individuals' effects. \cr
#' If using This method alone to model the data, it will very likely perform poorly.
#' Especially when there are only few observation per person, it becomes reliable when each individual has many observations. \cr
#' Regardless, this approach will still produce usable initial values for the next steps.
#'
#' @seealso [lctmc_2x2()]; [fmt_rowwise_trans()]; [gen_inits02_lctmc_2x2()]; [indiv_ctmc_2x2()]
#'
#' @importFrom nnet multinom
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom foreach registerDoSEQ
#' @importFrom doParallel registerDoParallel
#'
#' @example inst/examples/ex_running_lctmc.R
NULL

#' @rdname gen_inits01_lctmc
gen_inits01_lctmc_2x2 = function(df,
                                 df_Xmat,
                                 df_Wmat,
                                 df_dt,
                                 K,
                                 par_constraint,
                                 parallel_optim,
                                 N_sub,
                                 pct_keep,
                                 parallelize) {
  ### check specifications
  if (!all(pct_keep >= 0 & pct_keep <= 1)) {
    stop("`pct_keep` must be between 0 and 1, inclusive")
  }

  ### if "1" is not in `pct_keep`, append it
  pct_keep = append(pct_keep, 1)
  pct_keep = unique(pct_keep)

  ### list of unique subject
  subject_list = unique(df$id)

  ### for Step 1, uses at max `N_sub` subjects
  set.seed(123)
  step1_N = min(N_sub, length(subject_list))
  randID_select = sample(seq_along(subject_list), size = step1_N, replace = FALSE)
  subject_list = subject_list[randID_select]

  ### STEP 1  ~~>  start
  cat(" * using `N_sub` = ", step1_N, "/", length(unique(df$id)), " subjects for Step 1\n", sep = "")

  ### STEP 1  ~~>  estimate individual-wise CTMC
  if (!parallel_optim$run || !parallelize) {
    ## msg
    cat(" * starting non-parallelized individual fitting \n", sep = "")

    ## fit
    person_logQ.orig = indiv_ctmc_2x2(
      subject_list = subject_list,
      df = df,
      Xmat = df_Xmat,
      dt = df_dt,
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
      logQ = indiv_ctmc_2x2(
        subject_list = subject_list.sub,
        df = df,
        Xmat = df_Xmat,
        dt = df_dt,
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
    for (c in 1:6) {
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
    km = stats::kmeans(person_logQ[, 1:6], centers = K)

    ## re-arrange output such that most popular cluster always "1"
    km.cl = km$cluster
    km.order = order(-table(km.cl))
    for (k in 1:K) {
      km$cluster[km.cl == km.order[k]] = k
    }
    km$centers = km$centers[km.order, ]
    rownames(km$centers) = 1:K

    ## `Wmat_sub` is used for estimating alpha parameters
    Wmat_sub = df_Wmat[randID_select, ]
    Wmat_sub = data.frame(cbind(Wmat_sub[keep_obs, ], cl = km$cluster))

    ## if `par_constraint` is not NULL, then apply constrain, otherwise skip this step
    if (!is.null(par_constraint)) {
      # if alpha1.1 is constrained, then alpha1 for all classes will be constrained
      pc.par_names = names(par_constraint)[grepl(pattern = "\\.1$", names(par_constraint))]

      # loop through each ~ the sub(...) command extracts the covariate number e.g., alpha2.1 --> "2"
      for (pc in pc.par_names) {
        pc.data_name = paste("w", sub("(^alpha)(\\d)(\\.)(.+$)", "\\2", pc), sep = "")
        Wmat_sub[, pc.data_name] = par_constraint[pc] * Wmat_sub[, pc.data_name]
      }
    }

    ## multinomial model to determine alphas (`trace = FALSE` to suppress messages)
    multi_logistic = nnet::multinom(cl ~ w1 + w2, data = Wmat_sub, trace = FALSE)
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
    step1_out[[out_name]] = list(theta = unlist(theta), kmeans = km)
  }

  ### STEP 1  ~~>  extract theta for the iteration where `pct_keep == 1`
  step1_out.full = step1_out[[length(step1_out)]]$theta
  step1_out.full.kmeans = step1_out[[length(step1_out)]]$kmeans

  ### STEP 1  ~~>  determine the best option from Step 1
  constraint_index = names(par_constraint)[names(par_constraint) %in% names(step1_out[[1]]$theta)]
  theta.names.bik = gen_theta_names(K = K, type = "2x2", purpose = "bik")
  step1_out.log_PY = sapply(
    X = step1_out[-length(step1_out)],
    FUN = function(x) {
      ## apply constraint (NULLs will work)
      x.theta = x$theta
      x.theta[constraint_index] = par_constraint
      ## compute log(P(Y))
      y = bik_all_2x2(
        theta = x.theta,
        data = df,
        Xmat = df_Xmat,
        Wmat = df_Wmat,
        dt = df_dt,
        K = K,
        P.rs = FALSE,
        theta.names = theta.names.bik
      )
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
  step1_out.full.df = align_MLE_2x2(mle = step1_out.full, K = K)
  step1_out.full = step1_out.full.df$mle_theta
  names(step1_out.full) = step1_out.full.df$names

  ### STEP 1 ~~> print concluding messages
  if (!is.null(par_constraint)) {
    cat(" * ", "constrained parameters based on `par_constraint` \n", sep = "")
  } else {
    cat(" * ", "no parameter constraints provided, `par_constraint` is NULL \n", sep = "")
  }
  cat(" * best set of theta is when `pct_keep` = ", pct_keep[step1_out.best_index], "\n", sep = "")
  cat(" * log(P(Y)) = ", step1_out.log_PY[step1_out.best_index], " when evaluated at step1's best theta (`pct_keep`=", pct_keep[step1_out.best_index], ") \n", sep = "")

  ### STEP 1  ~~>  exits
  out = list(
    # step1 = step1_out, <-- takes up a lot of memory to save
    step1_full = list(theta = step1_out.full, kmeans = step1_out.full.kmeans),
    step1_best = step1_out.best
  )
  class(out) = append("lctmc_2x2.inits01", class(out))
  return(out)
}

#' @rdname gen_inits01_lctmc
gen_inits01_lctmc_3x3 = function(df,
                                 df_Xmat,
                                 df_Wmat,
                                 df_dt,
                                 K,
                                 par_constraint,
                                 parallel_optim,
                                 N_sub,
                                 pct_keep,
                                 parallelize) {
  ### check specifications
  if (!all(pct_keep >= 0 & pct_keep <= 1)) {
    stop("`pct_keep` must be between 0 and 1, inclusive")
  }

  ### if "1" is not in `pct_keep`, append it
  pct_keep = append(pct_keep, 1)
  pct_keep = unique(pct_keep)

  ### list of unique subject
  subject_list = unique(df$id)

  ### for Step 1, uses at max `N_sub` subjects
  set.seed(123)
  step1_N = min(N_sub, length(subject_list))
  randID_select = sample(seq_along(subject_list), size = step1_N, replace = FALSE)
  subject_list = subject_list[randID_select]

  ### STEP 1  ~~>  start
  cat(" * using `N_sub` = ", step1_N, "/", length(unique(df$id)), " subjects for Step 1\n", sep = "")

  ### STEP 1  ~~>  estimate individual-wise CTMC
  if (!parallel_optim$run || !parallelize) {
    ## msg
    cat(" * starting non-parallelized individual fitting \n", sep = "")

    ## fit
    person_logQ.orig = indiv_ctmc_3x3(
      subject_list = subject_list,
      df = df,
      Xmat = df_Xmat,
      dt = df_dt,
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
        Xmat = df_Xmat,
        dt = df_dt,
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

    ## `Wmat_sub` is used for estimating alpha parameters
    Wmat_sub = df_Wmat[randID_select, ]
    Wmat_sub = data.frame(cbind(Wmat_sub[keep_obs, ], cl = km$cluster))

    ## if `par_constraint` is not NULL, then apply constrain, otherwise skip this step
    if (!is.null(par_constraint)) {
      # if alpha1.1 is constrained, then alpha1 for all classes will be constrained
      pc.par_names = names(par_constraint)[grepl(pattern = "\\.1$", names(par_constraint))]

      # loop through each ~ the sub(...) command extracts the covariate number e.g., alpha2.1 --> "2"
      for (pc in pc.par_names) {
        pc.data_name = paste("w", sub("(^alpha)(\\d)(\\.)(.+$)", "\\2", pc), sep = "")
        Wmat_sub[, pc.data_name] = par_constraint[pc] * Wmat_sub[, pc.data_name]
      }
    }

    ## multinomial model to determine alphas (`trace = FALSE` to suppress messages)
    multi_logistic = nnet::multinom(cl ~ w1 + w2, data = Wmat_sub, trace = FALSE)
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
    step1_out[[out_name]] = list(theta = unlist(theta), kmeans = km)
  }

  ### STEP 1  ~~>  extract theta for the iteration where `pct_keep == 1`
  step1_out.full = step1_out[[length(step1_out)]]$theta
  step1_out.full.kmeans = step1_out[[length(step1_out)]]$kmeans

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
      y = bik_all_3x3(
        theta = x.theta,
        data = df,
        Xmat = df_Xmat,
        Wmat = df_Wmat,
        dt = df_dt,
        K = K,
        P.rs = FALSE,
        theta.names = theta.names.bik
      )
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
  step1_out.full.df = align_MLE_3x3(mle = step1_out.full, K = K)
  step1_out.full = step1_out.full.df$mle_theta
  names(step1_out.full) = step1_out.full.df$names

  ### STEP 1 ~~> print concluding messages
  if (!is.null(par_constraint)) {
    cat(" * ", "constrained parameters based on `par_constraint` \n", sep = "")
  } else {
    cat(" * ", "no parameter constraints provided, `par_constraint` is NULL \n", sep = "")
  }
  cat(" * best set of theta is when `pct_keep` = ", pct_keep[step1_out.best_index], "\n", sep = "")
  cat(" * log(P(Y)) = ", step1_out.log_PY[step1_out.best_index], " when evaluated at step1's best theta (`pct_keep`=", pct_keep[step1_out.best_index], ") \n\n", sep = "")

  ### STEP 1  ~~>  exits
  out = list(
    # step1 = step1_out, <-- takes up a lot of memory to save
    step1_full = list(theta = step1_out.full, kmeans = step1_out.full.kmeans),
    step1_best = step1_out.best
  )
  class(out) = append("lctmc_3x3.inits01", class(out))
  return(out)
}
